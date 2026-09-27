#pragma once
// Device solve of the H2 Color factorization (GPU backend).
//
// The factors the host solve reads (T, X_SR, X_NR, the LU of X_RR and its
// pivots, index lists) are uploaded once, at the first solve after a
// factorization, into the device heap, where they stay for later solves.
// Whether the device solve runs is decided collectively: every rank must be
// able to hold its factors (and every level must be a Color level with an
// LU-factored X_RR), else all ranks keep the host solve.
//
// A solve runs each level's Color waves on the device:
//   forward   x_R -= T^T x_S;  x_S += X_SR x_R;  neighbors += X_NR x_R
//             (neighbor updates of a wave summed in box order and applied
//             after it; those of other ranks' boxes sent after each wave)
//   diagonal  x_R = X_RR^{-1} x_R, right after the level's forward waves
//             (it touches only redundant entries, final by then)
//   backward  x_R += X_SR^T x_S + X_NR^T x_N;  x_S -= T x_R, waves reversed,
//             other ranks' neighbor vectors refreshed before each wave
// The hand-off between levels (gather_skeleton_to_parent,
// scatter_solution_to_children, with the process reductions) and the root
// solve stay on the host: a level's vectors go to the device and back once
// per sweep.  The arithmetic matches the host solve up to rounding (its
// sums of neighbor updates follow the host's thread chunks).
//
// H2_GPU_SOLVE=0 keeps the host solve.

#ifdef H2_HAVE_GPU

#include "device_heap.hpp"
#include "gpu_runtime.hpp"
#include "solve_kernels.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace fmm {
namespace gpu {

// H2_GPU_SOLVE_KEEP=0: never keep factors on the device during the
// factorization (all are uploaded at the first solve).
inline bool device_solve_keep() {
    static const bool keep = [] {
        const char* v = std::getenv("H2_GPU_SOLVE_KEEP");
        return v == nullptr || std::atoi(v) != 0;
    }();
    return keep;
}
// A level keeps its factors while the heap stays under this fraction of its
// capacity (H2_GPU_SOLVE_KEEP_FRACTION, default 0.6), which leaves room for
// the level's own peak.
inline double device_solve_keep_fraction() {
    static const double fraction = [] {
        const char* v = std::getenv("H2_GPU_SOLVE_KEEP_FRACTION");
        const double f = v != nullptr ? std::atof(v) : 0.6;
        return f > 0.0 && f <= 1.0 ? f : 0.6;
    }();
    return fraction;
}

// H2_GPU_SOLVE_CHECK=1: also run the host solve and print the difference.
inline bool device_solve_check() {
    static const bool check = [] {
        const char* v = std::getenv("H2_GPU_SOLVE_CHECK");
        return v != nullptr && std::atoi(v) != 0;
    }();
    return check;
}
// The host solve of a check runs with the device solve suspended.
inline bool& device_solve_suspended() {
    static bool suspended = false;
    return suspended;
}

inline bool device_solve_enabled() {
    static const bool enabled = [] {
        const char* v = std::getenv("H2_GPU_SOLVE");
        return v == nullptr || std::atoi(v) != 0;
    }();
    return enabled && color_gpu_enabled();
}

// One level of the device solve (on this rank): factors, and the launch
// and message tables of its sweeps (built once, at the first solve).
struct DeviceSolveLevel {
    bool active = false;
    std::vector<char*> blocks;          // factors and tables (heap)
    std::vector<SolveBox> boxes;        // local boxes, device pointers
    std::vector<SolveSlot> slots;       // their one-hop slots
    std::vector<int64_t> slot_morton;   // the neighbor of each slot
    int64_t points = 0;                 // level vector length (per right-hand side)
    int64_t work_points = 0;            // largest wave's X_NR rows
    int max_n = 0, max_r = 0;
    // other ranks' boxes one hop from local ones (multi-rank levels)
    std::vector<int> peers;
    struct Ghost { int64_t morton; int owner; int n; int64_t offset; };
    std::vector<Ghost> ghosts;
    std::unordered_map<int64_t, int> ghost_of;
    int64_t ghost_points = 0;
    std::vector<std::vector<int>> readers;  // per local box: peers (indices) that read its vector
    // A Color wave.  Message sizes are per peer, in points.
    struct Wave {
        int order0 = 0, count = 0;           // its boxes (forward order; reversed for the backward sweep)
        int max_ntot = 0;                    // largest X_NR row count of its boxes
        int item0 = 0, nitems = 0, max_count = 0;        // forward accumulation (local, outbox)
        int in_item0 = 0, in_nitems = 0, in_max_count = 0;  // updates from other ranks
        std::vector<int64_t> out, in;        // forward messages
        // other ranks' vectors refreshed before this wave: ref[0] for waves in
        // reverse (the solve's backward sweep), ref[1] in order (the multiply's
        // forward sweep)
        struct Refresh {
            int send0 = 0, nsend = 0, send_max_n = 0;
            int recv0 = 0, nrecv = 0, recv_max_n = 0;
            std::vector<int64_t> send, recv;
        };
        Refresh ref[2];
    };
    std::vector<Wave> waves;
    int ndiag = 0;
    int64_t out_max = 0, in_max = 0, send_max = 0, recv_max = 0;  // largest messages (points, all peers)
    // device tables
    const SolveBox* d_boxes = nullptr;
    const SolveSlot* d_slots = nullptr;
    const int* d_order = nullptr;       // forward order, wave after wave
    const int* d_rorder = nullptr;      // backward order (each wave reversed)
    const int* d_diag = nullptr;
    const SolveAccum* d_items = nullptr;
    const SolvePart* d_parts = nullptr;
    const SolveCopy* d_copies = nullptr;
};

// Factors of a level kept on the device during the factorization (all ranks
// of the level agreed), consumed by the first solve's preparation.
struct KeptSolveBox {
    const void* T = nullptr;  // elements of the factorization's type
    const void* xsr = nullptr;
    const void* xnr = nullptr;
    const void* lu = nullptr;
    const int* ipiv = nullptr;
    int k = 0, r = 0, ntot = 0;
};
struct KeptSolveLevel {
    std::vector<char*> blocks;
    std::unordered_map<int64_t, KeptSolveBox> boxes;  // Morton index -> factors
    double bytes = 0.0;
};

struct DeviceSolveStore {
    std::map<int, KeptSolveLevel> kept;  // level -> factors kept during the factorization
    const void* tree = nullptr;
    bool decided = false;
    bool usable = false;
    std::vector<DeviceSolveLevel> levels;
    double bytes = 0.0;
    PinnedBuffer host_vec, host_out, host_in;  // staging of the solves (kept: pinning is slow)
    bool empty() const { return levels.empty() && kept.empty(); }

    // A preparation's tables and uploaded factors (the kept factors stay).
    void release_prepared() {
        DeviceHeap& heap = DeviceHeap::instance();
        for (auto& lv : levels) {
            for (char* b : lv.blocks) heap.free(b);
        }
        levels.clear();
        tree = nullptr;
        decided = false;
        usable = false;
        bytes = 0.0;
    }
    void release() {
        release_prepared();
        DeviceHeap& heap = DeviceHeap::instance();
        for (auto& [level, kl] : kept) {
            for (char* b : kl.blocks) heap.free(b);
        }
        kept.clear();
        tree = nullptr;
        decided = false;
        usable = false;
        bytes = 0.0;
    }
};

inline DeviceSolveStore& device_solve_store() {
    static DeviceSolveStore store;
    return store;
}

// A factorization replaces the factors: drop the device copies.
inline void invalidate_device_solve() { device_solve_store().release(); }

// ---------------------------------------------------------------------------
// Messages with the one-hop neighbor ranks: a header of int64 and a payload
// of doubles per peer, in two rounds (sizes, then all messages at once).
struct PeerMessage {
    std::vector<int64_t> head;
    std::vector<double> data;
};

inline double exchange_peer_messages(MPI_Comm comm, const std::vector<int>& peers, const std::vector<PeerMessage>& out,
                                     std::vector<PeerMessage>& in, int tag) {
    const auto t0 = std::chrono::steady_clock::now();
    const size_t np = peers.size();
    in.assign(np, PeerMessage{});
    std::vector<uint64_t> sizes_out(2 * np), sizes_in(2 * np, 0);
    for (size_t i = 0; i < np; ++i) {
        sizes_out[2 * i] = out[i].head.size();
        sizes_out[2 * i + 1] = out[i].data.size();
    }
    std::vector<MPI_Request> reqs;
    reqs.reserve(4 * np);
    for (size_t i = 0; i < np; ++i) {
        reqs.emplace_back();
        MPI_Irecv(&sizes_in[2 * i], 2, MPI_UINT64_T, peers[i], tag, comm, &reqs.back());
    }
    for (size_t i = 0; i < np; ++i) {
        reqs.emplace_back();
        MPI_Isend(&sizes_out[2 * i], 2, MPI_UINT64_T, peers[i], tag, comm, &reqs.back());
    }
    MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
    reqs.clear();
    for (size_t i = 0; i < np; ++i) {
        in[i].head.resize(sizes_in[2 * i]);
        in[i].data.resize(sizes_in[2 * i + 1]);
        if (!in[i].head.empty()) {
            reqs.emplace_back();
            MPI_Irecv(in[i].head.data(), static_cast<int>(in[i].head.size()), MPI_INT64_T, peers[i], tag + 1, comm,
                      &reqs.back());
        }
        if (!in[i].data.empty()) {
            reqs.emplace_back();
            MPI_Irecv(in[i].data.data(), static_cast<int>(in[i].data.size()), MPI_DOUBLE, peers[i], tag + 2, comm,
                      &reqs.back());
        }
    }
    for (size_t i = 0; i < np; ++i) {
        if (!out[i].head.empty()) {
            reqs.emplace_back();
            MPI_Isend(out[i].head.data(), static_cast<int>(out[i].head.size()), MPI_INT64_T, peers[i], tag + 1, comm,
                      &reqs.back());
        }
        if (!out[i].data.empty()) {
            reqs.emplace_back();
            MPI_Isend(out[i].data.data(), static_cast<int>(out[i].data.size()), MPI_DOUBLE, peers[i], tag + 2, comm,
                      &reqs.back());
        }
    }
    MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// ---------------------------------------------------------------------------
// Preparation (collective over the tree's ranks)
// ---------------------------------------------------------------------------
namespace solve_detail {

template<typename CoordType, typename DataType>
int owner_of(ParallelTree<CoordType, DataType>* tree, TreeLevel<CoordType, DataType>& lvl, int level, int64_t morton) {
    const std::vector<uint64_t> single = {static_cast<uint64_t>(morton)};
    const std::vector<uint32_t> region =
        morton::assign_to_processes_nd(tree->dimension, single, lvl.num_active_processes, 1u << level);
    return lvl.morton_to_rank.at(static_cast<int>(region[0]));
}

// Color waves as the host solve bins them: boundary boxes by color, then
// interior boxes by color, Morton order within a wave.
template<typename CoordType, typename DataType>
std::vector<std::vector<int>> color_waves(const TreeLevel<CoordType, DataType>& lvl, int dimension) {
    const int num_colors = 1 << dimension;
    std::vector<std::vector<int>> waves(static_cast<size_t>(2 * num_colors));
    for (size_t b = 0; b < lvl.local_boxes.size(); ++b) {
        const int64_t morton = lvl.local_morton_start + static_cast<int64_t>(b);
        const int color = static_cast<int>(morton & (num_colors - 1));
        waves[static_cast<size_t>(color + (lvl.local_boxes[b].on_boundary ? 0 : num_colors))].push_back(
            static_cast<int>(b));
    }
    return waves;
}

// Per-box byte layout of the factors in a device block.
struct BoxLayout {
    size_t T = 0, xsr = 0, xnr = 0, lu = 0, ipiv = 0, skel = 0, red = 0, bytes = 0;
    bool has_T = false, has_xsr = false;
};

}  // namespace solve_detail

// The launch and message tables of a level's sweeps.  Each wave's
// accumulations (local neighbors, and per-peer blocks of the outbox) and
// vector refreshes are fixed by the factorization, so they are built once;
// one exchange with the neighbor ranks gives each rank the layout of what
// it will receive.  Returns the bytes of the tables on the device.
template<typename CoordType, typename DataType>
size_t build_solve_tables(ParallelTree<CoordType, DataType>* tree, TreeLevel<CoordType, DataType>& lvl, int level,
                          DeviceSolveLevel& L, const std::vector<int>& remote_lists,
                          const std::unordered_map<int64_t, size_t>& remote_list_at, double& comm_seconds) {
    const size_t I = sizeof(int);
    const size_t nb = L.boxes.size();
    const size_t np = L.peers.size();
    auto is_local = [&](int64_t m) {
        return m >= lvl.local_morton_start && m < lvl.local_morton_start + static_cast<int64_t>(nb);
    };
    std::unordered_map<int, size_t> peer_index;
    for (size_t i = 0; i < np; ++i) peer_index.emplace(L.peers[i], i);
    auto local_box = [&](int64_t m) -> const SolveBox& { return L.boxes[static_cast<size_t>(m - lvl.local_morton_start)]; };

    // waves (eliminated boxes only), work offsets
    const auto bins = solve_detail::color_waves(lvl, tree->dimension);
    std::vector<int> order, rorder, diag;
    L.waves.assign(bins.size(), DeviceSolveLevel::Wave{});
    for (size_t w = 0; w < bins.size(); ++w) {
        auto& W = L.waves[w];
        W.order0 = static_cast<int>(order.size());
        int64_t work = 0;
        for (int b : bins[w]) {
            SolveBox& sb = L.boxes[static_cast<size_t>(b)];
            if (sb.k == 0 || sb.r == 0) continue;
            sb.work = work;
            work += sb.ntot;
            W.max_ntot = std::max(W.max_ntot, sb.ntot);
            order.push_back(b);
        }
        W.count = static_cast<int>(order.size()) - W.order0;
        rorder.insert(rorder.end(), order.rbegin(), order.rbegin() + W.count);
        L.work_points = std::max(L.work_points, work);
    }
    for (size_t b = 0; b < nb; ++b) {
        if (L.boxes[b].k > 0 && L.boxes[b].r > 0) diag.push_back(static_cast<int>(b));
    }
    L.ndiag = static_cast<int>(diag.size());

    // forward: each wave's targets, in first-contribution order
    std::vector<SolveAccum> items;
    std::vector<SolvePart> parts;
    std::vector<std::vector<int64_t>> head(np);  // to each peer: (wave, morton, full, count)... then (wave, morton)...
    for (size_t w = 0; w < L.waves.size(); ++w) {
        auto& W = L.waves[w];
        struct Target { int64_t morton; int count; int full; std::vector<SolvePart> parts; };
        std::vector<Target> targets;
        std::unordered_map<int64_t, size_t> target_of;
        for (int i = W.order0; i < W.order0 + W.count; ++i) {
            const SolveBox& sb = L.boxes[static_cast<size_t>(order[static_cast<size_t>(i)])];
            for (int a = 0; a < sb.nslots; ++a) {
                const size_t s = static_cast<size_t>(sb.slot0 + a);
                const SolveSlot& sl = L.slots[s];
                const int64_t m = L.slot_morton[s];
                auto it = target_of.find(m);
                if (it == target_of.end()) {
                    it = target_of.emplace(m, targets.size()).first;
                    targets.push_back(Target{m, sl.count, sl.full, {}});
                }
                Target& t = targets[it->second];
                if (t.count != sl.count || t.full != sl.full) throw std::runtime_error("device solve: inconsistent updates");
                t.parts.push_back(SolvePart{sb.work, sl.row0, sb.ntot, 0});
            }
        }
        W.out.assign(np, 0);
        std::vector<std::pair<size_t, int64_t>> remote(targets.size(), {0, -1});  // peer, offset in its block
        for (size_t t = 0; t < targets.size(); ++t) {
            if (is_local(targets[t].morton)) continue;
            auto pit = peer_index.find(solve_detail::owner_of(tree, lvl, level, targets[t].morton));
            if (pit == peer_index.end()) throw std::runtime_error("device solve: update for a non-neighbor rank");
            remote[t] = {pit->second, W.out[pit->second]};
            W.out[pit->second] += targets[t].count;
            auto& h = head[pit->second];
            h.insert(h.end(), {static_cast<int64_t>(w), targets[t].morton, targets[t].full, targets[t].count});
        }
        std::vector<int64_t> prefix(np + 1, 0);
        for (size_t q = 0; q < np; ++q) prefix[q + 1] = prefix[q] + W.out[q];
        L.out_max = std::max(L.out_max, prefix[np]);
        W.item0 = static_cast<int>(items.size());
        for (size_t t = 0; t < targets.size(); ++t) {
            const Target& tg = targets[t];
            SolveAccum it{};
            it.count = tg.count;
            it.nparts = static_cast<int>(tg.parts.size());
            it.part0 = static_cast<int64_t>(parts.size());
            parts.insert(parts.end(), tg.parts.begin(), tg.parts.end());
            if (remote[t].second < 0) {
                const SolveBox& target = local_box(tg.morton);
                it.dst = target.vec;
                it.ld = target.n;
                it.rows = tg.full ? nullptr : target.skel;
                it.outbox = 0;
            } else {
                it.dst = prefix[remote[t].first] + remote[t].second;
                it.ld = tg.count;
                it.rows = nullptr;
                it.outbox = 1;
            }
            W.max_count = std::max(W.max_count, tg.count);
            items.push_back(it);
        }
        W.nitems = static_cast<int>(items.size()) - W.item0;
    }

    // Before each wave of a sweep that reads neighbors, the vectors other
    // ranks read that changed since the last refresh (all of them before the
    // sweep's first wave): waves in reverse (kind 0, the solve's backward
    // sweep) or in order (kind 1, the multiply's forward sweep).
    std::vector<SolveCopy> copies;
    const int nw = static_cast<int>(L.waves.size());
    std::vector<std::vector<int64_t>> head_ref[2] = {std::vector<std::vector<int64_t>>(np),
                                                     std::vector<std::vector<int64_t>>(np)};
    for (int kind = 0; kind < 2; ++kind) {
        for (int step = 0; step < nw; ++step) {
            const int w = kind == 0 ? nw - 1 - step : step;
            auto& R = L.waves[static_cast<size_t>(w)].ref[kind];
            R.send.assign(np, 0);
            std::vector<int> changed;
            if (step == 0) {
                for (size_t b = 0; b < nb; ++b) changed.push_back(static_cast<int>(b));
            } else {
                const auto& P = L.waves[static_cast<size_t>(kind == 0 ? w + 1 : w - 1)];
                changed.assign(order.begin() + P.order0, order.begin() + P.order0 + P.count);
                std::sort(changed.begin(), changed.end());
            }
            std::vector<std::vector<int>> per_peer(np);
            for (int b : changed) {
                for (int q : L.readers[static_cast<size_t>(b)]) per_peer[static_cast<size_t>(q)].push_back(b);
            }
            std::vector<int64_t> prefix(np + 1, 0);
            for (size_t q = 0; q < np; ++q) {
                for (int b : per_peer[q]) R.send[q] += L.boxes[static_cast<size_t>(b)].n;
                prefix[q + 1] = prefix[q] + R.send[q];
            }
            L.send_max = std::max(L.send_max, prefix[np]);
            R.send0 = static_cast<int>(copies.size());
            for (size_t q = 0; q < np; ++q) {
                int64_t at = prefix[q];
                for (int b : per_peer[q]) {
                    const SolveBox& sb = L.boxes[static_cast<size_t>(b)];
                    copies.push_back(SolveCopy{sb.vec, at, sb.n});
                    at += sb.n;
                    R.send_max_n = std::max(R.send_max_n, sb.n);
                    head_ref[kind][q].insert(head_ref[kind][q].end(), {static_cast<int64_t>(w), lvl.local_morton_start + b});
                }
            }
            R.nsend = static_cast<int>(copies.size()) - R.send0;
        }
    }

    // what the neighbor ranks send us, and where it goes
    if (np > 0) {
        std::vector<PeerMessage> out(np), in;
        for (size_t q = 0; q < np; ++q) {
            auto& h = out[q].head;
            h.push_back(static_cast<int64_t>(head[q].size() / 4));
            h.insert(h.end(), head[q].begin(), head[q].end());
            for (int kind = 0; kind < 2; ++kind) {
                h.push_back(static_cast<int64_t>(head_ref[kind][q].size() / 2));
                h.insert(h.end(), head_ref[kind][q].begin(), head_ref[kind][q].end());
            }
        }
        comm_seconds += exchange_peer_messages(tree->comm, L.peers, out, in, 746);
        struct Fwd { int64_t morton; int full, count; };
        std::vector<std::vector<std::vector<Fwd>>> fwd(L.waves.size(), std::vector<std::vector<Fwd>>(np));
        std::vector<std::vector<std::vector<int64_t>>> bwd[2];
        for (auto& v : bwd) v.assign(L.waves.size(), std::vector<std::vector<int64_t>>(np));
        for (size_t q = 0; q < np; ++q) {
            const auto& h = in[q].head;
            size_t at = 0;
            const int64_t nf = h.at(at++);
            for (int64_t e = 0; e < nf; ++e, at += 4) {
                fwd.at(static_cast<size_t>(h.at(at)))[q].push_back(
                    Fwd{h.at(at + 1), static_cast<int>(h.at(at + 2)), static_cast<int>(h.at(at + 3))});
            }
            for (int kind = 0; kind < 2; ++kind) {
                const int64_t nbk = h.at(at++);
                for (int64_t e = 0; e < nbk; ++e, at += 2) {
                    bwd[kind].at(static_cast<size_t>(h.at(at)))[q].push_back(h.at(at + 1));
                }
            }
            if (at != h.size()) throw std::runtime_error("device solve: malformed layout message");
        }
        for (size_t w = 0; w < L.waves.size(); ++w) {
            auto& W = L.waves[w];
            // forward: updates of local boxes, parts in peer order
            W.in.assign(np, 0);
            for (size_t q = 0; q < np; ++q) {
                for (const Fwd& f : fwd[w][q]) W.in[q] += f.count;
            }
            std::vector<int64_t> prefix(np + 1, 0);
            for (size_t q = 0; q < np; ++q) prefix[q + 1] = prefix[q] + W.in[q];
            L.in_max = std::max(L.in_max, prefix[np]);
            struct Target { int64_t morton; int count; int full; std::vector<SolvePart> parts; };
            std::vector<Target> targets;
            std::unordered_map<int64_t, size_t> target_of;
            for (size_t q = 0; q < np; ++q) {
                int64_t at = prefix[q];
                for (const Fwd& f : fwd[w][q]) {
                    if (!is_local(f.morton)) throw std::runtime_error("device solve: update for a box of another rank");
                    auto it = target_of.find(f.morton);
                    if (it == target_of.end()) {
                        it = target_of.emplace(f.morton, targets.size()).first;
                        targets.push_back(Target{f.morton, f.count, f.full, {}});
                    }
                    targets[it->second].parts.push_back(SolvePart{at, 0, f.count, 1});
                    at += f.count;
                }
            }
            W.in_item0 = static_cast<int>(items.size());
            for (const Target& t : targets) {
                const SolveBox& target = local_box(t.morton);
                if (t.full ? t.count != target.n : t.count != target.k) {
                    throw std::runtime_error("device solve: update of box " + std::to_string(t.morton) + " has another size");
                }
                SolveAccum it{};
                it.dst = target.vec;
                it.ld = target.n;
                it.count = t.count;
                it.rows = t.full ? nullptr : target.skel;
                it.outbox = 0;
                it.nparts = static_cast<int>(t.parts.size());
                it.part0 = static_cast<int64_t>(parts.size());
                parts.insert(parts.end(), t.parts.begin(), t.parts.end());
                W.in_max_count = std::max(W.in_max_count, t.count);
                items.push_back(it);
            }
            W.in_nitems = static_cast<int>(items.size()) - W.in_item0;
            // refreshed vectors into the ghost area, for either order
            for (int kind = 0; kind < 2; ++kind) {
                auto& R = W.ref[kind];
                R.recv.assign(np, 0);
                R.recv0 = static_cast<int>(copies.size());
                for (size_t q = 0; q < np; ++q) {
                    for (int64_t m : bwd[kind][w][q]) R.recv[q] += L.ghosts[static_cast<size_t>(L.ghost_of.at(m))].n;
                }
                prefix.assign(np + 1, 0);
                for (size_t q = 0; q < np; ++q) prefix[q + 1] = prefix[q] + R.recv[q];
                L.recv_max = std::max(L.recv_max, prefix[np]);
                for (size_t q = 0; q < np; ++q) {
                    int64_t at = prefix[q];
                    for (int64_t m : bwd[kind][w][q]) {
                        const auto& g = L.ghosts[static_cast<size_t>(L.ghost_of.at(m))];
                        copies.push_back(SolveCopy{at, g.offset, g.n});
                        at += g.n;
                        R.recv_max_n = std::max(R.recv_max_n, g.n);
                    }
                }
                R.nrecv = static_cast<int>(copies.size()) - R.recv0;
            }
        }
    }

    // ---- tables on the device
    const size_t t_lists = 0;
    const size_t t_slots = align_up(t_lists + remote_lists.size() * I);
    const size_t t_boxes = align_up(t_slots + L.slots.size() * sizeof(SolveSlot));
    const size_t t_order = align_up(t_boxes + nb * sizeof(SolveBox));
    const size_t t_rorder = align_up(t_order + order.size() * I);
    const size_t t_diag = align_up(t_rorder + rorder.size() * I);
    const size_t t_items = align_up(t_diag + diag.size() * I);
    const size_t t_parts = align_up(t_items + items.size() * sizeof(SolveAccum));
    const size_t t_copies = align_up(t_parts + parts.size() * sizeof(SolvePart));
    const size_t t_end = align_up(t_copies + copies.size() * sizeof(SolveCopy));
    char* tables = DeviceHeap::instance().alloc_resident(std::max<size_t>(t_end, 1));
    L.blocks.push_back(tables);
    for (size_t s = 0; s < L.slots.size(); ++s) {
        SolveSlot& sl = L.slots[s];
        if (sl.ghost && !sl.full) {
            sl.skel = reinterpret_cast<const int*>(tables + t_lists) + remote_list_at.at(L.slot_morton[s]);
        }
    }
    std::vector<char> image(t_end, 0);
    auto put = [&](size_t at, const void* data, size_t bytes) {
        if (bytes > 0) std::memcpy(image.data() + at, data, bytes);
    };
    put(t_lists, remote_lists.data(), remote_lists.size() * I);
    put(t_slots, L.slots.data(), L.slots.size() * sizeof(SolveSlot));
    put(t_boxes, L.boxes.data(), nb * sizeof(SolveBox));
    put(t_order, order.data(), order.size() * I);
    put(t_rorder, rorder.data(), rorder.size() * I);
    put(t_diag, diag.data(), diag.size() * I);
    put(t_items, items.data(), items.size() * sizeof(SolveAccum));
    put(t_parts, parts.data(), parts.size() * sizeof(SolvePart));
    put(t_copies, copies.data(), copies.size() * sizeof(SolveCopy));
    check_cuda(cudaMemcpy(tables, image.data(), t_end, cudaMemcpyHostToDevice), "solve tables");
    L.d_slots = reinterpret_cast<const SolveSlot*>(tables + t_slots);
    L.d_boxes = reinterpret_cast<const SolveBox*>(tables + t_boxes);
    L.d_order = reinterpret_cast<const int*>(tables + t_order);
    L.d_rorder = reinterpret_cast<const int*>(tables + t_rorder);
    L.d_diag = reinterpret_cast<const int*>(tables + t_diag);
    L.d_items = reinterpret_cast<const SolveAccum*>(tables + t_items);
    L.d_parts = reinterpret_cast<const SolvePart*>(tables + t_parts);
    L.d_copies = reinterpret_cast<const SolveCopy*>(tables + t_copies);
    return t_end;
}

template<typename CoordType, typename DataType>
bool prepare_device_solve(ParallelTree<CoordType, DataType>* tree, int verbosity) {
    DeviceSolveStore& store = device_solve_store();
    if (store.decided && store.tree == tree) return store.usable;
    store.release_prepared();
    store.tree = tree;
    store.decided = true;
    {
        using clock = std::chrono::steady_clock;
        const auto t0 = clock::now();
        const size_t D = sizeof(DataType), I = sizeof(int);
        const int leaf = tree->num_levels - 1;
        int rank = 0;
        MPI_Comm_rank(tree->comm, &rank);

        // ---- eligibility and size, on every rank
        std::string reason;
        double need = 0.0;
        for (int level = 2; level <= leaf && reason.empty(); ++level) {
            auto& lvl = tree->levels[static_cast<size_t>(level)];
            if (!lvl.is_process_active) continue;
            if (tree->level_uses_CA(level)) {
                reason = "CA levels";
                break;
            }
            const bool kept_level = store.kept.count(level) != 0;
            for (const auto& box : lvl.local_boxes) {
                const int64_t n = box.num_points, k = static_cast<int64_t>(box.skeleton_indices.size()),
                              r = static_cast<int64_t>(box.redundant_indices.size());
                if (static_cast<size_t>(n) * D > kSolveSharedLimit) {
                    reason = "boxes too large";
                    break;
                }
                if (r > 0 && k > 0 && box.X_RR.format != MatrixStorage<DataType>::LU_FACTORED) {
                    reason = "X_RR not LU-factored (use H2_XRR_factor=1)";
                    break;
                }
                const int64_t ntot = box.X_NR.is_allocated() ? box.X_NR.rows : 0;
                need += static_cast<double>((kept_level ? 0 : (2 * k * r + ntot * r + r * r) * D) + (r + k + r) * I + 1024);
            }
        }
        DeviceHeap& heap = DeviceHeap::instance();
        heap.ensure_initialized();
        const double room = static_cast<double>(heap.capacity() - heap.used());
        int ok = reason.empty() && need * 1.05 + (size_t{256} << 20) < room ? 1 : 0;
        if (reason.empty() && !ok) reason = "not enough device memory";
        int all_ok = 0;
        MPI_Allreduce(&ok, &all_ok, 1, MPI_INT, MPI_MIN, tree->comm);
        double need_max = need;
        MPI_Allreduce(MPI_IN_PLACE, &need_max, 1, MPI_DOUBLE, MPI_MAX, tree->comm);
        if (!all_ok) {
            store.release();  // the kept factors too
            store.tree = tree;
            store.decided = true;
            if (verbosity >= 0 && rank == 0) {
                std::printf("GPU solve: off (%s on some rank); host solve\n", reason.empty() ? "another rank" : reason.c_str());
                std::fflush(stdout);
            }
            return false;
        }

        Context& ctx = Context::instance();
        ctx.activate();
        cudaStream_t stream = ctx.stream();
        store.levels.resize(static_cast<size_t>(leaf + 1));
        double upload_seconds = 0.0, register_seconds = 0.0, pack_seconds = 0.0, wait_seconds = 0.0;
        double kept_bytes = 0.0;
        for (const auto& [level, kl] : store.kept) kept_bytes += kl.bytes;
        cudaEvent_t staged[2] = {nullptr, nullptr};
        int turn = 0;

        for (int level = 2; level <= leaf; ++level) {
            auto& lvl = tree->levels[static_cast<size_t>(level)];
            DeviceSolveLevel& L = store.levels[static_cast<size_t>(level)];
            if (!lvl.is_process_active) continue;
            L.active = true;
            const size_t nb = lvl.local_boxes.size();
            auto is_local = [&](int64_t m) {
                return m >= lvl.local_morton_start && m < lvl.local_morton_start + static_cast<int64_t>(nb);
            };

            // ---- other ranks' neighbors: sizes and skeletons, and who reads ours
            const auto tr = clock::now();
            L.readers.assign(nb, {});
            std::map<int64_t, std::vector<int64_t>> remote_skel;  // morton -> skeleton
            if (lvl.num_active_processes > 1) {
                L.peers = compute_one_hop_neighbor_ranks(tree, lvl, level);
                std::unordered_map<int, size_t> peer_index;
                for (size_t i = 0; i < L.peers.size(); ++i) peer_index.emplace(L.peers[i], i);
                std::vector<PeerMessage> out(L.peers.size()), in;
                std::map<int64_t, int> wanted;
                for (const auto& box : lvl.local_boxes) {
                    for (int64_t m : box.one_hop) {
                        if (!is_local(m)) wanted.emplace(m, solve_detail::owner_of(tree, lvl, level, m));
                    }
                }
                for (const auto& [m, owner] : wanted) {
                    auto it = peer_index.find(owner);
                    if (it == peer_index.end()) throw std::runtime_error("device solve: neighbor on a non-neighbor rank");
                    out[it->second].head.push_back(m);
                }
                register_seconds += exchange_peer_messages(tree->comm, L.peers, out, in, 740);
                // replies: n, k, skeleton per requested box
                std::vector<PeerMessage> reply(L.peers.size()), got;
                for (size_t i = 0; i < L.peers.size(); ++i) {
                    for (int64_t m : in[i].head) {
                        if (!is_local(m)) throw std::runtime_error("device solve: request for a box of another rank");
                        const auto& box = lvl.local_boxes[static_cast<size_t>(m - lvl.local_morton_start)];
                        L.readers[static_cast<size_t>(m - lvl.local_morton_start)].push_back(static_cast<int>(i));
                        reply[i].head.push_back(m);
                        reply[i].head.push_back(box.num_points);
                        reply[i].head.push_back(static_cast<int64_t>(box.skeleton_indices.size()));
                        reply[i].head.insert(reply[i].head.end(), box.skeleton_indices.begin(), box.skeleton_indices.end());
                    }
                }
                register_seconds += exchange_peer_messages(tree->comm, L.peers, reply, got, 743);
                for (size_t i = 0; i < L.peers.size(); ++i) {
                    const auto& h = got[i].head;
                    for (size_t at = 0; at < h.size();) {
                        const int64_t m = h[at], n = h[at + 1], k = h[at + 2];
                        L.ghost_of.emplace(m, static_cast<int>(L.ghosts.size()));
                        L.ghosts.push_back({m, L.peers[i], static_cast<int>(n), L.ghost_points});
                        L.ghost_points += n;
                        remote_skel[m].assign(h.begin() + static_cast<int64_t>(at) + 3,
                                              h.begin() + static_cast<int64_t>(at) + 3 + k);
                        at += 3 + static_cast<size_t>(k);
                    }
                }
            }
            register_seconds += std::chrono::duration<double>(clock::now() - tr).count();

            // ---- layout: per box [T | X_SR | X_NR | LU] doubles, [ipiv | skel | red] ints
            // (lists only for a level whose factors stayed on the device)
            const KeptSolveLevel* kept = nullptr;
            if (auto kit = store.kept.find(level); kit != store.kept.end()) kept = &kit->second;
            std::vector<const KeptSolveBox*> kept_box(nb, nullptr);
            std::vector<int64_t> base(nb);
            std::vector<solve_detail::BoxLayout> lay(nb);
            std::vector<int> ntot(nb, 0);
            std::vector<size_t> chunk_of(nb), offset_in(nb);
            std::vector<size_t> chunk_bytes;
            constexpr size_t kChunk = size_t{256} << 20;
            for (size_t b = 0; b < nb; ++b) {
                const auto& box = lvl.local_boxes[b];
                base[b] = L.points;
                L.points += box.num_points;
                const size_t k = box.skeleton_indices.size(), r = box.redundant_indices.size();
                const bool eliminated = k > 0 && r > 0;
                solve_detail::BoxLayout& l = lay[b];
                ntot[b] = eliminated && box.X_NR.is_allocated() && !box.one_hop.empty() ? static_cast<int>(box.X_NR.rows) : 0;
                l.has_T = eliminated && box.interpolation_matrix.is_allocated();
                l.has_xsr = eliminated && box.X_SR.is_allocated();
                if (kept != nullptr && eliminated) {
                    auto it = kept->boxes.find(box.morton_index);
                    if (it == kept->boxes.end() || it->second.k != static_cast<int>(k) || it->second.r != static_cast<int>(r) ||
                        it->second.ntot != ntot[b] || (l.has_T && it->second.T == nullptr)) {
                        throw std::runtime_error("device solve: kept factors of box " + std::to_string(box.morton_index) +
                                                 " do not match");
                    }
                    kept_box[b] = &it->second;
                }
                size_t at = 0;
                if (eliminated && kept_box[b] == nullptr) {
                    l.T = at;    at = align_up(at + (l.has_T ? k * r * D : 0));
                    l.xsr = at;  at = align_up(at + (l.has_xsr ? k * r * D : 0));
                    l.xnr = at;  at = align_up(at + static_cast<size_t>(ntot[b]) * r * D);
                    l.lu = at;   at = align_up(at + r * r * D);
                    l.ipiv = at; at = align_up(at + r * I);
                }
                l.skel = at; at = align_up(at + k * I);
                l.red = at;  at = align_up(at + r * I);
                l.bytes = at;
                if (chunk_bytes.empty() || (chunk_bytes.back() > 0 && chunk_bytes.back() + at > kChunk)) chunk_bytes.push_back(0);
                chunk_of[b] = chunk_bytes.size() - 1;
                offset_in[b] = chunk_bytes.back();
                chunk_bytes.back() += at;
                L.max_n = std::max(L.max_n, static_cast<int>(box.num_points));
                L.max_r = std::max(L.max_r, static_cast<int>(r));
            }
            std::vector<size_t> chunk_first(chunk_bytes.size() + 1, nb);  // boxes of chunk c: [first[c], first[c + 1])
            for (size_t b = nb; b-- > 0;) chunk_first[chunk_of[b]] = b;
            std::vector<char*> chunks(chunk_bytes.size());
            for (size_t c = 0; c < chunks.size(); ++c) {
                chunks[c] = heap.alloc_resident(std::max<size_t>(chunk_bytes[c], 1));
                L.blocks.push_back(chunks[c]);
            }
            // ---- upload, chunk by chunk through two pinned staging buffers
            const auto tu = clock::now();
            for (size_t c = 0; c < chunks.size(); ++c) {
                if (staged[turn] == nullptr) {
                    check_cuda(cudaEventCreateWithFlags(&staged[turn], cudaEventDisableTiming), "cudaEventCreate");
                } else {
                    const auto tw = clock::now();
                    check_cuda(cudaEventSynchronize(staged[turn]), "solve upload");
                    wait_seconds += std::chrono::duration<double>(clock::now() - tw).count();
                }
                char* h = static_cast<char*>(pinned_pool().staging[turn].reserve(std::max<size_t>(chunk_bytes[c], 1)));
                const auto tp = clock::now();
                // pieces of about 1 MB (a matrix's column range, or a box's lists), so all threads copy
                struct Piece { size_t box; int what; size_t col0, col1; };  // what: 0 T, 1 X_SR, 2 X_NR, 3 LU, 4 lists
                std::vector<Piece> pieces;
                for (size_t b = chunk_first[c]; b < chunk_first[c + 1]; ++b) {
                    const auto& box = lvl.local_boxes[b];
                    const size_t k = box.skeleton_indices.size(), r = box.redundant_indices.size();
                    pieces.push_back(Piece{b, 4, 0, 0});
                    if (k == 0 || r == 0 || kept_box[b] != nullptr) continue;
                    auto split = [&](int what, size_t rows) {
                        const size_t step = std::max<size_t>(1, (size_t{1} << 20) / std::max<size_t>(rows * D, 1));
                        for (size_t j = 0; j < r; j += step) pieces.push_back(Piece{b, what, j, std::min(r, j + step)});
                    };
                    if (lay[b].has_T) split(0, k);
                    if (lay[b].has_xsr) split(1, k);
                    if (ntot[b] > 0) split(2, static_cast<size_t>(ntot[b]));
                    split(3, r);
                }
                #pragma omp parallel for schedule(dynamic, 1)
                for (int64_t pi = 0; pi < static_cast<int64_t>(pieces.size()); ++pi) {
                    const Piece& pc = pieces[static_cast<size_t>(pi)];
                    const auto& box = lvl.local_boxes[pc.box];
                    const solve_detail::BoxLayout& l = lay[pc.box];
                    char* p = h + offset_in[pc.box];
                    const size_t k = box.skeleton_indices.size(), r = box.redundant_indices.size();
                    auto put = [&](size_t at, const MatrixStorage<DataType>& m, size_t rows) {
                        DataType* d = reinterpret_cast<DataType*>(p + at);
                        for (size_t j = pc.col0; j < pc.col1; ++j)
                            std::memcpy(d + j * rows, m.data.data() + j * static_cast<size_t>(m.lda), rows * D);
                    };
                    switch (pc.what) {
                        case 0: put(l.T, box.interpolation_matrix, k); break;
                        case 1: put(l.xsr, box.X_SR, k); break;
                        case 2: put(l.xnr, box.X_NR, static_cast<size_t>(ntot[pc.box])); break;
                        case 3: put(l.lu, box.X_RR, r); break;
                        default: {
                            if (k > 0 && r > 0 && kept_box[pc.box] == nullptr) {
                                std::memcpy(p + l.ipiv, box.X_RR_pivots.data(), r * I);
                            }
                            int* sk = reinterpret_cast<int*>(p + l.skel);
                            for (size_t i = 0; i < k; ++i) sk[i] = static_cast<int>(box.skeleton_indices[i]);
                            int* rd = reinterpret_cast<int*>(p + l.red);
                            for (size_t i = 0; i < r; ++i) rd[i] = static_cast<int>(box.redundant_indices[i]);
                        }
                    }
                }
                pack_seconds += std::chrono::duration<double>(clock::now() - tp).count();
                check_cuda(cudaMemcpyAsync(chunks[c], h, chunk_bytes[c], cudaMemcpyHostToDevice, stream), "solve upload");
                check_cuda(cudaEventRecord(staged[turn], stream), "cudaEventRecord");
                turn ^= 1;
                store.bytes += static_cast<double>(chunk_bytes[c]);
            }
            upload_seconds += std::chrono::duration<double>(clock::now() - tu).count();

            // ---- boxes and slots
            L.boxes.resize(nb);
            std::vector<int> remote_lists;  // skeletons of other ranks' neighbors
            std::unordered_map<int64_t, size_t> remote_list_at;
            for (const auto& [m, sk] : remote_skel) {
                remote_list_at.emplace(m, remote_lists.size());
                for (int64_t x : sk) remote_lists.push_back(static_cast<int>(x));
            }
            for (size_t b = 0; b < nb; ++b) {
                const auto& box = lvl.local_boxes[b];
                const solve_detail::BoxLayout& l = lay[b];
                const char* p = chunks[chunk_of[b]] + offset_in[b];
                SolveBox& sb = L.boxes[b];
                sb.vec = base[b];
                sb.n = static_cast<int>(box.num_points);
                sb.k = static_cast<int>(box.skeleton_indices.size());
                sb.r = static_cast<int>(box.redundant_indices.size());
                sb.ntot = ntot[b];
                sb.skel = reinterpret_cast<const int*>(p + l.skel);
                sb.red = reinterpret_cast<const int*>(p + l.red);
                const bool eliminated = sb.k > 0 && sb.r > 0;
                if (const KeptSolveBox* kb = kept_box[b]) {
                    sb.T = l.has_T ? kb->T : nullptr;
                    sb.xsr = l.has_xsr ? kb->xsr : nullptr;
                    sb.xnr = kb->xnr;
                    sb.lu = kb->lu;
                    sb.ipiv = kb->ipiv;
                } else {
                    sb.T = eliminated && l.has_T ? p + l.T : nullptr;
                    sb.xsr = eliminated && l.has_xsr ? p + l.xsr : nullptr;
                    sb.xnr = eliminated ? p + l.xnr : nullptr;
                    sb.lu = eliminated ? p + l.lu : nullptr;
                    sb.ipiv = eliminated ? reinterpret_cast<const int*>(p + l.ipiv) : nullptr;
                }
                sb.slot0 = static_cast<int>(L.slots.size());
                sb.nslots = 0;
                sb.work = 0;
                if (ntot[b] == 0) continue;  // no neighbor terms (or no elimination: the vector only)
                if (box.use_full_set.size() != box.one_hop.size()) {
                    throw std::runtime_error("device solve: use_full_set missing for box " + std::to_string(box.morton_index));
                }
                // neighbor counts: checked against those cached at factorization
                static const std::vector<int64_t> none;
                const auto& cached = b < lvl.solve_neighbor_size.size() ? lvl.solve_neighbor_size[b] : none;
                int row = 0;
                for (size_t a = 0; a < box.one_hop.size(); ++a) {
                    const int64_t m = box.one_hop[a];
                    const bool full = box.use_full_set[a] == 1;
                    SolveSlot sl{};
                    sl.full = full ? 1 : 0;
                    sl.row0 = row;
                    if (is_local(m)) {
                        const size_t nbx = static_cast<size_t>(m - lvl.local_morton_start);
                        const auto& nbox = lvl.local_boxes[nbx];
                        sl.ghost = 0;
                        sl.vec = base[nbx];
                        sl.n = static_cast<int>(nbox.num_points);
                        sl.count = full ? sl.n : static_cast<int>(nbox.skeleton_indices.size());
                        sl.skel = reinterpret_cast<const int*>(chunks[chunk_of[nbx]] + offset_in[nbx] + lay[nbx].skel);
                    } else {
                        auto it = L.ghost_of.find(m);
                        if (it == L.ghost_of.end()) throw std::runtime_error("device solve: neighbor " + std::to_string(m) + " unknown");
                        const auto& g = L.ghosts[static_cast<size_t>(it->second)];
                        sl.ghost = 1;
                        sl.vec = g.offset;
                        sl.n = g.n;
                        sl.count = full ? g.n : static_cast<int>(remote_skel.at(m).size());
                        sl.skel = nullptr;  // set once the remote lists are on the device
                    }
                    if (cached.size() == box.one_hop.size() && cached[a] != sl.count) {
                        throw std::runtime_error("device solve: neighbor size mismatch for box " +
                                                 std::to_string(box.morton_index));
                    }
                    row += sl.count;
                    L.slots.push_back(sl);
                    L.slot_morton.push_back(m);
                    ++sb.nslots;
                }
                if (row != ntot[b]) throw std::runtime_error("device solve: X_NR rows do not match the neighbors");
            }
            // ---- the sweeps' tables
            store.bytes += static_cast<double>(
                build_solve_tables(tree, lvl, level, L, remote_lists, remote_list_at, register_seconds));
        }
        check_cuda(cudaStreamSynchronize(stream), "solve upload");
        for (cudaEvent_t e : staged) {
            if (e != nullptr) cudaEventDestroy(e);
        }
        store.usable = true;
        double seconds = std::chrono::duration<double>(clock::now() - t0).count();
        double bytes = store.bytes;
        MPI_Allreduce(MPI_IN_PLACE, &seconds, 1, MPI_DOUBLE, MPI_MAX, tree->comm);
        MPI_Allreduce(MPI_IN_PLACE, &bytes, 1, MPI_DOUBLE, MPI_MAX, tree->comm);
        MPI_Allreduce(MPI_IN_PLACE, &kept_bytes, 1, MPI_DOUBLE, MPI_MAX, tree->comm);
        std::string kept_levels;
        for (const auto& [level, kl] : store.kept) kept_levels += (kept_levels.empty() ? "" : ",") + std::to_string(level);
        (void)verbosity;
        if (rank == 0) {  // once per factorization, also for a silent first solve
            std::printf("GPU solve: set up in %.2f s; factors kept on the device during the factorization: %s "
                        "(up to %.2f GB per rank), uploaded now: up to %.2f GB per rank (upload %.2f s [host packing "
                        "%.2f, waiting for the copies %.2f], neighbor registration %.2f s on rank 0)\n",
                        seconds, kept_levels.empty() ? "none" : ("levels " + kept_levels).c_str(), kept_bytes / 1e9,
                        bytes / 1e9, upload_seconds, pack_seconds, wait_seconds, register_seconds);
            std::fflush(stdout);
        }
        return true;
    }
}

// ---------------------------------------------------------------------------
// The sweeps
// ---------------------------------------------------------------------------
struct DeviceSolveTimes {
    double forward = 0.0, diagonal = 0.0, backward = 0.0, transfer = 0.0, comm = 0.0;
    int64_t direct_messages = 0, staged_messages = 0;  // exchanges moving data: in device memory, through the host
};

template<typename CoordType, typename DataType>
class DeviceSolveRun {
    using S = typename DeviceScalar<DataType>::type;  // vector elements on the device
    static constexpr int kWords = static_cast<int>(sizeof(S) / sizeof(double));  // MPI_DOUBLE per element

public:
    DeviceSolveRun(ParallelTree<CoordType, DataType>* tree, int nrhs) : tree_(tree), nrhs_(nrhs) {}

    // a level's vectors: host solve data -> device
    S* upload(DeviceSolveLevel& L, std::vector<SolveDataRequest<CoordType, DataType>>& data) {
        const auto t0 = std::chrono::steady_clock::now();
        const size_t len = static_cast<size_t>(L.points) * nrhs_;
        S* d = heap_.alloc<S>(std::max<size_t>(len, 1) * sizeof(S));
        S* h = static_cast<S*>(host_vec_.reserve(std::max<size_t>(len, 1) * sizeof(S)));
        for (size_t b = 0; b < L.boxes.size(); ++b) {
            const auto& v = data[b].left_side;
            if (v.size() != static_cast<size_t>(L.boxes[b].n) * nrhs_) throw std::runtime_error("device solve: vector size");
            std::memcpy(h + L.boxes[b].vec * nrhs_, v.data(), v.size() * sizeof(DataType));
        }
        check_cuda(cudaMemcpyAsync(d, h, len * sizeof(S), cudaMemcpyHostToDevice, stream_), "solve upload");
        check_cuda(cudaStreamSynchronize(stream_), "solve upload");
        times.transfer += seconds_since(t0);
        return d;
    }
    void download(DeviceSolveLevel& L, S* d, std::vector<SolveDataRequest<CoordType, DataType>>& data) {
        const auto t0 = std::chrono::steady_clock::now();
        const size_t len = static_cast<size_t>(L.points) * nrhs_;
        S* h = static_cast<S*>(host_vec_.reserve(std::max<size_t>(len, 1) * sizeof(S)));
        check_cuda(cudaMemcpyAsync(h, d, len * sizeof(S), cudaMemcpyDeviceToHost, stream_), "solve download");
        check_cuda(cudaStreamSynchronize(stream_), "solve download");
        for (size_t b = 0; b < L.boxes.size(); ++b) {
            auto& v = data[b].left_side;
            std::memcpy(v.data(), h + L.boxes[b].vec * nrhs_, v.size() * sizeof(DataType));
        }
        heap_.free(d);
        times.transfer += seconds_since(t0);
    }

    // The forward waves of a level, then its diagonal solves.
    void forward(DeviceSolveLevel& L, S* vec) {
        const auto t0 = std::chrono::steady_clock::now();
        S* work = alloc_points(L.work_points);
        S* outbox = alloc_message(L.out_max);
        S* inbox = alloc_message(L.in_max);
        for (const auto& W : L.waves) {
            launch_solve_forward(L.d_boxes, L.d_order + W.order0, W.count, vec, work, nrhs_, L.max_n, stream_);
            launch_solve_accum(L.d_items + W.item0, L.d_parts, W.nitems, W.max_count, nrhs_, vec, work, outbox, inbox,
                               stream_);
            if (L.peers.empty()) continue;
            exchange(L, W.out, W.in, outbox, inbox, 750);
            launch_solve_accum(L.d_items + W.in_item0, L.d_parts, W.in_nitems, W.in_max_count, nrhs_, vec, work, outbox,
                               inbox, stream_);
        }
        check_cuda(cudaStreamSynchronize(stream_), "forward sweep");
        free_points(work);
        free_message(outbox);
        free_message(inbox);
        times.forward += seconds_since(t0);
        const auto t1 = std::chrono::steady_clock::now();
        launch_solve_diagonal(L.d_boxes, L.d_diag, L.ndiag, vec, nrhs_, L.max_r, stream_);
        check_cuda(cudaStreamSynchronize(stream_), "diagonal solves");
        times.diagonal += seconds_since(t1);
    }

    // The backward waves of a level, in reverse, other ranks' neighbor
    // vectors refreshed before each.
    void backward(DeviceSolveLevel& L, S* vec) {
        const auto t0 = std::chrono::steady_clock::now();
        S* work = alloc_points(L.work_points);
        S* ghost = alloc_points(L.ghost_points);
        S* sendbox = alloc_message(L.send_max);
        S* inbox = alloc_message(L.recv_max);
        for (int w = static_cast<int>(L.waves.size()) - 1; w >= 0; --w) {
            const auto& W = L.waves[static_cast<size_t>(w)];
            if (!L.peers.empty()) refresh(L, W.ref[0], vec, ghost, sendbox, inbox);
            launch_solve_backward(L.d_boxes, L.d_slots, L.d_rorder + W.order0, W.count, vec, ghost, work, nrhs_, L.max_n,
                                  W.max_ntot, stream_);
        }
        check_cuda(cudaStreamSynchronize(stream_), "backward sweep");
        free_points(work);
        free_points(ghost);
        free_message(sendbox);
        free_message(inbox);
        times.backward += seconds_since(t0);
    }

    // The multiply F x (color_CA/apply_mul.hpp): a level's forward W waves
    // in order (reading neighbors, other ranks' vectors refreshed before
    // each), then its diagonal multiplies (x_R final by then).
    void mul_forward(DeviceSolveLevel& L, S* vec) {
        const auto t0 = std::chrono::steady_clock::now();
        S* work = alloc_points(L.work_points);
        S* ghost = alloc_points(L.ghost_points);
        S* sendbox = alloc_message(L.send_max);
        S* inbox = alloc_message(L.recv_max);
        for (const auto& W : L.waves) {
            if (!L.peers.empty()) refresh(L, W.ref[1], vec, ghost, sendbox, inbox);
            launch_mul_forward(L.d_boxes, L.d_slots, L.d_order + W.order0, W.count, vec, ghost, work, nrhs_, L.max_n,
                               W.max_ntot, stream_);
        }
        check_cuda(cudaStreamSynchronize(stream_), "multiply forward");
        free_points(work);
        free_points(ghost);
        free_message(sendbox);
        free_message(inbox);
        times.forward += seconds_since(t0);
        const auto t1 = std::chrono::steady_clock::now();
        launch_mul_diagonal(L.d_boxes, L.d_diag, L.ndiag, vec, nrhs_, L.max_r, stream_);
        check_cuda(cudaStreamSynchronize(stream_), "diagonal multiplies");
        times.diagonal += seconds_since(t1);
    }
    // A level's backward V waves in reverse (reverse Morton order within a
    // wave), their neighbor updates applied after each as in the solve's
    // forward sweep.
    void mul_backward(DeviceSolveLevel& L, S* vec) {
        const auto t0 = std::chrono::steady_clock::now();
        S* work = alloc_points(L.work_points);
        S* outbox = alloc_message(L.out_max);
        S* inbox = alloc_message(L.in_max);
        for (int w = static_cast<int>(L.waves.size()) - 1; w >= 0; --w) {
            const auto& W = L.waves[static_cast<size_t>(w)];
            launch_mul_backward(L.d_boxes, L.d_rorder + W.order0, W.count, vec, work, nrhs_, L.max_n, stream_);
            launch_solve_accum(L.d_items + W.item0, L.d_parts, W.nitems, W.max_count, nrhs_, vec, work, outbox, inbox,
                               stream_);
            if (L.peers.empty()) continue;
            exchange(L, W.out, W.in, outbox, inbox, 770);
            launch_solve_accum(L.d_items + W.in_item0, L.d_parts, W.in_nitems, W.in_max_count, nrhs_, vec, work, outbox,
                               inbox, stream_);
        }
        check_cuda(cudaStreamSynchronize(stream_), "multiply backward");
        free_points(work);
        free_message(outbox);
        free_message(inbox);
        times.backward += seconds_since(t0);
    }

    DeviceSolveTimes times;

private:
    static double seconds_since(std::chrono::steady_clock::time_point t) {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
    }
    S* alloc_points(int64_t points) {
        return points > 0 ? heap_.alloc<S>(static_cast<size_t>(points) * nrhs_ * sizeof(S)) : nullptr;
    }
    void free_points(S* p) {
        if (p != nullptr) heap_.free(p);
    }
    // Message buffers: in the MPI exchange arena when MPI is CUDA-aware
    // (sent and received in device memory), else staged through the host.
    S* alloc_message(int64_t points) {
        if (points <= 0) return nullptr;
        const size_t bytes = static_cast<size_t>(points) * nrhs_ * sizeof(S);
        static const bool direct = [] {  // H2_GPU_SOLVE_DIRECT=0: stage through the host
            const char* v = std::getenv("H2_GPU_SOLVE_DIRECT");
            return v == nullptr || std::atoi(v) != 0;
        }();
        if (direct && device_exchange_enabled()) {
            if (char* p = DeviceHeap::exchange_arena().try_alloc(bytes)) return reinterpret_cast<S*>(p);
        }
        return heap_.alloc<S>(bytes);
    }
    void free_message(S* p) {
        if (p == nullptr) return;
        DeviceHeap& arena = DeviceHeap::exchange_arena();
        if (arena.owns(p)) {
            arena.free(p);
        } else {
            heap_.free(p);
        }
    }
    // Other ranks' vectors that changed, into the ghost area.
    void refresh(DeviceSolveLevel& L, const DeviceSolveLevel::Wave::Refresh& R, S* vec, S* ghost,
                 S* sendbox, S* inbox) {
        launch_solve_copy(L.d_copies + R.send0, R.nsend, R.send_max_n, nrhs_, vec, sendbox, stream_);
        exchange(L, R.send, R.recv, sendbox, inbox, 760);
        launch_solve_copy(L.d_copies + R.recv0, R.nrecv, R.recv_max_n, nrhs_, inbox, ghost, stream_);
    }

    // Messages of known size with the neighbor ranks: out[q] and in[q] points
    // to and from peer q, in device memory when both buffers are in the MPI
    // exchange arena, else staged through the host.
    void exchange(DeviceSolveLevel& L, const std::vector<int64_t>& out, const std::vector<int64_t>& in,
                  const S* d_out, S* d_in, int tag) {
        const size_t np = L.peers.size();
        int64_t total_out = 0, total_in = 0;
        for (size_t q = 0; q < np; ++q) {
            total_out += out[q];
            total_in += in[q];
        }
        check_cuda(cudaStreamSynchronize(stream_), "solve messages");
        const auto t0 = std::chrono::steady_clock::now();
        const DeviceHeap& arena = DeviceHeap::exchange_arena();
        if ((total_out == 0 || arena.owns(d_out)) && (total_in == 0 || arena.owns(d_in))) {
            // CUDA-aware MPI: straight from and into device memory
            std::vector<MPI_Request> reqs;
            int64_t at = 0;
            for (size_t q = 0; q < np; ++q) {
                if (in[q] > 0) {
                    reqs.emplace_back();
                    MPI_Irecv(d_in + at * nrhs_, static_cast<int>(in[q] * nrhs_ * kWords), MPI_DOUBLE, L.peers[q], tag, tree_->comm,
                              &reqs.back());
                }
                at += in[q];
            }
            at = 0;
            for (size_t q = 0; q < np; ++q) {
                if (out[q] > 0) {
                    reqs.emplace_back();
                    MPI_Isend(d_out + at * nrhs_, static_cast<int>(out[q] * nrhs_ * kWords), MPI_DOUBLE, L.peers[q], tag,
                              tree_->comm, &reqs.back());
                }
                at += out[q];
            }
            MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
            if (total_out + total_in > 0) ++times.direct_messages;
            times.comm += seconds_since(t0);
            return;
        }
        if (total_out + total_in > 0) ++times.staged_messages;
        S* h_out = static_cast<S*>(host_out_.reserve(std::max<int64_t>(total_out, 1) * nrhs_ * sizeof(S)));
        S* h_in = static_cast<S*>(host_in_.reserve(std::max<int64_t>(total_in, 1) * nrhs_ * sizeof(S)));
        if (total_out > 0) {
            check_cuda(cudaMemcpy(h_out, d_out, static_cast<size_t>(total_out) * nrhs_ * sizeof(S),
                                  cudaMemcpyDeviceToHost), "solve messages");
        }
        std::vector<MPI_Request> reqs;
        int64_t at = 0;
        for (size_t q = 0; q < np; ++q) {
            if (in[q] > 0) {
                reqs.emplace_back();
                MPI_Irecv(h_in + at * nrhs_, static_cast<int>(in[q] * nrhs_ * kWords), MPI_DOUBLE, L.peers[q], tag, tree_->comm,
                          &reqs.back());
            }
            at += in[q];
        }
        at = 0;
        for (size_t q = 0; q < np; ++q) {
            if (out[q] > 0) {
                reqs.emplace_back();
                MPI_Isend(h_out + at * nrhs_, static_cast<int>(out[q] * nrhs_ * kWords), MPI_DOUBLE, L.peers[q], tag, tree_->comm,
                          &reqs.back());
            }
            at += out[q];
        }
        MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
        if (total_in > 0) {
            check_cuda(cudaMemcpyAsync(d_in, h_in, static_cast<size_t>(total_in) * nrhs_ * sizeof(S),
                                       cudaMemcpyHostToDevice, stream_), "solve messages");
        }
        times.comm += seconds_since(t0);
    }

    ParallelTree<CoordType, DataType>* tree_;
    int nrhs_;
    DeviceHeap& heap_ = DeviceHeap::instance();
    cudaStream_t stream_ = Context::instance().stream();
    PinnedBuffer& host_vec_ = device_solve_store().host_vec;
    PinnedBuffer& host_out_ = device_solve_store().host_out;
    PinnedBuffer& host_in_ = device_solve_store().host_in;
};

// The solve sweeps on the device (solve data initialized as by the host
// solve).  Collective over the tree's ranks.
template<typename CoordType, typename DataType>
void device_solve_sweeps(ParallelTree<CoordType, DataType>* tree,
                         std::vector<std::vector<SolveDataRequest<CoordType, DataType>>>& solve_data, int nrhs,
                         int verbosity) {
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    Context::instance().activate();
    DeviceSolveStore& store = device_solve_store();
    const int leaf = tree->num_levels - 1;
    int rank = 0;
    MPI_Comm_rank(tree->comm, &rank);
    DeviceSolveRun<CoordType, DataType> run(tree, nrhs);
    double host_seconds = 0.0;
    // forward (with each level's diagonal solves)
    for (int level = leaf; level >= 1; --level) {
        auto& lvl = tree->levels[static_cast<size_t>(level)];
        if (level >= 2 && lvl.is_process_active) {
            DeviceSolveLevel& L = store.levels[static_cast<size_t>(level)];
            auto* vec = run.upload(L, solve_data[static_cast<size_t>(level)]);
            run.forward(L, vec);
            run.download(L, vec, solve_data[static_cast<size_t>(level)]);
        }
        const auto th = clock::now();
        gather_skeleton_to_parent(lvl, tree->levels[static_cast<size_t>(level - 1)], solve_data[static_cast<size_t>(level)],
                                  solve_data[static_cast<size_t>(level - 1)], tree->dimension, tree->comm);
        host_seconds += std::chrono::duration<double>(clock::now() - th).count();
    }
    const auto th = clock::now();
    auto& root = tree->levels[0];
    if (root.is_process_active && !root.local_boxes.empty()) apply_diagonal_solve(root, solve_data[0][0], false);
    host_seconds += std::chrono::duration<double>(clock::now() - th).count();
    // backward
    for (int level = 1; level <= leaf; ++level) {
        auto& lvl = tree->levels[static_cast<size_t>(level)];
        auto& parent = tree->levels[static_cast<size_t>(level - 1)];
        const auto ts = clock::now();
        if (lvl.is_process_active || parent.is_process_active) {
            scatter_solution_to_children(lvl, parent, solve_data[static_cast<size_t>(level)],
                                         solve_data[static_cast<size_t>(level - 1)], tree->dimension, tree->comm);
        }
        host_seconds += std::chrono::duration<double>(clock::now() - ts).count();
        if (level >= 2 && lvl.is_process_active) {
            DeviceSolveLevel& L = store.levels[static_cast<size_t>(level)];
            auto* vec = run.upload(L, solve_data[static_cast<size_t>(level)]);
            run.backward(L, vec);
            run.download(L, vec, solve_data[static_cast<size_t>(level)]);
        }
    }
    const double total = std::chrono::duration<double>(clock::now() - t0).count();
    if (verbosity >= 0 && rank == 0) {
        const auto& t = run.times;
        std::printf("GPU solve: %.3f s (rank 0: forward %.3f, diagonal %.3f, backward %.3f, of which MPI %.3f%s; "
                    "vector transfers %.3f; host level hand-off and root %.3f)\n",
                    total, t.forward, t.diagonal, t.backward, t.comm,
                    t.direct_messages > 0 ? (t.staged_messages > 0 ? " partly in device memory" : " in device memory")
                                          : (t.staged_messages > 0 ? " through the host" : ""),
                    t.transfer, host_seconds);
        std::fflush(stdout);
    }
}

// The multiply F x on the device (data initialized as by the host multiply;
// hierarchical_mul_parallel).  Collective over the tree's ranks.
template<typename CoordType, typename DataType>
void device_mul_sweeps(ParallelTree<CoordType, DataType>* tree,
                       std::vector<std::vector<SolveDataRequest<CoordType, DataType>>>& data, int nrhs, bool verbose) {
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    Context::instance().activate();
    DeviceSolveStore& store = device_solve_store();
    const int leaf = tree->num_levels - 1;
    int rank = 0;
    MPI_Comm_rank(tree->comm, &rank);
    DeviceSolveRun<CoordType, DataType> run(tree, nrhs);
    double host_seconds = 0.0;
    // forward W (with each level's diagonal multiplies)
    for (int level = leaf; level >= 1; --level) {
        auto& lvl = tree->levels[static_cast<size_t>(level)];
        if (level >= 2 && lvl.is_process_active) {
            DeviceSolveLevel& L = store.levels[static_cast<size_t>(level)];
            auto* vec = run.upload(L, data[static_cast<size_t>(level)]);
            run.mul_forward(L, vec);
            run.download(L, vec, data[static_cast<size_t>(level)]);
        }
        const auto th = clock::now();
        gather_skeleton_to_parent(lvl, tree->levels[static_cast<size_t>(level - 1)], data[static_cast<size_t>(level)],
                                  data[static_cast<size_t>(level - 1)], tree->dimension, tree->comm);
        host_seconds += std::chrono::duration<double>(clock::now() - th).count();
    }
    const auto th = clock::now();
    auto& root = tree->levels[0];
    if (root.is_process_active && !root.local_boxes.empty()) apply_diagonal_multiply(root, data[0][0], false);
    host_seconds += std::chrono::duration<double>(clock::now() - th).count();
    // backward V
    for (int level = 1; level <= leaf; ++level) {
        auto& lvl = tree->levels[static_cast<size_t>(level)];
        auto& parent = tree->levels[static_cast<size_t>(level - 1)];
        const auto ts = clock::now();
        if (lvl.is_process_active || parent.is_process_active) {
            scatter_solution_to_children(lvl, parent, data[static_cast<size_t>(level)],
                                         data[static_cast<size_t>(level - 1)], tree->dimension, tree->comm);
        }
        host_seconds += std::chrono::duration<double>(clock::now() - ts).count();
        if (level >= 2 && lvl.is_process_active) {
            DeviceSolveLevel& L = store.levels[static_cast<size_t>(level)];
            auto* vec = run.upload(L, data[static_cast<size_t>(level)]);
            run.mul_backward(L, vec);
            run.download(L, vec, data[static_cast<size_t>(level)]);
        }
    }
    const double total = std::chrono::duration<double>(clock::now() - t0).count();
    if (verbose && rank == 0) {
        const auto& t = run.times;
        std::printf("GPU multiply: %.3f s (rank 0: forward %.3f, diagonal %.3f, backward %.3f, of which MPI %.3f%s; "
                    "vector transfers %.3f; host level hand-off and root %.3f)\n",
                    total, t.forward, t.diagonal, t.backward, t.comm,
                    t.direct_messages > 0 ? (t.staged_messages > 0 ? " partly in device memory" : " in device memory")
                                          : (t.staged_messages > 0 ? " through the host" : ""),
                    t.transfer, host_seconds);
        std::fflush(stdout);
    }
}

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
