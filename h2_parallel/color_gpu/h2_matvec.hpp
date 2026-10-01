#pragma once
// Device matvec of the compression-only H2 (hierarchical_h2_mul_parallel
// with H2_use_gpu), with the blocks the compression kept on the device
// (h2_matvec_store.hpp).  The whole matvec runs on the device: the input
// is uploaded and the output downloaded once, the level hand-offs
// (gather_skeleton_to_parent, scatter_solution_to_children) run on the
// device, and each level exchanges the vectors other ranks read in one
// round of messages whose layout the compression fixed (commit_device_matvec
// plans them once).  Only a hand-off across a process reduction goes
// through the host functions.  The arithmetic matches the host matvec up to
// rounding.

#ifdef H2_HAVE_GPU

#include "device_heap.hpp"
#include "device_solve.hpp"
#include "gpu_runtime.hpp"
#include "h2_matvec_kernels.hpp"
#include "h2_matvec_store.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <map>
#include <memory>
#include <unordered_map>
#include <vector>

namespace fmm {
namespace gpu {

// Totals of the device matvecs since the compression (seconds, this rank).
struct MatvecStats {
    int64_t calls = 0;
    double total = 0.0;
    double upward = 0.0;        // upward passes and device hand-offs
    double coupling = 0.0;      // skeleton vectors, interactions, downward passes, device hand-offs
    double near = 0.0;          // leaf near blocks
    double mpi = 0.0;           // messages between ranks (within coupling and near)
    double host_handoff = 0.0;  // hand-offs across a process reduction (host)
    double transfer = 0.0;      // input upload, output download
};
inline MatvecStats& matvec_stats() {
    static MatvecStats stats;
    return stats;
}

// The host matvec runs instead while a check computes its reference.
inline bool& device_matvec_suspended() {
    static bool suspended = false;
    return suspended;
}

inline bool device_matvec_usable(const void* tree) {
    const DeviceMatvecStore& store = device_matvec_store();
    return store.usable && store.tree == tree && !device_matvec_suspended();
}

namespace matvec_detail {

constexpr int kPlanTag = 810;

// Layout of the vectors one rank receives from its peers: per peer a
// segment of the listed boxes in list order.
struct Where {
    int64_t base, stride, off;
};

// One kind of message of a level (skeleton vectors, or the leaf's input
// vectors): the requests of this rank, the layouts, and the send spans.
template<typename CoordType, typename DataType>
void plan_exchange(ParallelTree<CoordType, DataType>* tree, TreeLevel<CoordType, DataType>& lvl, int level,
                   const std::vector<int>& peers, const std::map<int64_t, int>& remote,
                   const std::vector<std::vector<std::pair<int64_t, int>>>& requested,  // per peer: (morton, size)
                   const std::unordered_map<int64_t, size_t>& local_index, const H2MatvecLevel& lv, bool near,
                   H2Exchange& ex, std::unordered_map<int64_t, Where>& where, std::vector<H2Span>& spans) {
    const size_t np = peers.size();
    ex.peers = peers;
    ex.recv_base.assign(np, 0);
    ex.recv_size.assign(np, 0);
    ex.send_base.assign(np, 0);
    ex.send_size.assign(np, 0);
    std::unordered_map<int, size_t> peer_index;
    for (size_t i = 0; i < np; ++i) peer_index.emplace(peers[i], i);
    // receive: segments in peer order, boxes in Morton order
    std::vector<std::vector<std::pair<int64_t, int>>> mine(np);
    for (const auto& [morton, size] : remote) {
        auto it = peer_index.find(solve_detail::owner_of(tree, lvl, level, morton));
        if (it == peer_index.end()) throw std::runtime_error("device matvec: source on a non-neighbor rank");
        mine[it->second].emplace_back(morton, size);
    }
    for (size_t p = 0; p < np; ++p) {
        ex.recv_base[p] = ex.recv_points;
        int64_t off = 0;
        for (const auto& [morton, size] : mine[p]) off += size;
        ex.recv_size[p] = off;
        off = 0;
        for (const auto& [morton, size] : mine[p]) {
            where[morton] = Where{ex.recv_base[p], ex.recv_size[p], off};
            off += size;
        }
        ex.recv_points += ex.recv_size[p];
    }
    // send: the boxes each peer asked for, in its order
    for (size_t p = 0; p < np; ++p) {
        ex.send_base[p] = ex.send_points;
        int64_t seg = 0;
        for (const auto& [morton, size] : requested[p]) seg += size;
        ex.send_size[p] = seg;
        int64_t off = 0;
        for (const auto& [morton, size] : requested[p]) {
            auto it = local_index.find(morton);
            if (it == local_index.end()) throw std::runtime_error("device matvec: request for a box of another rank");
            const H2Box& hb = lv.boxes[it->second];
            if (size != (near ? hb.n : hb.k)) {
                throw std::runtime_error("device matvec: box " + std::to_string(morton) +
                                         " has another size on its owner");
            }
            if (size > 0) {
                spans.push_back(H2Span{near ? hb.vec : hb.q, near ? lv.points : lv.q_points, ex.send_base[p], seg, off,
                                       size});
                ex.max_len = std::max(ex.max_len, size);
            }
            off += size;
        }
        ex.send_points += seg;
    }
}

}  // namespace matvec_detail

// End of a compression: the device matvec runs only if every rank kept all
// of its levels; then its messages and hand-offs are planned.  Collective
// over the tree's ranks.
template<typename CoordType, typename DataType>
void commit_device_matvec(ParallelTree<CoordType, DataType>* tree, bool verbose) {
    DeviceMatvecStore& store = device_matvec_store();
    const int leaf = tree->num_levels - 1;
    int ok = store.building && !store.failed && leaf >= 2 ? 1 : 0;
    for (int level = 2; ok && level <= leaf; ++level) {
        const auto& lvl = tree->levels[static_cast<size_t>(level)];
        if (lvl.is_process_active && (static_cast<size_t>(level) >= store.levels.size() ||
                                      !store.levels[static_cast<size_t>(level)].active)) {
            ok = 0;
        }
    }
    int all_ok = 0;
    MPI_Allreduce(&ok, &all_ok, 1, MPI_INT, MPI_MIN, tree->comm);
    store.building = false;
    matvec_stats() = MatvecStats{};
    int rank = 0;
    MPI_Comm_rank(tree->comm, &rank);
    if (!all_ok) {
        store.release();
        if (verbose && rank == 0) {
            std::printf("  GPU matvec: off (%s); host matvec\n",
                        !device_matvec_enabled() ? "H2_GPU_MATVEC=0" : "the blocks do not fit on every rank's device");
            std::fflush(stdout);
        }
        return;
    }
    Context::instance().activate();
    cudaStream_t stream = Context::instance().stream();
    const int nc = morton::children_per_box(tree->dimension);
    double comm_seconds = 0.0;
    for (int level = 2; level <= leaf; ++level) {
        auto& lvl = tree->levels[static_cast<size_t>(level)];
        if (!lvl.is_process_active) continue;
        H2MatvecLevel& lv = store.levels[static_cast<size_t>(level)];
        std::unordered_map<int64_t, size_t> local_index;
        for (size_t b = 0; b < lvl.local_boxes.size(); ++b) local_index.emplace(lvl.local_boxes[b].morton_index, b);
        // remote sources of the blocks: skeleton vectors, and at the leaf the input vectors
        std::map<int64_t, int> remote_q, remote_x;
        for (size_t i = 0; i < lv.block_K.size(); ++i) {
            if (local_index.count(lv.block_source[i])) continue;
            (static_cast<int>(i) < lv.ncoupling ? remote_q : remote_x)[lv.block_source[i]] = lv.block_cols[i];
        }
        std::vector<int> peers;
        if (lvl.num_active_processes > 1) peers = compute_one_hop_neighbor_ranks(tree, lvl, level);
        const size_t np = peers.size();
        // requests: per peer, (morton, size) of the skeleton vectors, then of the input vectors
        std::vector<std::vector<std::pair<int64_t, int>>> want_q(np), want_x(np), asked_q(np), asked_x(np);
        if (np > 0) {
            std::unordered_map<int, size_t> peer_index;
            for (size_t i = 0; i < np; ++i) peer_index.emplace(peers[i], i);
            auto group = [&](const std::map<int64_t, int>& remote, std::vector<std::vector<std::pair<int64_t, int>>>& out) {
                for (const auto& [morton, size] : remote) {
                    auto it = peer_index.find(solve_detail::owner_of(tree, lvl, level, morton));
                    if (it == peer_index.end()) throw std::runtime_error("device matvec: source on a non-neighbor rank");
                    out[it->second].emplace_back(morton, size);
                }
            };
            group(remote_q, want_q);
            group(remote_x, want_x);
            std::vector<PeerMessage> out(np), in;
            for (size_t p = 0; p < np; ++p) {
                auto& h = out[p].head;
                for (const auto* list : {&want_q[p], &want_x[p]}) {
                    h.push_back(static_cast<int64_t>(list->size()));
                    for (const auto& [morton, size] : *list) {
                        h.push_back(morton);
                        h.push_back(size);
                    }
                }
            }
            comm_seconds += exchange_peer_messages(tree->comm, peers, out, in, matvec_detail::kPlanTag);
            for (size_t p = 0; p < np; ++p) {
                const auto& h = in[p].head;
                size_t at = 0;
                for (auto* list : {&asked_q[p], &asked_x[p]}) {
                    const int64_t count = h.at(at++);
                    for (int64_t e = 0; e < count; ++e, at += 2) {
                        list->emplace_back(h.at(at), static_cast<int>(h.at(at + 1)));
                    }
                }
                if (at != h.size()) throw std::runtime_error("device matvec: malformed plan message");
            }
        }
        std::unordered_map<int64_t, matvec_detail::Where> where_q, where_x;
        std::vector<H2Span> spans_q, spans_x;
        matvec_detail::plan_exchange(tree, lvl, level, peers, remote_q, asked_q, local_index, lv, false, lv.coupling,
                                     where_q, spans_q);
        matvec_detail::plan_exchange(tree, lvl, level, peers, remote_x, asked_x, local_index, lv, true, lv.near,
                                     where_x, spans_x);
        // block addressing
        std::vector<H2BlockRef> refs(lv.block_K.size());
        for (size_t i = 0; i < refs.size(); ++i) {
            const bool near = static_cast<int>(i) >= lv.ncoupling;
            const int64_t sm = lv.block_source[i];
            H2BlockRef& ref = refs[i];
            ref.K = lv.block_K[i];
            ref.cols = lv.block_cols[i];
            auto it = local_index.find(sm);
            if (lv.block_partial[i] >= 0) {  // a pair's partial vector
                ref.K = nullptr;
                ref.ghost = 2;
                ref.base = 0;
                ref.stride = lv.partial_points;
                ref.off = lv.block_partial[i];
            } else if (it != local_index.end()) {
                const H2Box& sb = lv.boxes[it->second];
                ref.ghost = 0;
                ref.base = 0;
                ref.stride = near ? lv.points : lv.q_points;
                ref.off = near ? sb.vec : sb.q;
            } else {
                const auto& w = (near ? where_x : where_q).at(sm);
                ref.ghost = 1;
                ref.base = w.base;
                ref.stride = w.stride;
                ref.off = w.off;
            }
        }
        // hand-off to the parent level on the device: no process reduction
        std::vector<H2Handoff> handoff;
        if (level > 2) {
            const auto& parent = tree->levels[static_cast<size_t>(level - 1)];
            lv.handoff_device = parent.num_active_processes == lvl.num_active_processes && parent.is_process_active;
            if (lv.handoff_device) {
                const H2MatvecLevel& pv = store.levels[static_cast<size_t>(level - 1)];
                if (lvl.local_boxes.size() != pv.boxes.size() * static_cast<size_t>(nc)) {
                    throw std::runtime_error("device matvec: parent and child levels do not match");
                }
                for (size_t b = 0; b < lvl.local_boxes.size(); ++b) {
                    const size_t p = b / static_cast<size_t>(nc);
                    int64_t off = 0;
                    for (size_t c = p * nc; c < b; ++c) off += lv.boxes[c].k;
                    handoff.push_back(H2Handoff{lv.boxes[b].vec, pv.boxes[p].vec + off, lv.boxes[b].skel, lv.boxes[b].k});
                    if (b % nc == static_cast<size_t>(nc) - 1 && off + lv.boxes[b].k != pv.boxes[p].n) {
                        throw std::runtime_error("device matvec: parent points do not match the children's skeletons");
                    }
                }
            }
        }
        // tables to the device
        size_t bytes = 0;
        const size_t o_refs = bytes;    bytes = align_up(bytes + refs.size() * sizeof(H2BlockRef));
        const size_t o_q = bytes;       bytes = align_up(bytes + spans_q.size() * sizeof(H2Span));
        const size_t o_x = bytes;       bytes = align_up(bytes + spans_x.size() * sizeof(H2Span));
        const size_t o_hand = bytes;    bytes = align_up(bytes + handoff.size() * sizeof(H2Handoff));
        const size_t o_pairs = bytes;   bytes = align_up(bytes + lv.pairs.size() * sizeof(H2Pair));
        char* d = DeviceHeap::instance().alloc_resident(std::max<size_t>(bytes, 1));
        lv.allocations.push_back(d);
        std::vector<char> image(std::max<size_t>(bytes, 1), 0);
        auto put = [&](size_t at, const void* data, size_t n) {
            if (n > 0) std::memcpy(image.data() + at, data, n);
        };
        put(o_refs, refs.data(), refs.size() * sizeof(H2BlockRef));
        put(o_q, spans_q.data(), spans_q.size() * sizeof(H2Span));
        put(o_x, spans_x.data(), spans_x.size() * sizeof(H2Span));
        put(o_hand, handoff.data(), handoff.size() * sizeof(H2Handoff));
        put(o_pairs, lv.pairs.data(), lv.pairs.size() * sizeof(H2Pair));
        check_cuda(cudaMemcpyAsync(d, image.data(), bytes, cudaMemcpyHostToDevice, stream), "matvec plan");
        check_cuda(cudaStreamSynchronize(stream), "matvec plan");
        lv.d_blocks = reinterpret_cast<const H2BlockRef*>(d + o_refs);
        lv.coupling.d_spans = reinterpret_cast<const H2Span*>(d + o_q);
        lv.coupling.nspans = static_cast<int>(spans_q.size());
        lv.near.d_spans = reinterpret_cast<const H2Span*>(d + o_x);
        lv.near.nspans = static_cast<int>(spans_x.size());
        lv.d_handoff = reinterpret_cast<const H2Handoff*>(d + o_hand);
        lv.nhandoff = static_cast<int>(handoff.size());
        lv.d_pairs = reinterpret_cast<const H2Pair*>(d + o_pairs);
        lv.bytes += static_cast<double>(bytes);
    }
    double bytes = 0.0;
    for (const auto& lv : store.levels) bytes += lv.bytes;
    double max_bytes = bytes;
    MPI_Allreduce(MPI_IN_PLACE, &max_bytes, 1, MPI_DOUBLE, MPI_MAX, tree->comm);
    store.usable = true;
    store.tree = tree;
    store.bytes = bytes;
    if (verbose && rank == 0) {
        std::printf("  GPU matvec: on (blocks kept on the device: up to %.2f GB per rank; messages %s; "
                    "planned with %.3f s of messages on rank 0)\n",
                    max_bytes / 1e9,
                    device_matvec_direct() && device_exchange_enabled() ? "in device memory" : "through the host",
                    comm_seconds);
        std::fflush(stdout);
    }
}

// The matvec on the device.  Collective over the tree's ranks (the level
// exchanges and the hand-offs across process reductions).
template<typename CoordType, typename DataType>
void device_h2_mul(ParallelTree<CoordType, DataType>* tree, const std::vector<DataType>& input,
                   std::vector<DataType>& output, int nrhs) {
    using S = typename DeviceScalar<DataType>::type;
    using clock = std::chrono::steady_clock;
    auto seconds_since = [](clock::time_point t) { return std::chrono::duration<double>(clock::now() - t).count(); };
    constexpr int kWords = static_cast<int>(sizeof(S) / sizeof(double));
    const auto t_start = clock::now();
    auto& stats = matvec_stats();
    ++stats.calls;
    Context& ctx = Context::instance();
    ctx.activate();
    cudaStream_t stream = ctx.stream();
    DeviceHeap& heap = DeviceHeap::instance();
    DeviceMatvecStore& store = device_matvec_store();
    const int leaf = tree->num_levels - 1;
    const int nlev = tree->num_levels;
    int rank = 0;
    MPI_Comm_rank(tree->comm, &rank);
    static PinnedBuffer* send_host = new PinnedBuffer;  // staging (never freed: see PinnedPool)
    static PinnedBuffer* recv_host = new PinnedBuffer;

    std::vector<S*> src(static_cast<size_t>(nlev), nullptr), tgt(static_cast<size_t>(nlev), nullptr);
    std::vector<S*> owned;
    auto alloc = [&](int64_t points) -> S* {
        if (points <= 0) return nullptr;
        S* p = heap.alloc<S>(static_cast<size_t>(points) * nrhs * sizeof(S));
        owned.push_back(p);
        return p;
    };
    auto zeroed = [&](int64_t points) -> S* {
        S* p = alloc(points);
        if (p != nullptr) {
            check_cuda(cudaMemsetAsync(p, 0, static_cast<size_t>(points) * nrhs * sizeof(S), stream), "matvec zero");
        }
        return p;
    };
    auto lv_of = [&](int level) -> H2MatvecLevel& { return store.levels[static_cast<size_t>(level)]; };
    auto active = [&](int level) { return tree->levels[static_cast<size_t>(level)].is_process_active; };
    // one round of messages with the level's peers; returns the ghost area
    auto exchange = [&](const H2Exchange& ex, const S* from, int tag) -> S* {
        if (ex.peers.empty()) return nullptr;
        S* send = alloc(ex.send_points);
        S* recv = alloc(ex.recv_points);
        launch_h2_pack<S>(ex.d_spans, ex.nspans, from, send, nrhs, ex.max_len, stream);
        check_cuda(cudaStreamSynchronize(stream), "matvec messages");
        const auto t0 = clock::now();
        const bool direct = device_matvec_direct() && device_exchange_enabled();
        S* sbuf = send;
        S* rbuf = recv;
        if (!direct) {
            sbuf = static_cast<S*>(send_host->reserve(std::max<int64_t>(ex.send_points, 1) * nrhs * sizeof(S)));
            rbuf = static_cast<S*>(recv_host->reserve(std::max<int64_t>(ex.recv_points, 1) * nrhs * sizeof(S)));
            if (ex.send_points > 0) {
                check_cuda(cudaMemcpy(sbuf, send, static_cast<size_t>(ex.send_points) * nrhs * sizeof(S),
                                      cudaMemcpyDeviceToHost), "matvec messages");
            }
        }
        std::vector<MPI_Request> reqs;
        for (size_t p = 0; p < ex.peers.size(); ++p) {
            if (ex.recv_size[p] > 0) {
                reqs.emplace_back();
                MPI_Irecv(rbuf + ex.recv_base[p] * nrhs, static_cast<int>(ex.recv_size[p] * nrhs * kWords), MPI_DOUBLE,
                          ex.peers[p], tag, tree->comm, &reqs.back());
            }
        }
        for (size_t p = 0; p < ex.peers.size(); ++p) {
            if (ex.send_size[p] > 0) {
                reqs.emplace_back();
                MPI_Isend(sbuf + ex.send_base[p] * nrhs, static_cast<int>(ex.send_size[p] * nrhs * kWords), MPI_DOUBLE,
                          ex.peers[p], tag, tree->comm, &reqs.back());
            }
        }
        MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
        if (!direct && ex.recv_points > 0) {
            check_cuda(cudaMemcpyAsync(recv, rbuf, static_cast<size_t>(ex.recv_points) * nrhs * sizeof(S),
                                       cudaMemcpyHostToDevice, stream), "matvec messages");
        }
        stats.mpi += seconds_since(t0);
        return recv;
    };
    // a level's device vector <-> host solve data (hand-offs across a reduction)
    auto to_host = [&](int level, const S* d, std::vector<SolveDataRequest<CoordType, DataType>>& data) {
        auto& lvl = tree->levels[static_cast<size_t>(level)];
        const H2MatvecLevel& L = lv_of(level);
        std::vector<DataType> h(static_cast<size_t>(L.points) * nrhs);
        if (!h.empty()) {
            check_cuda(cudaMemcpyAsync(h.data(), d, h.size() * sizeof(S), cudaMemcpyDeviceToHost, stream), "matvec");
            check_cuda(cudaStreamSynchronize(stream), "matvec");
        }
        data.resize(lvl.local_boxes.size());
        for (size_t b = 0; b < lvl.local_boxes.size(); ++b) {
            const auto& box = lvl.local_boxes[b];
            auto& r = data[b];
            r.initialize(box.morton_index, rank, box.num_points, nrhs);
            r.skeleton_indices = box.skeleton_indices;
            r.redundant_indices = box.redundant_indices;
            const H2Box& hb = L.boxes[b];
            for (int c = 0; c < nrhs; ++c) {
                for (int i = 0; i < hb.n; ++i) {
                    r.left_side[static_cast<size_t>(i + c * hb.n)] =
                        h[static_cast<size_t>(c * L.points + hb.vec + i)];
                }
            }
            r.right_side = r.left_side;
        }
    };
    auto to_device = [&](int level, const std::vector<SolveDataRequest<CoordType, DataType>>& data) -> S* {
        auto& lvl = tree->levels[static_cast<size_t>(level)];
        const H2MatvecLevel& L = lv_of(level);
        if (data.size() != lvl.local_boxes.size()) throw std::runtime_error("device matvec: hand-off data size");
        std::vector<DataType> h(static_cast<size_t>(L.points) * nrhs);
        for (size_t b = 0; b < data.size(); ++b) {
            const H2Box& hb = L.boxes[b];
            if (data[b].morton_index != lvl.local_boxes[b].morton_index ||
                data[b].left_side.size() != static_cast<size_t>(hb.n) * nrhs) {
                throw std::runtime_error("device matvec: hand-off data does not match the level's boxes");
            }
            for (int c = 0; c < nrhs; ++c) {
                for (int i = 0; i < hb.n; ++i) {
                    h[static_cast<size_t>(c * L.points + hb.vec + i)] =
                        data[b].left_side[static_cast<size_t>(i + c * hb.n)];
                }
            }
        }
        S* d = alloc(L.points);
        if (!h.empty()) {
            check_cuda(cudaMemcpyAsync(d, h.data(), h.size() * sizeof(S), cudaMemcpyHostToDevice, stream), "matvec");
            check_cuda(cudaStreamSynchronize(stream), "matvec");
        }
        return d;
    };

    // ---- the input (the leaf level's layout is the host's input layout)
    auto t0 = clock::now();
    S* orig = nullptr;
    if (active(leaf)) {
        const H2MatvecLevel& L = lv_of(leaf);
        if (input.size() != static_cast<size_t>(L.points) * nrhs) throw std::runtime_error("device matvec: input size");
        orig = alloc(L.points);
        src[static_cast<size_t>(leaf)] = alloc(L.points);
        if (L.points > 0) {
            check_cuda(cudaMemcpyAsync(orig, input.data(), input.size() * sizeof(S), cudaMemcpyHostToDevice, stream),
                       "matvec input");
            check_cuda(cudaMemcpyAsync(src[static_cast<size_t>(leaf)], orig, input.size() * sizeof(S),
                                       cudaMemcpyDeviceToDevice, stream), "matvec input");
        }
    }
    stats.transfer += seconds_since(t0);

    // ---- upward: q_B = x_B[S] + T_B x_B[R], child skeletons to the parents
    t0 = clock::now();
    for (int level = leaf; level >= 2; --level) {
        if (active(level)) {
            const H2MatvecLevel& L = lv_of(level);
            launch_h2_upward<S>(L.d_boxes, static_cast<int>(L.boxes.size()), src[static_cast<size_t>(level)], L.points,
                                nrhs, L.max_r, stream);
        }
        if (level == 2) continue;
        if (active(level) && lv_of(level).handoff_device) {
            const H2MatvecLevel& L = lv_of(level);
            const H2MatvecLevel& P = lv_of(level - 1);
            src[static_cast<size_t>(level - 1)] = alloc(P.points);
            launch_h2_handoff<S>(L.d_handoff, L.nhandoff, src[static_cast<size_t>(level)], L.points,
                                 src[static_cast<size_t>(level - 1)], P.points, nrhs, true, stream);
        } else {
            const auto th = clock::now();
            std::vector<SolveDataRequest<CoordType, DataType>> child, parent;
            if (active(level)) to_host(level, src[static_cast<size_t>(level)], child);
            gather_skeleton_to_parent(tree->levels[static_cast<size_t>(level)],
                                      tree->levels[static_cast<size_t>(level - 1)], child, parent, tree->dimension,
                                      tree->comm);
            if (active(level - 1)) src[static_cast<size_t>(level - 1)] = to_device(level - 1, parent);
            stats.host_handoff += seconds_since(th);
        }
    }
    check_cuda(cudaStreamSynchronize(stream), "matvec upward");
    stats.upward += seconds_since(t0);

    // ---- interactions and downward: y_B[S] += sum K q, y_B[R] += T^T y_B[S],
    // parents to the children's skeletons
    t0 = clock::now();
    for (int level = 2; level <= leaf; ++level) {
        if (active(level)) {
            const H2MatvecLevel& L = lv_of(level);
            if (tgt[static_cast<size_t>(level)] == nullptr) tgt[static_cast<size_t>(level)] = zeroed(L.points);
            S* q = alloc(L.q_points);
            launch_h2_gather_q<S>(L.d_boxes, static_cast<int>(L.boxes.size()), src[static_cast<size_t>(level)], L.points,
                                  q, L.q_points, nrhs, stream);
            S* ghost = exchange(L.coupling, q, 900 + 2 * level);
            launch_h2_accumulate<S>(L.d_boxes, static_cast<int>(L.boxes.size()), L.d_blocks,
                                    tgt[static_cast<size_t>(level)], L.points, q, ghost, nullptr, nrhs, false, L.max_k,
                                    L.max_q_cols, stream);
            launch_h2_downward<S>(L.d_boxes, static_cast<int>(L.boxes.size()), tgt[static_cast<size_t>(level)], L.points,
                                  nrhs, L.max_k, stream);
        }
        if (level == leaf) continue;
        if (active(level + 1) && lv_of(level + 1).handoff_device) {
            const H2MatvecLevel& C = lv_of(level + 1);
            const H2MatvecLevel& L = lv_of(level);
            tgt[static_cast<size_t>(level + 1)] = zeroed(C.points);
            launch_h2_handoff<S>(C.d_handoff, C.nhandoff, tgt[static_cast<size_t>(level + 1)], C.points,
                                 tgt[static_cast<size_t>(level)], L.points, nrhs, false, stream);
        } else {
            const auto th = clock::now();
            std::vector<SolveDataRequest<CoordType, DataType>> parent, child;
            if (active(level)) to_host(level, tgt[static_cast<size_t>(level)], parent);
            if (active(level + 1)) {
                auto& clvl = tree->levels[static_cast<size_t>(level + 1)];
                child.resize(clvl.local_boxes.size());
                for (size_t b = 0; b < child.size(); ++b) {
                    const auto& box = clvl.local_boxes[b];
                    child[b].initialize(box.morton_index, rank, box.num_points, nrhs);
                    child[b].skeleton_indices = box.skeleton_indices;
                    child[b].redundant_indices = box.redundant_indices;
                }
            }
            scatter_solution_to_children(tree->levels[static_cast<size_t>(level + 1)],
                                         tree->levels[static_cast<size_t>(level)], child, parent, tree->dimension,
                                         tree->comm);
            if (active(level + 1)) tgt[static_cast<size_t>(level + 1)] = to_device(level + 1, child);
            stats.host_handoff += seconds_since(th);
        }
    }
    check_cuda(cudaStreamSynchronize(stream), "matvec downward");
    stats.coupling += seconds_since(t0);

    // ---- leaf near blocks: y += sum K x
    t0 = clock::now();
    if (active(leaf)) {
        const H2MatvecLevel& L = lv_of(leaf);
        S* ghost = exchange(L.near, orig, 990);
        S* partial = alloc(L.partial_points);  // the pairs' products, each block read once
        launch_h2_pairs<S>(L.d_pairs, static_cast<int>(L.pairs.size()), orig, L.points, partial, L.partial_points, nrhs,
                           L.max_pair_cols, stream);
        launch_h2_accumulate<S>(L.d_boxes, static_cast<int>(L.boxes.size()), L.d_blocks, tgt[static_cast<size_t>(leaf)],
                                L.points, orig, ghost, partial, nrhs, true, L.max_n, L.max_near_cols, stream);
    }
    check_cuda(cudaStreamSynchronize(stream), "matvec near");
    stats.near += seconds_since(t0);

    // ---- the output
    t0 = clock::now();
    output.assign(input.size(), DataType{0});
    if (active(leaf) && !output.empty()) {
        check_cuda(cudaMemcpy(output.data(), tgt[static_cast<size_t>(leaf)], output.size() * sizeof(S),
                              cudaMemcpyDeviceToHost), "matvec output");
    }
    for (S* p : owned) heap.free(p);
    stats.transfer += seconds_since(t0);
    stats.total += seconds_since(t_start);
}

// The device matvec when the compression kept the blocks on every rank's
// device; false (nothing done) otherwise.
template<typename CoordType, typename DataType>
bool run_device_h2_mul(ParallelTree<CoordType, DataType>* tree, const std::vector<DataType>& input,
                       std::vector<DataType>& output, int nrhs) {
    if constexpr (gpu_data_type<DataType>) {
        activate_operator(tree);
        if (!device_matvec_usable(tree)) return false;
        device_h2_mul(tree, input, output, nrhs);
        return true;
    } else {
        (void)tree;
        (void)input;
        (void)output;
        (void)nrhs;
        return false;
    }
}

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
