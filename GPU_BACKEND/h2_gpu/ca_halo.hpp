#pragma once
// Device halo of a replicated CA level whose blocks stayed on the device
// (h2_gpu/GPU_CA_PLAN.md, M2).
//
// When a device CA level follows another device level, the transition keeps
// the parent blocks on the device, so the level's host boxes carry no blocks
// and the host halo gather (gather_CA_factorization_data) moves only their
// metadata: points, index and neighbor lists, assisting boxes.  The ghosts'
// blocks then move here, between the devices: every rank sends, for each
// ghost G a peer requested from it (the gather's request plan), G's Schur
// block and its pair blocks (G, X), X one hop from G, except the pairs the
// peer already has or gets with another ghost:
//   - X is the peer's own box: the peer's transition built its copy of the
//     pair, bitwise equal to this one (h2_gpu/level_eliminator.hpp, the
//     transition's two-hop fills on CA levels);
//   - X is also a ghost of the peer and X < G: the pair comes with X.
// Both sides derive the same list, in the same order (G in the request
// order, its Schur block, then X ascending), from data they share: the
// request sets (exchanged by the host gather's helpers), the one-hop lists
// and the box sizes.  Blocks keep the device layout (rows = the higher box).
// The ghosts' blocks are allocated first (added to the level's blocks for
// the eliminator to adopt); the messages then move in rounds of one chunk
// per message, packed straight from the blocks into slots (the exchange
// arena with CUDA-aware MPI, else small heap slots staged through the host)
// and unpacked straight into the new blocks: no message-sized buffers, which
// a fragmented heap may not hold.

#ifdef H2_HAVE_GPU

#include "level_eliminator.hpp"

#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <climits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace fmm {
namespace gpu {

struct CaHaloStats {
    double total = 0.0, plan = 0.0, pack = 0.0, mpi = 0.0, unpack = 0.0;
    double bytes_sent = 0.0, bytes_received = 0.0;
    int64_t blocks_received = 0;
    bool staged = false;  // through host buffers (MPI without device memory)
};

namespace ca_halo_detail {

// A block of a ghost's payload: its Schur block (lo == hi) or the pair
// (lo, hi), rows x cols as stored on the device (rows = hi's points).
struct Piece {
    int64_t lo = 0, hi = 0;
    int rows = 0, cols = 0;
    size_t elems() const { return static_cast<size_t>(rows) * static_cast<size_t>(cols); }
};

// Messages in chunks under MPI's int count.
constexpr size_t kChunk = size_t{1} << 30;

inline int count_of(size_t bytes) { return static_cast<int>(std::min(bytes, kChunk)); }

}  // namespace ca_halo_detail

template<typename CoordType, typename DataType>
CaHaloStats exchange_ca_ghost_blocks(ParallelTree<CoordType, DataType>* tree, int lvl,
                                     DeviceLevelBlocks<DataType>& blocks, MPI_Comm level_comm) {
  if constexpr (!gpu_data_type<DataType>) {
    (void)tree; (void)lvl; (void)blocks; (void)level_comm;
    throw std::runtime_error("exchange_ca_ghost_blocks: data type not run on the GPU");
  } else {
    using S = typename DeviceScalar<DataType>::type;
    using ca_halo_detail::Piece;
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    CaHaloStats stats;
    auto& level = tree->levels[static_cast<size_t>(lvl)];
    Context& ctx = Context::instance();
    ctx.activate();
    cudaStream_t stream = ctx.stream();
    DeviceHeap& heap = DeviceHeap::instance();
    const size_t D = sizeof(S);

    // ---- who sends which ghosts to whom (the host gather's plan), and the
    //      peers' request sets
    const auto plan = make_CA_request_plan(tree, lvl, level.ghost_id, level_comm);
    int comm_size = 0, comm_rank = 0;
    MPI_Comm_size(level_comm, &comm_size);
    MPI_Comm_rank(level_comm, &comm_rank);
    const auto peer_sets = exchange_CA_full_requested_ghost_sets(plan, level.ghost_id, level_comm);
    std::unordered_map<int, int> comm_of_global;
    for (int c = 0; c < comm_size; ++c) comm_of_global.emplace(plan.global_rank_by_comm[static_cast<size_t>(c)], c);
    std::unordered_map<int64_t, int> owner_cache;
    auto owner_comm = [&](int64_t m) {
        auto it = owner_cache.find(m);
        if (it != owner_cache.end()) return it->second;
        const std::vector<uint64_t> one = {static_cast<uint64_t>(m)};
        const auto region = morton::assign_to_processes_nd(tree->dimension, one, level.num_active_processes,
                                                           uint32_t{1} << lvl);
        const int owner = comm_of_global.at(level.morton_to_rank.at(static_cast<int>(region.front())));
        owner_cache.emplace(m, owner);
        return owner;
    };
    const std::unordered_set<int64_t> my_requests(level.ghost_id.begin(), level.ghost_id.end());

    // the payload of ghost g for receiver `peer`, from the neighbors' sizes
    auto payload = [](int64_t g, int n_g, std::vector<int64_t> hop, auto&& peer_has, auto&& also_requested,
                      auto&& size_of, std::vector<Piece>& out) {
        out.push_back(Piece{g, g, n_g, n_g});
        std::sort(hop.begin(), hop.end());
        hop.erase(std::unique(hop.begin(), hop.end()), hop.end());
        for (int64_t x : hop) {
            if (x == g || peer_has(x) || (also_requested(x) && x < g)) continue;
            const int n_x = size_of(x);
            const int64_t lo = std::min(g, x), hi = std::max(g, x);
            out.push_back(Piece{lo, hi, hi == g ? n_g : n_x, hi == g ? n_x : n_g});
        }
    };

    // sent: the ghosts peers requested from this rank (their shapes from the
    // device blocks)
    std::vector<std::vector<Piece>> sends(static_cast<size_t>(comm_size));
    std::vector<size_t> send_elems(static_cast<size_t>(comm_size), 0);
    for (int peer = 0; peer < comm_size; ++peer) {
        const auto& requested = plan.incoming[static_cast<size_t>(peer)];
        if (requested.empty()) continue;
        const auto& peer_set = peer_sets[static_cast<size_t>(peer)];
        for (int64_t g : requested) {
            const auto* box = level.find_local_box(g);
            if (box == nullptr) throw std::runtime_error("CA device halo: requested ghost " + std::to_string(g) + " is not local");
            payload(
                g, static_cast<int>(box->num_points), box->one_hop,
                [&](int64_t x) { return owner_comm(x) == peer; },
                [&](int64_t x) { return peer_set.count(x) != 0; },
                [&](int64_t x) -> int {
                    const auto it = blocks.edges.find((static_cast<uint64_t>(std::min(g, x)) << 32) |
                                                      static_cast<uint64_t>(std::max(g, x)));
                    if (it == blocks.edges.end()) {
                        throw std::runtime_error("CA device halo: pair (" + std::to_string(g) + ", " +
                                                 std::to_string(x) + ") is not on the device");
                    }
                    return x > g ? it->second.rows : it->second.cols;
                },
                sends[static_cast<size_t>(peer)]);
        }
        for (const Piece& p : sends[static_cast<size_t>(peer)]) send_elems[static_cast<size_t>(peer)] += p.elems();
    }
    // received: this rank's ghosts, from their owners (shapes from the
    // gathered metadata)
    auto size_on_this_rank = [&](int64_t x) -> int {
        if (const auto* gb = level.find_ghost_box(x)) return static_cast<int>(gb->num_points);
        const auto it = level.assisting_box_points_for_kernel_evaluation.find(x);
        if (it == level.assisting_box_points_for_kernel_evaluation.end()) {
            throw std::runtime_error("CA device halo: neighbor " + std::to_string(x) + " of a ghost is unknown here");
        }
        return static_cast<int>(level.assisting_boxes[static_cast<size_t>(it->second)].indices.size());
    };
    std::vector<std::vector<Piece>> recvs(static_cast<size_t>(comm_size));
    std::vector<size_t> recv_elems(static_cast<size_t>(comm_size), 0);
    for (int peer = 0; peer < comm_size; ++peer) {
        const auto& requested = plan.outgoing[static_cast<size_t>(peer)];
        for (int64_t g : requested) {
            const auto* gb = level.find_ghost_box(g);
            if (gb == nullptr) throw std::runtime_error("CA device halo: ghost " + std::to_string(g) + " has no metadata");
            payload(
                g, static_cast<int>(gb->num_points), gb->one_hop,
                [&](int64_t x) { return level.find_local_box(x) != nullptr; },
                [&](int64_t x) { return my_requests.count(x) != 0; }, size_on_this_rank,
                recvs[static_cast<size_t>(peer)]);
        }
        for (const Piece& p : recvs[static_cast<size_t>(peer)]) recv_elems[static_cast<size_t>(peer)] += p.elems();
    }
    // both sides must agree on every message's size
    {
        std::vector<int64_t> out(send_elems.begin(), send_elems.end()), in(static_cast<size_t>(comm_size), 0);
        MPI_Alltoall(out.data(), 1, MPI_INT64_T, in.data(), 1, MPI_INT64_T, level_comm);
        for (int peer = 0; peer < comm_size; ++peer) {
            if (static_cast<size_t>(in[static_cast<size_t>(peer)]) != recv_elems[static_cast<size_t>(peer)]) {
                throw std::runtime_error("CA device halo: rank " + std::to_string(comm_rank) + " expects " +
                                         std::to_string(recv_elems[static_cast<size_t>(peer)]) +
                                         " elements from rank " + std::to_string(peer) + ", which sends " +
                                         std::to_string(in[static_cast<size_t>(peer)]));
            }
        }
    }
    stats.plan = std::chrono::duration<double>(clock::now() - t0).count();

    // ---- the device blocks of the messages; the ghosts' new blocks
    auto t = clock::now();
    std::vector<std::vector<S*>> send_src(static_cast<size_t>(comm_size)), recv_dst(static_cast<size_t>(comm_size));
    for (int peer = 0; peer < comm_size; ++peer) {
        for (const Piece& p : sends[static_cast<size_t>(peer)]) {
            const DeviceMatrixT<S>* src = nullptr;
            if (p.lo == p.hi) {
                const auto it = blocks.schur.find(p.lo);
                if (it == blocks.schur.end()) throw std::runtime_error("CA device halo: Schur block missing");
                src = &it->second;
            } else {
                src = &blocks.edges.at((static_cast<uint64_t>(p.lo) << 32) | static_cast<uint64_t>(p.hi));
            }
            if (src->rows != p.rows || src->cols != p.cols) throw std::runtime_error("CA device halo: block shape");
            send_src[static_cast<size_t>(peer)].push_back(src->ptr);
            stats.bytes_sent += static_cast<double>(p.elems() * D);
        }
    }
    for (int peer = 0; peer < comm_size; ++peer) {
        for (const Piece& p : recvs[static_cast<size_t>(peer)]) {
            DeviceMatrixT<S> m{heap.alloc_resident<S>(p.elems() * D), p.rows, p.cols};
            const bool fresh = p.lo == p.hi
                ? blocks.schur.emplace(p.lo, m).second
                : blocks.edges.emplace((static_cast<uint64_t>(p.lo) << 32) | static_cast<uint64_t>(p.hi), m).second;
            if (!fresh) {
                throw std::runtime_error("CA device halo: block (" + std::to_string(p.lo) + ", " +
                                         std::to_string(p.hi) + ") received twice");
            }
            recv_dst[static_cast<size_t>(peer)].push_back(m.ptr);
            stats.bytes_received += static_cast<double>(p.elems() * D);
            ++stats.blocks_received;
        }
    }
    stats.plan += std::chrono::duration<double>(clock::now() - t).count();

    // ---- rounds: each moves one chunk of every message (a message is its
    //      blocks back to back), packed straight from the blocks into a slot
    //      and unpacked from the slot into the ghosts' blocks.  Slots: the
    //      exchange arena with CUDA-aware MPI (registered once for MPI;
    //      buffers of the main heap made the transfers ~30x slower), else
    //      small heap slots staged through the host.  Every rank runs the
    //      same number of rounds, so each send meets its receive.
    DeviceHeap& arena = DeviceHeap::exchange_arena();
    const bool device_mpi = device_exchange_enabled() && arena.initialized();
    stats.staged = !device_mpi;
    int messages = 0;
    for (int peer = 0; peer < comm_size; ++peer) {
        messages += (send_elems[static_cast<size_t>(peer)] > 0) + (recv_elems[static_cast<size_t>(peer)] > 0);
    }
    // (one chunk size on every rank: a message's pieces must match on both ends)
    unsigned long long chunk_all = size_t{64} << 20;
    if (device_mpi) {
        chunk_all = std::min((arena.capacity() - arena.used()) / static_cast<size_t>(std::max(messages, 1)),
                             ca_halo_detail::kChunk) / 256 * 256;
    }
    MPI_Allreduce(MPI_IN_PLACE, &chunk_all, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, level_comm);
    const size_t chunk = static_cast<size_t>(chunk_all);
    if (chunk < (size_t{1} << 20)) {
        throw std::runtime_error("CA device halo: the exchange arena has no room (BPACK_GPU_EXCHANGE_MB)");
    }
    const size_t chunk_elems = chunk / D;
    int rounds = 0;
    for (int peer = 0; peer < comm_size; ++peer) {
        const size_t most = std::max(send_elems[static_cast<size_t>(peer)], recv_elems[static_cast<size_t>(peer)]);
        rounds = std::max(rounds, static_cast<int>((most + chunk_elems - 1) / chunk_elems));
    }
    MPI_Allreduce(MPI_IN_PLACE, &rounds, 1, MPI_INT, MPI_MAX, level_comm);
    std::vector<S*> send_slot(static_cast<size_t>(comm_size), nullptr), recv_slot(static_cast<size_t>(comm_size), nullptr);
    std::vector<std::vector<char>> host_send(static_cast<size_t>(comm_size)), host_recv(static_cast<size_t>(comm_size));
    for (int peer = 0; peer < comm_size; ++peer) {
        auto slot = [&]() -> S* { return reinterpret_cast<S*>(device_mpi ? arena.alloc(chunk) : heap.alloc(chunk)); };
        if (send_elems[static_cast<size_t>(peer)] > 0) {
            send_slot[static_cast<size_t>(peer)] = slot();
            if (!device_mpi) host_send[static_cast<size_t>(peer)].resize(chunk);
        }
        if (recv_elems[static_cast<size_t>(peer)] > 0) {
            recv_slot[static_cast<size_t>(peer)] = slot();
            if (!device_mpi) host_recv[static_cast<size_t>(peer)].resize(chunk);
        }
    }
    // Copies between the elements [e0, e1) of a message and its slot: per
    // block, its partial first column, whole columns, and partial last column
    // (items no wider than the block, as the gather kernel's grid needs).
    MetaBuilder& meta = pinned_pool().meta;  // (no eliminator runs now)
    DeviceBuffer meta_device;
    std::vector<GatherItemT<S>> copies;
    int max_m = 0, max_n = 0;
    auto add_copies = [&](const std::vector<Piece>& pieces, const std::vector<S*>& ptrs, size_t e0, size_t e1, S* slot,
                          bool to_slot) {
        size_t off = 0;
        for (size_t i = 0; i < pieces.size() && off < e1; ++i) {
            const size_t e = pieces[i].elems();
            const size_t u0 = std::max(off, e0), u1 = std::min(off + e, e1);
            if (u0 < u1) {
                const int R = pieces[i].rows;
                const size_t l0 = u0 - off, l1 = u1 - off;  // within the block
                auto part = [&](size_t first, int m, int n) {
                    S* block_at = ptrs[i] + first;
                    S* slot_at = slot + (first + off - e0);
                    copies.push_back(to_slot ? GatherItemT<S>{slot_at, R, m, n, block_at, 1, R, IndexList{}, IndexList{}}
                                             : GatherItemT<S>{block_at, R, m, n, slot_at, 1, R, IndexList{}, IndexList{}});
                    max_m = std::max(max_m, m);
                    max_n = std::max(max_n, n);
                };
                const size_t c0 = l0 / static_cast<size_t>(R), r0 = l0 % static_cast<size_t>(R);
                const size_t c1 = l1 / static_cast<size_t>(R), r1 = l1 % static_cast<size_t>(R);
                if (c0 == c1) {
                    part(l0, static_cast<int>(r1 - r0), 1);
                } else {
                    size_t c = c0;
                    if (r0 > 0) {
                        part(l0, static_cast<int>(R - static_cast<int>(r0)), 1);
                        ++c;
                    }
                    if (c1 > c) part(c * static_cast<size_t>(R), R, static_cast<int>(c1 - c));
                    if (r1 > 0) part(c1 * static_cast<size_t>(R), static_cast<int>(r1), 1);
                }
            }
            off += e;
        }
    };
    auto run_copies = [&] {
        if (!copies.empty()) {
            meta.clear();
            const size_t o = meta.append(copies);
            char* md = meta.upload(meta_device, stream);
            launch_gather(reinterpret_cast<const GatherItemT<S>*>(md + o), static_cast<int>(copies.size()), max_m,
                          max_n, md, stream);
        }
        check_cuda(cudaStreamSynchronize(stream), "CA device halo copies");
        copies.clear();
        max_m = max_n = 0;
    };
    auto span = [&](size_t total, int round) -> std::pair<size_t, size_t> {
        const size_t e0 = std::min(total, static_cast<size_t>(round) * chunk_elems);
        return {e0, std::min(total, e0 + chunk_elems)};
    };
    constexpr int kTag = 4700;
    for (int round = 0; round < rounds; ++round) {
        t = clock::now();
        for (int peer = 0; peer < comm_size; ++peer) {
            const auto [e0, e1] = span(send_elems[static_cast<size_t>(peer)], round);
            if (e0 < e1) add_copies(sends[static_cast<size_t>(peer)], send_src[static_cast<size_t>(peer)], e0, e1,
                                    send_slot[static_cast<size_t>(peer)], true);
        }
        run_copies();
        stats.pack += std::chrono::duration<double>(clock::now() - t).count();
        t = clock::now();
        std::vector<MPI_Request> requests;
        for (int peer = 0; peer < comm_size; ++peer) {
            const auto [e0, e1] = span(recv_elems[static_cast<size_t>(peer)], round);
            if (e0 == e1) continue;
            void* at = device_mpi ? static_cast<void*>(recv_slot[static_cast<size_t>(peer)])
                                  : static_cast<void*>(host_recv[static_cast<size_t>(peer)].data());
            requests.emplace_back();
            MPI_Irecv(at, static_cast<int>((e1 - e0) * D), MPI_BYTE, peer, kTag, level_comm, &requests.back());
        }
        for (int peer = 0; peer < comm_size; ++peer) {
            const auto [e0, e1] = span(send_elems[static_cast<size_t>(peer)], round);
            if (e0 == e1) continue;
            void* at = send_slot[static_cast<size_t>(peer)];
            if (!device_mpi) {
                check_cuda(cudaMemcpy(host_send[static_cast<size_t>(peer)].data(), at, (e1 - e0) * D,
                                      cudaMemcpyDeviceToHost),
                           "CA device halo staging");
                at = host_send[static_cast<size_t>(peer)].data();
            }
            requests.emplace_back();
            MPI_Isend(at, static_cast<int>((e1 - e0) * D), MPI_BYTE, peer, kTag, level_comm, &requests.back());
        }
        MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
        if (!device_mpi) {
            for (int peer = 0; peer < comm_size; ++peer) {
                const auto [e0, e1] = span(recv_elems[static_cast<size_t>(peer)], round);
                if (e0 == e1) continue;
                check_cuda(cudaMemcpy(recv_slot[static_cast<size_t>(peer)], host_recv[static_cast<size_t>(peer)].data(),
                                      (e1 - e0) * D, cudaMemcpyHostToDevice),
                           "CA device halo staging");
            }
        }
        stats.mpi += std::chrono::duration<double>(clock::now() - t).count();
        t = clock::now();
        for (int peer = 0; peer < comm_size; ++peer) {
            const auto [e0, e1] = span(recv_elems[static_cast<size_t>(peer)], round);
            if (e0 < e1) add_copies(recvs[static_cast<size_t>(peer)], recv_dst[static_cast<size_t>(peer)], e0, e1,
                                    recv_slot[static_cast<size_t>(peer)], false);
        }
        run_copies();
        stats.unpack += std::chrono::duration<double>(clock::now() - t).count();
    }
    meta.clear();
    for (S* slot : send_slot) {
        if (slot != nullptr) (device_mpi ? arena : heap).free(slot);
    }
    for (S* slot : recv_slot) {
        if (slot != nullptr) (device_mpi ? arena : heap).free(slot);
    }
    stats.total = std::chrono::duration<double>(clock::now() - t0).count();
    return stats;
  }
}

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
