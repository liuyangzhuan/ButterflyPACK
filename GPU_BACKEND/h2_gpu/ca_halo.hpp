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
// ghost sets of the regions (from the level's layout, halo_geometry.hpp),
// the one-hop lists and the box sizes.  Blocks keep the device layout (rows = the higher box).
// The ghosts' blocks are allocated first (added to the level's blocks for
// the eliminator to adopt); the messages then move in rounds of one chunk
// per message, packed straight from the blocks into slots (the exchange
// arena with CUDA-aware MPI, else small heap slots staged through the host)
// and unpacked straight into the new blocks (device_pieces.hpp): no
// message-sized buffers, which a fragmented heap may not hold.

#ifdef H2_HAVE_GPU

#include "device_pieces.hpp"

#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <climits>
#include <map>
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

    // ---- who sends which ghosts to whom, from the level's region layout
    //      (halo_geometry.hpp, as the host gather at level entry), and the
    //      peers' ghost sets; messages over the tree's communicator
    const RegionLayout L = RegionLayout::of(tree->dimension, lvl, level.num_active_processes);
    if (!L.uniform || level.my_morton_id < 0) {
        throw std::runtime_error("CA device halo: the level's regions are not uniform cubes");
    }
    const uint32_t my_region = static_cast<uint32_t>(level.my_morton_id);
    MPI_Comm comm = tree->comm;
    const int comm_rank = tree->mpi_rank;
    auto rank_of_region = [&](uint32_t q) { return level.morton_to_rank.at(static_cast<int>(q)); };
    auto owner_rank = [&](int64_t m) { return rank_of_region(L.region_of(m)); };
    // the peers: ranks this one sends ghosts to (their ghost sets within
    // this region) or receives its ghosts from
    struct Peer {
        int rank = -1;
        std::vector<int64_t> requested;          // ghosts of this rank it needs
        std::unordered_set<int64_t> ghost_set;   // its whole ghost set
        std::vector<int64_t> mine;               // this rank's ghosts it owns
    };
    std::map<int, Peer> peer_of;
    for_each_region_within(L, my_region, halo_region_radius(L), [&](uint32_t q) {
        std::vector<int64_t> requested = ghosts_of_region_in(L, q, my_region);
        if (requested.empty()) return;
        Peer& p = peer_of[rank_of_region(q)];
        p.rank = rank_of_region(q);
        p.requested = std::move(requested);
        const RegionGhosts rg = ghosts_of_region(L, q);
        p.ghost_set.insert(rg.ghosts.begin(), rg.ghosts.end());
    });
    for (int64_t g : level.ghost_id) {
        Peer& p = peer_of[owner_rank(g)];
        p.rank = owner_rank(g);
        p.mine.push_back(g);
    }
    std::vector<Peer*> peers;
    for (auto& kv : peer_of) peers.push_back(&kv.second);
    const int np = static_cast<int>(peers.size());
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

    // sent: the ghosts the peers need from this rank (their shapes from the
    // device blocks)
    std::vector<std::vector<Piece>> sends(static_cast<size_t>(np));
    std::vector<size_t> send_elems(static_cast<size_t>(np), 0);
    for (int peer = 0; peer < np; ++peer) {
        const Peer& p = *peers[static_cast<size_t>(peer)];
        if (p.requested.empty()) continue;
        for (int64_t g : p.requested) {
            const auto* box = level.find_local_box(g);
            if (box == nullptr) throw std::runtime_error("CA device halo: requested ghost " + std::to_string(g) + " is not local");
            payload(
                g, static_cast<int>(box->num_points), box->one_hop,
                [&](int64_t x) { return owner_rank(x) == p.rank; },
                [&](int64_t x) { return p.ghost_set.count(x) != 0; },
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
        for (const Piece& p2 : sends[static_cast<size_t>(peer)]) send_elems[static_cast<size_t>(peer)] += p2.elems();
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
    std::vector<std::vector<Piece>> recvs(static_cast<size_t>(np));
    std::vector<size_t> recv_elems(static_cast<size_t>(np), 0);
    for (int peer = 0; peer < np; ++peer) {
        const Peer& p = *peers[static_cast<size_t>(peer)];
        for (int64_t g : p.mine) {
            const auto* gb = level.find_ghost_box(g);
            if (gb == nullptr) throw std::runtime_error("CA device halo: ghost " + std::to_string(g) + " has no metadata");
            payload(
                g, static_cast<int>(gb->num_points), gb->one_hop,
                [&](int64_t x) { return level.find_local_box(x) != nullptr; },
                [&](int64_t x) { return my_requests.count(x) != 0; }, size_on_this_rank,
                recvs[static_cast<size_t>(peer)]);
        }
        for (const Piece& p2 : recvs[static_cast<size_t>(peer)]) recv_elems[static_cast<size_t>(peer)] += p2.elems();
    }
    // both sides must agree on every message's size (point to point with
    // the peers: a message's receiver expects it exactly when its sender
    // computed it, from the same layout)
    {
        std::vector<int64_t> out(send_elems.begin(), send_elems.end()), in(static_cast<size_t>(np), 0);
        std::vector<MPI_Request> requests;
        for (int peer = 0; peer < np; ++peer) {
            if (recv_elems[static_cast<size_t>(peer)] == 0) continue;
            requests.emplace_back();
            MPI_Irecv(&in[static_cast<size_t>(peer)], 1, MPI_INT64_T, peers[static_cast<size_t>(peer)]->rank, 4699, comm,
                      &requests.back());
        }
        for (int peer = 0; peer < np; ++peer) {
            if (send_elems[static_cast<size_t>(peer)] == 0) continue;
            requests.emplace_back();
            MPI_Isend(&out[static_cast<size_t>(peer)], 1, MPI_INT64_T, peers[static_cast<size_t>(peer)]->rank, 4699, comm,
                      &requests.back());
        }
        MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
        for (int peer = 0; peer < np; ++peer) {
            if (static_cast<size_t>(in[static_cast<size_t>(peer)]) != recv_elems[static_cast<size_t>(peer)]) {
                throw std::runtime_error("CA device halo: rank " + std::to_string(comm_rank) + " expects " +
                                         std::to_string(recv_elems[static_cast<size_t>(peer)]) +
                                         " elements from rank " + std::to_string(peers[static_cast<size_t>(peer)]->rank) +
                                         ", which sends " + std::to_string(in[static_cast<size_t>(peer)]));
            }
        }
    }
    stats.plan = std::chrono::duration<double>(clock::now() - t0).count();

    // ---- the device blocks of the messages; the ghosts' new blocks
    auto t = clock::now();
    std::vector<std::vector<S*>> send_src(static_cast<size_t>(np)), recv_dst(static_cast<size_t>(np));
    for (int peer = 0; peer < np; ++peer) {
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
    for (int peer = 0; peer < np; ++peer) {
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

    // ---- the payloads: rounds of one chunk per message through the
    //      exchange arena (device_pieces.hpp)
    std::vector<int> peer_ranks;
    std::vector<std::vector<DevicePiece>> send_shapes(static_cast<size_t>(np)), recv_shapes(static_cast<size_t>(np));
    for (int peer = 0; peer < np; ++peer) {
        peer_ranks.push_back(peers[static_cast<size_t>(peer)]->rank);
        for (const Piece& p : sends[static_cast<size_t>(peer)]) send_shapes[static_cast<size_t>(peer)].push_back({p.rows, p.cols});
        for (const Piece& p : recvs[static_cast<size_t>(peer)]) recv_shapes[static_cast<size_t>(peer)].push_back({p.rows, p.cols});
    }
    const PieceExchangeStats x = exchange_device_pieces<S>(comm, level_comm, peer_ranks, send_shapes, send_src,
                                                          recv_shapes, recv_dst, 4700, stream, "CA device halo");
    stats.pack = x.pack;
    stats.mpi = x.mpi;
    stats.unpack = x.unpack;
    stats.staged = x.staged;
    stats.total = std::chrono::duration<double>(clock::now() - t0).count();
    return stats;
  }
}

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
