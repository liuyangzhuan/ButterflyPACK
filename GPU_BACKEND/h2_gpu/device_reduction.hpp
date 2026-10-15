#pragma once
// Device process reduction (h2_gpu/GPU_CA_PLAN.md, M7c).
//
// At a process reduction (the parent level has fewer active ranks) the
// device transition of every child rank has built its parents' blocks on
// the device: the Schur block of each of its parents, the near-field blocks
// among them, and its copies of the pairs with other ranks' parents.  Instead
// of copying them to the host for the host reduction (serialize, send,
// deserialize, re-upload at the next level's start), the blocks move between
// the devices to the parents' new owner (`level.parent_level_owner`): a host
// header per sender (kind, key, rows, cols per block), then the payloads
// through the exchange arena in chunked rounds (device_pieces.hpp).  A pair
// whose two parents' children live on two child ranks of the same owner is
// held by both; the rank holding the lower parent contributes it (the owner
// drops its own copy when it holds the higher one).  A pair across two owners
// goes to both, as the parent level's layout has it.  The owner merges the
// blocks into its own DeviceLevelBlocks, which the next level's eliminator
// adopts; the senders free theirs.  The host reduction keeps carrying the
// boxes' metadata (points, index and neighbor lists), now without block
// payloads.

#ifdef H2_HAVE_GPU

#include "device_pieces.hpp"

#include <mpi.h>

#include <chrono>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace fmm {
namespace gpu {

struct ReductionStats {
    double total = 0.0, plan = 0.0, pack = 0.0, mpi = 0.0, unpack = 0.0;
    double bytes_sent = 0.0, bytes_received = 0.0;
    int64_t blocks_sent = 0, blocks_received = 0, blocks_dropped = 0;  // dropped: copies another rank contributes
    int rounds = 0;
    bool staged = false;
};

// Collective over `level_comm` (the child level's active ranks: every one
// of them sends to its owner or is one).  `blocks` is the rank's parent
// blocks on entry; on return the owner's merged blocks, or null on a sender.
template<typename CoordType, typename DataType>
ReductionStats reduce_parent_blocks(ParallelTree<CoordType, DataType>* tree,
                                    const TreeLevel<CoordType, DataType>& level,
                                    const TreeLevel<CoordType, DataType>& parent_level,
                                    std::unique_ptr<DeviceLevelBlocks<DataType>>& blocks, MPI_Comm level_comm) {
  if constexpr (!gpu_data_type<DataType>) {
    (void)tree; (void)level; (void)parent_level; (void)blocks; (void)level_comm;
    throw std::runtime_error("reduce_parent_blocks: data type not run on the GPU");
  } else {
    using S = typename DeviceScalar<DataType>::type;
    using clock = std::chrono::steady_clock;
    const auto t0 = clock::now();
    ReductionStats stats;
    Context& ctx = Context::instance();
    ctx.activate();
    cudaStream_t stream = ctx.stream();
    DeviceHeap& heap = DeviceHeap::instance();
    const size_t D = sizeof(S);
    const int rank = tree->mpi_rank;
    MPI_Comm comm = tree->comm;
    const int dim = tree->dimension;
    const int nc = morton::children_per_box(dim);
    const int owner = level.parent_level_owner;
    const bool sender = level.is_process_active && owner != rank;
    const bool receiver = parent_level.is_process_active;
    if (level.is_process_active && !sender && !receiver) {
        throw std::runtime_error("reduce_parent_blocks: an active rank neither sends nor receives");
    }
    // the rank holding a parent's children (all in one region of the child
    // level), and a parent's owner at the parent level
    const RegionLayout Lc = RegionLayout::of(dim, level.level, level.num_active_processes);
    const RegionLayout Lp = RegionLayout::of(dim, level.level - 1, parent_level.num_active_processes);
    if (!Lc.uniform || !Lp.uniform || !level.regions_cubic || !parent_level.regions_cubic) {
        throw std::runtime_error("reduce_parent_blocks: the levels' regions are not uniform cubes");
    }
    auto child_rank_of = [&](int64_t pm) { return level.morton_to_rank.at(static_cast<int>(Lc.region_of(pm * nc))); };
    auto owner_of = [&](int64_t pm) { return parent_level.morton_to_rank.at(static_cast<int>(Lp.region_of(pm))); };
    // whether this rank contributes its copy of pair (lo, hi): a local pair
    // always; a pair with another child rank's parent when that rank has
    // another owner, or when this rank holds the lower parent
    auto contributes = [&](int64_t lo, int64_t hi) {
        const int rl = child_rank_of(lo), rh = child_rank_of(hi);
        if (rl == rh) return true;
        const bool mine_lo = rl == rank;
        const int64_t other = mine_lo ? hi : lo;
        if (owner_of(other) != owner) return true;
        return mine_lo;
    };
    auto edge_key = [](int64_t lo, int64_t hi) { return (static_cast<uint64_t>(lo) << 32) | static_cast<uint64_t>(hi); };

    // ---- the sender's list (its blocks, kept in order), the owner's drops
    struct Head { int64_t kind, lo, hi, rows, cols; };  // kind 0: Schur block (lo == hi), 1: pair
    std::vector<Head> heads;
    std::vector<DevicePiece> send_pieces;
    std::vector<S*> send_ptrs;
    if (level.is_process_active && blocks) {
        std::vector<int64_t> schur_keys;
        for (const auto& kv : blocks->schur) schur_keys.push_back(kv.first);
        std::sort(schur_keys.begin(), schur_keys.end());
        std::vector<uint64_t> edge_keys;
        for (const auto& kv : blocks->edges) edge_keys.push_back(kv.first);
        std::sort(edge_keys.begin(), edge_keys.end());
        if (sender) {
            for (int64_t pm : schur_keys) {
                const DeviceMatrixT<S>& m = blocks->schur.at(pm);
                if (m.ptr == nullptr) continue;
                heads.push_back({0, pm, pm, m.rows, m.cols});
                send_pieces.push_back({m.rows, m.cols});
                send_ptrs.push_back(m.ptr);
            }
        }
        for (uint64_t key : edge_keys) {
            const int64_t lo = static_cast<int64_t>(key >> 32), hi = static_cast<int64_t>(key & 0xffffffffu);
            const DeviceMatrixT<S>& m = blocks->edges.at(key);
            if (m.ptr == nullptr) continue;
            if (!contributes(lo, hi)) {
                if (receiver) {  // (the owner's own copy: the lower parent's rank sends it)
                    heap.free(m.ptr);
                    blocks->edges.erase(key);
                    ++stats.blocks_dropped;
                }
                continue;
            }
            if (sender) {
                heads.push_back({1, lo, hi, m.rows, m.cols});
                send_pieces.push_back({m.rows, m.cols});
                send_ptrs.push_back(m.ptr);
            }
        }
    }
    // ---- the headers: a count, then the list, from each child rank to its owner
    std::vector<int> peers;
    std::vector<std::vector<DevicePiece>> sends, recvs;
    std::vector<std::vector<S*>> send_src, recv_dst;
    std::vector<std::vector<int64_t>> recv_heads;
    {
        std::vector<MPI_Request> requests;
        std::vector<int64_t> counts;
        if (receiver) {
            for (int child : parent_level.children_senders) {
                if (child == rank) continue;
                peers.push_back(child);
            }
            counts.assign(peers.size(), 0);
            for (size_t i = 0; i < peers.size(); ++i) {
                requests.emplace_back();
                MPI_Irecv(&counts[i], 1, MPI_INT64_T, peers[i], 4800, comm, &requests.back());
            }
        }
        int64_t my_count = static_cast<int64_t>(heads.size());
        if (sender) {
            peers.push_back(owner);
            requests.emplace_back();
            MPI_Isend(&my_count, 1, MPI_INT64_T, owner, 4800, comm, &requests.back());
        }
        MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
        requests.clear();
        recv_heads.resize(peers.size());
        std::vector<int64_t> my_heads;
        if (receiver) {
            for (size_t i = 0; i < counts.size(); ++i) {
                recv_heads[i].resize(static_cast<size_t>(counts[i]) * 5);
                if (recv_heads[i].empty()) continue;
                requests.emplace_back();
                MPI_Irecv(recv_heads[i].data(), static_cast<int>(recv_heads[i].size()), MPI_INT64_T, peers[i], 4801,
                          comm, &requests.back());
            }
        }
        if (sender) {
            for (const Head& h : heads) {
                my_heads.push_back(h.kind);
                my_heads.push_back(h.lo);
                my_heads.push_back(h.hi);
                my_heads.push_back(h.rows);
                my_heads.push_back(h.cols);
            }
            if (!my_heads.empty()) {
                requests.emplace_back();
                MPI_Isend(my_heads.data(), static_cast<int>(my_heads.size()), MPI_INT64_T, owner, 4801, comm,
                          &requests.back());
            }
        }
        MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
    }
    // ---- the owner's new blocks, in the senders' order
    sends.assign(peers.size(), {});
    send_src.assign(peers.size(), {});
    recvs.assign(peers.size(), {});
    recv_dst.assign(peers.size(), {});
    if (receiver) {
        if (!blocks) blocks = std::make_unique<DeviceLevelBlocks<DataType>>();
        for (size_t i = 0; i < peers.size(); ++i) {
            const std::vector<int64_t>& h = recv_heads[i];
            for (size_t j = 0; j + 5 <= h.size(); j += 5) {
                const int64_t kind = h[j], lo = h[j + 1], hi = h[j + 2];
                const int rows = static_cast<int>(h[j + 3]), cols = static_cast<int>(h[j + 4]);
                if (static_cast<double>(rows) * cols * D > 1.0e9) {   // (8x diagnostics)
                    std::printf("  [gpu] device reduction: rank %d receives a %d x %d block (%.2f GB, kind %lld, %lld-%lld) from rank %d\n",
                                rank, rows, cols, static_cast<double>(rows) * cols * D / 1e9, static_cast<long long>(kind),
                                static_cast<long long>(lo), static_cast<long long>(hi), peers[i]);
                    std::fflush(stdout);
                }
                DeviceMatrixT<S> m{heap.alloc_resident<S>(std::max<size_t>(static_cast<size_t>(rows) * cols, 1) * D, ResidentKind::Reduction),
                                   rows, cols};
                const bool fresh = kind == 0 ? blocks->schur.emplace(lo, m).second
                                             : blocks->edges.emplace(edge_key(lo, hi), m).second;
                if (!fresh) {
                    throw std::runtime_error("device reduction: block (" + std::to_string(lo) + ", " +
                                             std::to_string(hi) + ") from rank " + std::to_string(peers[i]) +
                                             " is held already");
                }
                recvs[i].push_back({rows, cols});
                recv_dst[i].push_back(m.ptr);
                stats.bytes_received += static_cast<double>(m.rows) * m.cols * D;
                ++stats.blocks_received;
            }
        }
    }
    if (sender) {
        sends.back() = send_pieces;
        send_src.back() = send_ptrs;
        for (const DevicePiece& p : send_pieces) stats.bytes_sent += static_cast<double>(p.elems() * D);
        stats.blocks_sent = static_cast<int64_t>(send_pieces.size());
    }
    stats.plan = std::chrono::duration<double>(clock::now() - t0).count();

    // ---- the payloads
    const PieceExchangeStats x = exchange_device_pieces<S>(comm, level_comm, peers, sends, send_src, recvs, recv_dst,
                                                          4802, stream, "device reduction");
    stats.pack = x.pack;
    stats.mpi = x.mpi;
    stats.unpack = x.unpack;
    stats.rounds = x.rounds;
    stats.staged = x.staged;
    // ---- the senders' blocks are the owner's now
    if (sender && blocks) {
        for (auto& kv : blocks->schur) heap.free(kv.second.ptr);
        for (auto& kv : blocks->edges) heap.free(kv.second.ptr);
        blocks.reset();
    }
    stats.total = std::chrono::duration<double>(clock::now() - t0).count();
    return stats;
  }
}

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
