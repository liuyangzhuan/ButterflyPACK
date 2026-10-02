#pragma once
// Device-to-device messages between ranks: a message is a list of device
// blocks (column-major, rows x cols, each contiguous) laid back to back,
// moved in rounds of one chunk per message, packed straight from the blocks
// into a slot and unpacked straight from the slot into the destination
// blocks: no message-sized buffers, which a fragmented heap may not hold.
// Slots: the exchange arena with CUDA-aware MPI (registered once for MPI;
// buffers of the main heap made the transfers ~30x slower), else small heap
// slots staged through the host.  Every rank of `sync_comm` runs the same
// number of rounds with the same chunk, so each send meets its receive.
// Shared by the device halo of a CA level (ca_halo.hpp) and the device
// process reduction (device_reduction.hpp).

#ifdef H2_HAVE_GPU

#include "level_eliminator.hpp"

#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fmm {
namespace gpu {

struct PieceExchangeStats {
    double pack = 0.0, mpi = 0.0, unpack = 0.0;
    int rounds = 0;
    bool staged = false;  // through host buffers (MPI without device memory)
};

// A block of a message: rows x cols elements as stored on the device.
struct DevicePiece {
    int rows = 0, cols = 0;
    size_t elems() const { return static_cast<size_t>(rows) * static_cast<size_t>(cols); }
};

namespace device_pieces_detail {
// Messages in chunks under MPI's int count.
constexpr size_t kChunk = size_t{1} << 30;
}  // namespace device_pieces_detail

// `peers[i]`: a rank of `comm`; `sends[i]` / `send_src[i]`: the pieces sent
// to it and their device blocks; `recvs[i]` / `recv_dst[i]`: the pieces
// expected from it and the blocks they land in (both sides must list the
// same shapes in the same order).  Collective over `sync_comm` (the chunk
// and the number of rounds).
template<typename S>
PieceExchangeStats exchange_device_pieces(MPI_Comm comm, MPI_Comm sync_comm, const std::vector<int>& peers,
                                          const std::vector<std::vector<DevicePiece>>& sends,
                                          const std::vector<std::vector<S*>>& send_src,
                                          const std::vector<std::vector<DevicePiece>>& recvs,
                                          const std::vector<std::vector<S*>>& recv_dst, int tag,
                                          cudaStream_t stream, const char* what) {
    using clock = std::chrono::steady_clock;
    PieceExchangeStats stats;
    const size_t D = sizeof(S);
    const int np = static_cast<int>(peers.size());
    DeviceHeap& heap = DeviceHeap::instance();
    DeviceHeap& arena = DeviceHeap::exchange_arena();
    std::vector<size_t> send_elems(static_cast<size_t>(np), 0), recv_elems(static_cast<size_t>(np), 0);
    int messages = 0;
    for (int peer = 0; peer < np; ++peer) {
        for (const DevicePiece& p : sends[static_cast<size_t>(peer)]) send_elems[static_cast<size_t>(peer)] += p.elems();
        for (const DevicePiece& p : recvs[static_cast<size_t>(peer)]) recv_elems[static_cast<size_t>(peer)] += p.elems();
        messages += (send_elems[static_cast<size_t>(peer)] > 0) + (recv_elems[static_cast<size_t>(peer)] > 0);
    }
    const bool device_mpi = device_exchange_enabled() && arena.initialized();
    stats.staged = !device_mpi;
    // (one chunk size on every rank: a message's pieces must match on both ends)
    unsigned long long chunk_all = size_t{64} << 20;
    if (device_mpi) {
        chunk_all = std::min((arena.capacity() - arena.used()) / static_cast<size_t>(std::max(messages, 1)),
                             device_pieces_detail::kChunk) / 256 * 256;
    }
    MPI_Allreduce(MPI_IN_PLACE, &chunk_all, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, sync_comm);
    const size_t chunk = static_cast<size_t>(chunk_all);
    if (chunk < (size_t{1} << 20)) {
        throw std::runtime_error(std::string(what) + ": the exchange arena has no room (BPACK_GPU_EXCHANGE_MB)");
    }
    const size_t chunk_elems = chunk / D;
    int rounds = 0;
    for (int peer = 0; peer < np; ++peer) {
        const size_t most = std::max(send_elems[static_cast<size_t>(peer)], recv_elems[static_cast<size_t>(peer)]);
        rounds = std::max(rounds, static_cast<int>((most + chunk_elems - 1) / chunk_elems));
    }
    MPI_Allreduce(MPI_IN_PLACE, &rounds, 1, MPI_INT, MPI_MAX, sync_comm);
    stats.rounds = rounds;
    std::vector<S*> send_slot(static_cast<size_t>(np), nullptr), recv_slot(static_cast<size_t>(np), nullptr);
    std::vector<std::vector<char>> host_send(static_cast<size_t>(np)), host_recv(static_cast<size_t>(np));
    for (int peer = 0; peer < np; ++peer) {
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
    auto add_copies = [&](const std::vector<DevicePiece>& pieces, const std::vector<S*>& ptrs, size_t e0, size_t e1,
                          S* slot, bool to_slot) {
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
        check_cuda(cudaStreamSynchronize(stream), what);
        copies.clear();
        max_m = max_n = 0;
    };
    auto span = [&](size_t total, int round) -> std::pair<size_t, size_t> {
        const size_t e0 = std::min(total, static_cast<size_t>(round) * chunk_elems);
        return {e0, std::min(total, e0 + chunk_elems)};
    };
    for (int round = 0; round < rounds; ++round) {
        auto t = clock::now();
        for (int peer = 0; peer < np; ++peer) {
            const auto [e0, e1] = span(send_elems[static_cast<size_t>(peer)], round);
            if (e0 < e1) add_copies(sends[static_cast<size_t>(peer)], send_src[static_cast<size_t>(peer)], e0, e1,
                                    send_slot[static_cast<size_t>(peer)], true);
        }
        run_copies();
        stats.pack += std::chrono::duration<double>(clock::now() - t).count();
        t = clock::now();
        std::vector<MPI_Request> requests;
        for (int peer = 0; peer < np; ++peer) {
            const auto [e0, e1] = span(recv_elems[static_cast<size_t>(peer)], round);
            if (e0 == e1) continue;
            void* at = device_mpi ? static_cast<void*>(recv_slot[static_cast<size_t>(peer)])
                                  : static_cast<void*>(host_recv[static_cast<size_t>(peer)].data());
            requests.emplace_back();
            MPI_Irecv(at, static_cast<int>((e1 - e0) * D), MPI_BYTE, peers[static_cast<size_t>(peer)], tag, comm,
                      &requests.back());
        }
        for (int peer = 0; peer < np; ++peer) {
            const auto [e0, e1] = span(send_elems[static_cast<size_t>(peer)], round);
            if (e0 == e1) continue;
            void* at = send_slot[static_cast<size_t>(peer)];
            if (!device_mpi) {
                check_cuda(cudaMemcpy(host_send[static_cast<size_t>(peer)].data(), at, (e1 - e0) * D,
                                      cudaMemcpyDeviceToHost),
                           what);
                at = host_send[static_cast<size_t>(peer)].data();
            }
            requests.emplace_back();
            MPI_Isend(at, static_cast<int>((e1 - e0) * D), MPI_BYTE, peers[static_cast<size_t>(peer)], tag, comm,
                      &requests.back());
        }
        MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
        if (!device_mpi) {
            for (int peer = 0; peer < np; ++peer) {
                const auto [e0, e1] = span(recv_elems[static_cast<size_t>(peer)], round);
                if (e0 == e1) continue;
                check_cuda(cudaMemcpy(recv_slot[static_cast<size_t>(peer)], host_recv[static_cast<size_t>(peer)].data(),
                                      (e1 - e0) * D, cudaMemcpyHostToDevice),
                           what);
            }
        }
        stats.mpi += std::chrono::duration<double>(clock::now() - t).count();
        t = clock::now();
        for (int peer = 0; peer < np; ++peer) {
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
    return stats;
}

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
