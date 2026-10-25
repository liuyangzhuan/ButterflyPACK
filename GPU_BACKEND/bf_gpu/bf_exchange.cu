// The butterfly sweep's inter-level exchange on the device.  See
// bf_exchange.hpp for what this mirrors and for the two semantics that are
// easy to get wrong (accumulate, never zero; no self-messages).

#include "bf_exchange.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>

namespace bf {
namespace gpu {

using fmm::gpu::DeviceHeap;
using fmm::gpu::check_cuda;
using fmm::gpu::device_exchange_enabled;

namespace {

constexpr int kThreads = 128;

void check_launch(const char* what) {
    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA launch failed in ") + what + ": " + cudaGetErrorString(status));
    }
}

// Where a record sits in its peer's buffer.  The header is kHeader values and
// the payload rows * nrhs, so the ordinal and the row sum before it place it
// without an exchange of offsets.
__host__ __device__ __forceinline__ int64_t slot_of(const ExchangeRecord& r, const int64_t* displ, int nrhs) {
    return displ[r.peer] + static_cast<int64_t>(kHeader) * r.ordinal + r.cumrows * nrhs;
}

// One block per record: thread 0 writes the header, the block copies the
// values.  A level's block is contiguous (rows x nrhs, lda = rows), so the
// copy is flat.
__global__ void __launch_bounds__(kThreads) pack_kernel(const DataType* level, const ExchangeRecord* recs,
                                                        const int64_t* displ, int nrhs, DataType* buf) {
    const ExchangeRecord r = recs[blockIdx.x];
    DataType* out = buf + slot_of(r, displ, nrhs);
    if (threadIdx.x == 0) {
        out[0] = static_cast<DataType>(r.index_i);
        out[1] = static_cast<DataType>(r.index_j);
        out[2] = static_cast<DataType>(r.rows);
        out[3] = static_cast<DataType>(nrhs);
    }
    const DataType* src = level + r.row_offset * nrhs;
    const int64_t n = static_cast<int64_t>(r.rows) * nrhs;
    for (int64_t i = threadIdx.x; i < n; i += blockDim.x) out[kHeader + i] = src[i];
}

// The mirror, adding into the target - the CPU accumulates too, and a rank
// may already hold its own contribution to the block.
__global__ void __launch_bounds__(kThreads) unpack_kernel(DataType* level, const ExchangeRecord* recs,
                                                          const int64_t* displ, int nrhs, const DataType* buf) {
    const ExchangeRecord r = recs[blockIdx.x];
    const DataType* in = buf + slot_of(r, displ, nrhs);
    DataType* dst = level + r.row_offset * nrhs;
    const int64_t n = static_cast<int64_t>(r.rows) * nrhs;
    for (int64_t i = threadIdx.x; i < n; i += blockDim.x) dst[i] += in[kHeader + i];
}

template<typename T>
T* device_copy(const std::vector<T>& host) {
    if (host.empty()) return nullptr;
    const size_t bytes = host.size() * sizeof(T);
    T* d = reinterpret_cast<T*>(DeviceHeap::instance().alloc_resident<char>(bytes, fmm::gpu::ResidentKind::Matvec));
    check_cuda(cudaMemcpy(d, host.data(), bytes, cudaMemcpyHostToDevice), "bf exchange plan upload");
    return d;
}

}  // namespace

void ExchangePlan::displacements(int nrhs) {
    send_displ.assign(static_cast<size_t>(nproc) + 1, 0);
    recv_displ.assign(static_cast<size_t>(nproc) + 1, 0);
    for (int p = 0; p < nproc; ++p) {
        const size_t u = static_cast<size_t>(p);
        send_displ[u + 1] = send_displ[u] + static_cast<int64_t>(kHeader) * send_recs[u] + send_rows[u] * nrhs;
        recv_displ[u + 1] = recv_displ[u] + static_cast<int64_t>(kHeader) * recv_recs[u] + recv_rows[u] * nrhs;
    }
    send_total = send_displ[static_cast<size_t>(nproc)];
    recv_total = recv_displ[static_cast<size_t>(nproc)];
}

namespace {

// Build the records and exchange the descriptors.  `dest[b]` is the
// destination rank of local block b of `src`, or -1 when this rank does not
// send it, as the Fortran decided.
ExchangePlan assemble(MPI_Comm comm, int nproc, const std::vector<int>& dest, const VecLevel& src,
                      VecLevel& dst) {
    ExchangePlan plan;
    plan.comm = comm;
    plan.nproc = nproc;
    MPI_Comm_rank(comm, &plan.myid);
    plan.send_recs.assign(static_cast<size_t>(nproc), 0);
    plan.recv_recs.assign(static_cast<size_t>(nproc), 0);
    plan.send_rows.assign(static_cast<size_t>(nproc), 0);
    plan.recv_rows.assign(static_cast<size_t>(nproc), 0);

    // The source block of each record is kept beside it: the pair-wise
    // exchange passes the same level as source and target, so a target slot
    // sized by an arrival below re-packs the level and moves the sources too.
    std::vector<std::pair<ExchangeRecord, size_t>> sends;
    for (int i = 0; i < src.nr; ++i) {
        for (int j = 0; j < src.nc; ++j) {
            const size_t b = static_cast<size_t>(i) * src.nc + j;
            if (src.rows[b] <= 0 || dest[b] < 0) continue;
            ExchangeRecord r;
            r.peer = dest[b];
            r.index_i = i * src.inc_r + src.idx_r;
            r.index_j = j * src.inc_c + src.idx_c;
            r.rows = src.rows[b];
            sends.emplace_back(r, b);
        }
    }
    // by peer, so each peer's records are one contiguous message
    std::stable_sort(sends.begin(), sends.end(),
                     [](const std::pair<ExchangeRecord, size_t>& a,
                        const std::pair<ExchangeRecord, size_t>& b) { return a.first.peer < b.first.peer; });
    std::vector<size_t> send_block;
    send_block.reserve(sends.size());
    plan.send.reserve(sends.size());
    for (auto& e : sends) {
        plan.send.push_back(e.first);
        send_block.push_back(e.second);
    }
    for (ExchangeRecord& r : plan.send) {
        const size_t p = static_cast<size_t>(r.peer);
        r.ordinal = plan.send_recs[p];
        r.cumrows = plan.send_rows[p];
        plan.send_recs[p] += 1;
        plan.send_rows[p] += r.rows;
    }

    // The receiver needs each arrival's (index_i, index_j, rows) to place it;
    // one exchange of descriptors here, none per matvec.
    MPI_Alltoall(plan.send_recs.data(), 1, MPI_INT, plan.recv_recs.data(), 1, MPI_INT, comm);
    std::vector<int> scount(static_cast<size_t>(nproc), 0), rcount(static_cast<size_t>(nproc), 0);
    std::vector<int> sdispl(static_cast<size_t>(nproc), 0), rdispl(static_cast<size_t>(nproc), 0);
    int stotal = 0, rtotal = 0;
    for (int p = 0; p < nproc; ++p) {
        const size_t u = static_cast<size_t>(p);
        scount[u] = 3 * plan.send_recs[u];
        rcount[u] = 3 * plan.recv_recs[u];
        sdispl[u] = stotal;
        rdispl[u] = rtotal;
        stotal += scount[u];
        rtotal += rcount[u];
    }
    std::vector<int> sdesc(static_cast<size_t>(std::max(stotal, 1))), rdesc(static_cast<size_t>(std::max(rtotal, 1)));
    {
        std::vector<int> at = sdispl;
        for (const ExchangeRecord& r : plan.send) {
            int& k = at[static_cast<size_t>(r.peer)];
            sdesc[static_cast<size_t>(k)] = r.index_i;
            sdesc[static_cast<size_t>(k) + 1] = r.index_j;
            sdesc[static_cast<size_t>(k) + 2] = r.rows;
            k += 3;
        }
    }
    MPI_Alltoallv(sdesc.data(), scount.data(), sdispl.data(), MPI_INT, rdesc.data(), rcount.data(), rdispl.data(),
                  MPI_INT, comm);

    for (int p = 0; p < nproc; ++p) {
        for (int k = 0; k < plan.recv_recs[static_cast<size_t>(p)]; ++k) {
            const size_t at = static_cast<size_t>(rdispl[static_cast<size_t>(p)] + 3 * k);
            ExchangeRecord r;
            r.peer = p;
            r.index_i = rdesc[at];
            r.index_j = rdesc[at + 1];
            r.rows = rdesc[at + 2];
            r.ordinal = k;
            r.cumrows = plan.recv_rows[static_cast<size_t>(p)];
            plan.recv_rows[static_cast<size_t>(p)] += r.rows;
            plan.recv.push_back(r);
        }
    }

    // the target is local, so placing an arrival needs no further exchange
    auto block_of = [&](const ExchangeRecord& r) {
        const int i = (r.index_i - dst.idx_r) / std::max(dst.inc_r, 1);
        const int j = (r.index_j - dst.idx_c) / std::max(dst.inc_c, 1);
        if (i < 0 || i >= dst.nr || j < 0 || j >= dst.nc) {
            throw std::runtime_error("bf::gpu::assemble: arrival (" + std::to_string(r.index_i) + ", " +
                                     std::to_string(r.index_j) + ") is outside the target grid");
        }
        return static_cast<size_t>(i) * dst.nc + j;
    };
    // A slot this exchange is the first to fill is sized by the arrival: the
    // pair-wise broadcast exists precisely to deliver the partner block that
    // no local product writes, and the redistribution delivers the whole
    // level.  A slot this rank also computes must agree with the arrival
    // instead - that is where a mismatch means a real disagreement.
    bool grew = false;
    for (const ExchangeRecord& r : plan.recv) {
        const size_t b = block_of(r);
        if (dst.rows[b] == 0) {
            dst.rows[b] = r.rows;
            grew = true;
        } else if (dst.rows[b] != r.rows) {
            throw std::runtime_error("bf::gpu::assemble: arrival of " + std::to_string(r.rows) +
                                     " rows into a block of " + std::to_string(dst.rows[b]));
        }
    }
    if (grew) dst.pack(src.nrhs);
    // resolved last, so both sides read the final layout
    for (size_t k = 0; k < plan.send.size(); ++k) plan.send[k].row_offset = src.row_offset[send_block[k]];
    for (ExchangeRecord& r : plan.recv) r.row_offset = dst.row_offset[block_of(r)];

    for (int p = 0; p < nproc; ++p) {
        if (plan.send_recs[static_cast<size_t>(p)] > 0 || plan.recv_recs[static_cast<size_t>(p)] > 0) {
            plan.peers.push_back(p);
        }
    }
    plan.d_send = device_copy(plan.send);
    plan.d_recv = device_copy(plan.recv);
    return plan;
}

}  // namespace

ExchangePlan build_plan(const BfBridge& bridge, int level, char mode, char collect,
                        const VecLevel& level_in, VecLevel& level_out) {
    if (bridge.nproc <= 1) {
        ExchangePlan plan;
        plan.comm = bridge.comm;
        plan.nproc = 1;
        return plan;
    }
    // The Fortran decides, for the whole local grid at once, which block goes
    // where; -1 means this rank does not send it.
    std::vector<int> dest(static_cast<size_t>(std::max(level_in.nr * level_in.nc, 1)), -1);
    bridge.exchange_dest(level, mode, collect, level_in, dest.data());
    return assemble(bridge.comm, bridge.nproc, dest, level_in, level_out);
}

ExchangePlan build_all2all_plan(const BfBridge& bridge, int level, char mode_new, const VecLevel& src,
                                VecLevel& dst) {
    if (bridge.nproc <= 1) {
        // one rank: the orderings coincide, so the level is already in place
        dst = src;
        ExchangePlan plan;
        plan.comm = bridge.comm;
        plan.nproc = 1;
        return plan;
    }
    // the level's grid after the redistribution, which is also what the
    // second half of the sweep reads it as
    const Range2 r = bridge.all2all_range(level, mode_new);
    dst.idx_r = r.idx_r; dst.inc_r = r.inc_r; dst.nr = r.nr;
    dst.idx_c = r.idx_c; dst.inc_c = r.inc_c; dst.nc = r.nc;
    dst.rows.assign(static_cast<size_t>(std::max(r.nr, 0)) * static_cast<size_t>(std::max(r.nc, 0)), 0);
    dst.pack(src.nrhs);
    std::vector<int> dest(static_cast<size_t>(std::max(src.nr * src.nc, 1)), -1);
    bridge.all2all_dest(level, mode_new, src, dest.data());
    return assemble(bridge.comm, bridge.nproc, dest, src, dst);
}

void release_plan(ExchangePlan& plan) {
    DeviceHeap& heap = DeviceHeap::instance();
    for (void* p : {static_cast<void*>(plan.d_send), static_cast<void*>(plan.d_recv),
                    static_cast<void*>(plan.d_send_displ), static_cast<void*>(plan.d_recv_displ)}) {
        if (p != nullptr) heap.free(reinterpret_cast<char*>(p));
    }
    plan.d_send = nullptr;
    plan.d_recv = nullptr;
    plan.d_send_displ = nullptr;
    plan.d_recv_displ = nullptr;
    plan.displ_nrhs = 0;
}

void launch_pack(const DataType* level, const ExchangeRecord* recs, int nrecs, const int64_t* displ, int nrhs,
                 DataType* buf, cudaStream_t stream) {
    if (nrecs <= 0) return;
    pack_kernel<<<nrecs, kThreads, 0, stream>>>(level, recs, displ, nrhs, buf);
    check_launch("bf exchange pack");
}

void launch_unpack(DataType* level, const ExchangeRecord* recs, int nrecs, const int64_t* displ, int nrhs,
                   const DataType* buf, cudaStream_t stream) {
    if (nrecs <= 0) return;
    unpack_kernel<<<nrecs, kThreads, 0, stream>>>(level, recs, displ, nrhs, buf);
    check_launch("bf exchange unpack");
}

void exchange_level(ExchangePlan& plan, const VecLevel& src, VecLevel& dst, int nrhs, int tag,
                    cudaStream_t stream) {
    if (plan.empty()) return;
    plan.displacements(nrhs);
    const size_t displ_bytes = (static_cast<size_t>(plan.nproc) + 1) * sizeof(int64_t);
    if (plan.d_send_displ == nullptr) {
        DeviceHeap& heap = DeviceHeap::instance();
        plan.d_send_displ = reinterpret_cast<int64_t*>(
            heap.alloc_resident<char>(displ_bytes, fmm::gpu::ResidentKind::Matvec));
        plan.d_recv_displ = reinterpret_cast<int64_t*>(
            heap.alloc_resident<char>(displ_bytes, fmm::gpu::ResidentKind::Matvec));
        plan.displ_nrhs = 0;
    }
    if (plan.displ_nrhs != nrhs) {
        check_cuda(cudaMemcpyAsync(plan.d_send_displ, plan.send_displ.data(), displ_bytes, cudaMemcpyHostToDevice,
                                   stream), "bf exchange displacements");
        check_cuda(cudaMemcpyAsync(plan.d_recv_displ, plan.recv_displ.data(), displ_bytes, cudaMemcpyHostToDevice,
                                   stream), "bf exchange displacements");
        plan.displ_nrhs = nrhs;
    }

    // Message buffers: in the exchange arena with GPU-aware MPI (registered
    // once; sent and received in device memory), else from the heap and
    // staged through pinned host memory.
    DeviceHeap& arena = DeviceHeap::exchange_arena();
    DeviceHeap& heap = DeviceHeap::instance();
    const size_t sbytes = static_cast<size_t>(std::max<int64_t>(plan.send_total, 1)) * sizeof(DataType);
    const size_t rbytes = static_cast<size_t>(std::max<int64_t>(plan.recv_total, 1)) * sizeof(DataType);
    DataType* sdev = nullptr;
    DataType* rdev = nullptr;
    bool from_arena = false;
    if (device_exchange_enabled() && arena.initialized()) {
        char* s = arena.try_alloc(sbytes);
        char* r = s != nullptr ? arena.try_alloc(rbytes) : nullptr;
        if (r != nullptr) {
            sdev = reinterpret_cast<DataType*>(s);
            rdev = reinterpret_cast<DataType*>(r);
            from_arena = true;
        } else if (s != nullptr) {
            arena.free(s);
        }
    }
    if (!from_arena) {
        sdev = reinterpret_cast<DataType*>(heap.alloc<char>(sbytes));
        rdev = reinterpret_cast<DataType*>(heap.alloc<char>(rbytes));
    }

    launch_pack(src.data, plan.d_send, static_cast<int>(plan.send.size()), plan.d_send_displ, nrhs, sdev, stream);

    // staging buffers kept for the run (never freed: see PinnedPool)
    static fmm::gpu::PinnedBuffer* send_host = new fmm::gpu::PinnedBuffer;
    static fmm::gpu::PinnedBuffer* recv_host = new fmm::gpu::PinnedBuffer;
    DataType* sbuf = sdev;
    DataType* rbuf = rdev;
    if (!from_arena) {
        sbuf = static_cast<DataType*>(send_host->reserve(sbytes));
        rbuf = static_cast<DataType*>(recv_host->reserve(rbytes));
    }
    check_cuda(cudaStreamSynchronize(stream), "bf exchange pack");
    if (!from_arena && plan.send_total > 0) {
        check_cuda(cudaMemcpy(sbuf, sdev, static_cast<size_t>(plan.send_total) * sizeof(DataType),
                              cudaMemcpyDeviceToHost), "bf exchange stage out");
    }

    constexpr int kWords = static_cast<int>(sizeof(DataType) / sizeof(double));
    std::vector<MPI_Request> reqs;
    reqs.reserve(2 * plan.peers.size());
    for (int p : plan.peers) {
        const size_t u = static_cast<size_t>(p);
        const int64_t count = plan.recv_displ[u + 1] - plan.recv_displ[u];
        if (count > 0) {
            reqs.emplace_back();
            MPI_Irecv(rbuf + plan.recv_displ[u], static_cast<int>(count * kWords), MPI_DOUBLE, p, tag, plan.comm,
                      &reqs.back());
        }
    }
    for (int p : plan.peers) {
        const size_t u = static_cast<size_t>(p);
        const int64_t count = plan.send_displ[u + 1] - plan.send_displ[u];
        if (count > 0) {
            reqs.emplace_back();
            MPI_Isend(sbuf + plan.send_displ[u], static_cast<int>(count * kWords), MPI_DOUBLE, p, tag, plan.comm,
                      &reqs.back());
        }
    }
    MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);

    if (!from_arena && plan.recv_total > 0) {
        check_cuda(cudaMemcpyAsync(rdev, rbuf, static_cast<size_t>(plan.recv_total) * sizeof(DataType),
                                   cudaMemcpyHostToDevice, stream), "bf exchange stage in");
    }
    launch_unpack(dst.data, plan.d_recv, static_cast<int>(plan.recv.size()), plan.d_recv_displ, nrhs, rdev, stream);
    check_cuda(cudaStreamSynchronize(stream), "bf exchange unpack");

    if (from_arena) {
        arena.free(reinterpret_cast<char*>(rdev));
        arena.free(reinterpret_cast<char*>(sdev));
    } else {
        heap.free(reinterpret_cast<char*>(rdev));
        heap.free(reinterpret_cast<char*>(sdev));
    }
}

}  // namespace gpu
}  // namespace bf
