#pragma once

#ifdef H2_HAVE_GPU

// The butterfly sweep's inter-level exchange, on the device
// (GPU_BACKEND/bf_gpu).  This is option (B) of the three we weighed: the
// Fortran keeps the decision - which block goes to which rank - and the data
// never leaves the device.
//
// BF_exchange_matvec (SRC/Bplus_utilities.f90) does three separable things:
//
//   1. decides a destination per local block, from an explicit rule plus
//      GetBlockPID, and sends only from the head of the sub-group holding the
//      block (its other ranks hold the same copy).  Reused here through
//      c_bf_get_block_pid and c_bf_get_sub_head, so there is one copy of the
//      rule, not two;
//   2. packs each block as [index_i, index_j, Nrow, Ncol, values...], the
//      header stored as DT values and read back with NINT;
//   3. exchanges the buffers and ACCUMULATES the arrivals into the target.
//
// Here (2) and (3) run on the device: a pack kernel gathers the blocks, MPI
// sends device pointers, and an unpack kernel adds.
//
// Two things the CPU code decides that are easy to get wrong here:
//
//   * it accumulates into whatever the target block already holds, and
//     zeroes a block only when the arrival is what first allocates it.  On
//     the device every block of a level has storage, so the exchange must
//     NOT zero: a rank that computed its own contribution to a target block
//     would lose it.  The precondition is instead that a level's buffer is
//     zero where no GEMM wrote it, which the sweep guarantees by zeroing at
//     allocation.
//   * a rank never messages itself; a block whose destination is this rank
//     is already in place.
//
// The receive buffers come from DeviceHeap::exchange_arena(), which exists
// because CUDA-aware MPI registers the whole allocation holding a buffer at
// a cost that grows with it; without GPU-aware MPI the buffers are staged
// through pinned host memory.  With NCCL (TPL_NCCL_LIBRARIES) the messages go
// through it instead, on the stream: measured on one Perlmutter node, Cray
// MPICH's GPU path turned each message into a pipeline of small copies and
// event polls (about 6,500 copies and 19,000 polls for one 256-column matvec
// of a 262144^2 butterfly on 4 GPUs, ~7 GB/s over NVLink).
//
// The plan is built once per (block, level, mode, collect) and reused: what
// depends on nrhs is only the per-peer sizes, which scale, so a matvec
// recomputes displacements on the host (microseconds) and runs no collective
// of its own.  A butterfly matvec is called tens of thousands of times in a
// STRUMPACK factorization, so a per-call collective here would cost more
// than the arithmetic.

#include <mpi.h>
#include <cstdint>
#include <vector>

#include "bf_ranges.hpp"
#include "bf_vec.hpp"
#include "../h2_gpu/device_heap.hpp"

namespace bf {
namespace gpu {

// Record header: index_i, index_j, Nrow, Ncol, as the CPU writes it.
constexpr int kHeader = 4;

// One block to send or to receive.  `row_offset` is into the level's device
// buffer in ROWS - the source for a send, the resolved target for a receive,
// which the receiver fills in at build time from its own offset table.  Rows
// rather than elements is what makes a record independent of nrhs.
// `ordinal` and `cumrows` place the record inside its peer's buffer:
//     slot = displ[peer] + kHeader * ordinal + cumrows * nrhs
// so the slots follow from nrhs without another exchange.
struct ExchangeRecord {
    int peer = 0;
    int index_i = 0, index_j = 0;
    int rows = 0;
    int ordinal = 0;
    int64_t cumrows = 0;
    int64_t row_offset = 0;
};

// A level's exchange.  Sizes are held per peer as a record count and a row
// sum, from which `displacements(nrhs)` gives the element counts.
struct ExchangePlan {
    int nproc = 1, myid = 0;
    MPI_Comm comm = MPI_COMM_NULL;

    std::vector<ExchangeRecord> send, recv;
    std::vector<int> send_recs, recv_recs;        // per peer, nproc entries
    std::vector<int64_t> send_rows, recv_rows;    // per peer
    std::vector<int> peers;                       // peers with traffic either way

    // device copies of the records, built once
    ExchangeRecord* d_send = nullptr;
    ExchangeRecord* d_recv = nullptr;
    int64_t* d_send_displ = nullptr;              // nproc + 1, refreshed per nrhs
    int64_t* d_recv_displ = nullptr;
    int displ_nrhs = 0;                           // the nrhs those displacements hold

    std::vector<int64_t> send_displ, recv_displ;  // host, nproc + 1
    int64_t send_total = 0, recv_total = 0;

    // the group's NCCL communicator (an ncclComm_t) when the library has
    // NCCL: the exchange then runs on the stream, with no host
    // synchronization, instead of through MPI.  Shared, never freed here.
    void* nccl = nullptr;

    bool empty() const { return peers.empty(); }

    // Element counts and displacements for this nrhs (host side).
    void displacements(int nrhs);
};

// Build the plan for one level from the Fortran's own answers, which arrive
// through `bridge`.  `level` is the sweep level, `mode` 'R' or 'C' and
// `collect` 'R' (reduction) or 'B' (broadcast), as BF_exchange_matvec takes
// them; `level_in`/`level_out` are the vector levels the records leave and
// enter.  Collective over the block's communicator.
ExchangePlan build_plan(const BfBridge& bridge, int level, char mode, char collect,
                        const VecLevel& level_in, VecLevel& level_out);

// The midpoint redistribution (BF_all2all_vec_n_ker with kerflag = 0), which
// moves a level between the row-wise and column-wise orderings.  `level` is
// the source level and `mode_new` the destination ordering, as that routine
// takes them.  `dst` is filled in here, grid and heights both: its grid is
// the target range, and the height of a block it receives is only known from
// the arrival, since no product on this rank writes it.  Collective.
ExchangePlan build_all2all_plan(const BfBridge& bridge, int level, char mode_new, const VecLevel& src,
                                VecLevel& dst);

void release_plan(ExchangePlan& plan);

// A vector block the ranks of a leaf's group hold together: a low-rank block
// (level_butterfly 0) on several ranks, each holding some rows of its one
// leaf.  (ButterflyPACK's trees give every rank at least two leaf boxes, so
// a butterfly proper never shares one.)  After the opening factor every rank
// of the group has its rows' part
// of the block, and the CPU sums them on all of them (MPI_ALLREDUCE); before
// the closing factor only the group's head has the block, and the CPU
// broadcasts it from there (MPI_Bcast from rank 0 of the group).  `vec` is
// the sweep's vector level holding the block (an index into SweepPlan::vec,
// -1 for its midpoint source) and `block` its index in that level's grid.
// The sum always goes through MPI, on the host: the same call as the CPU's,
// so a group of two ranks still matches it to the bit.  The broadcast goes
// through NCCL when the group has it.
struct GroupPlan {
    MPI_Comm comm = MPI_COMM_NULL;
    int nproc = 1;
    bool sum = false;
    int vec = -1;
    size_t block = 0;
    void* nccl = nullptr;
};

// Collective over the group (an NCCL communicator for a broadcast).
GroupPlan build_group_plan(MPI_Comm comm, bool sum, int vec, size_t block);

// The sum or the broadcast of n elements at data, on the device, ordered on
// the stream.
void group_level(const GroupPlan& g, DataType* data, int64_t n, cudaStream_t stream);

// Pack, exchange and accumulate, for `nrhs` right-hand sides.  `src` is the
// level the records come from and `dst` the one they go into; the CPU passes
// the same array for both, and so may a caller here.
void exchange_level(ExchangePlan& plan, const VecLevel& src, VecLevel& dst, int nrhs,
                    int tag, cudaStream_t stream);

// Kernels, implemented in bf_exchange.cu.
void launch_pack(const DataType* level, const ExchangeRecord* recs, int nrecs, const int64_t* displ,
                 int nrhs, DataType* buf, cudaStream_t stream);
void launch_unpack(DataType* level, const ExchangeRecord* recs, int nrecs, const int64_t* displ,
                   int nrhs, const DataType* buf, cudaStream_t stream);

}  // namespace gpu
}  // namespace bf

#endif  // H2_HAVE_GPU
