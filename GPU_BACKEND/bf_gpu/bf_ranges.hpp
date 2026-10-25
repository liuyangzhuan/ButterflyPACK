#pragma once

#ifdef H2_HAVE_GPU

// What bf_gpu needs to learn from ButterflyPACK about a butterfly's
// distribution (GPU_BACKEND/bf_gpu).
//
// At one process every level has idx = 1 and inc = 1, so a global block index
// is its local index and the grids are the plain powers of two.  At several
// processes neither holds: a rank owns a strided subset of each level's
// blocks, the vector grids differ from the kernel grids, and below the
// sweep's midpoint a rank holding a single block column (or row) widens its
// range to the pair so the pair-wise broadcast has somewhere to land.  All of
// that is decided by GetLocalBlockRange, GetBlockPID and the sweep's own
// conversions, so it is read from there rather than recomputed here - the
// arithmetic is easy to get subtly wrong, and a subtly wrong grid gives a
// plausible product.
//
// The answers arrive through BfBridge rather than through direct calls to the
// C entry points, because those carry a precision prefix (d_c_bf_...,
// z_c_bf_... - see PrecisionPreprocessing.sh) while GPU_BACKEND is compiled
// once.  The precision-copied wrapper supplies the closures; a test harness
// supplies its own.  They are called when a plan is built, never per matvec.

#include <functional>
#include <mpi.h>

#include "bf_device.hpp"
#include "bf_vec.hpp"

namespace bf {
namespace gpu {

// One level's grid: the row and column ranges together.
struct Range2 {
    int idx_r = 1, inc_r = 1, nr = 0;
    int idx_c = 1, inc_c = 1, nc = 0;
};

struct BfBridge {
    MPI_Comm comm = MPI_COMM_NULL;
    int nproc = 1;
    // BF_exchange_matvec's rule over a level's whole local grid: dest[b] is
    // the zero-based rank in `comm` that local block b goes to, or -1 when
    // this rank does not send it.
    std::function<void(int level, char mode, char collect, const VecLevel& src, int* dest)> exchange_dest;
    // BF_all2all_vec_n_ker's rule (kerflag = 0), and the grid it leaves the
    // level in.
    std::function<void(int level, char mode_new, const VecLevel& src, int* dest)> all2all_dest;
    std::function<Range2(int level, char mode_new)> all2all_range;
    // The sub-group that owns one block of a level: its rank count, its
    // communicator and whether this rank heads it.  The closing product of
    // the sweep broadcasts its single vector block inside that sub-group when
    // the sub-group spans more than one rank.
    struct Group {
        int nproc = 1;
        MPI_Comm comm = MPI_COMM_NULL;
        bool is_head = true;
    };
    std::function<Group(int level, int index_i, int index_j, char mode)> block_group;

    bool usable() const {
        return nproc <= 1 || (exchange_dest && all2all_dest && all2all_range && block_group);
    }
};

// Fill in a resident butterfly's grids from the four range sets the Fortran
// returned: GetLocalBlockRange in both orderings (levels 0 .. L+1) and the
// vector grids of each orientation (vector levels 0 .. L+2).
inline void set_ranges(BfDevice& out, const SweepRanges& kernel_r, const SweepRanges& kernel_c,
                       const SweepRanges& vec_n, const SweepRanges& vec_t, int nproc, int myid) {
    const int want_k = out.level_butterfly + 2;
    const int want_v = out.level_butterfly + 3;
    if (kernel_r.levels() < want_k || kernel_c.levels() < want_k || vec_n.levels() < want_v ||
        vec_t.levels() < want_v) {
        throw std::invalid_argument("bf::gpu::set_ranges: a range set is too short for this butterfly");
    }
    out.kernel_r = kernel_r;
    out.kernel_c = kernel_c;
    out.vec_n = vec_n;
    out.vec_t = vec_t;
    out.nproc = nproc;
    out.myid = myid;
}

}  // namespace gpu
}  // namespace bf

#endif  // H2_HAVE_GPU
