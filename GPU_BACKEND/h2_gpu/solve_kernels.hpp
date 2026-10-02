#pragma once
// Kernels of the device solve (H2 GPU backend): the per-box steps of the
// forward (V^{-1}) and backward (W^{-1}) Color sweeps, the diagonal X_RR
// solves, ordered accumulation of neighbor updates, and vector copies.  One
// thread block per box.  A level's vectors are one array holding each box's
// n x nrhs block (column-major) at nrhs * (its point offset); the work,
// ghost, and message buffers follow the same convention.  All offsets in the
// tables below are in points, so the tables do not depend on nrhs or on
// where a solve's buffers live.  The vectors and factors are real (double)
// or complex (dcomplex) as the launch's element type T; the factors of a box
// are untyped pointers in the tables.  See device_solve.hpp.

#ifdef H2_HAVE_GPU

#include <cuda_runtime.h>

#include "gpu_common/device_scalar.hpp"

#include <cstdint>

namespace fmm {
namespace gpu {

// A box of a level.  The factors are the host solve's: T (interpolation,
// k x r), X_SR (-X_SR X_RR^{-1}, k x r), X_NR (-X_NR X_RR^{-1}, ntot x r),
// the LU of X_RR (r x r) and its 1-based row interchanges, of element type T.
struct SolveBox {
    int64_t vec;          // its vector in the level vector
    int n, k, r, ntot;
    const int* skel;      // k positions in the box
    const int* red;       // r positions
    const void* T;        // ld k (null: no step)
    const void* xsr;      // ld k (null: no step)
    const void* xnr;      // ld ntot
    const void* lu;       // ld r
    const int* ipiv;
    int slot0, nslots;    // its one-hop slots
    int64_t work;         // its ntot rows in the work buffer of its wave
};

// One-hop neighbor of a box: rows [row0, row0 + count) of its X_NR.
struct SolveSlot {
    int row0, count;
    int full;             // all points of the neighbor (1), or its skeleton (0)
    int ghost;            // its vector is in the ghost area (another rank's box)
    int64_t vec;          // its vector (level vector or ghost area)
    int n;                // its points (column stride)
    const int* skel;      // its skeleton positions (full == 0)
};

// A term of an accumulation: count rows of work (row0 on in the block at
// base, column stride ld) or of the inbox (a block at base, ld = count).
struct SolvePart {
    int64_t base;
    int row0;
    int ld;
    int inbox;            // 1: from the inbox, 0: from the work buffer
};
// dst rows (count x nrhs) (+)= the sum of its parts, in order.
struct SolveAccum {
    int64_t dst;          // a vector (level vector), or a block of the outbox
    int ld;               // column stride of dst
    int count;
    const int* rows;      // dst row of each entry (a skeleton), or null for 0..count-1
    int outbox;           // 1: dst = sum into the outbox, 0: vector += sum
    int nparts;
    int64_t part0;        // its parts in the parts table
};

// A box vector (n x nrhs) copied between two buffers.
struct SolveCopy {
    int64_t src, dst;
    int n;
};

// Contiguous words (4 bytes) copied on the device.
struct CopySpan {
    const void* src;
    void* dst;
    int64_t words;
};

// Largest per-box dynamic shared memory the solve kernels may use.
constexpr size_t kSolveSharedLimit = 96 * 1024;

// forward: x_R -= T^T x_S;  x_S += X_SR x_R;  work = X_NR x_R
template<typename T>
void launch_solve_forward(const SolveBox* boxes, const int* wave, int count, T* vec, T* work, int nrhs, int max_n,
                          cudaStream_t stream);
// backward: x_R += X_SR^T x_S + X_NR^T x_N;  x_S -= T x_R (x_N gathered into work)
template<typename T>
void launch_solve_backward(const SolveBox* boxes, const SolveSlot* slots, const int* wave, int count, T* vec,
                           const T* ghost, T* work, int nrhs, int max_n, int max_ntot, cudaStream_t stream);
// x_R = X_RR^{-1} x_R (the LU)
template<typename T>
void launch_solve_diagonal(const SolveBox* boxes, const int* list, int count, T* vec, int nrhs, int max_r,
                           cudaStream_t stream);
template<typename T>
void launch_solve_accum(const SolveAccum* items, const SolvePart* parts, int count, int max_count, int nrhs, T* vec,
                        const T* work, T* outbox, const T* inbox, cudaStream_t stream);
template<typename T>
void launch_solve_copy(const SolveCopy* items, int count, int max_n, int nrhs, const T* src, T* dst,
                       cudaStream_t stream);
void launch_copy_spans(const CopySpan* items, int count, int64_t max_words, cudaStream_t stream);

// The multiply F x (the factorization applied), core/apply_mul.hpp:
// forward W:  x_S += T x_R;  x_R -= X_SR^T x_S + X_NR^T x_N  (x_N gathered into work)
template<typename T>
void launch_mul_forward(const SolveBox* boxes, const SolveSlot* slots, const int* wave, int count, T* vec,
                        const T* ghost, T* work, int nrhs, int max_n, int max_ntot, cudaStream_t stream);
// x_R = X_RR x_R (P L U)
template<typename T>
void launch_mul_diagonal(const SolveBox* boxes, const int* list, int count, T* vec, int nrhs, int max_r,
                         cudaStream_t stream);
// backward V:  x_S -= X_SR x_R;  work = -X_NR x_R;  x_R += T^T x_S
template<typename T>
void launch_mul_backward(const SolveBox* boxes, const int* wave, int count, T* vec, T* work, int nrhs, int max_n,
                         cudaStream_t stream);

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
