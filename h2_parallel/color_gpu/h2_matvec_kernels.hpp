#pragma once
// Kernels of the device matvec of the compression-only H2 (h2_matvec.hpp),
// the steps of hierarchical_h2_mul_parallel on one level.  One thread block
// per box.  A level's vector holds column c of box b at c * ld + b.vec
// (ld: the level's points), as the host matvec's input and output arrays
// hold the leaf level; compact skeleton vectors likewise at c * ldq + b.q.
// Factors and blocks are untyped pointers of the launch's element type T
// (double or dcomplex).

#ifdef H2_HAVE_GPU

#include <cuda_runtime.h>

#include "device_scalar.hpp"

#include <cstdint>

namespace fmm {
namespace gpu {

// A box of one level.
struct H2Box {
    int64_t vec;          // its points in the level vector
    int64_t q;            // its compact skeleton vector
    int n, k, r;
    const int* skel;      // k positions
    const int* red;       // r positions
    const void* T;        // interpolation (k x r, ld k), null when r == 0
    int block0, nblocks;  // its coupling blocks, in the host's order
    int near0, nnear;     // its near blocks (leaf level)
};

// A block of a target: K (rows x cols, ld rows; rows = the target's k for
// a coupling block, n for a near block) times the source vector whose
// column c is at base * nrhs + c * stride + off of the local (source 0),
// ghost (1) or partial (2) array; K null: that vector (rows entries) itself.
struct H2BlockRef {
    const void* K;
    int cols;
    int ghost;
    int64_t base;
    int64_t stride;
    int64_t off;
};

// A pair of local boxes t > s of a symmetric kernel, whose near blocks are
// K (rows t x cols s, ld rows) and its transpose: partial a = K x_s (rows
// entries at pa) and b = K^T x_t (cols entries at pb).
struct H2Pair {
    const void* K;
    int rows, cols;
    int64_t xt, xs;   // the boxes' vectors
    int64_t pa, pb;   // in the partial array
};

// A box's skeleton entries in the parent level: positions skel of the
// child's vector at `child` <-> the entries at `parent` of the parent's.
struct H2Handoff {
    int64_t child;
    int64_t parent;
    const int* skel;
    int k;
};

// dst[dst_base * nrhs + c * dst_stride + dst_off + i] = src[c * src_stride + src + i], i < len
struct H2Span {
    int64_t src;
    int64_t src_stride;
    int64_t dst_base;
    int64_t dst_stride;
    int64_t dst_off;
    int len;
};

// Largest per-box dynamic shared memory of these kernels.
constexpr size_t kH2MatvecSharedLimit = 96 * 1024;
constexpr int kH2AccumThreads = 256;

// x_S += T x_R (upward projection)
template<typename T>
void launch_h2_upward(const H2Box* boxes, int count, T* vec, int64_t ld, int nrhs, int max_r, cudaStream_t stream);
// q = x_S (compact skeleton vectors)
template<typename T>
void launch_h2_gather_q(const H2Box* boxes, int count, const T* vec, int64_t ld, T* q, int64_t ldq, int nrhs,
                        cudaStream_t stream);
// y[rows] += sum over the box's blocks (in order) of K x, rows = the
// skeleton (coupling blocks) or all points (near blocks, `near`)
template<typename T>
void launch_h2_accumulate(const H2Box* boxes, int count, const H2BlockRef* blocks, T* vec, int64_t ld, const T* local,
                          const T* ghost, const T* partial, int nrhs, bool near, int max_rows, int max_cols,
                          cudaStream_t stream);
// the partial vectors of the pairs (x: vectors at c * ldx + offset; partial
// array: c * ldp + offset)
template<typename T>
void launch_h2_pairs(const H2Pair* pairs, int count, const T* x, int64_t ldx, T* partial, int64_t ldp, int nrhs,
                     int max_cols, cudaStream_t stream);
// largest near-block width the pair kernel takes
constexpr int kH2PairMaxCols = 512;
// x_R += T^T x_S (downward interpolation)
template<typename T>
void launch_h2_downward(const H2Box* boxes, int count, T* vec, int64_t ld, int nrhs, int max_k, cudaStream_t stream);
// gather: parent entries = child skeleton entries; else scatter: child
// skeleton entries = parent entries
template<typename T>
void launch_h2_handoff(const H2Handoff* items, int count, T* child, int64_t ldc, T* parent, int64_t ldp, int nrhs,
                       bool gather, cudaStream_t stream);
template<typename T>
void launch_h2_pack(const H2Span* spans, int count, const T* src, T* dst, int nrhs, int max_len, cudaStream_t stream);

// Shared memory of launch_h2_accumulate for the largest block.
template<typename T>
inline size_t h2_accumulate_shared(int max_rows, int max_cols) {
    return (static_cast<size_t>(max_rows) + kH2AccumThreads + static_cast<size_t>(max_cols)) * sizeof(T);
}

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
