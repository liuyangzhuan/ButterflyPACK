#pragma once
// Batched element kernels of the H2 Color GPU backend (device_kernels.cu).
//
// Every launch processes a list of items that lives in device memory.  Items
// hold absolute pointers into device matrices; index lists are byte offsets
// into a per-wave metadata block, whose device address is passed with the
// launch.  All matrices are column-major.  Items and launches are templates
// on the element type T (double or dcomplex, device_scalar.hpp); the names
// without the T suffix are the real items.

#ifdef H2_HAVE_GPU

#include "gpu_common/device_scalar.hpp"
#include "../../GPU_INTERFACE/bpack_gpu_kernels.cuh"  // PointTable, IndexList, EvalItemT, OrderedSketchItemT

#include <cuda_runtime.h>

#include <cstdint>

namespace fmm {
namespace gpu {

// The entries of the matrix come from the application's evaluator
// (evaluator.hpp); its kernels, and the items they read (PointTable,
// IndexList, EvalItemT, OrderedSketchItemT), are in
// GPU_INTERFACE/bpack_gpu_kernels.cuh.
using EvalItem = EvalItemT<double>;

// out(i, j) = src[rows[i] * rs + cols[j] * cs]
template<typename T>
struct GatherItemT {
    T* out;
    int ld;
    int m;
    int n;
    const T* src;
    int64_t rs;
    int64_t cs;
    IndexList rows;
    IndexList cols;
};
using GatherItem = GatherItemT<double>;

// out[out_rows[i] * ors + j * ocs] = src[src_rows[i] * rs + j * cs], or += with
// add (i < m, j < n): rows gathered and scattered in one.  The items of a
// launch write distinct entries.
template<typename T>
struct IndexCopyItemT {
    T* out;
    int64_t ors;
    int64_t ocs;
    IndexList out_rows;
    const T* src;
    int64_t rs;
    int64_t cs;
    IndexList src_rows;
    int m;
    int n;
    int add;
};

// out[out_rows[i] * ors + out_cols[j] * ocs] = src[src_rows[i] * rs +
// src_cols[j] * cs] (i < m, j < n), or += with add: rows and columns indexed
// on both sides, so one launch covers a gather into a dense block, the
// scatter of a block into indexed rows and columns of an output, and the two
// in one.  An index list left default is the identity (index_at), and a
// transposed source is the same item with rs and cs swapped.  The items of a
// launch write distinct entries, except with add 2 (atomic adds).  row_w and
// col_w (byte offsets in meta, -1: none) weight entry (i, j) by row_w[i] and
// col_w[j] (doubles): a source row or column listed several times, each
// with its weight, sketches the block.
template<typename T>
struct BlockCopyItemT {
    T* out;
    int64_t ors;
    int64_t ocs;
    IndexList out_rows;
    IndexList out_cols;
    const T* src;
    int64_t rs;
    int64_t cs;
    IndexList src_rows;
    IndexList src_cols;
    int m;
    int n;
    int add;
    int64_t row_w = -1;
    int64_t col_w = -1;
};
using BlockCopyItem = BlockCopyItemT<double>;

// out(i, j) = a(rows[i], cols[j]) of a column-major matrix (leading dimension
// lda) held as chunks of `chunk` elements, chunk k at pool + slot[k] * chunk
// (the chunks a gather reads, decompressed: H2_partial's records)
template<typename T>
struct ChunkGatherItemT {
    T* out;
    int ld;
    int m;
    int n;
    const T* pool;
    IndexList slot;
    int64_t chunk;
    int64_t lda;
    IndexList rows;
    IndexList cols;
};

// out[i * ors + j * ocs] = a(i, j) + b(i, j)
template<typename T>
struct AddStoreItemT {
    T* out;
    int64_t ors;
    int64_t ocs;
    int m;
    int n;
    const T* a;
    int lda;
    const T* b;
    int ldb;
};
using AddStoreItem = AddStoreItemT<double>;

// x(i, j) += t(i, j) + t(j, i) on an n x n block
// a (n x n, leading dimension ld) = I
template<typename T>
struct IdentityItemT {
    T* a;
    int ld;
    int n;
};
using IdentityItem = IdentityItemT<double>;

template<typename T>
struct SymAddItemT {
    T* x;
    int ldx;
    const T* t;
    int ldt;
    int n;
};
using SymAddItem = SymAddItemT<double>;

// Undo LU row pivoting on the columns of b (rows x n): for i = n-1 .. 0 swap
// columns i and piv[i] - 1 (piv holds 1-based LAPACK pivots).
template<typename T>
struct ColumnSwapItemT {
    T* b;
    int ldb;
    int rows;
    const int* piv;
    int n;
};
using ColumnSwapItem = ColumnSwapItemT<double>;

// A = P L U, rebuilt from the LU that replaced it.  X_RR_full is the
// symmetrized X_RR as it stood before getrf factored it in place, so a box
// that keeps its solve factors need not store it: the product of the kept
// factors is the same matrix to rounding, which is ~1e-16 relative against
// the 1e-3 tolerances the fill sources are held to.
// The product itself is a batched GEMM (magma's, not a hand-rolled one: the
// first version of this was a dot product per thread with stride-ldlu loads
// and cost 2.1 s a run).  So L and U are split out first, multiplied, and the
// pivots undone on the rows of the result.
template<typename T>
struct LuSplitItemT {
    const T* lu;     // the factorization, r x r, leading dimension ldlu
    int ldlu;
    T* l;            // unit lower, r x r, leading dimension r
    T* u;            // upper, r x r, leading dimension r
    int r;
};
using LuSplitItem = LuSplitItemT<double>;

// For i = r-1 .. 0 swap rows i and piv[i] - 1 of a (r x cols), undoing the
// pivoting of getrf.
template<typename T>
struct RowSwapItemT {
    T* a;
    int lda;
    int cols;
    const int* piv;
    int r;
    int forward;  // 0: the pivots undone (P A, A = P L U rebuilt); 1: applied (P^T A, as getrs)
};
using RowSwapItem = RowSwapItemT<double>;

// ---- device construction of the sketch lists (device_sketch.cu)
//
// Sketched rows row_base .. row_base + count - 1 are the points
// slot_base + (list ? list[i] : i).
struct RowBlockDesc {
    int row_base;
    int count;
    int slot_base;
    const int* list;
};

// Sketched rows row_base + i (i < count) read source rows
// src_offset + (list ? list[i] : i) of a fill source.
struct RunDesc {
    int row_base;
    int count;
    int src_offset;
    const int* list;
};

constexpr int kMaxSourceRuns = 128;  // ring blocks listing one fill source

// (the packing of the sketch list entries, kEntryRowBits etc., and the
// owner warps: bpack_gpu_kernels.cuh)

// One box: draws of its rows (std::mt19937_64 from seed, sk per row; direct:
// row q of the target is row q of the sketch, id_sketch_shape), point slots
// of its rows, and the destination lists of its kernel rows.
struct SketchListsItem {
    uint64_t seed;
    int d;
    int sk;
    int direct;
    int rows;
    int num_blocks;
    int64_t blocks_offset;  // RowBlockDesc array in the metadata
    int* draws;             // rows * sk
    int* rows_out;          // rows
    int* ptr;               // kSketchOwners + 1
    int* entries;           // rows * sk
};

// One fill source of a box: destination lists of its stored rows.
struct SourceListsItem {
    int box;                // index of the box's SketchListsItem
    int num_runs;           // <= kMaxSourceRuns
    int64_t runs_offset;    // RunDesc array in the metadata
    int* ptr;               // kSketchOwners + 1
    int* entries;           // sum of the run counts * sk
    int* row_index;         // sum of the run counts: source row of each sketched row
};

using OrderedSketchItem = OrderedSketchItemT<double>;
// (the sketch of kernel rows: Evaluator::sketch_rows; of stored rows:
// launch_ordered_sketch_stored, evaluator_device.hpp)

template<typename T>
void launch_lu_split(const LuSplitItemT<T>* items, int count, int max_r, cudaStream_t stream);
template<typename T>
void launch_row_swaps(const RowSwapItemT<T>* items, int count, cudaStream_t stream);

void launch_sketch_lists(const SketchListsItem* boxes, int num_boxes, int max_blocks,
                         const SourceListsItem* sources, int num_sources, const char* meta, cudaStream_t stream);

// dst (n x m, leading dimension ld_dst) = src (m x n)^T
template<typename T>
struct TransposeItemT {
    const T* src;
    int ld_src;
    T* dst;
    int ld_dst;
    int m;
    int n;
};
using TransposeItem = TransposeItemT<double>;

// Interpolative decomposition of one sketch a (m x n, leading dimension lda),
// in place (device_id.cu).  On exit: rank K; jpvt, the n column indices with
// the skeleton first; T = R11^{-1} R12 in rows 0..K-1 of columns K..n-1 (for
// K = 0, row 0 of columns 1..n-1 is zero); norm, the Frobenius norm of the
// input; flag 1 for a non-finite input, 2 for a non-finite T.  tol >= 0: the
// item's own rank tolerance instead of the launch's.  The rank rule is
// compute_id_complex's: |R_kk| > max(tol |R_11|, abs_tol) on the input's scale
// (abs_tol, launch_qrcp's: the H2_s_op floor, 0 off).
template<typename T>
struct QrcpItemT {
    T* a;
    int m;
    int n;
    int lda;
    int* jpvt;
    int* rank;
    double* norm;
    int* flag;
    double tol = -1.0;
    // optional: the Householder scalars of the steps taken (at least min(m, n)
    // entries; the reflectors stay below the diagonal as LAPACK's geqp3 leaves
    // them), for a caller forming Q (the butterfly's range finder, bf_gpu)
    T* tau = nullptr;
};
using QrcpItem = QrcpItemT<double>;

// Device work (bytes, 16-byte aligned) that launch_qrcp needs for `count`
// items of at most max_n columns; 0 when it needs none.
template<typename T>
size_t qrcp_work_bytes(int count, int max_n);
// batch_independent: each box's result independent of its batch (one fixed
// configuration), for replicated CA levels.
template<typename T>
void launch_qrcp(const QrcpItemT<T>* items, int count, int max_n, double tol, double abs_tol, void* work,
                 cudaStream_t stream, bool batch_independent = false);
// The same ID by a blocked QRCP (LAPACK dlaqps: a GEMM update a panel), for
// tall items of many columns; max_m, the items' most rows.  Its results
// differ from launch_qrcp's in rounding; a box's are the same bit for bit
// whatever its batch.  The work block (bytes, 16-byte aligned) is always
// needed.
template<typename T>
size_t qrcp_blocked_work_bytes(int count, int max_m, int max_n);
template<typename T>
void launch_qrcp_blocked(const QrcpItemT<T>* items, int count, int max_m, int max_n, double tol, double abs_tol,
                         void* work, cudaStream_t stream);

// C_i = alpha op(A_i) op(B_i) + beta C_i on the FP64 tensor cores
// (device_gemm.cu); sizes, leading dimensions and pointers are device arrays
// as MAGMA's vbatched GEMM takes them.  One thread block per entry of blocks
// (device): {i, tile row << 16 | tile column} of a tile x tile output tile,
// tile 32, 48 or 64.  C_i is not read when beta = 0.
void launch_dgemm_vbatched_tc(bool trans_a, bool trans_b, const int* m, const int* n, const int* k, double alpha,
                              const double* const* a, const int* lda, const double* const* b, const int* ldb,
                              double beta, double* const* c, const int* ldc, const int2* blocks, int num_blocks,
                              int tile, cudaStream_t stream);
// The same on double complex data, C_i = op(A_i) B_i (+ C_i when accumulate),
// tile 32 or 64: each complex product as two real DMMA products of the
// operands' real views - the real part with A's imaginary entries negated,
// the imaginary part with B's (re, im) pairs swapped - so every element sums
// its k products in cuCfma's order, MAGMA's zgemm's, and the result is
// MAGMA's to the bit.  (Alpha 1 and beta 0 or 1 only: MAGMA's complex
// scaling is left to MAGMA.)
void launch_zgemm_vbatched_tc(bool trans_a, const int* m, const int* n, const int* k, const dcomplex* const* a,
                              const int* lda, const dcomplex* const* b, const int* ldb, bool accumulate,
                              dcomplex* const* c, const int* ldc, const int2* blocks, int num_blocks, int tile,
                              cudaStream_t stream);

template<typename T>
void launch_gather(const GatherItemT<T>* items, int count, int max_m, int max_n,
                   const char* meta, cudaStream_t stream);
template<typename T>
void launch_add_store(const AddStoreItemT<T>* items, int count, int max_m, int max_n,
                      cudaStream_t stream);
template<typename T>
void launch_index_copy(const IndexCopyItemT<T>* items, int count, int max_m, int max_n, const char* meta,
                       cudaStream_t stream);
template<typename T>
void launch_block_copy(const BlockCopyItemT<T>* items, int count, int max_m, int max_n, const char* meta,
                       cudaStream_t stream);
template<typename T>
void launch_chunk_gather(const ChunkGatherItemT<T>* items, int count, int max_m, int max_n, const char* meta,
                         cudaStream_t stream);
template<typename T>
void launch_sym_add(const SymAddItemT<T>* items, int count, int max_n, cudaStream_t stream);
// a := (a + a^T) / 2 of each item's n x n matrix (the symmetric part)
template<typename T>
void launch_symmetrize(const IdentityItemT<T>* items, int count, int max_n, cudaStream_t stream);
// ... of one n x n matrix (leading dimension ld)
template<typename T>
void launch_symmetrize(T* a, int ld, int n, cudaStream_t stream);

// target(i, j) += ((part_0(i, j) + part_1(i, j)) + ...) + part_{n-1}(i, j):
// parts summed in order first (as the host sums transported deltas), each
// rows x cols with leading dimension rows; their pointers are a
// meta-relative array.
template<typename T>
struct SumAddItemT {
    T* target;
    int ldt;
    int rows;
    int cols;
    int64_t parts_offset;
    int nparts;
};
using SumAddItem = SumAddItemT<double>;
template<typename T>
void launch_sum_add(const SumAddItemT<T>* items, int count, int max_rows, int max_cols, const char* meta,
                    cudaStream_t stream);
template<typename T>
void launch_identity(const IdentityItemT<T>* items, int count, int max_n, cudaStream_t stream);
template<typename T>
void launch_column_swaps(const ColumnSwapItemT<T>* items, int count, cudaStream_t stream);

// max |a[i]| over a contiguous block, one result per item.  Not templated:
// the element types are arrays of doubles in memory, so a complex block is
// viewed as 2n doubles, which is the component-wise bound the lossy
// compression's delta needs.
struct AbsMaxItem {
    const double* a;
    int64_t n;
    double* out;
};
void launch_abs_max(const AbsMaxItem* items, int count, cudaStream_t stream);
template<typename T>
void launch_transpose(const TransposeItemT<T>* items, int count, int max_m, int max_n, cudaStream_t stream);
// Copy `bytes` (a multiple of 16, both pointers 16-byte aligned) with SM
// stores.  Used to write small results straight into mapped pinned host
// memory, so they do not queue behind bulk copies on the copy engine.
void launch_copy_bytes(void* dst, const void* src, size_t bytes, cudaStream_t stream);

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
