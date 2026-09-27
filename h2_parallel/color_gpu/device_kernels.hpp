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

#include "device_scalar.hpp"

#include <cuda_runtime.h>

#include <cstdint>

namespace fmm {
namespace gpu {

// Device-evaluable kernel K(x, y) registered by the application.
//   kind 1 (real):    K = p[1] if the global ids match, else p[0] / |x - y|  (3D)
//   kind 2 (complex): K = (p[3], p[4]) if the global ids match, else
//                     (p[1], p[2]) e^{i p[0] r} / (p[5] r), r = |x - y|  (3D;
//                     symmetric Helmholtz, p[5] = 4 pi as the host forms it)
constexpr int kKernelParams = 8;
struct KernelSpec {
    int kind = 0;
    double p[kKernelParams] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
};

// Coordinates (xyz per slot) and global ids of the points of one level.
struct PointTable {
    const double* xyz = nullptr;
    const int64_t* ids = nullptr;
};

// Index i of a list: base plus i, the int32 entry i of the list at
// meta + offset (offset >= 0), or entry i of the device array ptr.
struct IndexList {
    int64_t offset = -1;
    int base = 0;
    const int* ptr = nullptr;
};

// out(i, j) = K(point rows[i], point cols[j])
template<typename T>
struct EvalItemT {
    T* out;
    int ld;
    int m;
    int n;
    IndexList rows;
    IndexList cols;
};
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

// Sketch list entries pack (row, destination, sign): row in the low
// kEntryRowBits bits, destination above, sign in bit 31 (set = negative).
// The draws of a box are partitioned by owner warp: destination mod
// kSketchOwners.
constexpr int kEntryRowBits = 20;
constexpr int kEntryRowMask = (1 << kEntryRowBits) - 1;
constexpr int kEntryDestMask = (1 << (31 - kEntryRowBits)) - 1;
constexpr int kSketchOwners = 16;

// One box: draws of its rows (std::mt19937_64 from seed, sk per row), point
// slots of its rows, and the destination lists of its kernel rows.
struct SketchListsItem {
    uint64_t seed;
    int d;
    int sk;
    int rows;
    int num_blocks;
    int64_t blocks_offset;  // RowBlockDesc array in the metadata
    int* draws;             // rows * sk
    int* rows_out;          // rows
    int* ptr;               // kSketchOwners + 1
    int* entries;           // rows * sk
    int* hist;              // 32 * kSketchOwners scratch
};

// One fill source of a box: destination lists of its stored rows.
struct SourceListsItem {
    int box;                // index of the box's SketchListsItem
    int num_runs;           // <= kMaxSourceRuns
    int64_t runs_offset;    // RunDesc array in the metadata
    int* ptr;               // kSketchOwners + 1
    int* entries;           // sum of the run counts * sk
    int* hist;              // 32 * kSketchOwners scratch
    int* row_index;         // sum of the run counts: source row of each sketched row
};

// out (d x ncols, leading dimension ldo) = sparse-sign sketch of `rows` rows:
// out(i, c) = sum over entries (row, i, sign), in list order, of
// sign * scale * value(row, c), where value is K(point row_slots[row],
// point col_base + c) for kernel rows, else src[c + row_index[row] *
// row_stride].  ptr / entries: the lists of the owner warps.
template<typename T>
struct OrderedSketchItemT {
    T* out;
    int ldo;
    int d;
    int ncols;
    int rows;
    const int* ptr;
    const int* entries;
    double scale;
    const int* row_slots;
    int col_base;
    const T* src;
    int row_stride;
    const int* row_index;
};
using OrderedSketchItem = OrderedSketchItemT<double>;

template<typename T>
void launch_ordered_sketch(const OrderedSketchItemT<T>* items, int count, int max_d, int max_cols, bool kernel_rows,
                           KernelSpec spec, PointTable points, cudaStream_t stream);

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
// input; flag 1 for a non-finite input, 2 for a non-finite T.
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
};
using QrcpItem = QrcpItemT<double>;

// Device work (bytes, 16-byte aligned) that launch_qrcp needs for `count`
// items of at most max_n columns; 0 when it needs none.
template<typename T>
size_t qrcp_work_bytes(int count, int max_n);
template<typename T>
void launch_qrcp(const QrcpItemT<T>* items, int count, int max_n, double tol, void* work, cudaStream_t stream);

// C_i = alpha op(A_i) op(B_i) + beta C_i on the FP64 tensor cores
// (device_gemm.cu); sizes, leading dimensions and pointers are device arrays
// as MAGMA's vbatched GEMM takes them.  One thread block per entry of blocks
// (device): {i, tile row << 16 | tile column} of a tile x tile output tile,
// tile 32, 48 or 64.  C_i is not read when beta = 0.
void launch_dgemm_vbatched_tc(bool trans_a, bool trans_b, const int* m, const int* n, const int* k, double alpha,
                              const double* const* a, const int* lda, const double* const* b, const int* ldb,
                              double beta, double* const* c, const int* ldc, const int2* blocks, int num_blocks,
                              int tile, cudaStream_t stream);

template<typename T>
void launch_eval(const EvalItemT<T>* items, int count, int max_m, int max_n,
                 const char* meta, KernelSpec spec, PointTable points, cudaStream_t stream);
template<typename T>
void launch_gather(const GatherItemT<T>* items, int count, int max_m, int max_n,
                   const char* meta, cudaStream_t stream);
template<typename T>
void launch_add_store(const AddStoreItemT<T>* items, int count, int max_m, int max_n,
                      cudaStream_t stream);
template<typename T>
void launch_sym_add(const SymAddItemT<T>* items, int count, int max_n, cudaStream_t stream);

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
template<typename T>
void launch_transpose(const TransposeItemT<T>* items, int count, int max_m, int max_n, cudaStream_t stream);
// Copy `bytes` (a multiple of 16, both pointers 16-byte aligned) with SM
// stores.  Used to write small results straight into mapped pinned host
// memory, so they do not queue behind bulk copies on the copy engine.
void launch_copy_bytes(void* dst, const void* src, size_t bytes, cudaStream_t stream);

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
