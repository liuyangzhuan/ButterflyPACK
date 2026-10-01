#pragma once
// ButterflyPACK GPU backends (H2 and HODLR): the kernels that evaluate the
// matrix's entries with a user-defined entry evaluator, instantiated for the
// evaluator in the user's code (bpack_gpu_entry.cuh) or compiled at run time
// with NVRTC (an entry given as CUDA source text, e.g. from Python).
// doc/gpu_kernels.md describes the interface; this file is part of it only
// through bpack_gpu_entry.cuh.  It must stay self-contained for NVRTC: no
// system headers when __CUDACC_RTC__ is defined.
//
// An entry evaluator is a trivially copyable functor (passed by value to the
// kernels, at most 4 KB) with
//
//   using value_type = double;                    // or bpack::gpu::dcomplex
//   __device__ value_type operator()(const double* x, int64_t i,
//                                    const double* y, int64_t j) const;
//
// returning K(i, j) for the points with 0-based global ids i and j, whose
// coordinates are x[0 .. dim - 1] and y[0 .. dim - 1].

#ifndef __CUDACC_RTC__
#include <cuda_runtime.h>
#include <cmath>
#include <cstdint>
#else
typedef long long int64_t;
typedef unsigned long long uint64_t;
#endif

// The version of this interface: user code compiled against these headers
// registers with the version it saw (bpack_gpu_entry.cuh), and the library
// refuses another.
#define BPACK_GPU_INTERFACE_VERSION 1

namespace fmm {
namespace gpu {

// ---- double complex values on the device (the layout of std::complex<double>)

struct __align__(16) dcomplex {
    double re;
    double im;
    dcomplex() = default;
    __host__ __device__ constexpr dcomplex(double r, double i = 0.0) : re(r), im(i) {}
};

__host__ __device__ inline dcomplex operator+(dcomplex a, dcomplex b) { return {a.re + b.re, a.im + b.im}; }
__host__ __device__ inline dcomplex operator-(dcomplex a, dcomplex b) { return {a.re - b.re, a.im - b.im}; }
__host__ __device__ inline dcomplex operator-(dcomplex a) { return {-a.re, -a.im}; }
__host__ __device__ inline dcomplex operator*(dcomplex a, dcomplex b) {
    return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
}
__host__ __device__ inline dcomplex operator*(double a, dcomplex b) { return {a * b.re, a * b.im}; }
__host__ __device__ inline dcomplex operator*(dcomplex a, double b) { return {a.re * b, a.im * b}; }
__host__ __device__ inline dcomplex operator/(dcomplex a, double b) { return {a.re / b, a.im / b}; }
__host__ __device__ inline dcomplex operator/(dcomplex a, dcomplex b) {
    if (fabs(b.re) >= fabs(b.im)) {
        const double r = b.im / b.re, d = b.re + b.im * r;
        return {(a.re + a.im * r) / d, (a.im - a.re * r) / d};
    }
    const double r = b.re / b.im, d = b.re * r + b.im;
    return {(a.re * r + a.im) / d, (a.im * r - a.re) / d};
}
__host__ __device__ inline dcomplex& operator+=(dcomplex& a, dcomplex b) { a.re += b.re; a.im += b.im; return a; }
__host__ __device__ inline dcomplex& operator-=(dcomplex& a, dcomplex b) { a.re -= b.re; a.im -= b.im; return a; }
__host__ __device__ inline dcomplex& operator*=(dcomplex& a, dcomplex b) { a = a * b; return a; }
__host__ __device__ inline dcomplex& operator/=(dcomplex& a, double b) { a.re /= b; a.im /= b; return a; }
__host__ __device__ inline bool operator==(dcomplex a, dcomplex b) { return a.re == b.re && a.im == b.im; }
__host__ __device__ inline bool operator!=(dcomplex a, dcomplex b) { return !(a == b); }
__host__ __device__ inline dcomplex conj(dcomplex a) { return {a.re, -a.im}; }
__host__ __device__ inline double conj(double a) { return a; }

// ---- what the kernels read

// Coordinates (dim per point) and 0-based global ids of the points of a
// level, by slot.
struct PointTable {
    const double* xyz = nullptr;
    const int64_t* ids = nullptr;
    int dim = 3;
};

// Index i of a list: base plus i, the int32 entry i of the list at
// meta + offset (offset >= 0), or entry i of the device array ptr.
struct IndexList {
    int64_t offset = -1;
    int base = 0;
    const int* ptr = nullptr;
};

// out(i, j) = K(point rows[i], point cols[j]), column-major
template<typename T>
struct EvalItemT {
    T* out;
    int ld;
    int m;
    int n;
    IndexList rows;
    IndexList cols;
};

// Sketch list entries pack (row, destination, sign): row in the low
// kEntryRowBits bits, destination above, sign in bit 31 (set = negative).
// The draws of a box are partitioned by owner warp: destination mod
// kSketchOwners.
constexpr int kEntryRowBits = 20;
constexpr int kEntryRowMask = (1 << kEntryRowBits) - 1;
constexpr int kEntryDestMask = (1 << (31 - kEntryRowBits)) - 1;
constexpr int kSketchOwners = 16;

// out (d x ncols, leading dimension ldo) = sparse-sign sketch of `rows` rows:
// out(i, c) = sum over entries (row, i, sign), in list order, of
// sign * scale * value(row, c), where value is K(point col_base + c, point
// row_slots[row]) for kernel rows, else src[c + row_index[row] *
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
    const int* row_index;  // null: row i is stored row i
};

// ---- the entry kernel: thread blocks of kEvalTile x kEvalTile entries

constexpr int kEvalTile = 32;
constexpr int kEvalTileRows = 8;  // thread rows; each thread covers kEvalTile / kEvalTileRows columns

#if defined(__CUDACC__) || defined(__CUDACC_RTC__)  // (the kernels: compiled by nvcc or NVRTC only)

__device__ __forceinline__ int index_at(const char* meta, const IndexList& list, int i) {
    if (list.ptr != nullptr) return list.base + list.ptr[i];
    return list.base + (list.offset < 0 ? i : reinterpret_cast<const int*>(meta + list.offset)[i]);
}

template<typename T, typename Entry>
__global__ void eval_kernel(const EvalItemT<T>* items, const char* meta, Entry entry, PointTable points) {
    const EvalItemT<T> item = items[blockIdx.x];
    const int i = blockIdx.y * kEvalTile + threadIdx.x;
    const int j0 = blockIdx.z * kEvalTile;
    if (i >= item.m || j0 >= item.n) return;
    const int si = index_at(meta, item.rows, i);
    const double* x = points.xyz + static_cast<int64_t>(points.dim) * si;
    const int64_t x_id = points.ids[si];
    for (int jj = threadIdx.y; jj < kEvalTile; jj += kEvalTileRows) {
        const int j = j0 + jj;
        if (j >= item.n) break;
        const int sj = index_at(meta, item.cols, j);
        item.out[i + static_cast<int64_t>(j) * item.ld] =
            entry(x, x_id, points.xyz + static_cast<int64_t>(points.dim) * sj, points.ids[sj]);
    }
}

#endif  // __CUDACC__

// ---- the ordered sketch of kernel rows (or stored rows)
//
// A block owns kSketchTileCols columns and the destinations [i0, i0 +
// dest_block); the output tile lives in shared memory and kSketchTileRows
// rows of values at a time are formed once (kernel entries, or source rows
// read once) and applied by the warps, warp w owning the destinations i = w
// mod kSketchOwners.  Every destination sums in the draw order, so the
// result matches the per-destination lists exactly.

constexpr int kSketchTileCols = 32;
constexpr int kSketchTileRows = 64;
constexpr int kSketchPitch = kSketchTileCols + 1;  // padded rows of the shared tiles
constexpr int kSketchRowsPerWarp = kSketchTileRows / kSketchOwners;
constexpr int kSketchCachedDim = 4;  // a column point of at most this many coordinates is held in registers

#if defined(__CUDACC__) || defined(__CUDACC_RTC__)

// A stand-in entry for sketches of stored rows only.
template<typename T>
struct NoEntry {
    using value_type = T;
    __device__ T operator()(const double*, int64_t, const double*, int64_t) const { return T(0.0); }
};

// Values of this thread's rows of the chunk starting at r0 (rows
// r0 + warp + q * kSketchOwners, its column c): the row indices first, then
// the values, so the loads of a thread are independent.
template<typename T, bool kKernelRows, typename Entry>
__device__ __forceinline__ void chunk_values(const OrderedSketchItemT<T>& item, const Entry& entry,
                                             const PointTable& points, const double* col_x, int64_t col_id,
                                             bool column, int c, int warp, int r0, T (&v)[kSketchRowsPerWarp]) {
    int src_row[kSketchRowsPerWarp];
#pragma unroll
    for (int q = 0; q < kSketchRowsPerWarp; ++q) {
        const int row = r0 + warp + q * kSketchOwners;
        src_row[q] = row < item.rows ? (kKernelRows ? item.row_slots[row] : (item.row_index ? item.row_index[row] : row))
                                     : -1;
    }
#pragma unroll
    for (int q = 0; q < kSketchRowsPerWarp; ++q) {
        v[q] = T(0.0);
        if (src_row[q] >= 0 && column) {
            if (kKernelRows) {
                v[q] = entry(col_x, col_id, points.xyz + static_cast<int64_t>(points.dim) * src_row[q],
                             points.ids[src_row[q]]);
            } else {
                v[q] = item.src[c + static_cast<int64_t>(src_row[q]) * item.row_stride];
            }
        }
    }
}

// y + s v for a real s: y + s v componentwise for complex values (the host's
// axpy with a real coefficient)
__device__ __forceinline__ double sketch_fma(double s, double v, double y) { return fma(s, v, y); }
__device__ __forceinline__ dcomplex sketch_fma(double s, dcomplex v, dcomplex y) {
    return dcomplex(fma(s, v.re, y.re), fma(s, v.im, y.im));
}

template<typename T, bool kKernelRows, typename Entry>
__global__ void __launch_bounds__(32 * kSketchOwners)
ordered_sketch_kernel(const OrderedSketchItemT<T>* items, Entry entry, PointTable points, int dest_block) {
    extern __shared__ __align__(16) unsigned char shared_bytes[];
    T* shared = reinterpret_cast<T*>(shared_bytes);
    const OrderedSketchItemT<T> item = items[blockIdx.x];
    const int c0 = blockIdx.y * kSketchTileCols;
    const int i0 = blockIdx.z * dest_block;
    if (c0 >= item.ncols || i0 >= item.d) return;
    const int i1 = min(item.d, i0 + dest_block);
    T* tile = shared;                                                            // (i1 - i0) x kSketchPitch
    T* values = shared + static_cast<int64_t>(dest_block) * kSketchPitch;        // kSketchTileRows x kSketchPitch
    const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5;
    for (int e = threadIdx.x; e < (i1 - i0) * kSketchPitch; e += blockDim.x) tile[e] = T(0.0);
    const int c = c0 + lane;
    const bool column = c < item.ncols;
    double col_cached[kSketchCachedDim] = {0.0, 0.0, 0.0, 0.0};
    const double* col_x = col_cached;
    int64_t col_id = 0;
    if (kKernelRows && column) {
        const int sc = item.col_base + c;
        const double* x = points.xyz + static_cast<int64_t>(points.dim) * sc;
        if (points.dim <= kSketchCachedDim) {
            for (int k = 0; k < points.dim; ++k) col_cached[k] = x[k];
        } else {
            col_x = x;
        }
        col_id = points.ids[sc];
    }
    // this warp's entries, 32 at a time: one coalesced load, then broadcasts
    const int* list = item.entries + item.ptr[warp];
    const int* end = item.entries + item.ptr[warp + 1];
    // +-scale by the entry's sign bit (the same fma as with a branch)
    const unsigned long long scale_bits = static_cast<unsigned long long>(__double_as_longlong(item.scale));
    auto signed_scale = [&](int e) {
        return __longlong_as_double(static_cast<long long>(
            scale_bits ^ (static_cast<unsigned long long>(static_cast<unsigned>(e) & 0x80000000u) << 32)));
    };
    const bool all_dests = i0 == 0 && i1 == item.d;  // no range test
    int batch = 0, batch_count = 0, batch_pos = 0;
    T next[kSketchRowsPerWarp];
    chunk_values<T, kKernelRows, Entry>(item, entry, points, col_x, col_id, column, c, warp, 0, next);
    for (int r0 = 0; r0 < item.rows; r0 += kSketchTileRows) {
        __syncthreads();  // the value tile is free
#pragma unroll
        for (int q = 0; q < kSketchRowsPerWarp; ++q) values[(warp + q * kSketchOwners) * kSketchPitch + lane] = next[q];
        __syncthreads();
        // the next chunk's loads are in flight while this one is applied
        if (r0 + kSketchTileRows < item.rows) {
            chunk_values<T, kKernelRows, Entry>(item, entry, points, col_x, col_id, column, c, warp,
                                                r0 + kSketchTileRows, next);
        }
        const int r1 = r0 + kSketchTileRows;
        const T* vrow = values + lane;
        T* ycol = tile + lane;
        while (true) {
            if (batch_pos == batch_count) {
                if (list >= end) break;
                batch_count = end - list < 32 ? static_cast<int>(end - list) : 32;
                batch = lane < batch_count ? list[lane] : 0;
                list += batch_count;
                batch_pos = 0;
            }
            // entries are in row order: those of this chunk are a run from batch_pos
            const unsigned here = __ballot_sync(0xffffffffu, lane >= batch_pos && lane < batch_count &&
                                                                 (batch & kEntryRowMask) < r1);
            const int stop = batch_pos + __popc(here);
            int j = batch_pos;
            if (all_dests) {
                // two entries at a time: different destinations update independently
                for (; j + 1 < stop; j += 2) {
                    const int e1 = __shfl_sync(0xffffffffu, batch, j);
                    const int e2 = __shfl_sync(0xffffffffu, batch, j + 1);
                    const int a1 = ((e1 >> kEntryRowBits) & kEntryDestMask) * kSketchPitch;
                    const int a2 = ((e2 >> kEntryRowBits) & kEntryDestMask) * kSketchPitch;
                    const T v1 = vrow[((e1 & kEntryRowMask) - r0) * kSketchPitch];
                    const T v2 = vrow[((e2 & kEntryRowMask) - r0) * kSketchPitch];
                    if (a1 != a2) {
                        const T y1 = ycol[a1], y2 = ycol[a2];
                        ycol[a1] = sketch_fma(signed_scale(e1), v1, y1);
                        ycol[a2] = sketch_fma(signed_scale(e2), v2, y2);
                    } else {
                        ycol[a1] = sketch_fma(signed_scale(e2), v2, sketch_fma(signed_scale(e1), v1, ycol[a1]));
                    }
                }
                if (j < stop) {
                    const int e = __shfl_sync(0xffffffffu, batch, j);
                    const int a = ((e >> kEntryRowBits) & kEntryDestMask) * kSketchPitch;
                    ycol[a] = sketch_fma(signed_scale(e), vrow[((e & kEntryRowMask) - r0) * kSketchPitch], ycol[a]);
                    ++j;
                }
            } else {
                for (; j < stop; ++j) {
                    const int e = __shfl_sync(0xffffffffu, batch, j);
                    const int i = (e >> kEntryRowBits) & kEntryDestMask;
                    if (i >= i0 && i < i1) {
                        T& y = ycol[(i - i0) * kSketchPitch];
                        y = sketch_fma(signed_scale(e), vrow[((e & kEntryRowMask) - r0) * kSketchPitch], y);
                    }
                }
            }
            batch_pos = stop;
            if (batch_pos < batch_count) break;  // the rest belongs to later chunks
        }
    }
    __syncthreads();
    for (int col = warp; col < kSketchTileCols; col += kSketchOwners) {
        if (c0 + col >= item.ncols) break;
        for (int i = lane; i < i1 - i0; i += 32) {
            item.out[(i0 + i) + static_cast<int64_t>(c0 + col) * item.ldo] = tile[i * kSketchPitch + col];
        }
    }
}

#endif  // __CUDACC__

#ifndef __CUDACC_RTC__

// ---- launch geometry (shared by the library's launches and the user's)

inline dim3 eval_grid(int count, int max_m, int max_n) {
    return dim3(static_cast<unsigned>(count), static_cast<unsigned>((max_m + kEvalTile - 1) / kEvalTile),
                static_cast<unsigned>((max_n + kEvalTile - 1) / kEvalTile));
}
inline dim3 eval_block() { return dim3(kEvalTile, kEvalTileRows); }

// All destinations in one block when the tile fits the shared memory of a
// block: a split recomputes every value and rescans every list per part.
struct SketchGeometry {
    dim3 grid;
    dim3 block;
    size_t shared;
    int dest_block;
};
inline SketchGeometry sketch_geometry(int count, int max_d, int max_cols, size_t scalar_bytes) {
    int device = 0, optin = 0;
    cudaGetDevice(&device);
    cudaDeviceGetAttribute(&optin, cudaDevAttrMaxSharedMemoryPerBlockOptin, device);
    const int max_block = optin / static_cast<int>(kSketchPitch * scalar_bytes) - kSketchTileRows;
    const int parts = (max_d + max_block - 1) / max_block;
    SketchGeometry g;
    g.dest_block = (max_d + parts - 1) / parts;
    g.shared = static_cast<size_t>(g.dest_block + kSketchTileRows) * kSketchPitch * scalar_bytes;
    g.grid = dim3(static_cast<unsigned>(count), static_cast<unsigned>((max_cols + kSketchTileCols - 1) / kSketchTileCols),
                  static_cast<unsigned>((max_d + g.dest_block - 1) / g.dest_block));
    g.block = dim3(32 * kSketchOwners);
    return g;
}

#endif  // __CUDACC_RTC__

}  // namespace gpu
}  // namespace fmm
