// Batched element kernels of the H2 Color GPU backend.  See device_kernels.hpp.

#include "device_kernels.hpp"
#include "kernel_eval.cuh"

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace fmm {
namespace gpu {

namespace {

constexpr int kTile = 32;      // rows and columns per thread block
constexpr int kTileRows = 8;   // thread rows; each thread covers kTile / kTileRows columns

void check_launch(const char* what) {
    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA launch failed in ") + what + ": " +
                                 cudaGetErrorString(status));
    }
}

__device__ __forceinline__ int index_at(const char* meta, const IndexList& list, int i) {
    if (list.ptr != nullptr) return list.base + list.ptr[i];
    return list.base + (list.offset < 0 ? i : reinterpret_cast<const int*>(meta + list.offset)[i]);
}

template<typename T, int Kind>
__global__ void eval_kernel(const EvalItemT<T>* items, const char* meta, KernelSpec spec, PointTable points) {
    const EvalItemT<T> item = items[blockIdx.x];
    const int i = blockIdx.y * kTile + threadIdx.x;
    const int j0 = blockIdx.z * kTile;
    if (i >= item.m || j0 >= item.n) return;
    const int si = index_at(meta, item.rows, i);
    const double* x = points.xyz + 3 * static_cast<int64_t>(si);
    const int64_t x_id = points.ids[si];
    for (int jj = threadIdx.y; jj < kTile; jj += kTileRows) {
        const int j = j0 + jj;
        if (j >= item.n) break;
        const int sj = index_at(meta, item.cols, j);
        item.out[i + static_cast<int64_t>(j) * item.ld] =
            kernel_value<T, Kind>(spec, x, x_id, points.xyz + 3 * static_cast<int64_t>(sj), points.ids[sj]);
    }
}

template<typename T>
__global__ void gather_kernel(const GatherItemT<T>* items, const char* meta) {
    const GatherItemT<T> item = items[blockIdx.x];
    const int i = blockIdx.y * kTile + threadIdx.x;
    const int j0 = blockIdx.z * kTile;
    if (i >= item.m || j0 >= item.n) return;
    const int64_t row = index_at(meta, item.rows, i);
    for (int jj = threadIdx.y; jj < kTile; jj += kTileRows) {
        const int j = j0 + jj;
        if (j >= item.n) break;
        const int64_t col = index_at(meta, item.cols, j);
        item.out[i + static_cast<int64_t>(j) * item.ld] = item.src[row * item.rs + col * item.cs];
    }
}

template<typename T>
__global__ void add_store_kernel(const AddStoreItemT<T>* items) {
    const AddStoreItemT<T> item = items[blockIdx.x];
    const int i = blockIdx.y * kTile + threadIdx.x;
    const int j0 = blockIdx.z * kTile;
    if (i >= item.m || j0 >= item.n) return;
    for (int jj = threadIdx.y; jj < kTile; jj += kTileRows) {
        const int j = j0 + jj;
        if (j >= item.n) break;
        item.out[i * item.ors + j * item.ocs] =
            item.a[i + static_cast<int64_t>(j) * item.lda] +
            item.b[i + static_cast<int64_t>(j) * item.ldb];
    }
}

template<typename T>
__global__ void sym_add_kernel(const SymAddItemT<T>* items) {
    const SymAddItemT<T> item = items[blockIdx.x];
    const int i = blockIdx.y * kTile + threadIdx.x;
    const int j0 = blockIdx.z * kTile;
    if (i >= item.n || j0 >= item.n) return;
    for (int jj = threadIdx.y; jj < kTile; jj += kTileRows) {
        const int j = j0 + jj;
        if (j >= item.n) break;
        T& x = item.x[i + static_cast<int64_t>(j) * item.ldx];
        x += item.t[i + static_cast<int64_t>(j) * item.ldt] + item.t[j + static_cast<int64_t>(i) * item.ldt];
    }
}

template<typename T>
__device__ void symmetrize_tile(const IdentityItemT<T>& item) {
    const int i = blockIdx.y * kTile + threadIdx.x;
    const int j0 = blockIdx.z * kTile;
    if (i >= item.n || j0 >= item.n) return;
    for (int jj = threadIdx.y; jj < kTile; jj += kTileRows) {
        const int j = j0 + jj;
        if (j >= item.n) break;
        if (i >= j) continue;  // each pair once, from its upper entry
        T& upper = item.a[i + static_cast<int64_t>(j) * item.ld];
        T& lower = item.a[j + static_cast<int64_t>(i) * item.ld];
        const T mean = (upper + lower) * 0.5;
        upper = mean;
        lower = mean;
    }
}

template<typename T>
__global__ void symmetrize_kernel(const IdentityItemT<T>* items) {
    symmetrize_tile(items[blockIdx.x]);
}

template<typename T>
__global__ void symmetrize_one_kernel(IdentityItemT<T> item) {
    symmetrize_tile(item);
}

template<typename T>
__global__ void sum_add_kernel(const SumAddItemT<T>* items, const char* meta) {
    const SumAddItemT<T> item = items[blockIdx.x];
    const int i = blockIdx.y * kTile + threadIdx.x;
    const int j0 = blockIdx.z * kTile;
    if (i >= item.rows || j0 >= item.cols) return;
    const T* const* parts = reinterpret_cast<const T* const*>(meta + item.parts_offset);
    for (int jj = threadIdx.y; jj < kTile; jj += kTileRows) {
        const int j = j0 + jj;
        if (j >= item.cols) break;
        const int64_t at = i + static_cast<int64_t>(j) * item.rows;
        T sum = parts[0][at];
        for (int p = 1; p < item.nparts; ++p) sum += parts[p][at];
        item.target[i + static_cast<int64_t>(j) * item.ldt] += sum;
    }
}

template<typename T>
__global__ void identity_kernel(const IdentityItemT<T>* items) {
    const IdentityItemT<T> item = items[blockIdx.x];
    const int i = blockIdx.y * kTile + threadIdx.x;
    const int j0 = blockIdx.z * kTile;
    if (i >= item.n || j0 >= item.n) return;
    for (int jj = threadIdx.y; jj < kTile; jj += kTileRows) {
        const int j = j0 + jj;
        if (j >= item.n) break;
        item.a[i + static_cast<int64_t>(j) * item.ld] = T(i == j ? 1.0 : 0.0);
    }
}

template<typename T>
__global__ void column_swap_kernel(const ColumnSwapItemT<T>* items) {
    const ColumnSwapItemT<T> item = items[blockIdx.x];
    // Each thread owns whole rows, so the sequential swaps need no barrier.
    for (int row = threadIdx.x; row < item.rows; row += blockDim.x) {
        for (int i = item.n - 1; i >= 0; --i) {
            const int p = item.piv[i] - 1;
            if (p == i) continue;
            T* ci = item.b + static_cast<int64_t>(i) * item.ldb + row;
            T* cp = item.b + static_cast<int64_t>(p) * item.ldb + row;
            const T t = *ci;
            *ci = *cp;
            *cp = t;
        }
    }
}

template<typename T>
__global__ void transpose_kernel(const TransposeItemT<T>* items) {
    __shared__ T tile[kTile][kTile + 1];
    const TransposeItemT<T> item = items[blockIdx.x];
    const int i0 = blockIdx.y * kTile;
    const int j0 = blockIdx.z * kTile;
    if (i0 >= item.m || j0 >= item.n) return;  // uniform over the block
    for (int jj = threadIdx.y; jj < kTile; jj += kTileRows) {
        const int i = i0 + threadIdx.x, j = j0 + jj;
        if (i < item.m && j < item.n) tile[jj][threadIdx.x] = item.src[i + static_cast<int64_t>(j) * item.ld_src];
    }
    __syncthreads();
    for (int ii = threadIdx.y; ii < kTile; ii += kTileRows) {
        const int j = j0 + threadIdx.x, i = i0 + ii;
        if (i < item.m && j < item.n) item.dst[j + static_cast<int64_t>(i) * item.ld_dst] = tile[threadIdx.x][ii];
    }
}

dim3 tile_grid(int count, int max_m, int max_n) {
    return dim3(static_cast<unsigned>(count),
                static_cast<unsigned>((max_m + kTile - 1) / kTile),
                static_cast<unsigned>((max_n + kTile - 1) / kTile));
}

}  // namespace

template<typename T>
void launch_eval(const EvalItemT<T>* items, int count, int max_m, int max_n,
                 const char* meta, KernelSpec spec, PointTable points, cudaStream_t stream) {
    if (count <= 0 || max_m <= 0 || max_n <= 0) return;
    const dim3 grid = tile_grid(count, max_m, max_n), block(kTile, kTileRows);
    if constexpr (is_complex_scalar<T>) {
        if (kernel_kind_of<T>(spec, "launch_eval") == 3) {
            eval_kernel<T, 3><<<grid, block, 0, stream>>>(items, meta, spec, points);
        } else {
            eval_kernel<T, 2><<<grid, block, 0, stream>>>(items, meta, spec, points);
        }
    } else {
        kernel_kind_of<T>(spec, "launch_eval");
        eval_kernel<T, 1><<<grid, block, 0, stream>>>(items, meta, spec, points);
    }
    check_launch("eval_kernel");
}

template<typename T>
void launch_gather(const GatherItemT<T>* items, int count, int max_m, int max_n,
                   const char* meta, cudaStream_t stream) {
    if (count <= 0 || max_m <= 0 || max_n <= 0) return;
    gather_kernel<T><<<tile_grid(count, max_m, max_n), dim3(kTile, kTileRows), 0, stream>>>(items, meta);
    check_launch("gather_kernel");
}

template<typename T>
void launch_add_store(const AddStoreItemT<T>* items, int count, int max_m, int max_n,
                      cudaStream_t stream) {
    if (count <= 0 || max_m <= 0 || max_n <= 0) return;
    add_store_kernel<T><<<tile_grid(count, max_m, max_n), dim3(kTile, kTileRows), 0, stream>>>(items);
    check_launch("add_store_kernel");
}

template<typename T>
void launch_sym_add(const SymAddItemT<T>* items, int count, int max_n, cudaStream_t stream) {
    if (count <= 0 || max_n <= 0) return;
    sym_add_kernel<T><<<tile_grid(count, max_n, max_n), dim3(kTile, kTileRows), 0, stream>>>(items);
    check_launch("sym_add_kernel");
}

template<typename T>
void launch_symmetrize(const IdentityItemT<T>* items, int count, int max_n, cudaStream_t stream) {
    if (count <= 0 || max_n <= 0) return;
    symmetrize_kernel<T><<<tile_grid(count, max_n, max_n), dim3(kTile, kTileRows), 0, stream>>>(items);
    check_launch("symmetrize_kernel");
}

template<typename T>
void launch_symmetrize(T* a, int ld, int n, cudaStream_t stream) {
    if (n <= 0) return;
    symmetrize_one_kernel<T><<<tile_grid(1, n, n), dim3(kTile, kTileRows), 0, stream>>>(IdentityItemT<T>{a, ld, n});
    check_launch("symmetrize_one_kernel");
}

template<typename T>
void launch_sum_add(const SumAddItemT<T>* items, int count, int max_rows, int max_cols, const char* meta,
                    cudaStream_t stream) {
    if (count <= 0 || max_rows <= 0 || max_cols <= 0) return;
    sum_add_kernel<T><<<tile_grid(count, max_rows, max_cols), dim3(kTile, kTileRows), 0, stream>>>(items, meta);
    check_launch("sum_add_kernel");
}

template<typename T>
void launch_identity(const IdentityItemT<T>* items, int count, int max_n, cudaStream_t stream) {
    if (count <= 0 || max_n <= 0) return;
    identity_kernel<T><<<tile_grid(count, max_n, max_n), dim3(kTile, kTileRows), 0, stream>>>(items);
    check_launch("identity_kernel");
}

template<typename T>
void launch_transpose(const TransposeItemT<T>* items, int count, int max_m, int max_n, cudaStream_t stream) {
    if (count <= 0 || max_m <= 0 || max_n <= 0) return;
    transpose_kernel<T><<<tile_grid(count, max_m, max_n), dim3(kTile, kTileRows), 0, stream>>>(items);
    check_launch("transpose_kernel");
}

template<typename T>
void launch_column_swaps(const ColumnSwapItemT<T>* items, int count, cudaStream_t stream) {
    if (count <= 0) return;
    column_swap_kernel<T><<<static_cast<unsigned>(count), 128, 0, stream>>>(items);
    check_launch("column_swap_kernel");
}

#define H2_ELEMENT_LAUNCHES(T)                                                                                     \
    template void launch_eval<T>(const EvalItemT<T>*, int, int, int, const char*, KernelSpec, PointTable,          \
                                 cudaStream_t);                                                                    \
    template void launch_gather<T>(const GatherItemT<T>*, int, int, int, const char*, cudaStream_t);               \
    template void launch_add_store<T>(const AddStoreItemT<T>*, int, int, int, cudaStream_t);                       \
    template void launch_sym_add<T>(const SymAddItemT<T>*, int, int, cudaStream_t);                               \
    template void launch_symmetrize<T>(const IdentityItemT<T>*, int, int, cudaStream_t);                          \
    template void launch_symmetrize<T>(T*, int, int, cudaStream_t);                                               \
    template void launch_sum_add<T>(const SumAddItemT<T>*, int, int, int, const char*, cudaStream_t);              \
    template void launch_identity<T>(const IdentityItemT<T>*, int, int, cudaStream_t);                             \
    template void launch_transpose<T>(const TransposeItemT<T>*, int, int, int, cudaStream_t);                      \
    template void launch_column_swaps<T>(const ColumnSwapItemT<T>*, int, cudaStream_t);
H2_ELEMENT_LAUNCHES(double)
H2_ELEMENT_LAUNCHES(dcomplex)
#undef H2_ELEMENT_LAUNCHES

__global__ void copy_bytes_kernel(int4* __restrict__ dst, const int4* __restrict__ src, size_t count) {
    for (size_t i = blockIdx.x * static_cast<size_t>(blockDim.x) + threadIdx.x; i < count;
         i += static_cast<size_t>(gridDim.x) * blockDim.x) {
        dst[i] = src[i];
    }
}

void launch_copy_bytes(void* dst, const void* src, size_t bytes, cudaStream_t stream) {
    if (bytes == 0) return;
    if (bytes % sizeof(int4) != 0 || reinterpret_cast<uintptr_t>(dst) % sizeof(int4) != 0 ||
        reinterpret_cast<uintptr_t>(src) % sizeof(int4) != 0) {
        throw std::invalid_argument("launch_copy_bytes: size and pointers must be 16-byte aligned");
    }
    const size_t count = bytes / sizeof(int4);
    const unsigned blocks = static_cast<unsigned>(std::min<size_t>((count + 255) / 256, 1024));
    copy_bytes_kernel<<<blocks, 256, 0, stream>>>(static_cast<int4*>(dst), static_cast<const int4*>(src), count);
    check_launch("copy_bytes_kernel");
}

}  // namespace gpu
}  // namespace fmm
