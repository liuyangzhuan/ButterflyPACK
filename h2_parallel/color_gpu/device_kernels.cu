// Batched element kernels of the H2 Color GPU backend.  See device_kernels.hpp.

#include "device_kernels.hpp"

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

__device__ __forceinline__ double kernel_value(const KernelSpec& spec,
                                               const double* x, int64_t x_id,
                                               const double* y, int64_t y_id) {
    // kind 1 (inverse distance); the only kind accepted on the host so far
    if (x_id == y_id) return spec.p[1];
    const double dx = x[0] - y[0];
    const double dy = x[1] - y[1];
    const double dz = x[2] - y[2];
    const double r = sqrt(dx * dx + dy * dy + dz * dz);
    return spec.p[0] / r;
}

__global__ void eval_kernel(const EvalItem* items, const char* meta, KernelSpec spec,
                            PointTable points) {
    const EvalItem item = items[blockIdx.x];
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
            kernel_value(spec, x, x_id, points.xyz + 3 * static_cast<int64_t>(sj), points.ids[sj]);
    }
}

__global__ void gather_kernel(const GatherItem* items, const char* meta) {
    const GatherItem item = items[blockIdx.x];
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

__global__ void add_store_kernel(const AddStoreItem* items) {
    const AddStoreItem item = items[blockIdx.x];
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

__global__ void sym_add_kernel(const SymAddItem* items) {
    const SymAddItem item = items[blockIdx.x];
    const int i = blockIdx.y * kTile + threadIdx.x;
    const int j0 = blockIdx.z * kTile;
    if (i >= item.n || j0 >= item.n) return;
    for (int jj = threadIdx.y; jj < kTile; jj += kTileRows) {
        const int j = j0 + jj;
        if (j >= item.n) break;
        double& x = item.x[i + static_cast<int64_t>(j) * item.ldx];
        x += item.t[i + static_cast<int64_t>(j) * item.ldt] + item.t[j + static_cast<int64_t>(i) * item.ldt];
    }
}

__global__ void sum_add_kernel(const SumAddItem* items, const char* meta) {
    const SumAddItem item = items[blockIdx.x];
    const int i = blockIdx.y * kTile + threadIdx.x;
    const int j0 = blockIdx.z * kTile;
    if (i >= item.rows || j0 >= item.cols) return;
    const double* const* parts = reinterpret_cast<const double* const*>(meta + item.parts_offset);
    for (int jj = threadIdx.y; jj < kTile; jj += kTileRows) {
        const int j = j0 + jj;
        if (j >= item.cols) break;
        const int64_t at = i + static_cast<int64_t>(j) * item.rows;
        double sum = parts[0][at];
        for (int p = 1; p < item.nparts; ++p) sum += parts[p][at];
        item.target[i + static_cast<int64_t>(j) * item.ldt] += sum;
    }
}

__global__ void identity_kernel(const IdentityItem* items) {
    const IdentityItem item = items[blockIdx.x];
    const int i = blockIdx.y * kTile + threadIdx.x;
    const int j0 = blockIdx.z * kTile;
    if (i >= item.n || j0 >= item.n) return;
    for (int jj = threadIdx.y; jj < kTile; jj += kTileRows) {
        const int j = j0 + jj;
        if (j >= item.n) break;
        item.a[i + static_cast<int64_t>(j) * item.ld] = i == j ? 1.0 : 0.0;
    }
}

__global__ void column_swap_kernel(const ColumnSwapItem* items) {
    const ColumnSwapItem item = items[blockIdx.x];
    // Each thread owns whole rows, so the sequential swaps need no barrier.
    for (int row = threadIdx.x; row < item.rows; row += blockDim.x) {
        for (int i = item.n - 1; i >= 0; --i) {
            const int p = item.piv[i] - 1;
            if (p == i) continue;
            double* ci = item.b + static_cast<int64_t>(i) * item.ldb + row;
            double* cp = item.b + static_cast<int64_t>(p) * item.ldb + row;
            const double t = *ci;
            *ci = *cp;
            *cp = t;
        }
    }
}

__global__ void transpose_kernel(const TransposeItem* items) {
    __shared__ double tile[kTile][kTile + 1];
    const TransposeItem item = items[blockIdx.x];
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

void launch_eval(const EvalItem* items, int count, int max_m, int max_n,
                 const char* meta, KernelSpec spec, PointTable points, cudaStream_t stream) {
    if (count <= 0 || max_m <= 0 || max_n <= 0) return;
    if (spec.kind != 1) throw std::runtime_error("launch_eval: unsupported device kernel kind");
    eval_kernel<<<tile_grid(count, max_m, max_n), dim3(kTile, kTileRows), 0, stream>>>(
        items, meta, spec, points);
    check_launch("eval_kernel");
}

void launch_gather(const GatherItem* items, int count, int max_m, int max_n,
                   const char* meta, cudaStream_t stream) {
    if (count <= 0 || max_m <= 0 || max_n <= 0) return;
    gather_kernel<<<tile_grid(count, max_m, max_n), dim3(kTile, kTileRows), 0, stream>>>(items, meta);
    check_launch("gather_kernel");
}

void launch_add_store(const AddStoreItem* items, int count, int max_m, int max_n,
                      cudaStream_t stream) {
    if (count <= 0 || max_m <= 0 || max_n <= 0) return;
    add_store_kernel<<<tile_grid(count, max_m, max_n), dim3(kTile, kTileRows), 0, stream>>>(items);
    check_launch("add_store_kernel");
}

void launch_sym_add(const SymAddItem* items, int count, int max_n, cudaStream_t stream) {
    if (count <= 0 || max_n <= 0) return;
    sym_add_kernel<<<tile_grid(count, max_n, max_n), dim3(kTile, kTileRows), 0, stream>>>(items);
    check_launch("sym_add_kernel");
}

void launch_sum_add(const SumAddItem* items, int count, int max_rows, int max_cols, const char* meta,
                    cudaStream_t stream) {
    if (count <= 0 || max_rows <= 0 || max_cols <= 0) return;
    sum_add_kernel<<<tile_grid(count, max_rows, max_cols), dim3(kTile, kTileRows), 0, stream>>>(items, meta);
    check_launch("sum_add_kernel");
}

void launch_identity(const IdentityItem* items, int count, int max_n, cudaStream_t stream) {
    if (count <= 0 || max_n <= 0) return;
    identity_kernel<<<tile_grid(count, max_n, max_n), dim3(kTile, kTileRows), 0, stream>>>(items);
    check_launch("identity_kernel");
}

void launch_transpose(const TransposeItem* items, int count, int max_m, int max_n, cudaStream_t stream) {
    if (count <= 0 || max_m <= 0 || max_n <= 0) return;
    transpose_kernel<<<tile_grid(count, max_m, max_n), dim3(kTile, kTileRows), 0, stream>>>(items);
    check_launch("transpose_kernel");
}

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

void launch_column_swaps(const ColumnSwapItem* items, int count, cudaStream_t stream) {
    if (count <= 0) return;
    column_swap_kernel<<<static_cast<unsigned>(count), 128, 0, stream>>>(items);
    check_launch("column_swap_kernel");
}

}  // namespace gpu
}  // namespace fmm
