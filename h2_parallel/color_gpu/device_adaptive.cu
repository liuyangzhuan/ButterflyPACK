// Device kernels of the adaptive ID training rows.  See device_adaptive.hpp.
// Compiled without FMA contraction (--fmad=false, CMakeLists.txt), as the
// device ID: on a replicated CA level every copy of a box must select the
// same rows.

#include "device_adaptive.hpp"

#include <stdexcept>
#include <string>

namespace fmm {
namespace gpu {

namespace {

constexpr int kThreads = 256;
constexpr int kWarps = kThreads / 32;

void check_launch(const char* what) {
    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA launch failed in ") + what + ": " + cudaGetErrorString(status));
    }
}

// Sum of v over the block, the same on every thread: lanes of a warp by
// shuffles, then the warps in order.
__device__ double block_sum(double v, double* scratch) {
    for (int o = 16; o > 0; o >>= 1) v += __shfl_xor_sync(0xffffffffu, v, o);
    __syncthreads();
    if ((threadIdx.x & 31) == 0) scratch[threadIdx.x >> 5] = v;
    __syncthreads();
    double sum = 0.0;
    for (int w = 0; w < kWarps; ++w) sum += scratch[w];
    return sum;
}

template<typename T>
__global__ void __launch_bounds__(kThreads) id_residual_kernel(const IdResidualItemT<T>* items, const char* meta) {
    __shared__ double scratch[kWarps];
    const IdResidualItemT<T> item = items[blockIdx.x];
    const int rows = item.rows, n = item.n, tid = threadIdx.x;
    const IdRefT<T>* ids = reinterpret_cast<const IdRefT<T>*>(meta + item.ids_offset);
    for (int e = 0; e < item.nids; ++e) {
        const IdRefT<T> id = ids[e];
        const int rank = id.rank_ptr != nullptr ? *id.rank_ptr : id.rank;  // (the same on every thread)
        if (rank <= 0) continue;
        if (rank < n) {
            // redundant columns: each entry by one thread, the skeleton's
            // terms in order (the skeleton columns are only read here)
            const int64_t entries = static_cast<int64_t>(rows) * (n - rank);
            for (int64_t p = tid; p < entries; p += kThreads) {
                const int row = static_cast<int>(p % rows), j = static_cast<int>(p / rows);
                const T* t = id.f + static_cast<int64_t>(rank + j) * id.ld;
                T sum = T(0.0);
                for (int i = 0; i < rank; ++i) {
                    sum += item.r[row + static_cast<int64_t>(id.jpvt[i]) * item.ld] * t[i];
                }
                T* x = item.r + row + static_cast<int64_t>(id.jpvt[rank + j]) * item.ld;
                *x = *x - sum;
            }
            __syncthreads();
        }
        const int skeleton = rank < n ? rank : n;  // rank n: every column, in the natural order
        const int64_t entries = static_cast<int64_t>(rows) * skeleton;
        for (int64_t p = tid; p < entries; p += kThreads) {
            const int row = static_cast<int>(p % rows), i = static_cast<int>(p / rows);
            item.r[row + static_cast<int64_t>(rank < n ? id.jpvt[i] : i) * item.ld] = T(0.0);
        }
        __syncthreads();
    }
    if (item.norm == nullptr) return;
    double ssq = 0.0;
    const int64_t entries = static_cast<int64_t>(rows) * n;
    for (int64_t p = tid; p < entries; p += kThreads) {
        ssq += abs2(item.r[p % rows + (p / rows) * item.ld]);
    }
    ssq = block_sum(ssq, scratch);
    if (tid == 0) *item.norm = sqrt(ssq);
}

template<typename T>
__global__ void __launch_bounds__(kThreads) append_rows_kernel(const AppendRowsItemT<T>* items) {
    const AppendRowsItemT<T> item = items[blockIdx.x];
    const int count = *item.rank;
    if (count <= 0) return;
    const bool natural = count >= item.s;
    const int64_t entries = static_cast<int64_t>(count) * item.n;
    for (int64_t p = threadIdx.x; p < entries; p += kThreads) {
        const int i = static_cast<int>(p % count);
        const int64_t col = p / count;
        item.dst[item.row0 + i + col * item.ldd] = item.src[(natural ? i : item.jpvt[i]) + col * item.lds];
    }
}

}  // namespace

template<typename T>
void launch_id_residual(const IdResidualItemT<T>* items, int count, const char* meta, cudaStream_t stream) {
    if (count <= 0) return;
    id_residual_kernel<T><<<static_cast<unsigned>(count), kThreads, 0, stream>>>(items, meta);
    check_launch("id_residual_kernel");
}

template<typename T>
void launch_append_rows(const AppendRowsItemT<T>* items, int count, cudaStream_t stream) {
    if (count <= 0) return;
    append_rows_kernel<T><<<static_cast<unsigned>(count), kThreads, 0, stream>>>(items);
    check_launch("append_rows_kernel");
}

template void launch_id_residual<double>(const IdResidualItemT<double>*, int, const char*, cudaStream_t);
template void launch_id_residual<dcomplex>(const IdResidualItemT<dcomplex>*, int, const char*, cudaStream_t);
template void launch_append_rows<double>(const AppendRowsItemT<double>*, int, cudaStream_t);
template void launch_append_rows<dcomplex>(const AppendRowsItemT<dcomplex>*, int, cudaStream_t);

}  // namespace gpu
}  // namespace fmm
