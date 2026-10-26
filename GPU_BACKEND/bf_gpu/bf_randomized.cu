// The kernels of a butterfly's construction from a black-box matvec on the
// device (bf_randomized.hpp): the random vectors, the panels the range finder
// factorizes, Q formed from launch_qrcp's reflectors and split into kernel
// blocks, and the least-squares solve of the level where the halves meet.

#include "bf_randomized.hpp"

#include <stdexcept>
#include <string>
#include <type_traits>

namespace bf {
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

__device__ __forceinline__ uint64_t mix64(uint64_t z) {
    z += 0x9e3779b97f4a7c15ull;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    return z ^ (z >> 31);
}

// x(row0 + r, col0 + c) uniform in [-1, 1) - a complex x's two parts each -
// from a counter-based hash of the row's global index (hash_row0 + row0 + r)
// and the column: the same numbers whatever the launch, and whichever rank
// holds the row
template<typename T>
__global__ void fill_uniform_kernel(T* x, int64_t ld, int64_t row0, int64_t nrows, int64_t col0, int64_t ncols,
                                    uint64_t seed, int64_t hash_row0) {
    const int64_t total = nrows * ncols;
    for (int64_t t = blockIdx.x * static_cast<int64_t>(blockDim.x) + threadIdx.x; t < total;
         t += static_cast<int64_t>(gridDim.x) * blockDim.x) {
        const int64_t r = t % nrows, c = t / nrows;
        const uint64_t h = mix64(seed ^ mix64(static_cast<uint64_t>(hash_row0 + row0 + r) * 0x100000001b3ull +
                                              static_cast<uint64_t>(col0 + c)));
        const double u = static_cast<double>(h >> 11) * (1.0 / 9007199254740992.0);  // [0, 1)
        if constexpr (std::is_same<T, double>::value) {
            x[(row0 + r) + (col0 + c) * ld] = 2.0 * u - 1.0;
        } else {
            const double v = static_cast<double>(mix64(h ^ 0x632be59bd9b4e019ull) >> 11) * (1.0 / 9007199254740992.0);
            x[(row0 + r) + (col0 + c) * ld] = T(2.0 * u - 1.0, 2.0 * v - 1.0);
        }
    }
}

template<typename DataType>
__global__ void fill_value_kernel(DataType* x, int64_t n, DataType v) {
    for (int64_t t = blockIdx.x * static_cast<int64_t>(blockDim.x) + threadIdx.x; t < n;
         t += static_cast<int64_t>(gridDim.x) * blockDim.x) {
        x[t] = v;
    }
}

__device__ __forceinline__ double warp_sum(double v) {
    for (int o = 16; o > 0; o >>= 1) v += __shfl_down_sync(0xffffffffu, v, o);
    return __shfl_sync(0xffffffffu, v, 0);
}
__device__ __forceinline__ fmm::gpu::dcomplex warp_sum(fmm::gpu::dcomplex v) {
    return fmm::gpu::dcomplex(warp_sum(v.re), warp_sum(v.im));
}

// One thread block per panel: rows m0 of s0 then m1 of s1, columns
// col0 .. col0 + ncols, into dst (or its transpose), and the largest
// magnitude when asked.
template<typename DataType>
__global__ void __launch_bounds__(kThreads) panel_kernel(const PanelItem<DataType>* items, double* amax) {
    const PanelItem<DataType> p = items[blockIdx.x];
    const int m = p.m0 + p.m1;
    const int64_t total = static_cast<int64_t>(m) * p.ncols;
    double big = 0.0;
    for (int64_t t = threadIdx.x; t < total; t += kThreads) {
        const int r = static_cast<int>(t % m), c = static_cast<int>(t / m);
        const DataType v0 = r < p.m0 ? p.s0[r + static_cast<int64_t>(p.col0 + c) * p.ld0]
                                     : p.s1[(r - p.m0) + static_cast<int64_t>(p.col0 + c) * p.ld1];
        const DataType v = p.conj ? fmm::gpu::conj(v0) : v0;
        if (p.transposed) {
            p.dst[c + static_cast<int64_t>(r) * p.ldd] = v;
        } else {
            p.dst[r + static_cast<int64_t>(c) * p.ldd] = v;
        }
        big = fmax(big, fmm::gpu::magnitude(v));
    }
    if (amax == nullptr) return;
    __shared__ double s[kWarps];
    for (int o = 16; o > 0; o >>= 1) big = fmax(big, __shfl_down_sync(0xffffffffu, big, o));
    if ((threadIdx.x & 31) == 0) s[threadIdx.x >> 5] = big;
    __syncthreads();
    if (threadIdx.x == 0) {
        double b = 0.0;
        for (int w = 0; w < kWarps; ++w) b = fmax(b, s[w]);
        amax[blockIdx.x] = b;
    }
}

// One thread block per panel: Q = H_0 ... H_{rank-1} [I; 0] (LAPACK dorg2r,
// zung2r: H_i = I - tau_i v_i v_i^H) from the reflectors launch_qrcp left
// below the diagonal of a, then Q split into the two outputs.  A warp per
// column of Q.
template<typename DataType>
__global__ void __launch_bounds__(kThreads) range_kernel(const RangeItem<DataType>* items) {
    using fmm::gpu::conj;
    const RangeItem<DataType> it = items[blockIdx.x];
    const int m = it.m, k = it.rank;
    const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5;
    DataType* Q = it.q;
    for (int64_t t = threadIdx.x; t < static_cast<int64_t>(m) * k; t += kThreads) {
        const int r = static_cast<int>(t % m), c = static_cast<int>(t / m);
        Q[t] = (!it.zero && r == c) ? DataType(1) : DataType(0);
    }
    __syncthreads();
    if (!it.zero) {
        for (int i = k - 1; i >= 0; --i) {
            const DataType tau = it.tau[i];
            const DataType* v = it.a + static_cast<int64_t>(i) * m;  // v_i = 1, then v[i+1 ..]
            for (int c = i + warp; c < k; c += kWarps) {
                DataType* q = Q + static_cast<int64_t>(c) * m;
                DataType s = DataType(0.0);
                for (int r = i + lane; r < m; r += 32) s += conj(r == i ? DataType(1.0) : v[r]) * q[r];
                s = warp_sum(s);
                const DataType w = tau * s;
                for (int r = i + lane; r < m; r += 32) q[r] -= w * (r == i ? DataType(1.0) : v[r]);
            }
            __syncthreads();
        }
    }
    // the outputs: rows 0 .. m0 and m0 .. m of Q, or of Q^T
    const int m1 = m - it.m0;
    for (int64_t t = threadIdx.x; t < static_cast<int64_t>(m) * k; t += kThreads) {
        const int r = static_cast<int>(t % m), c = static_cast<int>(t / m);
        const DataType v = Q[t];
        if (r < it.m0) {
            if (it.transposed) it.out0[c + static_cast<int64_t>(r) * k] = v;
            else it.out0[r + static_cast<int64_t>(c) * it.m0] = v;
        } else if (it.out1 != nullptr) {
            const int r1 = r - it.m0;
            if (it.transposed) it.out1[c + static_cast<int64_t>(r1) * k] = v;
            else it.out1[r1 + static_cast<int64_t>(c) * m1] = v;
        }
    }
}

// One thread block per solve: B := Q^H B with launch_qrcp's reflectors of A
// (each H_i^H = I - conj(tau_i) v_i v_i^H), X1 = R11^{-1} B(0:k, :) / scale
// (the factorization ran on A / scale), then X(jpvt[t], :) = X1(t, :),
// written transposed into the two outputs.  A warp per column of B.
template<typename DataType>
__global__ void __launch_bounds__(kThreads) solve_kernel(const SolveItem<DataType>* items) {
    using fmm::gpu::conj;
    const SolveItem<DataType> it = items[blockIdx.x];
    const int mm = it.mm, n = it.n;
    const int k = it.k != nullptr ? *it.k : 0;
    const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5;
    const int n1 = n - it.m0;
    for (int c = warp; c < n; c += kWarps) {
        DataType* b = it.b + static_cast<int64_t>(c) * mm;
        for (int i = 0; i < k; ++i) {
            const DataType* v = it.a + static_cast<int64_t>(i) * mm;
            DataType s = DataType(0.0);
            for (int r = i + lane; r < mm; r += 32) s += conj(r == i ? DataType(1.0) : v[r]) * b[r];
            s = warp_sum(s);
            const DataType w = conj(it.tau[i]) * s;
            for (int r = i + lane; r < mm; r += 32) b[r] -= w * (r == i ? DataType(1.0) : v[r]);
            __syncwarp();
        }
        // back substitution, R11 in the upper triangle of a
        for (int t = k - 1; t >= 0; --t) {
            DataType s = DataType(0.0);
            for (int u = t + 1 + lane; u < k; u += 32) s += it.a[t + static_cast<int64_t>(u) * mm] * b[u];
            s = warp_sum(s);
            if (lane == 0) b[t] = (b[t] - s) / it.a[t + static_cast<int64_t>(t) * mm];
            __syncwarp();
        }
        // X(jpvt[t], c) = b[t] / scale, into out0 (c < m0) or out1, transposed
        for (int t = lane; t < k; t += 32) {
            const int row = it.jpvt[t];
            const DataType x = b[t] / *it.scale;
            if (c < it.m0) it.out0[c + static_cast<int64_t>(row) * it.m0] = x;
            else it.out1[(c - it.m0) + static_cast<int64_t>(row) * n1] = x;
        }
    }
}

// One thread block per item.
template<typename DataType>
__global__ void __launch_bounds__(kThreads) copy_kernel(const CopyItem<DataType>* items) {
    const CopyItem<DataType> it = items[blockIdx.x];
    for (int64_t t = threadIdx.x; t < it.n; t += kThreads) it.dst[t] = it.src[t];
}

int grid_for(int64_t total) {
    const int64_t blocks = (total + kThreads - 1) / kThreads;
    return static_cast<int>(blocks < 4096 ? (blocks > 0 ? blocks : 1) : 4096);
}

}  // namespace

template<typename DataType>
void launch_fill_uniform(DataType* x, int64_t ld, int64_t row0, int64_t nrows, int64_t col0, int64_t ncols,
                         uint64_t seed, int64_t hash_row0, cudaStream_t stream) {
    if (nrows <= 0 || ncols <= 0) return;
    fill_uniform_kernel<DataType><<<grid_for(nrows * ncols), kThreads, 0, stream>>>(x, ld, row0, nrows, col0, ncols, seed,
                                                                          hash_row0);
    check_launch("launch_fill_uniform");
}

template<typename DataType>
void launch_fill_value(DataType* x, int64_t n, DataType v, cudaStream_t stream) {
    if (n <= 0) return;
    fill_value_kernel<DataType><<<grid_for(n), kThreads, 0, stream>>>(x, n, v);
    check_launch("launch_fill_value");
}

template<typename DataType>
void launch_panels(const PanelItem<DataType>* items, int count, double* amax, cudaStream_t stream) {
    if (count <= 0) return;
    panel_kernel<DataType><<<count, kThreads, 0, stream>>>(items, amax);
    check_launch("launch_panels");
}

template<typename DataType>
void launch_range(const RangeItem<DataType>* items, int count, cudaStream_t stream) {
    if (count <= 0) return;
    range_kernel<DataType><<<count, kThreads, 0, stream>>>(items);
    check_launch("launch_range");
}

template<typename DataType>
void launch_solve(const SolveItem<DataType>* items, int count, cudaStream_t stream) {
    if (count <= 0) return;
    solve_kernel<DataType><<<count, kThreads, 0, stream>>>(items);
    check_launch("launch_solve");
}

template<typename DataType>
void launch_copies(const CopyItem<DataType>* items, int count, cudaStream_t stream) {
    if (count <= 0) return;
    copy_kernel<DataType><<<count, kThreads, 0, stream>>>(items);
    check_launch("launch_copies");
}

// the double and the double complex libraries
#define BF_RANDOMIZED_INSTANTIATE(T)                                                                                  \
    template void launch_fill_uniform<T>(T*, int64_t, int64_t, int64_t, int64_t, int64_t, uint64_t, int64_t,          \
                                         cudaStream_t);                                                               \
    template void launch_fill_value<T>(T*, int64_t, T, cudaStream_t);                                                 \
    template void launch_panels<T>(const PanelItem<T>*, int, double*, cudaStream_t);                                   \
    template void launch_range<T>(const RangeItem<T>*, int, cudaStream_t);                                             \
    template void launch_solve<T>(const SolveItem<T>*, int, cudaStream_t);                                             \
    template void launch_copies<T>(const CopyItem<T>*, int, cudaStream_t);
BF_RANDOMIZED_INSTANTIATE(double)
BF_RANDOMIZED_INSTANTIATE(fmm::gpu::dcomplex)
#undef BF_RANDOMIZED_INSTANTIATE

}  // namespace gpu
}  // namespace bf
