// The kernels of a butterfly's construction from a black-box matvec on the
// device (bf_randomized.hpp): the random vectors, the panels the range finder
// factorizes, Q formed from launch_qrcp's reflectors and split into kernel
// blocks, and the least-squares solve of the level where the halves meet.

#include "bf_randomized.hpp"

#include <stdexcept>
#include <string>

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

// x(row0 + r, col0 + c) uniform in [-1, 1), from a counter-based hash of the
// row's global index (hash_row0 + row0 + r) and the column: the same numbers
// whatever the launch, and whichever rank holds the row
__global__ void fill_uniform_kernel(DataType* x, int64_t ld, int64_t row0, int64_t nrows, int64_t col0,
                                    int64_t ncols, uint64_t seed, int64_t hash_row0) {
    const int64_t total = nrows * ncols;
    for (int64_t t = blockIdx.x * static_cast<int64_t>(blockDim.x) + threadIdx.x; t < total;
         t += static_cast<int64_t>(gridDim.x) * blockDim.x) {
        const int64_t r = t % nrows, c = t / nrows;
        const uint64_t h = mix64(seed ^ mix64(static_cast<uint64_t>(hash_row0 + row0 + r) * 0x100000001b3ull +
                                              static_cast<uint64_t>(col0 + c)));
        const double u = static_cast<double>(h >> 11) * (1.0 / 9007199254740992.0);  // [0, 1)
        x[(row0 + r) + (col0 + c) * ld] = 2.0 * u - 1.0;
    }
}

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

// One thread block per panel: rows m0 of s0 then m1 of s1, columns
// col0 .. col0 + ncols, into dst (or its transpose), and the largest
// magnitude when asked.
__global__ void __launch_bounds__(kThreads) panel_kernel(const PanelItem* items, double* amax) {
    const PanelItem p = items[blockIdx.x];
    const int m = p.m0 + p.m1;
    const int64_t total = static_cast<int64_t>(m) * p.ncols;
    double big = 0.0;
    for (int64_t t = threadIdx.x; t < total; t += kThreads) {
        const int r = static_cast<int>(t % m), c = static_cast<int>(t / m);
        const DataType v = r < p.m0 ? p.s0[r + static_cast<int64_t>(p.col0 + c) * p.ld0]
                                    : p.s1[(r - p.m0) + static_cast<int64_t>(p.col0 + c) * p.ld1];
        if (p.transposed) {
            p.dst[c + static_cast<int64_t>(r) * p.ldd] = v;
        } else {
            p.dst[r + static_cast<int64_t>(c) * p.ldd] = v;
        }
        big = fmax(big, fabs(v));
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

// One thread block per panel: Q = H_0 ... H_{rank-1} [I; 0] (LAPACK dorg2r)
// from the reflectors launch_qrcp left below the diagonal of a, then Q split
// into the two outputs.  A warp per column of Q.
__global__ void __launch_bounds__(kThreads) range_kernel(const RangeItem* items) {
    const RangeItem it = items[blockIdx.x];
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
                double s = 0.0;
                for (int r = i + lane; r < m; r += 32) s += (r == i ? 1.0 : v[r]) * q[r];
                s = warp_sum(s);
                const double w = tau * s;
                for (int r = i + lane; r < m; r += 32) q[r] -= w * (r == i ? 1.0 : v[r]);
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

// One thread block per solve: B := Q^T B with launch_qrcp's reflectors of A,
// X1 = R11^{-1} B(0:k, :) / scale (the factorization ran on A / scale), then
// X(jpvt[t], :) = X1(t, :), written transposed into the two outputs.  A warp
// per column of B.
__global__ void __launch_bounds__(kThreads) solve_kernel(const SolveItem* items) {
    const SolveItem it = items[blockIdx.x];
    const int mm = it.mm, n = it.n;
    const int k = it.k != nullptr ? *it.k : 0;
    const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5;
    const int n1 = n - it.m0;
    for (int c = warp; c < n; c += kWarps) {
        DataType* b = it.b + static_cast<int64_t>(c) * mm;
        for (int i = 0; i < k; ++i) {
            const DataType* v = it.a + static_cast<int64_t>(i) * mm;
            double s = 0.0;
            for (int r = i + lane; r < mm; r += 32) s += (r == i ? 1.0 : v[r]) * b[r];
            s = warp_sum(s);
            const double w = it.tau[i] * s;
            for (int r = i + lane; r < mm; r += 32) b[r] -= w * (r == i ? 1.0 : v[r]);
            __syncwarp();
        }
        // back substitution, R11 in the upper triangle of a
        for (int t = k - 1; t >= 0; --t) {
            double s = 0.0;
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
__global__ void __launch_bounds__(kThreads) copy_kernel(const CopyItem* items) {
    const CopyItem it = items[blockIdx.x];
    for (int64_t t = threadIdx.x; t < it.n; t += kThreads) it.dst[t] = it.src[t];
}

int grid_for(int64_t total) {
    const int64_t blocks = (total + kThreads - 1) / kThreads;
    return static_cast<int>(blocks < 4096 ? (blocks > 0 ? blocks : 1) : 4096);
}

}  // namespace

void launch_fill_uniform(DataType* x, int64_t ld, int64_t row0, int64_t nrows, int64_t col0, int64_t ncols,
                         uint64_t seed, int64_t hash_row0, cudaStream_t stream) {
    if (nrows <= 0 || ncols <= 0) return;
    fill_uniform_kernel<<<grid_for(nrows * ncols), kThreads, 0, stream>>>(x, ld, row0, nrows, col0, ncols, seed,
                                                                          hash_row0);
    check_launch("launch_fill_uniform");
}

void launch_fill_value(DataType* x, int64_t n, DataType v, cudaStream_t stream) {
    if (n <= 0) return;
    fill_value_kernel<<<grid_for(n), kThreads, 0, stream>>>(x, n, v);
    check_launch("launch_fill_value");
}

void launch_panels(const PanelItem* items, int count, double* amax, cudaStream_t stream) {
    if (count <= 0) return;
    panel_kernel<<<count, kThreads, 0, stream>>>(items, amax);
    check_launch("launch_panels");
}

void launch_range(const RangeItem* items, int count, cudaStream_t stream) {
    if (count <= 0) return;
    range_kernel<<<count, kThreads, 0, stream>>>(items);
    check_launch("launch_range");
}

void launch_solve(const SolveItem* items, int count, cudaStream_t stream) {
    if (count <= 0) return;
    solve_kernel<<<count, kThreads, 0, stream>>>(items);
    check_launch("launch_solve");
}

void launch_copies(const CopyItem* items, int count, cudaStream_t stream) {
    if (count <= 0) return;
    copy_kernel<<<count, kThreads, 0, stream>>>(items);
    check_launch("launch_copies");
}

}  // namespace gpu
}  // namespace bf
