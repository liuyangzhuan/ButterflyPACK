// Small batched kernels of the HODLR GPU backend (see hodlr_kernels.hpp).

#include "hodlr_kernels.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>

namespace bpack {
namespace gpu {

using fmm::gpu::dcomplex;

namespace {

constexpr int kThreads = 256;

void check_launch(const char* what) {
    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA launch failed in ") + what + ": " + cudaGetErrorString(status));
    }
}

unsigned tiles(long long work, long long cap = 4096) {
    return static_cast<unsigned>(std::max<long long>(1, std::min<long long>((work + kThreads - 1) / kThreads, cap)));
}

template<typename T>
__global__ void row_swaps_kernel(const RowSwapItem<T>* items) {
    const RowSwapItem<T> it = items[blockIdx.x];
    const int j = blockIdx.y * blockDim.x + threadIdx.x;
    if (j >= it.n) return;
    T* col = it.b + static_cast<long long>(j) * it.ldb;
    for (int i = 0; i < it.m; ++i) {
        const int p = it.ipiv[i] - 1;
        if (p != i) {
            const T t = col[i];
            col[i] = col[p];
            col[p] = t;
        }
    }
}

// One block of kPivotThreads threads per item: the columns are spread over
// the threads, the reductions (argmax, Householder norm) go through shared
// memory in a fixed order.
constexpr int kPivotThreads = 1024;

__device__ __forceinline__ double sq_abs(double v) { return v * v; }
__device__ __forceinline__ double sq_abs(dcomplex v) { return v.re * v.re + v.im * v.im; }
__device__ __forceinline__ dcomplex conj_of(dcomplex v) { return {v.re, -v.im}; }
__device__ __forceinline__ double conj_of(double v) { return v; }
__device__ __forceinline__ double real_of(double v) { return v; }
__device__ __forceinline__ double real_of(dcomplex v) { return v.re; }
__device__ __forceinline__ double imag_of(double) { return 0.0; }
__device__ __forceinline__ double imag_of(dcomplex v) { return v.im; }

template<typename T>
__global__ void __launch_bounds__(kPivotThreads) qrcp_pivots_kernel(const PivotItem<T>* items, const char* meta) {
    const PivotItem<T> it = items[blockIdx.x];
    const int r = it.r, L = it.L;
    // w[l + i * L]: row i of the (permuted) column l, so that the threads,
    // one column each, read and write consecutive addresses
    T* w = it.work;
    auto at = [&](int i, int l) -> T& { return w[l + static_cast<long long>(i) * L]; };
    double* vn1 = it.norms;
    double* vn2 = it.norms + L;
    int* perm = it.perm;
    __shared__ double s_val[kPivotThreads];
    __shared__ int s_idx[kPivotThreads];
    __shared__ T s_tau;  // (the Householder vector overwrites rows k+1.. of column k, v_k = 1)
    const double tol3z = sqrt(1.1102230246251565e-16);  // sqrt(dlamch('Epsilon'))

    // copy, with the masked columns zeroed, and the column norms
    for (int l = threadIdx.x; l < L; l += blockDim.x) perm[l] = l;
    __syncthreads();
    const int* mask = it.mask_ptr != nullptr ? it.mask_ptr : reinterpret_cast<const int*>(meta + it.mask_offset);
    for (int m = threadIdx.x; m < it.nmask; m += blockDim.x) perm[mask[m]] = -1;  // (flag)
    __syncthreads();
    for (int l = threadIdx.x; l < L; l += blockDim.x) {
        const bool zero = perm[l] < 0;
        double s = 0.0;
        for (int i = 0; i < r; ++i) {
            const T v = zero ? T(0.0) : (it.transposed ? it.a[l + static_cast<long long>(i) * it.lda]
                                                       : it.a[i + static_cast<long long>(l) * it.lda]);
            at(i, l) = v;
            s += sq_abs(v);
        }
        vn1[l] = sqrt(s);
        vn2[l] = vn1[l];
        perm[l] = l;
    }
    __syncthreads();

    const int steps = min(it.p, min(r, L));
    for (int k = 0; k < steps; ++k) {
        // pivot: the first column of largest norm among k..L-1 (idamax)
        double best = -1.0;
        int bi = L;
        for (int l = k + threadIdx.x; l < L; l += blockDim.x) {
            if (vn1[l] > best) {
                best = vn1[l];
                bi = l;
            }
        }
        s_val[threadIdx.x] = best;
        s_idx[threadIdx.x] = bi;
        __syncthreads();
        for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
            if (threadIdx.x < stride) {
                const double v2 = s_val[threadIdx.x + stride];
                const int i2 = s_idx[threadIdx.x + stride];
                if (v2 > s_val[threadIdx.x] || (v2 == s_val[threadIdx.x] && i2 < s_idx[threadIdx.x])) {
                    s_val[threadIdx.x] = v2;
                    s_idx[threadIdx.x] = i2;
                }
            }
            __syncthreads();
        }
        const int pvt = s_idx[0];
        __syncthreads();
        if (pvt != k) {  // swap columns k and pvt (data, norms, permutation)
            for (int i = threadIdx.x; i < r; i += blockDim.x) {
                const T t = at(i, k);
                at(i, k) = at(i, pvt);
                at(i, pvt) = t;
            }
            if (threadIdx.x == 0) {
                const int tp = perm[k];
                perm[k] = perm[pvt];
                perm[pvt] = tp;
                vn1[pvt] = vn1[k];
                vn2[pvt] = vn2[k];
            }
        }
        __syncthreads();
        // Householder reflector of column k, rows k..r-1 (LAPACK ?larfg; the
        // real one leaves a single row alone, the complex one makes it real)
        if (threadIdx.x == 0) {
            double xnorm2 = 0.0;
            for (int i = k + 1; i < r; ++i) xnorm2 += sq_abs(at(i, k));
            const T alpha = at(k, k);
            const double ar = real_of(alpha), ai = imag_of(alpha);
            if (xnorm2 == 0.0 && ai == 0.0) {
                s_tau = T(0.0);
            } else {
                const double nrm = sqrt(ar * ar + ai * ai + xnorm2);
                const double beta = ar >= 0.0 ? -nrm : nrm;
                T tau;
                if constexpr (fmm::gpu::is_complex_scalar<T>) {
                    tau = T((beta - ar) / beta, -ai / beta);
                } else {
                    tau = T((beta - ar) / beta);
                }
                const T scale = T(1.0) / (alpha - T(beta));
                s_tau = tau;
                for (int i = k + 1; i < r; ++i) at(i, k) = at(i, k) * scale;
                at(k, k) = T(beta);
            }
        }
        __syncthreads();
        // apply H^H = I - conj(tau) v v^H to columns k+1..L-1, then downdate their norms (?laqp2)
        const T tau = s_tau;
        const T ctau = conj_of(tau);
        for (int l = k + 1 + threadIdx.x; l < L; l += blockDim.x) {
            if (tau != T(0.0)) {
                T dot = at(k, l);  // v = (1, w(k+1.., k))
                for (int i = k + 1; i < r; ++i) dot += conj_of(at(i, k)) * at(i, l);
                const T f = ctau * dot;
                at(k, l) -= f;
                for (int i = k + 1; i < r; ++i) at(i, l) -= f * at(i, k);
            }
            if (vn1[l] != 0.0) {
                double temp = sqrt(sq_abs(at(k, l))) / vn1[l];
                temp = 1.0 - temp * temp;
                temp = temp > 0.0 ? temp : 0.0;
                const double ratio = vn1[l] / vn2[l];
                const double temp2 = temp * ratio * ratio;
                if (temp2 <= tol3z) {
                    double s = 0.0;
                    for (int i = k + 1; i < r; ++i) s += sq_abs(at(i, l));
                    vn1[l] = sqrt(s);
                    vn2[l] = vn1[l];
                } else {
                    vn1[l] *= sqrt(temp);
                }
            }
        }
        __syncthreads();
    }
    for (int k = threadIdx.x; k < min(it.p, L); k += blockDim.x) it.out[k] = perm[k];
}

// ---- the column-pivoted QR pivots with the columns over thread blocks ----
// Per item, in norms + 2 L: the block maxima of the norms (nblk doubles,
// then nblk ints) and tau.  Each step: wide_pivot_kernel (one block per
// item) reduces the block maxima, swaps and forms the reflector;
// wide_apply_kernel applies it to the columns of its block, downdates their
// norms and leaves the block maximum for the next step.  The arithmetic per
// column is that of qrcp_pivots_kernel.
constexpr int kWideThreads = 256;

template<typename T>
struct WideLayout {
    int nblk;
    double* part_val;
    int* part_idx;
    T* tau;
    __device__ explicit WideLayout(const PivotItem<T>& it) {
        nblk = (it.L + kWideThreads - 1) / kWideThreads;
        part_val = it.norms + 2 * static_cast<long long>(it.L);
        part_idx = reinterpret_cast<int*>(part_val + nblk);
        tau = reinterpret_cast<T*>(part_val + 2 * nblk);
    }
};

// the largest (val, first idx) of the block into slot
__device__ void block_argmax(double val, int idx, double* s_val, int* s_idx, double* out_val, int* out_idx) {
    s_val[threadIdx.x] = val;
    s_idx[threadIdx.x] = idx;
    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (threadIdx.x < stride) {
            const double v2 = s_val[threadIdx.x + stride];
            const int i2 = s_idx[threadIdx.x + stride];
            if (v2 > s_val[threadIdx.x] || (v2 == s_val[threadIdx.x] && i2 < s_idx[threadIdx.x])) {
                s_val[threadIdx.x] = v2;
                s_idx[threadIdx.x] = i2;
            }
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        *out_val = s_val[0];
        *out_idx = s_idx[0];
    }
}

template<typename T>
__global__ void wide_init_kernel(const PivotItem<T>* items, const char* meta) {
    const PivotItem<T> it = items[blockIdx.x];
    const int l = blockIdx.y * blockDim.x + threadIdx.x;
    if (l < it.L) it.perm[l] = l;
}

template<typename T>
__global__ void wide_mask_kernel(const PivotItem<T>* items, const char* meta) {
    const PivotItem<T> it = items[blockIdx.x];
    const int* mask = it.mask_ptr != nullptr ? it.mask_ptr : reinterpret_cast<const int*>(meta + it.mask_offset);
    for (int m = threadIdx.x; m < it.nmask; m += blockDim.x) it.perm[mask[m]] = -1;
}

template<typename T>
__global__ void __launch_bounds__(kWideThreads) wide_copy_kernel(const PivotItem<T>* items) {
    const PivotItem<T> it = items[blockIdx.x];
    const WideLayout<T> lay(it);
    if (static_cast<int>(blockIdx.y) >= lay.nblk) return;
    __shared__ double s_val[kWideThreads];
    __shared__ int s_idx[kWideThreads];
    const int r = it.r, L = it.L;
    const int l = blockIdx.y * kWideThreads + threadIdx.x;
    double val = -1.0;
    int idx = L;
    if (l < L) {
        const bool zero = it.perm[l] < 0;
        double s = 0.0;
        for (int i = 0; i < r; ++i) {
            const T v = zero ? T(0.0) : (it.transposed ? it.a[l + static_cast<long long>(i) * it.lda]
                                                       : it.a[i + static_cast<long long>(l) * it.lda]);
            it.work[l + static_cast<long long>(i) * L] = v;
            s += sq_abs(v);
        }
        it.norms[l] = sqrt(s);
        it.norms[L + l] = it.norms[l];
        it.perm[l] = l;
        val = it.norms[l];
        idx = l;
    }
    block_argmax(val, idx, s_val, s_idx, lay.part_val + blockIdx.y, lay.part_idx + blockIdx.y);
}

template<typename T>
__global__ void __launch_bounds__(1024) wide_pivot_kernel(const PivotItem<T>* items, int k) {
    const PivotItem<T> it = items[blockIdx.x];
    const int r = it.r, L = it.L;
    if (k >= min(it.p, min(r, L))) return;
    const WideLayout<T> lay(it);
    T* w = it.work;
    auto at = [&](int i, int l) -> T& { return w[l + static_cast<long long>(i) * L]; };
    double* vn1 = it.norms;
    double* vn2 = it.norms + L;
    __shared__ double s_val[1024];
    __shared__ int s_idx[1024];
    __shared__ double r_val;
    __shared__ int r_idx;
    double best = -1.0;
    int bi = L;
    for (int b = threadIdx.x; b < lay.nblk; b += blockDim.x) {
        const double v = lay.part_val[b];
        const int i = lay.part_idx[b];
        if (v > best || (v == best && i < bi)) {
            best = v;
            bi = i;
        }
    }
    block_argmax(best, bi, s_val, s_idx, &r_val, &r_idx);
    __syncthreads();
    const int pvt = r_idx;
    if (pvt != k) {
        for (int i = threadIdx.x; i < r; i += blockDim.x) {
            const T t = at(i, k);
            at(i, k) = at(i, pvt);
            at(i, pvt) = t;
        }
        if (threadIdx.x == 0) {
            const int tp = it.perm[k];
            it.perm[k] = it.perm[pvt];
            it.perm[pvt] = tp;
            vn1[pvt] = vn1[k];
            vn2[pvt] = vn2[k];
        }
    }
    __syncthreads();
    if (threadIdx.x == 0) {
        double xnorm2 = 0.0;
        for (int i = k + 1; i < r; ++i) xnorm2 += sq_abs(at(i, k));
        const T alpha = at(k, k);
        const double ar = real_of(alpha), ai = imag_of(alpha);
        if (xnorm2 == 0.0 && ai == 0.0) {
            *lay.tau = T(0.0);
        } else {
            const double nrm = sqrt(ar * ar + ai * ai + xnorm2);
            const double beta = ar >= 0.0 ? -nrm : nrm;
            T tau;
            if constexpr (fmm::gpu::is_complex_scalar<T>) {
                tau = T((beta - ar) / beta, -ai / beta);
            } else {
                tau = T((beta - ar) / beta);
            }
            const T scale = T(1.0) / (alpha - T(beta));
            *lay.tau = tau;
            for (int i = k + 1; i < r; ++i) at(i, k) = at(i, k) * scale;
            at(k, k) = T(beta);
        }
    }
}

template<typename T>
__global__ void __launch_bounds__(kWideThreads) wide_apply_kernel(const PivotItem<T>* items, int k) {
    const PivotItem<T> it = items[blockIdx.x];
    const int r = it.r, L = it.L;
    if (k >= min(it.p, min(r, L))) return;
    const WideLayout<T> lay(it);
    if (static_cast<int>(blockIdx.y) >= lay.nblk) return;
    extern __shared__ unsigned char s_raw[];
    T* sv = reinterpret_cast<T*>(s_raw);  // the Householder vector, rows k+1 .. r-1
    __shared__ double s_val[kWideThreads];
    __shared__ int s_idx[kWideThreads];
    T* w = it.work;
    auto at = [&](int i, int l) -> T& { return w[l + static_cast<long long>(i) * L]; };
    double* vn1 = it.norms;
    double* vn2 = it.norms + L;
    const double tol3z = sqrt(1.1102230246251565e-16);
    for (int i = k + 1 + threadIdx.x; i < r; i += blockDim.x) sv[i] = at(i, k);
    const T tau = *lay.tau;
    __syncthreads();
    const T ctau = conj_of(tau);
    const int l = blockIdx.y * kWideThreads + threadIdx.x;
    double val = -1.0;
    int idx = L;
    if (l > k && l < L) {
        if (tau != T(0.0)) {
            T dot = at(k, l);
            for (int i = k + 1; i < r; ++i) dot += conj_of(sv[i]) * at(i, l);
            const T f = ctau * dot;
            at(k, l) -= f;
            for (int i = k + 1; i < r; ++i) at(i, l) -= f * sv[i];
        }
        if (vn1[l] != 0.0) {
            double temp = sqrt(sq_abs(at(k, l))) / vn1[l];
            temp = 1.0 - temp * temp;
            temp = temp > 0.0 ? temp : 0.0;
            const double ratio = vn1[l] / vn2[l];
            const double temp2 = temp * ratio * ratio;
            if (temp2 <= tol3z) {
                double s = 0.0;
                for (int i = k + 1; i < r; ++i) s += sq_abs(at(i, l));
                vn1[l] = sqrt(s);
                vn2[l] = vn1[l];
            } else {
                vn1[l] *= sqrt(temp);
            }
        }
        val = vn1[l];
        idx = l;
    }
    block_argmax(val, idx, s_val, s_idx, lay.part_val + blockIdx.y, lay.part_idx + blockIdx.y);
}

template<typename T>
__global__ void wide_out_kernel(const PivotItem<T>* items) {
    const PivotItem<T> it = items[blockIdx.x];
    const int k = blockIdx.y * blockDim.x + threadIdx.x;
    if (k < min(it.p, it.L)) it.out[k] = it.perm[k];
}

template<typename T>
__global__ void axpby_kernel(const AxpbyItem<T>* items) {
    const AxpbyItem<T> it = items[blockIdx.x];
    const long long total = static_cast<long long>(it.m) * it.n;
    for (long long e = static_cast<long long>(blockIdx.y) * blockDim.x + threadIdx.x; e < total;
         e += static_cast<long long>(gridDim.y) * blockDim.x) {
        const int i = static_cast<int>(e % it.m);
        const long long j = e / it.m;
        T v = it.alpha * it.a[i + j * it.lda];
        if (it.b != nullptr) v = v + it.beta * it.b[i + j * it.ldb];
        it.c[i + j * it.ldc] = v;
    }
}

template<typename T>
__global__ void unit_lower_kernel(const UnitLowerItem<T>* items) {
    const UnitLowerItem<T> it = items[blockIdx.x];
    const int j = blockIdx.y;
    if (j >= it.n) return;
    T* col = it.a + static_cast<long long>(j) * it.lda;
    for (int i = threadIdx.x; i <= j; i += blockDim.x) col[i] = i == j ? T(1.0) : T(0.0);
}

template<typename T>
__global__ void add_diagonal_kernel(const DiagItem<T>* items) {
    const DiagItem<T> it = items[blockIdx.x];
    const int i = blockIdx.y * blockDim.x + threadIdx.x;
    if (i < it.n) it.a[i + static_cast<long long>(i) * it.lda] += it.value;
}

template<typename T>
__global__ void symmetrize_kernel(const SquareItem<T>* items) {
    const SquareItem<T> it = items[blockIdx.x];
    const long long total = static_cast<long long>(it.n) * it.n;
    for (long long e = static_cast<long long>(blockIdx.y) * blockDim.x + threadIdx.x; e < total;
         e += static_cast<long long>(gridDim.y) * blockDim.x) {
        const int i = static_cast<int>(e % it.n);
        const int j = static_cast<int>(e / it.n);
        if (i >= j) continue;  // each pair once, by its upper entry
        T& upper = it.a[i + static_cast<long long>(j) * it.lda];
        T& lower = it.a[j + static_cast<long long>(i) * it.lda];
        const T v = (upper + lower) * 0.5;
        upper = v;
        lower = v;
    }
}

template<typename T>
__global__ void copy_diagonal_kernel(const DiagCopyItem<T>* items) {
    const DiagCopyItem<T> it = items[blockIdx.x];
    const int i = blockIdx.y * blockDim.x + threadIdx.x;
    if (i < it.n) it.out[i] = it.a[i + static_cast<long long>(i) * it.lda];
}

template<typename T>
__global__ void fnorm_kernel(const NormItem<T>* items) {
    __shared__ double partial[kThreads];
    const NormItem<T> it = items[blockIdx.x];
    const long long total = static_cast<long long>(it.m) * it.n;
    double s = 0.0;
    for (long long e = threadIdx.x; e < total; e += blockDim.x) {
        const int i = static_cast<int>(e % it.m);
        const long long j = e / it.m;
        s += fmm::gpu::abs2(it.a[i + j * it.lda]);
    }
    partial[threadIdx.x] = s;
    __syncthreads();
    for (int stride = kThreads / 2; stride > 0; stride >>= 1) {
        if (threadIdx.x < stride) partial[threadIdx.x] += partial[threadIdx.x + stride];
        __syncthreads();
    }
    if (threadIdx.x == 0) *it.norm = sqrt(partial[0]);
}

}  // namespace

template<typename T>
void launch_row_swaps(const RowSwapItem<T>* items, int count, int max_n, cudaStream_t stream) {
    if (count <= 0 || max_n <= 0) return;
    const dim3 grid(static_cast<unsigned>(count), static_cast<unsigned>((max_n + kThreads - 1) / kThreads));
    row_swaps_kernel<T><<<grid, kThreads, 0, stream>>>(items);
    check_launch("row_swaps_kernel");
}

template<typename T>
void launch_qrcp_pivots(const PivotItem<T>* items, int count, const char* meta, cudaStream_t stream) {
    if (count <= 0) return;
    qrcp_pivots_kernel<T><<<static_cast<unsigned>(count), kPivotThreads, 0, stream>>>(items, meta);
    check_launch("qrcp_pivots_kernel");
}

template<typename T>
void launch_qrcp_pivots_wide(const PivotItem<T>* items, int count, int max_r, int max_L, int max_steps, int max_p,
                             const char* meta, cudaStream_t stream) {
    if (count <= 0 || max_L <= 0) return;
    const unsigned nblk = static_cast<unsigned>((max_L + kWideThreads - 1) / kWideThreads);
    const dim3 cols(static_cast<unsigned>(count), nblk);
    wide_init_kernel<T><<<cols, kWideThreads, 0, stream>>>(items, meta);
    wide_mask_kernel<T><<<static_cast<unsigned>(count), kWideThreads, 0, stream>>>(items, meta);
    wide_copy_kernel<T><<<cols, kWideThreads, 0, stream>>>(items);
    check_launch("wide pivots (setup)");
    const size_t smem = static_cast<size_t>(std::max(max_r, 1)) * sizeof(T);
    for (int k = 0; k < max_steps; ++k) {
        wide_pivot_kernel<T><<<static_cast<unsigned>(count), 1024, 0, stream>>>(items, k);
        wide_apply_kernel<T><<<cols, kWideThreads, smem, stream>>>(items, k);
    }
    check_launch("wide pivots (steps)");
    wide_out_kernel<T><<<dim3(static_cast<unsigned>(count), static_cast<unsigned>((std::max(max_p, 1) + kWideThreads - 1) /
                                                                                    kWideThreads)),
                         kWideThreads, 0, stream>>>(items);
    check_launch("wide pivots (output)");
}

template<typename T>
void launch_unit_lower(const UnitLowerItem<T>* items, int count, int max_n, cudaStream_t stream) {
    if (count <= 0 || max_n <= 0) return;
    unit_lower_kernel<T><<<dim3(static_cast<unsigned>(count), static_cast<unsigned>(max_n)), kThreads, 0, stream>>>(items);
    check_launch("unit_lower_kernel");
}

template<typename T>
void launch_axpby(const AxpbyItem<T>* items, int count, int max_m, int max_n, cudaStream_t stream) {
    if (count <= 0 || max_m <= 0 || max_n <= 0) return;
    const dim3 grid(static_cast<unsigned>(count), tiles(static_cast<long long>(max_m) * max_n));
    axpby_kernel<T><<<grid, kThreads, 0, stream>>>(items);
    check_launch("axpby_kernel");
}

template<typename T>
void launch_add_diagonal(const DiagItem<T>* items, int count, int max_n, cudaStream_t stream) {
    if (count <= 0 || max_n <= 0) return;
    const dim3 grid(static_cast<unsigned>(count), static_cast<unsigned>((max_n + kThreads - 1) / kThreads));
    add_diagonal_kernel<T><<<grid, kThreads, 0, stream>>>(items);
    check_launch("add_diagonal_kernel");
}

template<typename T>
void launch_symmetrize(const SquareItem<T>* items, int count, int max_n, cudaStream_t stream) {
    if (count <= 0 || max_n <= 0) return;
    const dim3 grid(static_cast<unsigned>(count), tiles(static_cast<long long>(max_n) * max_n));
    symmetrize_kernel<T><<<grid, kThreads, 0, stream>>>(items);
    check_launch("symmetrize_kernel");
}

template<typename T>
void launch_copy_diagonal(const DiagCopyItem<T>* items, int count, int max_n, cudaStream_t stream) {
    if (count <= 0 || max_n <= 0) return;
    const dim3 grid(static_cast<unsigned>(count), static_cast<unsigned>((max_n + kThreads - 1) / kThreads));
    copy_diagonal_kernel<T><<<grid, kThreads, 0, stream>>>(items);
    check_launch("copy_diagonal_kernel");
}

template<typename T>
void launch_fnorm(const NormItem<T>* items, int count, cudaStream_t stream) {
    if (count <= 0) return;
    fnorm_kernel<T><<<static_cast<unsigned>(count), kThreads, 0, stream>>>(items);
    check_launch("fnorm_kernel");
}

template<typename T>
__global__ void zero_lower_kernel(const TriuItem<T>* items) {
    const TriuItem<T> it = items[blockIdx.x];
    const int j = blockIdx.y;
    if (j >= it.n) return;
    T* col = it.a + static_cast<long long>(j) * it.lda;
    for (int i = j + 1 + threadIdx.x; i < it.n; i += blockDim.x) col[i] = T(0.0);
}

template<typename T>
void launch_zero_lower(const TriuItem<T>* items, int count, int max_n, cudaStream_t stream) {
    if (count <= 0 || max_n <= 0) return;
    zero_lower_kernel<T><<<dim3(static_cast<unsigned>(count), static_cast<unsigned>(max_n)), kThreads, 0, stream>>>(items);
    check_launch("zero_lower_kernel");
}

template<typename T>
__global__ void col_scale_kernel(const ColScaleItem<T>* items) {
    const ColScaleItem<T> it = items[blockIdx.x];
    const long long total = static_cast<long long>(it.m) * it.n;
    for (long long e = static_cast<long long>(blockIdx.y) * blockDim.x + threadIdx.x; e < total;
         e += static_cast<long long>(gridDim.y) * blockDim.x) {
        const int i = static_cast<int>(e % it.m);
        const long long j = e / it.m;
        T v = it.a[i + j * it.lda];
        if (it.conj) v = fmm::gpu::conj(v);
        if (it.s != nullptr) v = v * it.s[j];
        it.c[i + j * it.ldc] = v;
    }
}

template<typename T>
void launch_col_scale(const ColScaleItem<T>* items, int count, int max_m, int max_n, cudaStream_t stream) {
    if (count <= 0 || max_m <= 0 || max_n <= 0) return;
    const dim3 grid(static_cast<unsigned>(count), tiles(static_cast<long long>(max_m) * max_n));
    col_scale_kernel<T><<<grid, kThreads, 0, stream>>>(items);
    check_launch("col_scale_kernel");
}

template<typename T>
__global__ void tinv_kernel(const TinvItem<T>* items) {
    const TinvItem<T> it = items[blockIdx.x];
    const int j = blockIdx.y;
    if (j >= it.k) return;
    T* col = it.t + static_cast<long long>(j) * it.k;
    const T* g = it.gram + static_cast<long long>(j) * it.ldg;
    for (int i = threadIdx.x; i < it.k; i += blockDim.x) {
        col[i] = i < j ? g[i] : (i == j ? T(1.0) / it.tau[j] : T(0.0));
    }
}

template<typename T>
void launch_tinv(const TinvItem<T>* items, int count, int max_k, cudaStream_t stream) {
    if (count <= 0 || max_k <= 0) return;
    tinv_kernel<T><<<dim3(static_cast<unsigned>(count), static_cast<unsigned>(max_k)), kThreads, 0, stream>>>(items);
    check_launch("tinv_kernel");
}

template<typename T>
__global__ void gather_rows_kernel(const T* a, int lda, int k, const int* rows, int n, T* out) {
    const long long total = static_cast<long long>(n) * k;
    for (long long e = static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x; e < total;
         e += static_cast<long long>(gridDim.x) * blockDim.x) {
        const int t = static_cast<int>(e % n);
        const long long c = e / n;
        out[e] = a[rows[t] + c * lda];
    }
}

template<typename T>
void launch_gather_rows(const T* a, int lda, int k, const int* rows, int n, T* out, cudaStream_t stream) {
    if (n <= 0 || k <= 0) return;
    gather_rows_kernel<T><<<tiles(static_cast<long long>(n) * k), kThreads, 0, stream>>>(a, lda, k, rows, n, out);
    check_launch("gather_rows_kernel");
}

#define BPACK_HODLR_KERNELS(T)                                                                          \
    template void launch_row_swaps<T>(const RowSwapItem<T>*, int, int, cudaStream_t);                  \
    template void launch_qrcp_pivots<T>(const PivotItem<T>*, int, const char*, cudaStream_t);           \
    template void launch_qrcp_pivots_wide<T>(const PivotItem<T>*, int, int, int, int, int, const char*,      \
                                             cudaStream_t);                                                \
    template void launch_axpby<T>(const AxpbyItem<T>*, int, int, int, cudaStream_t);                   \
    template void launch_unit_lower<T>(const UnitLowerItem<T>*, int, int, cudaStream_t);               \
    template void launch_add_diagonal<T>(const DiagItem<T>*, int, int, cudaStream_t);                  \
    template void launch_symmetrize<T>(const SquareItem<T>*, int, int, cudaStream_t);                  \
    template void launch_copy_diagonal<T>(const DiagCopyItem<T>*, int, int, cudaStream_t);             \
    template void launch_fnorm<T>(const NormItem<T>*, int, cudaStream_t);                             \
    template void launch_zero_lower<T>(const TriuItem<T>*, int, int, cudaStream_t);                   \
    template void launch_col_scale<T>(const ColScaleItem<T>*, int, int, int, cudaStream_t);           \
    template void launch_tinv<T>(const TinvItem<T>*, int, int, cudaStream_t);                          \
    template void launch_gather_rows<T>(const T*, int, int, const int*, int, T*, cudaStream_t);
BPACK_HODLR_KERNELS(double)
BPACK_HODLR_KERNELS(dcomplex)
#undef BPACK_HODLR_KERNELS

}  // namespace gpu
}  // namespace bpack
