// Batched interpolative decomposition of the H2 Color GPU backend: column-
// pivoted Householder QR with the host's rank rule, then T = R11^{-1} R12.
// Real (dgeqp3) or complex (zgeqp3) data.  See QrcpItemT in device_kernels.hpp.
//
// Every launch variant (one block of 256 or 1024 threads per box, or several
// blocks per box) gives a box the same pivots, rank, R, T and traced norm,
// bit for bit, whatever its batch: each column is worked on by one warp with
// the same code, block-wide results are exact (maxima) or summed in column
// order, and this file is compiled without FMA contraction (--fmad=false,
// CMakeLists.txt), which could otherwise differ between the template
// instances.  Replicated CA levels rely on it (tests/batch_determinism.cpp).

#include "device_kernels.hpp"

#include <cooperative_groups.h>

#include <algorithm>
#include <climits>
#include <stdexcept>
#include <string>

namespace fmm {
namespace gpu {

namespace {

constexpr int kMaxWarps = 32;

__device__ __forceinline__ double warp_sum(double v) {
    for (int o = 16; o > 0; o >>= 1) v += __shfl_xor_sync(0xffffffffu, v, o);
    return v;
}
__device__ __forceinline__ dcomplex warp_sum(dcomplex v) { return dcomplex(warp_sum(v.re), warp_sum(v.im)); }

template<int kWarps>
__device__ double block_max(double v, double* scratch) {
    for (int o = 16; o > 0; o >>= 1) v = fmax(v, __shfl_xor_sync(0xffffffffu, v, o));
    __syncthreads();
    if ((threadIdx.x & 31) == 0) scratch[threadIdx.x >> 5] = v;
    __syncthreads();
    double best = 0.0;
    for (int w = 0; w < kWarps; ++w) best = fmax(best, scratch[w]);
    return best;
}

// First index of the largest value (LAPACK idamax on nonnegative values).
template<int kWarps>
__device__ int block_argmax(double v, int idx, double* vscratch, int* iscratch) {
    for (int o = 16; o > 0; o >>= 1) {
        const double ov = __shfl_xor_sync(0xffffffffu, v, o);
        const int oi = __shfl_xor_sync(0xffffffffu, idx, o);
        if (ov > v || (ov == v && oi < idx)) {
            v = ov;
            idx = oi;
        }
    }
    __syncthreads();
    if ((threadIdx.x & 31) == 0) {
        vscratch[threadIdx.x >> 5] = v;
        iscratch[threadIdx.x >> 5] = idx;
    }
    __syncthreads();
    double best = vscratch[0];
    int best_idx = iscratch[0];
    for (int w = 1; w < kWarps; ++w) {
        if (vscratch[w] > best || (vscratch[w] == best && iscratch[w] < best_idx)) {
            best = vscratch[w];
            best_idx = iscratch[w];
        }
    }
    return best_idx;
}

// Loads of data other blocks write between barriers bypass L1 (it is not
// coherent across multiprocessors).
template<typename T>
__device__ __forceinline__ T ldg2(const T* p) { return __ldcg(p); }
__device__ __forceinline__ dcomplex ldg2(const dcomplex* p) {
    const double2 v = __ldcg(reinterpret_cast<const double2*>(p));
    return dcomplex(v.x, v.y);
}
template<bool kUncached, typename T>
__device__ __forceinline__ T load(const T* p) {
    if constexpr (kUncached) {
        return ldg2(p);
    } else {
        return *p;
    }
}

// The loops over a column below load kAhead rows a lane ahead of their sums
// and stores, which keep the order of the one-row-at-a-time loop and so its
// result bit for bit: a warp otherwise waited on one load at a time (the
// cooperative QRCP of level 2 at 64^3 stalled 60% of its cycles there, at
// 37% of the DRAM bandwidth; with 4 rows ahead 52%, its launches 27% faster,
// and 8 gained nothing more).
constexpr int kAhead = 4;

// x[u] = col[r + 32 u], u < kAhead
template<bool kUncached, typename T>
__device__ __forceinline__ void load_ahead(const T* col, int r, T (&x)[kAhead]) {
#pragma unroll
    for (int u = 0; u < kAhead; ++u) x[u] = load<kUncached>(col + r + 32 * u);
}

// Householder reflector of col[i:m) by one warp (LAPACK dlarfg): col[i] =
// beta, col[i+1:m) = v; returns tau and beta.
template<bool kUncached>
__device__ __forceinline__ void make_reflector(double* col, int i, int m, int lane, double& tau, double& beta) {
    const double alpha = load<kUncached>(col + i);
    double s = 0.0;
    int r = i + 1 + lane;
    for (; r + 32 * (kAhead - 1) < m; r += 32 * kAhead) {
        double x[kAhead];
        load_ahead<kUncached>(col, r, x);
#pragma unroll
        for (int u = 0; u < kAhead; ++u) s += x[u] * x[u];
    }
    for (; r < m; r += 32) {
        const double x = load<kUncached>(col + r);
        s += x * x;
    }
    const double xnorm = sqrt(warp_sum(s));
    tau = 0.0;
    beta = alpha;
    if (xnorm != 0.0) {
        beta = -copysign(sqrt(alpha * alpha + xnorm * xnorm), alpha);
        tau = (beta - alpha) / beta;
        const double scal = 1.0 / (alpha - beta);
        r = i + 1 + lane;
        for (; r + 32 * (kAhead - 1) < m; r += 32 * kAhead) {
            double x[kAhead];
            load_ahead<kUncached>(col, r, x);
#pragma unroll
            for (int u = 0; u < kAhead; ++u) col[r + 32 * u] = x[u] * scal;
        }
        for (; r < m; r += 32) col[r] = load<kUncached>(col + r) * scal;
        __syncwarp();
        if (lane == 0) col[i] = beta;
    }
}
// Complex (LAPACK zlarfg): beta is real, so the reflector also acts on a
// complex alpha with x = 0 (the last row), and the diagonal of R is real.
template<bool kUncached>
__device__ __forceinline__ void make_reflector(dcomplex* col, int i, int m, int lane, dcomplex& tau, double& beta) {
    const dcomplex alpha = load<kUncached>(col + i);
    double s = 0.0;
    int r = i + 1 + lane;
    for (; r + 32 * (kAhead - 1) < m; r += 32 * kAhead) {
        dcomplex x[kAhead];
        load_ahead<kUncached>(col, r, x);
#pragma unroll
        for (int u = 0; u < kAhead; ++u) s += abs2(x[u]);
    }
    for (; r < m; r += 32) s += abs2(load<kUncached>(col + r));
    const double xnorm = sqrt(warp_sum(s));
    tau = dcomplex(0.0);
    beta = alpha.re;
    if (xnorm != 0.0 || alpha.im != 0.0) {
        beta = -copysign(norm3d(alpha.re, alpha.im, xnorm), alpha.re);  // dlapy3
        tau = dcomplex((beta - alpha.re) / beta, -alpha.im / beta);
        const dcomplex scal = dcomplex(1.0) / (alpha - dcomplex(beta));
        r = i + 1 + lane;
        for (; r + 32 * (kAhead - 1) < m; r += 32 * kAhead) {
            dcomplex x[kAhead];
            load_ahead<kUncached>(col, r, x);
#pragma unroll
            for (int u = 0; u < kAhead; ++u) col[r + 32 * u] = x[u] * scal;
        }
        for (; r < m; r += 32) col[r] = load<kUncached>(col + r) * scal;
        __syncwarp();
        if (lane == 0) col[i] = dcomplex(beta);
    }
}

// H^H = I - conj(tau) v v^H applied to col by one warp (v_i = 1 implicit),
// v's rows from vat(r).
template<bool kUncached, typename T, typename V>
__device__ __forceinline__ void apply_reflector_with(V&& vat, T tau, T* col, int i, int m, int lane) {
    T s = lane == 0 ? load<kUncached>(col + i) : T(0.0);
    int r = i + 1 + lane;
    for (; r + 32 * (kAhead - 1) < m; r += 32 * kAhead) {
        T c[kAhead], v[kAhead];
        load_ahead<kUncached>(col, r, c);
#pragma unroll
        for (int u = 0; u < kAhead; ++u) v[u] = vat(r + 32 * u);
#pragma unroll
        for (int u = 0; u < kAhead; ++u) s += conj(v[u]) * c[u];
    }
    for (; r < m; r += 32) s += conj(vat(r)) * load<kUncached>(col + r);
    const T w = conj(tau) * warp_sum(s);
    if (w != T(0.0)) {
        if (lane == 0) col[i] = load<kUncached>(col + i) - w;
        r = i + 1 + lane;
        for (; r + 32 * (kAhead - 1) < m; r += 32 * kAhead) {
            T c[kAhead], v[kAhead];
            load_ahead<kUncached>(col, r, c);
#pragma unroll
            for (int u = 0; u < kAhead; ++u) v[u] = vat(r + 32 * u);
#pragma unroll
            for (int u = 0; u < kAhead; ++u) col[r + 32 * u] = c[u] - w * v[u];
        }
        for (; r < m; r += 32) col[r] = load<kUncached>(col + r) - w * vat(r);
    }
}
template<bool kUncached, typename T>
__device__ __forceinline__ void apply_reflector(const T* v, T tau, T* col, int i, int m, int lane) {
    apply_reflector_with<kUncached>([&](int r) { return load<kUncached>(v + r); }, tau, col, i, m, lane);
}
// ... with v's rows i+1 .. i+staged in v_sh (the block's copy, the rest
// from v): the same values, so the same result bit for bit.
template<typename T>
__device__ __forceinline__ void apply_reflector_staged(const T* v_sh, int staged, const T* v, T tau, T* col, int i, int m,
                                                       int lane) {
    const int split = i + 1 + staged;
    apply_reflector_with<true>([&](int r) { return r < split ? v_sh[r - i - 1] : ldg2(v + r); }, tau, col, i, m, lane);
}

// LAPACK dlaqp2's partial norm downdate of a trailing column after step i.
template<bool kUncached, typename T>
__device__ __forceinline__ void downdate_norm(const T* col, int i, int m, int lane, double tol3z, double* vn1j,
                                              double* vn2j) {
    const double n1 = load<kUncached>(vn1j);
    if (n1 != 0.0) {
        double temp = magnitude(load<kUncached>(col + i)) / n1;
        temp = fmax(1.0 - temp * temp, 0.0);
        const double ratio = n1 / load<kUncached>(vn2j);
        if (temp * ratio * ratio <= tol3z) {
            double s2 = 0.0;
            for (int r = i + 1 + lane; r < m; r += 32) s2 += abs2(load<kUncached>(col + r));
            const double fresh = i < m - 1 ? sqrt(warp_sum(s2)) : 0.0;
            if (lane == 0) *vn1j = *vn2j = fresh;
        } else if (lane == 0) {
            *vn1j = n1 * sqrt(temp);
        }
    }
}

// Column j of T = R11^{-1} R12 in place by one warp (no conjugates: ztrtrs
// 'N'); returns false for a non-finite entry.
template<bool kUncached, typename T>
__device__ __forceinline__ bool solve_r11(const T* A, int lda, int rank, T* col, int lane) {
    bool ok = true;
    for (int p = rank - 1; p >= 0; --p) {
        T s = T(0.0);
        for (int q = p + 1 + lane; q < rank; q += 32) {
            s += load<kUncached>(A + p + static_cast<int64_t>(q) * lda) * load<kUncached>(col + q);
        }
        s = warp_sum(s);
        const T x = (load<kUncached>(col + p) - s) / load<kUncached>(A + p + static_cast<int64_t>(p) * lda);
        __syncwarp();
        if (lane == 0) col[p] = x;
        __syncwarp();
        if (!is_finite(x)) ok = false;
    }
    return ok;
}

// One matrix per block.  Mirrors fmm::compute_id_complex: normalize by the
// largest magnitude, factor A P = Q R with column pivoting (LAPACK dlaqp2's
// pivot choice and norm downdating), stop at the first
// |R_ii| <= max(tol |R_00|, abs_tol / amax) (abs_tol, the H2_s_op floor, in
// the input's units; 0: off), and solve R11 T = R12 in place.  The factorization stops at the rank, so
// the redundant columns keep the order of that step (the host's full
// factorization permutes them further; the skeleton is the same).
template<typename T, int kThreads>
__global__ void __launch_bounds__(kThreads) qrcp_kernel(const QrcpItemT<T>* items, double tol, double abs_tol) {
    constexpr int kWarps = kThreads / 32;
    extern __shared__ double smem[];
    const QrcpItemT<T> item = items[blockIdx.x];
    const int m = item.m, n = item.n, lda = item.lda;
    if (item.tol >= 0.0) tol = item.tol;
    T* A = item.a;
    double* vn1 = smem;
    double* vn2 = smem + n;
    double* scratch = smem + 2 * n;                             // kMaxWarps doubles
    int* iscratch = reinterpret_cast<int*>(scratch + kMaxWarps);  // kMaxWarps ints
    int* jpvt = iscratch + kMaxWarps;                           // n ints
    __shared__ T s_tau;
    __shared__ double s_beta;
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5;

    // ---- norm (for traces), finiteness, normalization.  The norm sums the
    // columns' sums in column order, as the cooperative kernel does: every
    // launch variant traces the same value.
    double amax = 0.0, bad = 0.0;
    for (int j = warp; j < n; j += kWarps) {
        const T* col = A + static_cast<int64_t>(j) * lda;
        double s = 0.0;
        for (int r = lane; r < m; r += 32) {
            const T x = col[r];
            if (!is_finite(x)) bad = 1.0;
            amax = fmax(amax, magnitude(x));
            s += abs2(x);
        }
        s = warp_sum(s);
        if (lane == 0) vn2[j] = s;
    }
    bad = block_max<kWarps>(bad, scratch);
    amax = block_max<kWarps>(amax, scratch);
    if (tid == 0) {
        double ssq = 0.0;
        for (int j = 0; j < n; ++j) ssq += vn2[j];
        *item.norm = sqrt(ssq);
        *item.flag = bad != 0.0 ? 1 : 0;
    }
    for (int j = tid; j < n; j += kThreads) jpvt[j] = j;
    if (bad != 0.0 || amax == 0.0) {
        // zero target: rank 0 (the host keeps column 0 with a zero row of T)
        __syncthreads();
        for (int j = tid; j < n; j += kThreads) {
            item.jpvt[j] = j;
            if (j >= 1 && m > 0) A[static_cast<int64_t>(j) * lda] = T(0.0);
        }
        if (tid == 0) *item.rank = 0;
        return;
    }
    for (int j = warp; j < n; j += kWarps) {
        T* col = A + static_cast<int64_t>(j) * lda;
        for (int r = lane; r < m; r += 32) col[r] /= amax;
    }
    const double abs_floor = abs_tol / amax;  // (compute_id_complex: abs_tolerance / scale)
    __syncthreads();

    // ---- initial column norms
    for (int j = warp; j < n; j += kWarps) {
        const T* col = A + static_cast<int64_t>(j) * lda;
        double s = 0.0;
        for (int r = lane; r < m; r += 32) s += abs2(col[r]);
        s = warp_sum(s);
        if (lane == 0) vn1[j] = vn2[j] = sqrt(s);
    }
    __syncthreads();

    const double tol3z = sqrt(1.1102230246251565e-16);  // sqrt(dlamch('Epsilon'))
    const int kmax = m < n ? m : n;
    int rank = 0;
    double r00 = 0.0;
    for (int i = 0; i < kmax; ++i) {
        // pivot: largest remaining column norm, first index on ties
        double best = -1.0;
        int best_idx = 0x7fffffff;
        for (int j = i + tid; j < n; j += kThreads) {
            if (vn1[j] > best) {
                best = vn1[j];
                best_idx = j;
            }
        }
        const int pvt = block_argmax<kWarps>(best, best_idx, scratch, iscratch);
        if (pvt != i) {
            T* ci = A + static_cast<int64_t>(i) * lda;
            T* cp = A + static_cast<int64_t>(pvt) * lda;
            for (int r = tid; r < m; r += kThreads) {
                const T t = ci[r];
                ci[r] = cp[r];
                cp[r] = t;
            }
            if (tid == 0) {
                const int t = jpvt[i];
                jpvt[i] = jpvt[pvt];
                jpvt[pvt] = t;
                vn1[pvt] = vn1[i];
                vn2[pvt] = vn2[i];
            }
        }
        __syncthreads();

        // Householder reflector of column i (LAPACK dlarfg / zlarfg)
        if (warp == 0) {
            T tau;
            double beta;
            make_reflector<false>(A + static_cast<int64_t>(i) * lda, i, m, lane, tau, beta);
            if (lane == 0) {
                s_tau = tau;
                s_beta = beta;
                if (item.tau != nullptr) item.tau[i] = tau;
            }
        }
        __syncthreads();
        const T tau = s_tau;
        const double beta = s_beta;
        if (i == 0) r00 = fabs(beta);
        if (!(fabs(beta) > fmax(tol * r00, abs_floor))) break;  // rank rule of compute_id_complex
        rank = i + 1;

        // apply the reflector to the trailing columns and downdate their norms
        const T* v = A + static_cast<int64_t>(i) * lda;
        for (int j = i + 1 + warp; j < n; j += kWarps) {
            T* col = A + static_cast<int64_t>(j) * lda;
            apply_reflector<false>(v, tau, col, i, m, lane);
            __syncwarp();
            downdate_norm<false>(col, i, m, lane, tol3z, vn1 + j, vn2 + j);
            __syncwarp();
        }
        __syncthreads();
    }

    // ---- T = R11^{-1} R12, one warp per column, in place
    double t_bad = 0.0;
    if (rank == 0) {
        for (int j = 1 + tid; j < n; j += kThreads) A[static_cast<int64_t>(j) * lda] = T(0.0);
    } else {
        for (int j = rank + warp; j < n; j += kWarps) {
            if (!solve_r11<false>(A, lda, rank, A + static_cast<int64_t>(j) * lda, lane)) t_bad = 1.0;
        }
    }
    t_bad = block_max<kWarps>(t_bad, scratch);
    // rank 0 keeps column 0 (the host's convention), so report the natural order
    for (int j = tid; j < n; j += kThreads) item.jpvt[j] = rank == 0 ? j : jpvt[j];
    if (tid == 0) {
        *item.rank = rank;
        if (t_bad != 0.0) *item.flag = 2;
    }
}

// ---------------------------------------------------------------------------
// Few boxes: each box on `cpb` blocks of one cooperative launch.  Block
// `part` owns the column positions part, part + cpb, ...; the matrix stays in
// global memory and a grid barrier follows each pivot search and each
// reflector.  Every column is processed by one warp exactly as in
// qrcp_kernel, so pivots, R, T and the traced norm are the same.
// ---------------------------------------------------------------------------
constexpr int kCoopThreads = 256;

__host__ __device__ inline size_t coop_box_doubles(int max_n, int cpb) {
    // vn1, vn2, cand, partials (amax, bad), tau (2)/beta/r00, the columns' sums (traced norm)
    return 3 * static_cast<size_t>(max_n) + 3 * static_cast<size_t>(cpb) + 4;
}

// tau in box_d[0] (real) or box_d[0..1] (complex)
__device__ __forceinline__ void put_tau(double* p, double v) { p[0] = v; }
__device__ __forceinline__ void put_tau(double* p, dcomplex v) {
    p[0] = v.re;
    p[1] = v.im;
}
template<typename T>
__device__ __forceinline__ T get_tau(const double* p) {
    return make_scalar<T>(ldg2(p), is_complex_scalar<T> ? ldg2(p + 1) : 0.0);
}
__host__ __device__ inline size_t coop_box_ints(int cpb) { return static_cast<size_t>(cpb) + 2; }  // cand, stopped, rank

// smallest t with part + cpb * t >= from
__device__ __forceinline__ int first_owned(int from, int part, int cpb) {
    return from <= part ? 0 : (from - part + cpb - 1) / cpb;
}

// largest value, first index on ties, over the block; all threads get both
template<int kWarps>
__device__ void block_argmax_pair(double& v, int& idx, double* vscratch, int* iscratch) {
    for (int o = 16; o > 0; o >>= 1) {
        const double ov = __shfl_xor_sync(0xffffffffu, v, o);
        const int oi = __shfl_xor_sync(0xffffffffu, idx, o);
        if (ov > v || (ov == v && oi < idx)) {
            v = ov;
            idx = oi;
        }
    }
    __syncthreads();
    if ((threadIdx.x & 31) == 0) {
        vscratch[threadIdx.x >> 5] = v;
        iscratch[threadIdx.x >> 5] = idx;
    }
    __syncthreads();
    v = vscratch[0];
    idx = iscratch[0];
    for (int w = 1; w < kWarps; ++w) {
        if (vscratch[w] > v || (vscratch[w] == v && iscratch[w] < idx)) {
            v = vscratch[w];
            idx = iscratch[w];
        }
    }
}

template<typename T>
__global__ void __launch_bounds__(kCoopThreads)
qrcp_coop_kernel(const QrcpItemT<T>* items, int count, int cpb, int max_n, double tol, double abs_tol, double* dwork,
                 int* iwork) {
    namespace cg = cooperative_groups;
    cg::grid_group grid = cg::this_grid();
    constexpr int kWarps = kCoopThreads / 32;
    // the step's reflector, staged once a block: every owned column reads it
    // twice, and otherwise from L2 (the loads bypass L1)
    constexpr int kStaged = (40 << 10) / static_cast<int>(sizeof(T));
    __shared__ double v_raw[(40 << 10) / sizeof(double)];
    T* v_sh = reinterpret_cast<T*>(v_raw);
    __shared__ double scratch[kMaxWarps];
    __shared__ int iscratch[kMaxWarps];
    const int box = blockIdx.x / cpb, part = blockIdx.x % cpb;
    const QrcpItemT<T> item = items[box];
    const int m = item.m, n = item.n, lda = item.lda;
    if (item.tol >= 0.0) tol = item.tol;
    T* A = item.a;
    int* jpvt = item.jpvt;  // the working permutation
    double* vn1 = dwork + static_cast<size_t>(box) * coop_box_doubles(max_n, cpb);
    double* vn2 = vn1 + max_n;
    double* cand_v = vn2 + max_n;
    double* part_amax = cand_v + cpb;
    double* part_bad = part_amax + cpb;
    double* box_d = part_bad + cpb;  // tau (2), beta, r00
    double* col_ssq = box_d + 4;     // the columns' sums of squares (traced norm)
    int* cand_i = iwork + static_cast<size_t>(box) * coop_box_ints(cpb);
    int* box_i = cand_i + cpb;       // stopped, rank
    int* stopped_count = iwork + static_cast<size_t>(count) * coop_box_ints(cpb);
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5;
    const int stride = cpb;  // between owned positions

    // ---- norm, finiteness and largest magnitude of the owned columns
    double amax = 0.0, bad = 0.0;
    for (int t = warp;; t += kWarps) {
        const int j = part + stride * t;
        if (j >= n) break;
        const T* col = A + static_cast<int64_t>(j) * lda;
        double s = 0.0;
        for (int r = lane; r < m; r += 32) {
            const T x = ldg2(col + r);
            if (!is_finite(x)) bad = 1.0;
            amax = fmax(amax, magnitude(x));
            s += abs2(x);
        }
        s = warp_sum(s);
        if (lane == 0) col_ssq[j] = s;
    }
    bad = block_max<kWarps>(bad, scratch);
    amax = block_max<kWarps>(amax, scratch);
    if (tid == 0) {
        part_amax[part] = amax;
        part_bad[part] = bad;
        if (part == 0) {
            box_i[0] = 0;
            box_i[1] = 0;
        }
        if (blockIdx.x == 0) *stopped_count = 0;
    }
    for (int t = tid;; t += kCoopThreads) {
        const int j = part + stride * t;
        if (j >= n) break;
        jpvt[j] = j;
    }
    grid.sync();
    amax = 0.0;
    bad = 0.0;
    for (int c = 0; c < cpb; ++c) {
        amax = fmax(amax, ldg2(part_amax + c));
        bad = fmax(bad, ldg2(part_bad + c));
    }
    const bool zero = bad != 0.0 || amax == 0.0;
    const double abs_floor = zero ? 0.0 : abs_tol / amax;  // as qrcp_kernel
    if (part == 0 && tid == 0) {
        double ssq = 0.0;
        for (int j = 0; j < n; ++j) ssq += ldg2(col_ssq + j);
        *item.norm = sqrt(ssq);
        *item.flag = bad != 0.0 ? 1 : 0;
        if (zero) {
            box_i[0] = 1;  // rank 0
            atomicAdd(stopped_count, 1);
        }
    }
    if (!zero) {
        for (int t = warp;; t += kWarps) {
            const int j = part + stride * t;
            if (j >= n) break;
            T* col = A + static_cast<int64_t>(j) * lda;
            for (int r = lane; r < m; r += 32) col[r] = ldg2(col + r) / amax;
            __syncwarp();
            double s = 0.0;
            for (int r = lane; r < m; r += 32) s += abs2(ldg2(col + r));
            s = warp_sum(s);
            if (lane == 0) vn1[j] = vn2[j] = sqrt(s);
        }
    }
    grid.sync();

    const double tol3z = sqrt(1.1102230246251565e-16);
    const int kmax = m < n ? m : n;
    bool active = !zero;
    for (int i = 0;; ++i) {
        __syncthreads();  // the norms of the previous step's update
        // ---- candidates: largest remaining norm among the owned positions
        if (active && i < kmax) {
            double best = -1.0;
            int best_idx = INT_MAX;
            for (int t = first_owned(i, part, stride) + tid;; t += kCoopThreads) {
                const int j = part + stride * t;
                if (j >= n) break;
                const double nj = ldg2(vn1 + j);
                if (nj > best) {
                    best = nj;
                    best_idx = j;
                }
            }
            block_argmax_pair<kWarps>(best, best_idx, scratch, iscratch);
            if (tid == 0) {
                cand_v[part] = best;
                cand_i[part] = best_idx;
            }
        }
        grid.sync();
        // ---- the owner of position i: pivot, swap, reflector, rank rule
        if (active && part == (i < n ? i % cpb : 0)) {
            if (i >= kmax) {
                if (tid == 0) {
                    box_i[0] = 1;
                    box_i[1] = kmax;
                    atomicAdd(stopped_count, 1);
                }
            } else {
                double best = ldg2(cand_v);
                int pvt = ldg2(cand_i);
                for (int c = 1; c < cpb; ++c) {
                    const double cv = ldg2(cand_v + c);
                    const int ci = ldg2(cand_i + c);
                    if (cv > best || (cv == best && ci < pvt)) {
                        best = cv;
                        pvt = ci;
                    }
                }
                if (pvt != i) {
                    T* ci = A + static_cast<int64_t>(i) * lda;
                    T* cp = A + static_cast<int64_t>(pvt) * lda;
                    for (int r = tid; r < m; r += kCoopThreads) {
                        const T t = ldg2(ci + r);
                        ci[r] = ldg2(cp + r);
                        cp[r] = t;
                    }
                    if (tid == 0) {
                        const int t = ldg2(jpvt + i);
                        jpvt[i] = ldg2(jpvt + pvt);
                        jpvt[pvt] = t;
                        vn1[pvt] = ldg2(vn1 + i);
                        vn2[pvt] = ldg2(vn2 + i);
                    }
                }
                __syncthreads();
                if (warp == 0) {
                    T tau;
                    double beta;
                    make_reflector<true>(A + static_cast<int64_t>(i) * lda, i, m, lane, tau, beta);
                    if (lane == 0) {
                        if (i == 0) box_d[3] = fabs(beta);
                        const double r00 = i == 0 ? fabs(beta) : ldg2(box_d + 3);
                        put_tau(box_d, tau);
                        if (item.tau != nullptr) item.tau[i] = tau;
                        box_d[2] = beta;
                        if (!(fabs(beta) > fmax(tol * r00, abs_floor))) {  // rank rule of compute_id_complex
                            box_i[0] = 1;
                            box_i[1] = i;
                            atomicAdd(stopped_count, 1);
                        } else {
                            box_i[1] = i + 1;
                        }
                    }
                }
            }
        }
        grid.sync();
        if (*reinterpret_cast<volatile int*>(stopped_count) == count) break;
        if (!active) continue;
        if (ldg2(box_i) != 0) {
            active = false;
            continue;
        }
        // ---- the reflector on the owned trailing columns, norm downdates
        const T tau = get_tau<T>(box_d);
        const T* v = A + static_cast<int64_t>(i) * lda;
        const int staged = min(m - i - 1, kStaged);
        for (int r = tid; r < staged; r += kCoopThreads) v_sh[r] = ldg2(v + i + 1 + r);
        __syncthreads();
        for (int t = first_owned(i + 1, part, stride) + warp;; t += kWarps) {
            const int j = part + stride * t;
            if (j >= n) break;
            T* col = A + static_cast<int64_t>(j) * lda;
            apply_reflector_staged(v_sh, staged, v, tau, col, i, m, lane);
            __syncwarp();
            downdate_norm<true>(col, i, m, lane, tol3z, vn1 + j, vn2 + j);
            __syncwarp();
        }
    }

    // ---- T = R11^{-1} R12 on the owned columns (R11 is final: every block
    // passed the barrier after the last reflector)
    const int rank = ldg2(box_i + 1);
    double t_bad = 0.0;
    if (rank == 0) {
        for (int t = tid;; t += kCoopThreads) {
            const int j = part + stride * t;
            if (j >= n) break;
            if (j >= 1 && m > 0) A[static_cast<int64_t>(j) * lda] = T(0.0);
            jpvt[j] = j;  // rank 0 keeps column 0 (the host's convention)
        }
    } else {
        for (int t = first_owned(rank, part, stride) + warp;; t += kWarps) {
            const int j = part + stride * t;
            if (j >= n) break;
            if (!solve_r11<true>(A, lda, rank, A + static_cast<int64_t>(j) * lda, lane)) t_bad = 1.0;
        }
    }
    t_bad = block_max<kWarps>(t_bad, scratch);
    if (tid == 0) {
        if (part == 0) *item.rank = rank;
        if (t_bad != 0.0) atomicMax(item.flag, 2);
    }
}

// ---------------------------------------------------------------------------
// Blocked (LAPACK dlaqps / zlaqps): a step reads the trailing matrix once,
// F(:, k) = tau A^H v over the rows the panel has not updated, and a panel
// of kNb steps updates it once, A -= V F^H; the unblocked kernels above read
// it twice and write it once a step.  A step's own rows (the pivot column's
// update by the panel's reflectors, the reflector's norm, V^H x) are spread
// over the blocks in chunks of kChunk rows, each chunk's sums in a fixed
// order, and every block combines the chunks' sums in chunk order and forms
// the reflector itself: two barriers a step, as the cooperative kernel, and
// a third at a panel's end.  A cancelled norm downdate (dlaqps's test) does
// not end the panel, as in dlaqps: the column's norm is replaced by an upper
// bound until the panel's end, and the panel ends early only if a bound
// would be the next pivot (on the ID targets of 64^3 a downdate cancelled
// every fourth or fifth step, also under a test of eps^(3/4)).  F's rows and
// the column norms are kept by original column, so the pivot swap moves only
// A's columns and jpvt.  The variants (several blocks a box, or one of 256
// or 1024 threads) give a box the same results bit for bit, as above; they
// differ from the unblocked kernels' in rounding only.
// ---------------------------------------------------------------------------
constexpr int kNb = 32;     // panel width: F's row is a warp's lanes
constexpr int kChunk = 256;  // rows a chunk of the row phase
constexpr int kChunkWarps = kChunk / 32;

__host__ __device__ inline int blocked_chunks(int max_m) { return (max_m + kChunk - 1) / kChunk; }
// per box, in elements of T: F (by original column), the chunks' sums
__host__ __device__ inline size_t blocked_box_scalars(int max_m, int max_n) {
    return static_cast<size_t>(max_n) * kNb + static_cast<size_t>(blocked_chunks(max_m)) * (kNb + 1);
}
// vn1, vn2 (by original column), the columns' sums (traced norm); two sets of cand, partials (amax, bad)
__host__ __device__ inline size_t blocked_box_doubles(int max_n, int cpb) {
    return 3 * static_cast<size_t>(max_n) + 4 * static_cast<size_t>(cpb);
}
// two sets of cand position, original column and whether its norm is a bound
__host__ __device__ inline size_t blocked_box_ints(int cpb) { return 6 * static_cast<size_t>(cpb); }

__device__ __forceinline__ double shfl_from(double v, int src) { return __shfl_sync(0xffffffffu, v, src); }
__device__ __forceinline__ dcomplex shfl_from(dcomplex v, int src) {
    return dcomplex(shfl_from(v.re, src), shfl_from(v.im, src));
}

// The reflector of a column from its leading entry alpha and the norm of
// the rest (make_reflector's formulas): beta, tau and the scale of v.
__device__ __forceinline__ void reflector_scalars(double alpha, double xnorm, double& beta, double& tau, double& scal) {
    tau = 0.0;
    beta = alpha;
    scal = 1.0;
    if (xnorm != 0.0) {
        beta = -copysign(sqrt(alpha * alpha + xnorm * xnorm), alpha);
        tau = (beta - alpha) / beta;
        scal = 1.0 / (alpha - beta);
    }
}
__device__ __forceinline__ void reflector_scalars(dcomplex alpha, double xnorm, double& beta, dcomplex& tau,
                                                  dcomplex& scal) {
    tau = dcomplex(0.0);
    beta = alpha.re;
    scal = dcomplex(1.0);
    if (xnorm != 0.0 || alpha.im != 0.0) {
        beta = -copysign(norm3d(alpha.re, alpha.im, xnorm), alpha.re);
        tau = dcomplex((beta - alpha.re) / beta, -alpha.im / beta);
        scal = dcomplex(1.0) / (alpha - dcomplex(beta));
    }
}

// Box `blockIdx.x / cpb` on block `part` (cpb 1: one block a box).  Steps:
// [P] the rows (chunks over the blocks): the pivot from the candidates, the
//     swap, the previous step's v and beta (deferred: other blocks read the
//     column until the barrier), the pivot column's update by the panel's
//     reflectors and the chunks' sums of |x|^2 and V^H x;
// [E] every block: the reflector, then its positions: F's column k, the
//     row i update, the norm downdates (a cancelled one: a bound), the
//     candidates of the next step;
// [G] a panel's end (kNb steps, or a bound the next pivot): the trailing
//     rows' update, the bounded norms afresh, the candidates again.
template<typename T, int kThreads, bool kCoop>
__global__ void __launch_bounds__(kThreads)
qrcp_blocked_kernel(const QrcpItemT<T>* items, int count, int cpb, int max_m, int max_n, double tol, double abs_tol,
                    void* work) {
    constexpr int kWarps = kThreads / 32;
    constexpr int kGroups = kThreads / kChunk;  // the row phase's chunks at a time
    constexpr int kStaged = (40 << 10) / static_cast<int>(sizeof(T));
    constexpr int kTile = kStaged / kNb / 32 * 32;  // rows a tile of the panel update
    constexpr int kPer = kTile / 32;
    __shared__ __align__(16) double buf_raw[(40 << 10) / sizeof(double)];  // v, a tile of V, or the chunks' x
    __shared__ __align__(16) double small_raw[4 * (kNb + 1) * 2];
    __shared__ double scratch[kMaxWarps];
    __shared__ int iscratch[kMaxWarps];
    T* buf = reinterpret_cast<T*>(buf_raw);
    T* s_fp = reinterpret_cast<T*>(small_raw);  // the pivot's F row
    T* s_comb = s_fp + (kNb + 1);               // |x|^2 and u = V^H x over the chunks
    T* s_vrow = s_comb + (kNb + 1);             // V(i, k0 + l)
    T* s_aux = s_vrow + (kNb + 1);              // -tau V^H v
    const int box = blockIdx.x / cpb, part = blockIdx.x % cpb;
    const QrcpItemT<T> item = items[box];
    const int m = item.m, n = item.n, lda = item.lda;
    if (item.tol >= 0.0) tol = item.tol;
    T* A = item.a;
    int* jpvt = item.jpvt;  // the working permutation
    const size_t box_scalars = blocked_box_scalars(max_m, max_n);
    T* F = static_cast<T*>(work) + static_cast<size_t>(box) * box_scalars;  // F(orig, l) at orig * kNb + l
    T* partial = F + static_cast<size_t>(max_n) * kNb;                     // chunk c: |x|^2, then u(l)
    double* dbase = reinterpret_cast<double*>(static_cast<T*>(work) + static_cast<size_t>(count) * box_scalars);
    double* vn1 = dbase + static_cast<size_t>(box) * blocked_box_doubles(max_n, cpb);
    double* vn2 = vn1 + max_n;
    double* col_ssq = vn2 + max_n;
    // the candidates: set 0 from [E] (and the start), set 1 from [G], which
    // runs while other blocks may still read set 0 for the panel's end
    double* cand_v = col_ssq + max_n;
    double* part_amax = cand_v + 2 * cpb;
    double* part_bad = part_amax + cpb;
    int* ibase = reinterpret_cast<int*>(dbase + static_cast<size_t>(count) * blocked_box_doubles(max_n, cpb));
    int* cand_i = ibase + static_cast<size_t>(box) * blocked_box_ints(cpb);
    int* cand_o = cand_i + 2 * cpb;
    int* cand_b = cand_o + 2 * cpb;  // 1: the candidate's norm is a bound
    int* stopped_count = ibase + static_cast<size_t>(count) * blocked_box_ints(cpb);
    int* events = stopped_count + 1;  // the launch's possible panel ends at even and odd steps
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5;
    auto sync = [] {
        if constexpr (kCoop) {
            cooperative_groups::this_grid().sync();
        } else {
            __syncthreads();
        }
    };
    // the block's candidate pivot of step i (into set s): the largest norm over its positions from i, first on
    // ties
    auto candidates = [&](int i, int s) {
        double best = -1.0;
        int best_idx = INT_MAX;
        for (int t = first_owned(i, part, cpb) + tid;; t += kThreads) {
            const int j = part + cpb * t;
            if (j >= n) break;
            const double nj = ldg2(vn1 + ldg2(jpvt + j));
            if (nj > best) {
                best = nj;
                best_idx = j;
            }
        }
        block_argmax_pair<kWarps>(best, best_idx, scratch, iscratch);
        bool bound = false;
        if (tid == 0) {
            const int o = best_idx < n ? ldg2(jpvt + best_idx) : -1;
            bound = o >= 0 && ldg2(vn2 + o) < 0.0;
            const int at = s * cpb + part;
            cand_v[at] = best;
            cand_i[at] = best_idx;
            cand_o[at] = o;
            cand_b[at] = bound ? 1 : 0;
        }
        return bound;  // (thread 0's)
    };
    // the pivot of the next step from the box's candidates (every block the same), and whether its norm is a bound
    int pvt = 0, opvt = 0;
    bool pvt_bound = false;
    auto winner = [&](int s) {
        const int s0 = s * cpb;
        double best = ldg2(cand_v + s0);
        pvt = ldg2(cand_i + s0);
        opvt = ldg2(cand_o + s0);
        pvt_bound = ldg2(cand_b + s0) != 0;
        for (int c = s0 + 1; c < s0 + cpb; ++c) {
            const double cv = ldg2(cand_v + c);
            const int cj = ldg2(cand_i + c);
            if (cv > best || (cv == best && cj < pvt)) {
                best = cv;
                pvt = cj;
                opvt = ldg2(cand_o + c);
                pvt_bound = ldg2(cand_b + c) != 0;
            }
        }
    };

    // ---- norm, finiteness and largest magnitude of the owned columns (as qrcp_coop_kernel)
    double amax = 0.0, bad = 0.0;
    for (int t = warp;; t += kWarps) {
        const int j = part + cpb * t;
        if (j >= n) break;
        const T* col = A + static_cast<int64_t>(j) * lda;
        double s = 0.0;
        for (int r = lane; r < m; r += 32) {
            const T x = ldg2(col + r);
            if (!is_finite(x)) bad = 1.0;
            amax = fmax(amax, magnitude(x));
            s += abs2(x);
        }
        s = warp_sum(s);
        if (lane == 0) col_ssq[j] = s;
    }
    bad = block_max<kWarps>(bad, scratch);
    amax = block_max<kWarps>(amax, scratch);
    if (tid == 0) {
        part_amax[part] = amax;
        part_bad[part] = bad;
        if (kCoop && blockIdx.x == 0) {
            *stopped_count = 0;
            events[0] = events[1] = 0;
        }
    }
    for (int t = tid;; t += kThreads) {
        const int j = part + cpb * t;
        if (j >= n) break;
        jpvt[j] = j;
    }
    sync();
    amax = 0.0;
    bad = 0.0;
    for (int c = 0; c < cpb; ++c) {
        amax = fmax(amax, ldg2(part_amax + c));
        bad = fmax(bad, ldg2(part_bad + c));
    }
    const bool zero = bad != 0.0 || amax == 0.0;
    const double abs_floor = zero ? 0.0 : abs_tol / amax;  // as qrcp_kernel
    if (part == 0 && tid == 0) {
        double ssq = 0.0;
        for (int j = 0; j < n; ++j) ssq += ldg2(col_ssq + j);
        *item.norm = sqrt(ssq);
        *item.flag = bad != 0.0 ? 1 : 0;
        if (kCoop && zero) atomicAdd(stopped_count, 1);  // rank 0
    }
    if (!zero) {
        for (int t = warp;; t += kWarps) {
            const int j = part + cpb * t;
            if (j >= n) break;
            T* col = A + static_cast<int64_t>(j) * lda;
            for (int r = lane; r < m; r += 32) col[r] = ldg2(col + r) / amax;
            __syncwarp();
            double s = 0.0;
            for (int r = lane; r < m; r += 32) s += abs2(ldg2(col + r));
            s = warp_sum(s);
            if (lane == 0) vn1[j] = vn2[j] = sqrt(s);
        }
        __syncthreads();
        candidates(0, 0);
    }
    sync();
    if (!zero) winner(0);

    const double tol3z = sqrt(1.1102230246251565e-16);
    // A cancelled downdate's bound: the norm below the row is at most
    // n1 sqrt(temp + 2 delta) for an estimate n1 within delta relative
    // (dlaqps's test keeps delta near sqrt(eps); 128 times that here).
    const double slack = 256.0 * tol3z;
    const int kmax = m < n ? m : n;
    const int nchunks = blocked_chunks(m);
    bool active = !zero, finished = false;
    int rank = 0, k0 = 0;
    double r00 = 0.0, beta = 0.0, prev_beta = 0.0;
    T tau = T(0.0), scal = T(1.0), prev_scal = T(1.0);
    for (int i = 0;; ++i) {
        const int kk = i - k0;  // the step's column in the panel
        T* ci = A + static_cast<int64_t>(i) * lda;
        // ---- [P] the rows
        if (active) {
            if (tid < kk) s_fp[tid] = ldg2(F + static_cast<size_t>(opvt) * kNb + tid);
            if (part == 0 && tid == 0 && pvt != i) {
                const int t = ldg2(jpvt + i);
                jpvt[i] = ldg2(jpvt + pvt);
                jpvt[pvt] = t;
            }
            __syncthreads();
            T* cp = A + static_cast<int64_t>(pvt) * lda;
            T* cprev = A + static_cast<int64_t>(i > 0 ? i - 1 : 0) * lda;  // (used for i >= 1)
            const int c_first = pvt != i ? 0 : (i > 0 ? (i - 1) / kChunk : 0);
            const int c_live = (i + 1) / kChunk;  // the first chunk with rows below i
            const int group = tid / kChunk, gt = tid % kChunk, wg = warp % kChunkWarps;
            T* xs = buf + group * kChunk;
            for (int cb = c_first + part * kGroups; cb < nchunks; cb += cpb * kGroups) {
                const int c = cb + group;
                const int r = c * kChunk + gt;
                T x = T(0.0);
                if (c < nchunks && r < m) {
                    // the previous step's v and beta
                    T vprev = T(0.0);
                    if (i > 0) {
                        if (r > i - 1) {
                            vprev = ldg2(cprev + r) * prev_scal;
                            cprev[r] = vprev;
                        } else if (r == i - 1) {
                            cprev[r] = T(prev_beta);
                        }
                    }
                    // the swap
                    T a = ldg2(ci + r);
                    if (pvt != i) {
                        const T b = ldg2(cp + r);
                        cp[r] = a;
                        a = b;
                    }
                    if (r >= i) {
                        // A(r, i) -= sum_l V(r, k0 + l) conj(F(pvt, l))
                        T s = T(0.0);
#pragma unroll 8
                        for (int l = 0; l < kk; ++l) {
                            const T vl = l == kk - 1 ? vprev : ldg2(A + r + static_cast<int64_t>(k0 + l) * lda);
                            s += vl * conj(s_fp[l]);
                        }
                        a = a - s;
                        if (r > i) x = a;
                    }
                    if (pvt != i || r >= i) ci[r] = a;
                }
                xs[gt] = x;
                __syncthreads();
                // the chunk's sums by its warps: q 0 |x|^2, 1 + l conj(V(r, k0 + l)) x, lanes over its rows in order
                if (c < nchunks && c >= c_live) {
                    for (int q = wg; q <= kk; q += kChunkWarps) {
                        T s = T(0.0);
#pragma unroll
                        for (int u = 0; u < kChunk / 32; ++u) {
                            const int rr = lane + 32 * u, rw = c * kChunk + rr;
                            if (rw > i && rw < m) {
                                const T xv = xs[rr];
                                s += q == 0 ? T(abs2(xv))
                                            : conj(ldg2(A + rw + static_cast<int64_t>(k0 + q - 1) * lda)) * xv;
                            }
                        }
                        s = warp_sum(s);
                        if (lane == 0) partial[static_cast<size_t>(c) * (kNb + 1) + q] = s;
                    }
                }
                __syncthreads();
            }
        }
        sync();
        // ---- [E] the reflector, F's column kk, row i, the norms, the next candidates
        if (kCoop && blockIdx.x == 0 && tid == 0) events[(i + 1) & 1] = 0;  // (last read before this step's [P] barrier)
        if (active) {
            const int c_live = (i + 1) / kChunk;
            if (tid <= kk) {
                T s = T(0.0);
                for (int c = c_live; c < nchunks; ++c) s += ldg2(partial + static_cast<size_t>(c) * (kNb + 1) + tid);
                s_comb[tid] = s;
            }
            __syncthreads();
            reflector_scalars(ldg2(ci + i), sqrt(real_part(s_comb[0])), beta, tau, scal);
            if (i == 0) r00 = fabs(beta);
            if (!(fabs(beta) > fmax(tol * r00, abs_floor))) {  // rank rule of compute_id_complex
                rank = i;
                active = false;
                if (kCoop && part == 0 && tid == 0) atomicAdd(stopped_count, 1);
            } else {
                rank = i + 1;
                if (tid < kk) {
                    // h(l) = V(:, l)^H v = conj(V(i, l)) + u(l) scal
                    const T vil = ldg2(A + i + static_cast<int64_t>(k0 + tid) * lda);
                    s_vrow[tid] = vil;
                    s_aux[tid] = -(tau * (conj(vil) + s_comb[1 + tid] * scal));
                }
                const int staged = min(m - i - 1, kStaged);
                for (int r = tid; r < staged; r += kThreads) buf[r] = ldg2(ci + i + 1 + r) * scal;
                __syncthreads();
                const int split = i + 1 + staged;
                const T sc = scal;
                auto vat = [&](int r) { return r < split ? buf[r - i - 1] : ldg2(ci + r) * sc; };
                for (int t = first_owned(i + 1, part, cpb) + warp;; t += kWarps) {
                    const int j = part + cpb * t;
                    if (j >= n) break;
                    T* col = A + static_cast<int64_t>(j) * lda;
                    const int orig = ldg2(jpvt + j);
                    T* Fj = F + static_cast<size_t>(orig) * kNb;
                    // g = A(i:m, j)^H v (v(i) = 1)
                    const T aij0 = ldg2(col + i);
                    T g = lane == 0 ? conj(aij0) : T(0.0);
                    int r = i + 1 + lane;
                    for (; r + 32 * (kAhead - 1) < m; r += 32 * kAhead) {
                        T c[kAhead], v[kAhead];
                        load_ahead<true>(col, r, c);
#pragma unroll
                        for (int u = 0; u < kAhead; ++u) v[u] = vat(r + 32 * u);
#pragma unroll
                        for (int u = 0; u < kAhead; ++u) g += conj(c[u]) * v[u];
                    }
                    for (; r < m; r += 32) g += conj(ldg2(col + r)) * vat(r);
                    g = warp_sum(g);
                    // F(j, kk) = tau g + sum_l F(j, l) aux(l)
                    const T f_l = lane < kk ? ldg2(Fj + lane) : T(0.0);
                    const T fk = tau * g + warp_sum(lane < kk ? f_l * s_aux[lane] : T(0.0));
                    // A(i, j) -= sum_{l <= kk} V(i, k0 + l) conj(F(j, l)), V(i, i) = 1
                    const T d = warp_sum(lane < kk ? s_vrow[lane] * conj(f_l) : (lane == kk ? conj(fk) : T(0.0)));
                    if (lane == 0) {
                        Fj[kk] = fk;
                        const T aij = aij0 - d;
                        col[i] = aij;
                        const double n1 = ldg2(vn1 + orig), n2 = ldg2(vn2 + orig);
                        if (i < kmax - 1 && n1 != 0.0 && n2 >= 0.0) {  // dlaqps's downdate (a bound stays)
                            double temp = magnitude(aij) / n1;
                            temp = fmax(0.0, (1.0 + temp) * (1.0 - temp));
                            const double ratio = n1 / n2;
                            if (temp * ratio * ratio <= tol3z) {
                                vn1[orig] = n1 * sqrt(temp + slack);  // cancelled: a bound, afresh at the panel's end
                                vn2[orig] = -1.0;
                            } else {
                                vn1[orig] = n1 * sqrt(temp);
                            }
                        }
                    }
                }
                __syncthreads();
                if (i + 1 == kmax) {
                    active = false;
                    finished = true;
                    if (kCoop && part == 0 && tid == 0) atomicAdd(stopped_count, 1);
                } else {
                    // (a panel may end: its last step, or a bound among the candidates)
                    const bool bound = candidates(i + 1, 0);
                    if (kCoop && tid == 0 && ((part == 0 && kk == kNb - 1) || bound)) atomicAdd(events + (i & 1), 1);
                }
            }
        }
        sync();
        if constexpr (kCoop) {
            if (*reinterpret_cast<volatile int*>(stopped_count) == count) break;
        } else {
            if (!active) break;
        }
        // ---- [G] a panel's end: its last step, or the next pivot's norm a bound
        if (active) winner(0);
        const bool box_end = active && (kk == kNb - 1 || pvt_bound);
        bool any_end = box_end;
        if constexpr (kCoop) any_end = ldg2(events + (i & 1)) > 0;
        if (any_end) {
            if (box_end) {
                // A(r, j) -= sum_{l < kb} V(r, k0 + l) conj(F(j, l)) for rows r > i and the positions j > i;
                // V's column i still holds x (times scal)
                const int kb = kk + 1;
                for (int r0 = i + 1; r0 < m; r0 += kTile) {
                    const int rows = min(kTile, m - r0);
                    __syncthreads();
                    for (int e = tid; e < kb * kTile; e += kThreads) {
                        const int l = e / kTile, rr = e % kTile;
                        T v = T(0.0);
                        if (rr < rows) {
                            v = ldg2(A + r0 + rr + static_cast<int64_t>(k0 + l) * lda);
                            if (l == kb - 1) v = v * scal;
                        }
                        buf[e] = v;
                    }
                    __syncthreads();
                    for (int t = first_owned(i + 1, part, cpb) + warp;; t += kWarps) {
                        const int j = part + cpb * t;
                        if (j >= n) break;
                        T* col = A + static_cast<int64_t>(j) * lda;
                        const int orig = ldg2(jpvt + j);
                        const T f = lane < kb ? conj(ldg2(F + static_cast<size_t>(orig) * kNb + lane)) : T(0.0);
                        T s[kPer];
#pragma unroll
                        for (int u = 0; u < kPer; ++u) s[u] = T(0.0);
                        for (int l = 0; l < kb; ++l) {
                            const T fl = shfl_from(f, l);
#pragma unroll
                            for (int u = 0; u < kPer; ++u) s[u] += buf[l * kTile + lane + 32 * u] * fl;
                        }
#pragma unroll
                        for (int u = 0; u < kPer; ++u) {
                            const int rr = lane + 32 * u;
                            if (rr < rows) col[r0 + rr] = ldg2(col + r0 + rr) - s[u];
                        }
                    }
                }
                __syncthreads();
                // the bounded norms afresh, rows below i
                for (int t = first_owned(i + 1, part, cpb) + warp;; t += kWarps) {
                    const int j = part + cpb * t;
                    if (j >= n) break;
                    const int orig = ldg2(jpvt + j);
                    if (ldg2(vn2 + orig) >= 0.0) continue;
                    const T* col = A + static_cast<int64_t>(j) * lda;
                    double s2 = 0.0;
                    for (int r = i + 1 + lane; r < m; r += 32) s2 += abs2(ldg2(col + r));
                    const double fresh = sqrt(warp_sum(s2));
                    if (lane == 0) vn1[orig] = vn2[orig] = fresh;
                }
                __syncthreads();
                candidates(i + 1, 1);
                k0 = i + 1;
            }
            sync();
            if (box_end) winner(1);
        }
        prev_scal = scal;
        prev_beta = beta;
    }

    // ---- the last reflector's diagonal (deferred) of a box that ran to kmax; T = R11^{-1} R12 on the owned
    // columns
    if (finished && part == 0 && tid == 0) A[static_cast<int64_t>(kmax - 1) * lda + kmax - 1] = T(beta);
    sync();
    double t_bad = 0.0;
    if (rank == 0) {
        for (int t = tid;; t += kThreads) {
            const int j = part + cpb * t;
            if (j >= n) break;
            if (j >= 1 && m > 0) A[static_cast<int64_t>(j) * lda] = T(0.0);
            jpvt[j] = j;  // rank 0 keeps column 0 (the host's convention)
        }
    } else {
        for (int t = first_owned(rank, part, cpb) + warp;; t += kWarps) {
            const int j = part + cpb * t;
            if (j >= n) break;
            if (!solve_r11<true>(A, lda, rank, A + static_cast<int64_t>(j) * lda, lane)) t_bad = 1.0;
        }
    }
    t_bad = block_max<kWarps>(t_bad, scratch);
    if (tid == 0) {
        if (part == 0) *item.rank = rank;
        if (t_bad != 0.0) atomicMax(item.flag, 2);
    }
}

// Blocks per box of the cooperative path (0: one block per box instead),
// of qrcp_coop_kernel or (kBlocked) of qrcp_blocked_kernel.
template<typename T, bool kBlocked = false>
int coop_blocks_per_box(int count, int max_n) {
    static const int sms = [] {
        int device = 0, value = 0;
        cudaGetDevice(&device);
        cudaDeviceGetAttribute(&value, cudaDevAttrMultiProcessorCount, device);
        return value;
    }();
    static const int resident = [] {
        int per_sm = 0;
        if constexpr (kBlocked) {
            cudaOccupancyMaxActiveBlocksPerMultiprocessor(&per_sm, qrcp_blocked_kernel<T, kCoopThreads, true>,
                                                          kCoopThreads, 0);
        } else {
            cudaOccupancyMaxActiveBlocksPerMultiprocessor(&per_sm, qrcp_coop_kernel<T>, kCoopThreads, 0);
        }
        return per_sm * sms;
    }();
    if (count <= 0 || count > sms) return 0;
    const int by_columns = std::max(1, (max_n + 15) / 16);  // at least 16 columns per block
    const int cpb = std::min(std::min(2 * sms / count, by_columns), resident / count);
    // With fewer blocks per box the grid barriers and the uncached loads
    // cost more than one large block per box (measured: 27 boxes, 8 blocks
    // each, is slower; 7 boxes, 30 blocks each, is 2.5x faster).
    constexpr int kMinBlocksPerBox = 12;
    return cpb >= kMinBlocksPerBox ? cpb : 0;
}

template<typename T, int kThreads>
void launch_qrcp_with(const QrcpItemT<T>* items, int count, size_t shared, double tol, double abs_tol,
                      cudaStream_t stream) {
    if (shared > 48 * 1024) {
        const cudaError_t attr = cudaFuncSetAttribute(qrcp_kernel<T, kThreads>,
                                                      cudaFuncAttributeMaxDynamicSharedMemorySize,
                                                      static_cast<int>(shared));
        if (attr != cudaSuccess) throw std::runtime_error("launch_qrcp: box too large for shared memory");
    }
    qrcp_kernel<T, kThreads><<<static_cast<unsigned>(count), kThreads, shared, stream>>>(items, tol, abs_tol);
}

}  // namespace

template<typename T>
size_t qrcp_work_bytes(int count, int max_n) {
    const int cpb = coop_blocks_per_box<T>(count, max_n);
    if (cpb == 0) return 0;
    return static_cast<size_t>(count) * coop_box_doubles(max_n, cpb) * sizeof(double) +
           (static_cast<size_t>(count) * coop_box_ints(cpb) + 1) * sizeof(int);
}

// Few boxes (at most 2 x multiprocessors / kMinBlocksPerBox): the cooperative
// path, several blocks per box.  Otherwise one block per box: large blocks for
// few boxes, small ones (several per multiprocessor) for many.  The variants
// agree bitwise (see the top of the file); batch_independent always takes one
// block of 256 threads per box (kept for tests/batch_determinism.cpp).
template<typename T>
void launch_qrcp(const QrcpItemT<T>* items, int count, int max_n, double tol, double abs_tol, void* work,
                 cudaStream_t stream, bool batch_independent) {
    if (count <= 0) return;
    int cpb = batch_independent ? 0 : coop_blocks_per_box<T>(count, max_n);
    if (cpb > 0) {
        if (work == nullptr) throw std::runtime_error("launch_qrcp: the cooperative path needs its work block");
        double* dwork = static_cast<double*>(work);
        int* iwork = reinterpret_cast<int*>(dwork + static_cast<size_t>(count) * coop_box_doubles(max_n, cpb));
        void* args[] = {const_cast<QrcpItemT<T>**>(&items), &count, &cpb, &max_n, &tol, &abs_tol, &dwork, &iwork};
        const cudaError_t status = cudaLaunchCooperativeKernel(reinterpret_cast<const void*>(qrcp_coop_kernel<T>),
                                                               dim3(static_cast<unsigned>(count * cpb)),
                                                               dim3(kCoopThreads), args, 0, stream);
        if (status != cudaSuccess) {
            throw std::runtime_error(std::string("CUDA launch failed in qrcp_coop_kernel: ") + cudaGetErrorString(status));
        }
        return;
    }
    const size_t shared = static_cast<size_t>(2 * max_n + kMaxWarps) * sizeof(double) +
                          static_cast<size_t>(kMaxWarps + max_n) * sizeof(int);
    if (count < 216 && !batch_independent) {
        launch_qrcp_with<T, 1024>(items, count, shared, tol, abs_tol, stream);
    } else {
        launch_qrcp_with<T, 256>(items, count, shared, tol, abs_tol, stream);
    }
    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA launch failed in qrcp_kernel: ") + cudaGetErrorString(status));
    }
}

template<typename T>
size_t qrcp_blocked_work_bytes(int count, int max_m, int max_n) {
    if (count <= 0) return 0;
    const int cpb = std::max(1, coop_blocks_per_box<T, true>(count, max_n));
    return static_cast<size_t>(count) * blocked_box_scalars(max_m, max_n) * sizeof(T) +
           static_cast<size_t>(count) * blocked_box_doubles(max_n, cpb) * sizeof(double) +
           (static_cast<size_t>(count) * blocked_box_ints(cpb) + 4) * sizeof(int);
}

// As launch_qrcp: few boxes on several blocks each (cooperative), otherwise
// one block a box; always with its work block (qrcp_blocked_work_bytes).
template<typename T>
void launch_qrcp_blocked(const QrcpItemT<T>* items, int count, int max_m, int max_n, double tol, double abs_tol,
                         void* work, cudaStream_t stream) {
    if (count <= 0) return;
    if (work == nullptr) throw std::runtime_error("launch_qrcp_blocked: no work block");
    int cpb = coop_blocks_per_box<T, true>(count, max_n);
    if (cpb > 0) {
        void* args[] = {const_cast<QrcpItemT<T>**>(&items), &count, &cpb, &max_m, &max_n, &tol, &abs_tol, &work};
        const cudaError_t status = cudaLaunchCooperativeKernel(
            reinterpret_cast<const void*>(qrcp_blocked_kernel<T, kCoopThreads, true>),
            dim3(static_cast<unsigned>(count * cpb)), dim3(kCoopThreads), args, 0, stream);
        if (status != cudaSuccess) {
            throw std::runtime_error(std::string("CUDA launch failed in qrcp_blocked_kernel: ") +
                                     cudaGetErrorString(status));
        }
        return;
    }
    if (count < 216) {
        qrcp_blocked_kernel<T, 1024, false>
            <<<static_cast<unsigned>(count), 1024, 0, stream>>>(items, count, 1, max_m, max_n, tol, abs_tol, work);
    } else {
        qrcp_blocked_kernel<T, 256, false>
            <<<static_cast<unsigned>(count), 256, 0, stream>>>(items, count, 1, max_m, max_n, tol, abs_tol, work);
    }
    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA launch failed in qrcp_blocked_kernel: ") + cudaGetErrorString(status));
    }
}

template size_t qrcp_work_bytes<double>(int, int);
template size_t qrcp_work_bytes<dcomplex>(int, int);
template void launch_qrcp<double>(const QrcpItemT<double>*, int, int, double, double, void*, cudaStream_t, bool);
template void launch_qrcp<dcomplex>(const QrcpItemT<dcomplex>*, int, int, double, double, void*, cudaStream_t, bool);
template size_t qrcp_blocked_work_bytes<double>(int, int, int);
template size_t qrcp_blocked_work_bytes<dcomplex>(int, int, int);
template void launch_qrcp_blocked<double>(const QrcpItemT<double>*, int, int, int, double, double, void*,
                                          cudaStream_t);
template void launch_qrcp_blocked<dcomplex>(const QrcpItemT<dcomplex>*, int, int, int, double, double, void*,
                                            cudaStream_t);

}  // namespace gpu
}  // namespace fmm
