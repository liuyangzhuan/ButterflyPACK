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

// Householder reflector of col[i:m) by one warp (LAPACK dlarfg): col[i] =
// beta, col[i+1:m) = v; returns tau and beta.
template<bool kUncached>
__device__ __forceinline__ void make_reflector(double* col, int i, int m, int lane, double& tau, double& beta) {
    const double alpha = load<kUncached>(col + i);
    double s = 0.0;
    for (int r = i + 1 + lane; r < m; r += 32) {
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
        for (int r = i + 1 + lane; r < m; r += 32) col[r] = load<kUncached>(col + r) * scal;
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
    for (int r = i + 1 + lane; r < m; r += 32) s += abs2(load<kUncached>(col + r));
    const double xnorm = sqrt(warp_sum(s));
    tau = dcomplex(0.0);
    beta = alpha.re;
    if (xnorm != 0.0 || alpha.im != 0.0) {
        beta = -copysign(norm3d(alpha.re, alpha.im, xnorm), alpha.re);  // dlapy3
        tau = dcomplex((beta - alpha.re) / beta, -alpha.im / beta);
        const dcomplex scal = dcomplex(1.0) / (alpha - dcomplex(beta));
        for (int r = i + 1 + lane; r < m; r += 32) col[r] = load<kUncached>(col + r) * scal;
        __syncwarp();
        if (lane == 0) col[i] = dcomplex(beta);
    }
}

// H^H = I - conj(tau) v v^H applied to col by one warp (v_i = 1 implicit).
template<bool kUncached, typename T>
__device__ __forceinline__ void apply_reflector(const T* v, T tau, T* col, int i, int m, int lane) {
    T s = lane == 0 ? load<kUncached>(col + i) : T(0.0);
    for (int r = i + 1 + lane; r < m; r += 32) s += conj(load<kUncached>(v + r)) * load<kUncached>(col + r);
    const T w = conj(tau) * warp_sum(s);
    if (w != T(0.0)) {
        if (lane == 0) col[i] = load<kUncached>(col + i) - w;
        for (int r = i + 1 + lane; r < m; r += 32) col[r] = load<kUncached>(col + r) - w * load<kUncached>(v + r);
    }
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
// pivot choice and norm downdating), stop at the first |R_ii| <= tol |R_00|,
// and solve R11 T = R12 in place.  The factorization stops at the rank, so
// the redundant columns keep the order of that step (the host's full
// factorization permutes them further; the skeleton is the same).
template<typename T, int kThreads>
__global__ void __launch_bounds__(kThreads) qrcp_kernel(const QrcpItemT<T>* items, double tol) {
    constexpr int kWarps = kThreads / 32;
    extern __shared__ double smem[];
    const QrcpItemT<T> item = items[blockIdx.x];
    const int m = item.m, n = item.n, lda = item.lda;
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
            }
        }
        __syncthreads();
        const T tau = s_tau;
        const double beta = s_beta;
        if (i == 0) r00 = fabs(beta);
        if (!(fabs(beta) > tol * r00)) break;  // rank rule of compute_id_complex
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
qrcp_coop_kernel(const QrcpItemT<T>* items, int count, int cpb, int max_n, double tol, double* dwork, int* iwork) {
    namespace cg = cooperative_groups;
    cg::grid_group grid = cg::this_grid();
    constexpr int kWarps = kCoopThreads / 32;
    __shared__ double scratch[kMaxWarps];
    __shared__ int iscratch[kMaxWarps];
    const int box = blockIdx.x / cpb, part = blockIdx.x % cpb;
    const QrcpItemT<T> item = items[box];
    const int m = item.m, n = item.n, lda = item.lda;
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
                        box_d[2] = beta;
                        if (!(fabs(beta) > tol * r00)) {  // rank rule of compute_id_complex
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
        for (int t = first_owned(i + 1, part, stride) + warp;; t += kWarps) {
            const int j = part + stride * t;
            if (j >= n) break;
            T* col = A + static_cast<int64_t>(j) * lda;
            apply_reflector<true>(v, tau, col, i, m, lane);
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

// Blocks per box of the cooperative path (0: one block per box instead).
template<typename T>
int coop_blocks_per_box(int count, int max_n) {
    static const int sms = [] {
        int device = 0, value = 0;
        cudaGetDevice(&device);
        cudaDeviceGetAttribute(&value, cudaDevAttrMultiProcessorCount, device);
        return value;
    }();
    static const int resident = [] {
        int per_sm = 0;
        cudaOccupancyMaxActiveBlocksPerMultiprocessor(&per_sm, qrcp_coop_kernel<T>, kCoopThreads, 0);
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
void launch_qrcp_with(const QrcpItemT<T>* items, int count, size_t shared, double tol, cudaStream_t stream) {
    if (shared > 48 * 1024) {
        const cudaError_t attr = cudaFuncSetAttribute(qrcp_kernel<T, kThreads>,
                                                      cudaFuncAttributeMaxDynamicSharedMemorySize,
                                                      static_cast<int>(shared));
        if (attr != cudaSuccess) throw std::runtime_error("launch_qrcp: box too large for shared memory");
    }
    qrcp_kernel<T, kThreads><<<static_cast<unsigned>(count), kThreads, shared, stream>>>(items, tol);
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
void launch_qrcp(const QrcpItemT<T>* items, int count, int max_n, double tol, void* work, cudaStream_t stream,
                 bool batch_independent) {
    if (count <= 0) return;
    int cpb = batch_independent ? 0 : coop_blocks_per_box<T>(count, max_n);
    if (cpb > 0) {
        if (work == nullptr) throw std::runtime_error("launch_qrcp: the cooperative path needs its work block");
        double* dwork = static_cast<double*>(work);
        int* iwork = reinterpret_cast<int*>(dwork + static_cast<size_t>(count) * coop_box_doubles(max_n, cpb));
        void* args[] = {const_cast<QrcpItemT<T>**>(&items), &count, &cpb, &max_n, &tol, &dwork, &iwork};
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
        launch_qrcp_with<T, 1024>(items, count, shared, tol, stream);
    } else {
        launch_qrcp_with<T, 256>(items, count, shared, tol, stream);
    }
    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA launch failed in qrcp_kernel: ") + cudaGetErrorString(status));
    }
}

template size_t qrcp_work_bytes<double>(int, int);
template size_t qrcp_work_bytes<dcomplex>(int, int);
template void launch_qrcp<double>(const QrcpItemT<double>*, int, int, double, void*, cudaStream_t, bool);
template void launch_qrcp<dcomplex>(const QrcpItemT<dcomplex>*, int, int, double, void*, cudaStream_t, bool);

}  // namespace gpu
}  // namespace fmm
