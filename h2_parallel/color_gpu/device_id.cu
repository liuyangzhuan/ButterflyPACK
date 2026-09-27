// Batched interpolative decomposition of the H2 Color GPU backend: column-
// pivoted Householder QR with the host's rank rule, then T = R11^{-1} R12.
// See QrcpItem in device_kernels.hpp.

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

// Sum over the block in a fixed order (deterministic); all threads get it.
template<int kWarps>
__device__ double block_sum(double v, double* scratch) {
    v = warp_sum(v);
    __syncthreads();
    if ((threadIdx.x & 31) == 0) scratch[threadIdx.x >> 5] = v;
    __syncthreads();
    double total = 0.0;
    for (int w = 0; w < kWarps; ++w) total += scratch[w];
    return total;
}

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

// One matrix per block.  Mirrors fmm::compute_id_complex: normalize by the
// largest magnitude, factor A P = Q R with column pivoting (LAPACK dlaqp2's
// pivot choice and norm downdating), stop at the first |R_ii| <= tol |R_00|,
// and solve R11 T = R12 in place.  The factorization stops at the rank, so
// the redundant columns keep the order of that step (the host's full
// factorization permutes them further; the skeleton is the same).
template<int kThreads>
__global__ void __launch_bounds__(kThreads) qrcp_kernel(const QrcpItem* items, double tol) {
    constexpr int kWarps = kThreads / 32;
    extern __shared__ double smem[];
    const QrcpItem item = items[blockIdx.x];
    const int m = item.m, n = item.n, lda = item.lda;
    double* A = item.a;
    double* vn1 = smem;
    double* vn2 = smem + n;
    double* scratch = smem + 2 * n;                             // kMaxWarps doubles
    int* iscratch = reinterpret_cast<int*>(scratch + kMaxWarps);  // kMaxWarps ints
    int* jpvt = iscratch + kMaxWarps;                           // n ints
    __shared__ double s_tau, s_beta;
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5;

    // ---- norm (for traces), finiteness, normalization
    double amax = 0.0, ssq = 0.0, bad = 0.0;
    for (int j = warp; j < n; j += kWarps) {
        const double* col = A + static_cast<int64_t>(j) * lda;
        for (int r = lane; r < m; r += 32) {
            const double x = col[r];
            if (!isfinite(x)) bad = 1.0;
            amax = fmax(amax, fabs(x));
            ssq += x * x;
        }
    }
    bad = block_max<kWarps>(bad, scratch);
    ssq = block_sum<kWarps>(ssq, scratch);
    amax = block_max<kWarps>(amax, scratch);
    if (tid == 0) {
        *item.norm = sqrt(ssq);
        *item.flag = bad != 0.0 ? 1 : 0;
    }
    for (int j = tid; j < n; j += kThreads) jpvt[j] = j;
    if (bad != 0.0 || amax == 0.0) {
        // zero target: rank 0 (the host keeps column 0 with a zero row of T)
        __syncthreads();
        for (int j = tid; j < n; j += kThreads) {
            item.jpvt[j] = j;
            if (j >= 1 && m > 0) A[static_cast<int64_t>(j) * lda] = 0.0;
        }
        if (tid == 0) *item.rank = 0;
        return;
    }
    for (int j = warp; j < n; j += kWarps) {
        double* col = A + static_cast<int64_t>(j) * lda;
        for (int r = lane; r < m; r += 32) col[r] /= amax;
    }
    __syncthreads();

    // ---- initial column norms
    for (int j = warp; j < n; j += kWarps) {
        const double* col = A + static_cast<int64_t>(j) * lda;
        double s = 0.0;
        for (int r = lane; r < m; r += 32) s += col[r] * col[r];
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
            double* ci = A + static_cast<int64_t>(i) * lda;
            double* cp = A + static_cast<int64_t>(pvt) * lda;
            for (int r = tid; r < m; r += kThreads) {
                const double t = ci[r];
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

        // Householder reflector of column i (LAPACK dlarfg)
        if (warp == 0) {
            double* col = A + static_cast<int64_t>(i) * lda;
            const double alpha = col[i];
            double s = 0.0;
            for (int r = i + 1 + lane; r < m; r += 32) s += col[r] * col[r];
            const double xnorm = sqrt(warp_sum(s));
            double tau = 0.0, beta = alpha;
            if (xnorm != 0.0) {
                beta = -copysign(sqrt(alpha * alpha + xnorm * xnorm), alpha);
                tau = (beta - alpha) / beta;
                const double scal = 1.0 / (alpha - beta);
                for (int r = i + 1 + lane; r < m; r += 32) col[r] *= scal;
                __syncwarp();
                if (lane == 0) col[i] = beta;
            }
            if (lane == 0) {
                s_tau = tau;
                s_beta = beta;
            }
        }
        __syncthreads();
        const double tau = s_tau, beta = s_beta;
        if (i == 0) r00 = fabs(beta);
        if (!(fabs(beta) > tol * r00)) break;  // rank rule of compute_id_complex
        rank = i + 1;

        // apply the reflector to the trailing columns and downdate their norms
        const double* v = A + static_cast<int64_t>(i) * lda;
        for (int j = i + 1 + warp; j < n; j += kWarps) {
            double* col = A + static_cast<int64_t>(j) * lda;
            double s = lane == 0 ? col[i] : 0.0;
            for (int r = i + 1 + lane; r < m; r += 32) s += v[r] * col[r];
            const double w = tau * warp_sum(s);
            if (w != 0.0) {
                if (lane == 0) col[i] -= w;
                for (int r = i + 1 + lane; r < m; r += 32) col[r] -= w * v[r];
            }
            __syncwarp();
            const double n1 = vn1[j];
            if (n1 != 0.0) {
                double temp = fabs(col[i]) / n1;
                temp = fmax(1.0 - temp * temp, 0.0);
                const double ratio = n1 / vn2[j];
                if (temp * ratio * ratio <= tol3z) {
                    double s2 = 0.0;
                    for (int r = i + 1 + lane; r < m; r += 32) s2 += col[r] * col[r];
                    const double fresh = i < m - 1 ? sqrt(warp_sum(s2)) : 0.0;
                    if (lane == 0) vn1[j] = vn2[j] = fresh;
                } else if (lane == 0) {
                    vn1[j] = n1 * sqrt(temp);
                }
            }
            __syncwarp();
        }
        __syncthreads();
    }

    // ---- T = R11^{-1} R12, one warp per column, in place
    double t_bad = 0.0;
    if (rank == 0) {
        for (int j = 1 + tid; j < n; j += kThreads) A[static_cast<int64_t>(j) * lda] = 0.0;
    } else {
        for (int j = rank + warp; j < n; j += kWarps) {
            double* col = A + static_cast<int64_t>(j) * lda;
            for (int p = rank - 1; p >= 0; --p) {
                double s = 0.0;
                for (int q = p + 1 + lane; q < rank; q += 32) s += A[p + static_cast<int64_t>(q) * lda] * col[q];
                s = warp_sum(s);
                const double x = (col[p] - s) / A[p + static_cast<int64_t>(p) * lda];
                __syncwarp();
                if (lane == 0) col[p] = x;
                __syncwarp();
                if (!isfinite(x)) t_bad = 1.0;
            }
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
// qrcp_kernel, so pivots, R and T are the same; only the Frobenius norm
// (traces) is summed in another order.
// ---------------------------------------------------------------------------
constexpr int kCoopThreads = 256;

__host__ __device__ inline size_t coop_box_doubles(int max_n, int cpb) {
    return 2 * static_cast<size_t>(max_n) + 4 * static_cast<size_t>(cpb) + 4;  // vn1, vn2, cand, partials, tau/beta/r00
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

// Loads of data other blocks write between barriers bypass L1 (it is not
// coherent across multiprocessors).
template<typename T>
__device__ __forceinline__ T ldg2(const T* p) { return __ldcg(p); }

__global__ void __launch_bounds__(kCoopThreads)
qrcp_coop_kernel(const QrcpItem* items, int count, int cpb, int max_n, double tol, double* dwork, int* iwork) {
    namespace cg = cooperative_groups;
    cg::grid_group grid = cg::this_grid();
    constexpr int kWarps = kCoopThreads / 32;
    __shared__ double scratch[kMaxWarps];
    __shared__ int iscratch[kMaxWarps];
    const int box = blockIdx.x / cpb, part = blockIdx.x % cpb;
    const QrcpItem item = items[box];
    const int m = item.m, n = item.n, lda = item.lda;
    double* A = item.a;
    int* jpvt = item.jpvt;  // the working permutation
    double* vn1 = dwork + static_cast<size_t>(box) * coop_box_doubles(max_n, cpb);
    double* vn2 = vn1 + max_n;
    double* cand_v = vn2 + max_n;
    double* part_amax = cand_v + cpb;
    double* part_ssq = part_amax + cpb;
    double* part_bad = part_ssq + cpb;
    double* box_d = part_bad + cpb;  // tau, beta, r00
    int* cand_i = iwork + static_cast<size_t>(box) * coop_box_ints(cpb);
    int* box_i = cand_i + cpb;       // stopped, rank
    int* stopped_count = iwork + static_cast<size_t>(count) * coop_box_ints(cpb);
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5;
    const int stride = cpb;  // between owned positions

    // ---- norm, finiteness and largest magnitude of the owned columns
    double amax = 0.0, ssq = 0.0, bad = 0.0;
    for (int t = warp;; t += kWarps) {
        const int j = part + stride * t;
        if (j >= n) break;
        const double* col = A + static_cast<int64_t>(j) * lda;
        for (int r = lane; r < m; r += 32) {
            const double x = ldg2(col + r);
            if (!isfinite(x)) bad = 1.0;
            amax = fmax(amax, fabs(x));
            ssq += x * x;
        }
    }
    bad = block_max<kWarps>(bad, scratch);
    ssq = block_sum<kWarps>(ssq, scratch);
    amax = block_max<kWarps>(amax, scratch);
    if (tid == 0) {
        part_amax[part] = amax;
        part_ssq[part] = ssq;
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
    ssq = 0.0;
    bad = 0.0;
    for (int c = 0; c < cpb; ++c) {
        amax = fmax(amax, ldg2(part_amax + c));
        ssq += ldg2(part_ssq + c);
        bad = fmax(bad, ldg2(part_bad + c));
    }
    const bool zero = bad != 0.0 || amax == 0.0;
    if (part == 0 && tid == 0) {
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
            double* col = A + static_cast<int64_t>(j) * lda;
            for (int r = lane; r < m; r += 32) col[r] = ldg2(col + r) / amax;
            __syncwarp();
            double s = 0.0;
            for (int r = lane; r < m; r += 32) {
                const double x = ldg2(col + r);
                s += x * x;
            }
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
                    double* ci = A + static_cast<int64_t>(i) * lda;
                    double* cp = A + static_cast<int64_t>(pvt) * lda;
                    for (int r = tid; r < m; r += kCoopThreads) {
                        const double t = ldg2(ci + r);
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
                    double* col = A + static_cast<int64_t>(i) * lda;
                    const double alpha = ldg2(col + i);
                    double s = 0.0;
                    for (int r = i + 1 + lane; r < m; r += 32) {
                        const double x = ldg2(col + r);
                        s += x * x;
                    }
                    const double xnorm = sqrt(warp_sum(s));
                    double tau = 0.0, beta = alpha;
                    if (xnorm != 0.0) {
                        beta = -copysign(sqrt(alpha * alpha + xnorm * xnorm), alpha);
                        tau = (beta - alpha) / beta;
                        const double scal = 1.0 / (alpha - beta);
                        for (int r = i + 1 + lane; r < m; r += 32) col[r] = ldg2(col + r) * scal;
                        __syncwarp();
                        if (lane == 0) col[i] = beta;
                    }
                    if (lane == 0) {
                        if (i == 0) box_d[2] = fabs(beta);
                        const double r00 = i == 0 ? fabs(beta) : ldg2(box_d + 2);
                        box_d[0] = tau;
                        box_d[1] = beta;
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
        const double tau = ldg2(box_d);
        const double* v = A + static_cast<int64_t>(i) * lda;
        for (int t = first_owned(i + 1, part, stride) + warp;; t += kWarps) {
            const int j = part + stride * t;
            if (j >= n) break;
            double* col = A + static_cast<int64_t>(j) * lda;
            double s = lane == 0 ? ldg2(col + i) : 0.0;
            for (int r = i + 1 + lane; r < m; r += 32) s += ldg2(v + r) * ldg2(col + r);
            const double w = tau * warp_sum(s);
            if (w != 0.0) {
                if (lane == 0) col[i] = ldg2(col + i) - w;
                for (int r = i + 1 + lane; r < m; r += 32) col[r] = ldg2(col + r) - w * ldg2(v + r);
            }
            __syncwarp();
            const double n1 = ldg2(vn1 + j);
            if (n1 != 0.0) {
                double temp = fabs(ldg2(col + i)) / n1;
                temp = fmax(1.0 - temp * temp, 0.0);
                const double ratio = n1 / ldg2(vn2 + j);
                if (temp * ratio * ratio <= tol3z) {
                    double s2 = 0.0;
                    for (int r = i + 1 + lane; r < m; r += 32) {
                        const double x = ldg2(col + r);
                        s2 += x * x;
                    }
                    const double fresh = i < m - 1 ? sqrt(warp_sum(s2)) : 0.0;
                    if (lane == 0) vn1[j] = vn2[j] = fresh;
                } else if (lane == 0) {
                    vn1[j] = n1 * sqrt(temp);
                }
            }
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
            if (j >= 1 && m > 0) A[static_cast<int64_t>(j) * lda] = 0.0;
            jpvt[j] = j;  // rank 0 keeps column 0 (the host's convention)
        }
    } else {
        for (int t = first_owned(rank, part, stride) + warp;; t += kWarps) {
            const int j = part + stride * t;
            if (j >= n) break;
            double* col = A + static_cast<int64_t>(j) * lda;
            for (int p = rank - 1; p >= 0; --p) {
                double s = 0.0;
                for (int q = p + 1 + lane; q < rank; q += 32) s += ldg2(A + p + static_cast<int64_t>(q) * lda) * ldg2(col + q);
                s = warp_sum(s);
                const double x = (ldg2(col + p) - s) / ldg2(A + p + static_cast<int64_t>(p) * lda);
                __syncwarp();
                if (lane == 0) col[p] = x;
                __syncwarp();
                if (!isfinite(x)) t_bad = 1.0;
            }
        }
    }
    t_bad = block_max<kWarps>(t_bad, scratch);
    if (tid == 0) {
        if (part == 0) *item.rank = rank;
        if (t_bad != 0.0) atomicMax(item.flag, 2);
    }
}

// Blocks per box of the cooperative path (0: one block per box instead).
int coop_blocks_per_box(int count, int max_n) {
    static const int sms = [] {
        int device = 0, value = 0;
        cudaGetDevice(&device);
        cudaDeviceGetAttribute(&value, cudaDevAttrMultiProcessorCount, device);
        return value;
    }();
    static const int resident = [] {
        int per_sm = 0;
        cudaOccupancyMaxActiveBlocksPerMultiprocessor(&per_sm, qrcp_coop_kernel, kCoopThreads, 0);
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

template<int kThreads>
void launch_qrcp_with(const QrcpItem* items, int count, size_t shared, double tol, cudaStream_t stream) {
    if (shared > 48 * 1024) {
        const cudaError_t attr = cudaFuncSetAttribute(qrcp_kernel<kThreads>, cudaFuncAttributeMaxDynamicSharedMemorySize,
                                                      static_cast<int>(shared));
        if (attr != cudaSuccess) throw std::runtime_error("launch_qrcp: box too large for shared memory");
    }
    qrcp_kernel<kThreads><<<static_cast<unsigned>(count), kThreads, shared, stream>>>(items, tol);
}

}  // namespace

size_t qrcp_work_bytes(int count, int max_n) {
    const int cpb = coop_blocks_per_box(count, max_n);
    if (cpb == 0) return 0;
    return static_cast<size_t>(count) * coop_box_doubles(max_n, cpb) * sizeof(double) +
           (static_cast<size_t>(count) * coop_box_ints(cpb) + 1) * sizeof(int);
}

// Few boxes (at most 2 x multiprocessors / kMinBlocksPerBox): the cooperative
// path, several blocks per box.  Otherwise one block per box: large blocks for
// few boxes, small ones (several per multiprocessor) for many.
void launch_qrcp(const QrcpItem* items, int count, int max_n, double tol, void* work, cudaStream_t stream) {
    if (count <= 0) return;
    int cpb = coop_blocks_per_box(count, max_n);
    if (cpb > 0) {
        if (work == nullptr) throw std::runtime_error("launch_qrcp: the cooperative path needs its work block");
        double* dwork = static_cast<double*>(work);
        int* iwork = reinterpret_cast<int*>(dwork + static_cast<size_t>(count) * coop_box_doubles(max_n, cpb));
        void* args[] = {const_cast<QrcpItem**>(&items), &count, &cpb, &max_n, &tol, &dwork, &iwork};
        const cudaError_t status = cudaLaunchCooperativeKernel(reinterpret_cast<const void*>(qrcp_coop_kernel),
                                                               dim3(static_cast<unsigned>(count * cpb)),
                                                               dim3(kCoopThreads), args, 0, stream);
        if (status != cudaSuccess) {
            throw std::runtime_error(std::string("CUDA launch failed in qrcp_coop_kernel: ") + cudaGetErrorString(status));
        }
        return;
    }
    const size_t shared = static_cast<size_t>(2 * max_n + kMaxWarps) * sizeof(double) +
                          static_cast<size_t>(kMaxWarps + max_n) * sizeof(int);
    if (count < 216) {
        launch_qrcp_with<1024>(items, count, shared, tol, stream);
    } else {
        launch_qrcp_with<256>(items, count, shared, tol, stream);
    }
    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA launch failed in qrcp_kernel: ") + cudaGetErrorString(status));
    }
}

}  // namespace gpu
}  // namespace fmm
