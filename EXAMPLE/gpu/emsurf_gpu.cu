// GPU evaluators of EMSURF_Driver.cpp (gpu_evaluators.h): the entry Z(m, n)
// of the RWG basis functions of the edges m and n of a triangle mesh, as
// Zelem_EMSURF (EMSURF_Module.f90, e^{-jkr} convention), in the host's order
// of operations (Gauss points, singular self terms by ianalytic and
// ianalytic2).
//
// Both are block evaluators.  An RWG edge lives on two triangles, so an
// entry is a sum over 4 triangle pairs, and a triangle pair serves the up to
// 9 edge pairs on it.  For a block of row edges by column edges, the
// quadrature of every pair of their triangles runs once, into sums per pair
// that do not depend on the edges; the entries then combine the sums of
// their 4 pairs with the edges' geometry.  The block's distinct triangles are
// found on the host (BPACK_GPU_HOST_IDS).
//
//  - EFIE (cfie_alpha == 1): 9 sums a pair; symmetric.
//  - CFIE: alpha Z_EFIE + (1 - alpha) eta0 Z_MFIE, the EFIE's 9 sums and 27
//    of the MFIE a pair; not symmetric.
//
// The mesh tables (emsurf_get_gpu_tables_c, EMSURF_C_Bindings.f90):
//   params: wavenumber k, frequency, eps0, pi, Gauss points per triangle (at
//           most kMaxGauss), vertices, edges, alpha
//   reals:  vertex xyz (3 per vertex), the Gauss rule (ng1, ng2, ng3, w per
//           point); CFIE: then eta0 and the unit normal of each triangle
//   ints:   per edge (6) its two vertices, its two triangles and their
//           vertices opposite the edge; then per triangle (3) its vertices
//           (all 0-based; a missing triangle is -1)
#include "gpu_evaluators.h"

#include "bpack_gpu_entry.cuh"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {

using bpack::gpu::dcomplex;

constexpr int kMaxGauss = 7;

struct Mesh {
    const double* xyz;
    const double* gauss;    // ng1, ng2, ng3, w per point
    const double* normals;  // CFIE: 3 per triangle
    const int* edge;
    const int* tri;
    int nq;
    double k, freq, eps0, pi, alpha, eta0;
};

__device__ __forceinline__ const double* vertex(const Mesh& m, int v) { return m.xyz + 3 * static_cast<int64_t>(v); }

// c = a x b (rrcurl)
__device__ __forceinline__ void cross(const double* a, const double* b, double* c) {
    c[0] = a[1] * b[2] - b[1] * a[2];
    c[1] = -a[0] * b[2] + b[0] * a[2];
    c[2] = a[0] * b[1] - b[0] * a[1];
}
__device__ __forceinline__ double dot(const double* a, const double* b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
__device__ __forceinline__ double dist(const double* a, const double* b) {
    const double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
    return sqrt(dx * dx + dy * dy + dz * dz);
}

// Gauss points of triangle t (gau_grobal): ng1 x1 + ng2 x2 + ng3 x3,
// rounded as the host forms them (no fused multiply-adds)
__device__ __forceinline__ void gauss_points(const Mesh& m, int t, double (&p)[3][kMaxGauss]) {
    const int* v = m.tri + 3 * static_cast<int64_t>(t);
    const double* a = vertex(m, v[0]);
    const double* b = vertex(m, v[1]);
    const double* c = vertex(m, v[2]);
    for (int q = 0; q < m.nq; ++q) {
        const double g1 = m.gauss[4 * q], g2 = m.gauss[4 * q + 1], g3 = m.gauss[4 * q + 2];
        for (int d = 0; d < 3; ++d) {
            p[d][q] = __dadd_rn(__dadd_rn(__dmul_rn(g1, a[d]), __dmul_rn(g2, b[d])), __dmul_rn(g3, c[d]));
        }
    }
}

// ianalytic, ianalytic2(1) and ianalytic2(2) of triangle t at the point x:
// the singular part of the integrals of 1 / R, and of ng1 / R and ng2 / R,
// over the triangle
__device__ __noinline__ void analytic(const Mesh& m, int t, const double* x, double& i0, double& i1, double& i2) {
    const int* tv = m.tri + 3 * static_cast<int64_t>(t);
    const double* n1 = vertex(m, tv[0]);
    const double* n2 = vertex(m, tv[1]);
    const double* n3 = vertex(m, tv[2]);
    double a[3], b[3], w[3], u[3], v[3], l[3];
    for (int i = 0; i < 3; ++i) {
        a[i] = n2[i] - n1[i];
        b[i] = n3[i] - n1[i];
    }
    cross(a, b, w);
    const double area = 0.5 * sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
    for (int i = 0; i < 3; ++i) w[i] = w[i] / 2. / area;
    l[0] = dist(n3, n2);
    l[1] = dist(n3, n1);
    l[2] = dist(n1, n2);
    for (int i = 0; i < 3; ++i) u[i] = a[i] / l[2];
    cross(w, u, v);
    const double u3 = dot(u, b);
    const double v3 = 2. * area / l[2];
    for (int i = 0; i < 3; ++i) b[i] = x[i] - n1[i];
    const double u0 = dot(u, b), v0 = dot(v, b), w0 = dot(w, b);

    double s1[3], s2[3], tm[3], t0[3], tp[3];
    s1[0] = -((l[2] - u0) * (l[2] - u3) + v0 * v3) / l[0];
    s2[0] = ((u3 - u0) * (u3 - l[2]) + v3 * (v3 - v0)) / l[0];
    s1[1] = -(u3 * (u3 - u0) + v3 * (v3 - v0)) / l[1];
    s2[1] = (u0 * u3 + v0 * v3) / l[1];
    s1[2] = -u0;
    s2[2] = l[2] - u0;
    t0[0] = (v0 * (u3 - l[2]) + v3 * (l[2] - u0)) / l[0];
    t0[1] = (u0 * v3 - v0 * u3) / l[1];
    t0[2] = v0;
    tm[0] = sqrt((l[2] - u0) * (l[2] - u0) + v0 * v0);
    tp[0] = sqrt((u3 - u0) * (u3 - u0) + (v3 - v0) * (v3 - v0));
    tp[1] = sqrt(u0 * u0 + v0 * v0);
    tp[2] = tm[0];
    tm[1] = tp[0];
    tm[2] = tp[1];
    double rm[3], r0[3], rp[3], f2[3];
    const double aw0 = fabs(w0);
    double temp = 0.;
    for (int i = 0; i < 3; ++i) {
        rm[i] = sqrt(tm[i] * tm[i] + w0 * w0);
        r0[i] = sqrt(t0[i] * t0[i] + w0 * w0);
        rp[i] = sqrt(tp[i] * tp[i] + w0 * w0);
    }
    for (int i = 0; i < 3; ++i) {
        f2[i] = log((rp[i] + s2[i]) / (rm[i] + s1[i]));
        const double beta = atan(t0[i] * s2[i] / (r0[i] * r0[i] + aw0 * rp[i])) -
                            atan(t0[i] * s1[i] / (r0[i] * r0[i] + aw0 * rm[i]));
        temp = temp + t0[i] * f2[i] - aw0 * beta;
    }
    i0 = temp / area / 2;

    double f3[3], e1[3], e2[3], e3[3], m1[3], m2[3], m3[3];
    for (int i = 0; i < 3; ++i) f3[i] = s2[i] * rp[i] - s1[i] * rm[i] + r0[i] * r0[i] * f2[i];
    for (int i = 0; i < 3; ++i) {
        e1[i] = (n3[i] - n2[i]) / l[0];
        e2[i] = (n1[i] - n3[i]) / l[1];
        e3[i] = (n2[i] - n1[i]) / l[2];
    }
    cross(e1, w, m1);
    cross(e2, w, m2);
    cross(e3, w, m3);
    const double iua = dot(u, m1) * f3[0] / 2 + dot(u, m2) * f3[1] / 2 + dot(u, m3) * f3[2] / 2;
    const double iva = dot(v, m1) * f3[0] / 2 + dot(v, m2) * f3[1] / 2 + dot(v, m3) * f3[2] / 2;
    const double n01 = 1. - u0 / l[2] + v0 * (u3 / l[2] - 1.) / v3;
    const double n02 = u0 / l[2] - u3 * v0 / l[2] / v3;
    i1 = (n01 * temp - iua / l[2] + iva * (u3 / l[2] - 1.) / v3) / 2 / area;
    i2 = (n02 * temp + iua / l[2] - iva * u3 / l[2] / v3) / 2 / area;
}

// lm ln j z / 8 / pi^2 / freq / eps0: the EFIE entry from the sum z of its
// triangle pairs
__device__ __forceinline__ dcomplex efie_scale(const Mesh& m, double lengths, dcomplex z) {
    const double pi2 = m.pi * m.pi;
    double re = -(lengths * z.im), im = lengths * z.re;
    re = re / 8. / pi2 / m.freq / m.eps0;
    im = im / 8. / pi2 / m.freq / m.eps0;
    return dcomplex(re, im);
}

// ---- EFIE and CFIE: block evaluators by triangle pairs
//
// The EFIE entry of test edge m and source edge n is
//   Z = lm ln j (c1 - c2) / 8 / pi^2 / freq / eps0,
//   c1 - c2 = sum over the triangles T of m, S of n (signs sT, sS) of
//             sT sS sum_i w_i [k^2 (x_i - o_m) . (P_i - o_n I0_i) - 4 I0_i],
// with, per test point x_i of T, the host's sums over the source points of S
//   I0_i = sum_j w_j g_ij, I1_i = sum_j b1_j w_j g_ij, I2_i = sum_j b2_j w_j g_ij
// (imp, imp1, imp2; g = e^{-jkR}/R, analytic singular terms when T == S),
// P_i = Y1 I1_i + Y2 I2_i + Y3 (I0_i - I1_i - I2_i) and o_m, o_n the edges'
// vertices opposite them.  Relative to the third vertices (x_i - X3 =
// a1_i f1 + a2_i f2 with f1 = X1 - X3, f2 = X2 - X3 on T, e1 = Y1 - Y3, e2 =
// Y2 - Y3 on S; u = o_m - X3, v = o_n - Y3):
//   sum_i w_i [...] = k^2 (D - v.B - u.C + (u.v) A) - 4 A,
//   A = sum w_i I0_i, B = f1 H1 + f2 H2, C = e1 J1 + e2 J2,
//   D = (f1.e1) K11 + (f1.e2) K12 + (f2.e1) K21 + (f2.e2) K22,
// where the 9 sums of the pair, independent of the edges, are
//   A, H1 = sum w a1 I0, H2 = sum w a2 I0, J1 = sum w I1, J2 = sum w I2,
//   K11 = sum w a1 I1, K12 = sum w a1 I2, K21 = sum w a2 I1, K22 = sum w a2 I2.
//
// The MFIE part value_m of Zelem_EMSURF (times lm ln) on the pair (T, S) is,
// with a_m = x_i - o_m, a_n = y_j - o_n and n_T the unit normal of T,
//   T == S:  sT sS sum_i w_i (a_m . (x_i - o_n)) / (4 area)   (the identity term)
//   T != S: -sT sS sum_ij w_i w_j (n_T x (a_n x dG_ij)) . a_m,
//            dG_ij = (x_i - y_j) (1 + j k d) e^{-jkd} / (4 pi d^3).
// As (n x (a x g)) . b = (b . a)(n . g) - (b . g)(n . a), and a_m = sum_p
// alpha_ip F_p, a_n = sum_q beta_jq E_q with alpha_i = (a1_i, a2_i, 1), F =
// (f1, f2, -u), beta_j = (b1_j, b2_j, 1) (the source point's ng1, ng2), E =
// (e1, e2, -v), the T != S term is
//   -sT sS sum_pq [(F_p . E_q)(n_T . G_pq) - (n_T . E_q)(F_p . G_pq)],
// where the 27 sums of the pair (9 vectors), independent of the edges, are
//   G_pq = sum_ij w_i w_j alpha_ip beta_jq dG_ij.

constexpr int kPairSums = 9;
constexpr int kMfieSums = 27;  // G_pq[d] at 9 p + 3 q + d

// One block of a batch on the device: its edges (the batch's device ids),
// each edge's two triangles as indices into the block's distinct triangles
// (-1: none), those triangles, and the scratch of the pairs' sums.
struct EfieItem {
    int nrow, ncol, ntr, ntc;
    const int64_t* row_edges;
    const int64_t* col_edges;
    const int* row_tri;
    const int* col_tri;
    const int* tr;
    const int* tc;
    dcomplex* sums;   // ntr * ntc * kPairSums
    dcomplex* msums;  // CFIE: ntr * ntc * kMfieSums (EFIE: null)
    dcomplex* out;    // out[r + c * ld]
    int ld;
};

// The 9 sums of each triangle pair (test triangle of the rows, source
// triangle of the columns): one thread per pair.
__global__ void __launch_bounds__(128) efie_pairs_kernel(const EfieItem* items, Mesh m) {
    const EfieItem item = items[blockIdx.y];
    const int64_t pair = static_cast<int64_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (pair >= static_cast<int64_t>(item.ntr) * item.ntc) return;
    const int t = item.tr[pair / item.ntc];
    const int s = item.tc[pair % item.ntc];
    const double k = m.k;
    const int nq = m.nq;
    double y[3][kMaxGauss];
    gauss_points(m, s, y);
    const int* tv = m.tri + 3 * static_cast<int64_t>(t);
    const double* x1 = vertex(m, tv[0]);
    const double* x2 = vertex(m, tv[1]);
    const double* x3 = vertex(m, tv[2]);
    dcomplex sum[kPairSums];
#pragma unroll
    for (int q = 0; q < kPairSums; ++q) sum[q] = dcomplex(0.0, 0.0);
    for (int i = 0; i < nq; ++i) {
        const double a1 = m.gauss[4 * i], a2 = m.gauss[4 * i + 1], a3 = m.gauss[4 * i + 2], wi = m.gauss[4 * i + 3];
        double xi[3];
        for (int d = 0; d < 3; ++d) {
            xi[d] = __dadd_rn(__dadd_rn(__dmul_rn(a1, x1[d]), __dmul_rn(a2, x2[d])), __dmul_rn(a3, x3[d]));
        }
        dcomplex i0(0.0, 0.0), i1(0.0, 0.0), i2(0.0, 0.0);
#pragma unroll
        for (int j = 0; j < kMaxGauss; ++j) {
            if (j >= nq) break;
            const double yj[3] = {y[0][j], y[1][j], y[2][j]};
            const double d = dist(xi, yj);
            const double wj = m.gauss[4 * j + 3];
            const double w1 = m.gauss[4 * j] * wj, w2 = m.gauss[4 * j + 1] * wj;
            if (t == s && d == 0.0) {
                // the same triangle, i == j: the host's analytic singular part
                double s0, s1, s2;
                analytic(m, s, yj, s0, s1, s2);
                i0 += dcomplex(s0, wj * -k);
                i1 += dcomplex(s1, w1 * -k);
                i2 += dcomplex(s2, w2 * -k);
            } else {
                double sn, cs;
                sincos(k * d, &sn, &cs);
                const double inv = 1.0 / d;
                const dcomplex g(cs * inv, -sn * inv);  // e^{-jkd} / d
                i0 += wj * g;
                i1 += w1 * g;
                i2 += w2 * g;
            }
        }
        const double wa1 = wi * a1, wa2 = wi * a2;
        sum[0] += wi * i0;
        sum[1] += wa1 * i0;
        sum[2] += wa2 * i0;
        sum[3] += wi * i1;
        sum[4] += wi * i2;
        sum[5] += wa1 * i1;
        sum[6] += wa1 * i2;
        sum[7] += wa2 * i1;
        sum[8] += wa2 * i2;
    }
    dcomplex* out = item.sums + pair * kPairSums;
#pragma unroll
    for (int q = 0; q < kPairSums; ++q) out[q] = sum[q];
}

__device__ __forceinline__ void sub3(const double* a, const double* b, double* c) {
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
}

// The entries from the sums of their 4 triangle pairs: one thread per
// entry, rows fastest.
__global__ void __launch_bounds__(256) efie_entries_kernel(const EfieItem* items, Mesh m) {
    const EfieItem item = items[blockIdx.y];
    const int64_t e = static_cast<int64_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (e >= static_cast<int64_t>(item.nrow) * item.ncol) return;
    const int r = static_cast<int>(e % item.nrow);
    const int c = static_cast<int>(e / item.nrow);
    const double k2 = m.k * m.k;
    const int* Er = m.edge + 6 * item.row_edges[r];
    const int* Ec = m.edge + 6 * item.col_edges[c];
    const int* row_tri = item.row_tri + 2 * static_cast<int64_t>(r);
    const int* col_tri = item.col_tri + 2 * static_cast<int64_t>(c);
    dcomplex total(0.0, 0.0);
    for (int a = 0; a < 2; ++a) {
        if (row_tri[a] < 0) continue;
        const int* tv = m.tri + 3 * static_cast<int64_t>(item.tr[row_tri[a]]);
        const double* x3 = vertex(m, tv[2]);
        double f1[3], f2[3], u[3];
        sub3(vertex(m, tv[0]), x3, f1);
        sub3(vertex(m, tv[1]), x3, f2);
        sub3(vertex(m, Er[4 + a]), x3, u);
        for (int b = 0; b < 2; ++b) {
            if (col_tri[b] < 0) continue;
            const int* sv = m.tri + 3 * static_cast<int64_t>(item.tc[col_tri[b]]);
            const double* y3 = vertex(m, sv[2]);
            double e1[3], e2[3], v[3];
            sub3(vertex(m, sv[0]), y3, e1);
            sub3(vertex(m, sv[1]), y3, e2);
            sub3(vertex(m, Ec[4 + b]), y3, v);
            const dcomplex* M = item.sums + (static_cast<int64_t>(row_tri[a]) * item.ntc + col_tri[b]) * kPairSums;
            dcomplex vB(0.0, 0.0), uC(0.0, 0.0);
            for (int d = 0; d < 3; ++d) {
                vB += v[d] * (f1[d] * M[1] + f2[d] * M[2]);
                uC += u[d] * (e1[d] * M[3] + e2[d] * M[4]);
            }
            const dcomplex D = dot(f1, e1) * M[5] + dot(f1, e2) * M[6] + dot(f2, e1) * M[7] + dot(f2, e2) * M[8];
            const dcomplex F = k2 * (D - vB - uC + dot(u, v) * M[0]) - 4. * M[0];
            total += ((a == 0 ? 1.0 : -1.0) * (b == 0 ? 1.0 : -1.0)) * F;
        }
    }
    const double lengths = dist(vertex(m, Ec[0]), vertex(m, Ec[1])) * dist(vertex(m, Er[0]), vertex(m, Er[1]));
    item.out[r + static_cast<int64_t>(c) * item.ld] = efie_scale(m, lengths, total);
}

// The MFIE sums G_pq of each pair of distinct triangles: one thread per pair
// (a triangle with itself has only the identity term, which the entries take
// from its geometry).
__global__ void __launch_bounds__(128) mfie_pairs_kernel(const EfieItem* items, Mesh m) {
    const EfieItem item = items[blockIdx.y];
    const int64_t pair = static_cast<int64_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (pair >= static_cast<int64_t>(item.ntr) * item.ntc) return;
    const int t = item.tr[pair / item.ntc];
    const int s = item.tc[pair % item.ntc];
    if (t == s) return;
    const double k = m.k;
    const int nq = m.nq;
    double y[3][kMaxGauss];
    gauss_points(m, s, y);
    const int* tv = m.tri + 3 * static_cast<int64_t>(t);
    const double* x1 = vertex(m, tv[0]);
    const double* x2 = vertex(m, tv[1]);
    const double* x3 = vertex(m, tv[2]);
    dcomplex g[kMfieSums];
#pragma unroll
    for (int q = 0; q < kMfieSums; ++q) g[q] = dcomplex(0.0, 0.0);
    for (int i = 0; i < nq; ++i) {
        const double a1 = m.gauss[4 * i], a2 = m.gauss[4 * i + 1], a3 = m.gauss[4 * i + 2], wi = m.gauss[4 * i + 3];
        double xi[3];
        for (int d = 0; d < 3; ++d) {
            xi[d] = __dadd_rn(__dadd_rn(__dmul_rn(a1, x1[d]), __dmul_rn(a2, x2[d])), __dmul_rn(a3, x3[d]));
        }
        dcomplex h[9];  // sum_j w_j beta_jq dG_ij, at 3 q + d
#pragma unroll
        for (int q = 0; q < 9; ++q) h[q] = dcomplex(0.0, 0.0);
#pragma unroll
        for (int j = 0; j < kMaxGauss; ++j) {
            if (j >= nq) break;
            const double r[3] = {xi[0] - y[0][j], xi[1] - y[1][j], xi[2] - y[2][j]};
            const double d = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
            double sn, cs;
            sincos(k * d, &sn, &cs);
            const double kd = k * d;
            const double den = 4 * m.pi * (d * d * d);
            const dcomplex c((cs + kd * sn) / den, (kd * cs - sn) / den);  // (1 + j k d) e^{-jkd} / (4 pi d^3)
            const double wj = m.gauss[4 * j + 3];
            const double wb[3] = {m.gauss[4 * j] * wj, m.gauss[4 * j + 1] * wj, wj};
#pragma unroll
            for (int dd = 0; dd < 3; ++dd) {
                const dcomplex rc = r[dd] * c;
#pragma unroll
                for (int q = 0; q < 3; ++q) h[3 * q + dd] += wb[q] * rc;
            }
        }
        const double wa[3] = {wi * a1, wi * a2, wi};
#pragma unroll
        for (int p = 0; p < 3; ++p) {
#pragma unroll
            for (int q = 0; q < 9; ++q) g[9 * p + q] += wa[p] * h[q];
        }
    }
    dcomplex* out = item.msums + pair * kMfieSums;
#pragma unroll
    for (int q = 0; q < kMfieSums; ++q) out[q] = g[q];
}

// The MFIE identity term of the test and source functions of a triangle t
// (vertices om, on opposite their edges): sum_i w_i 0.5 (x_i - om) . (x_i -
// on) / (2 area), as the host forms it
__device__ __noinline__ double mfie_identity(const Mesh& m, int t, const double* om, const double* on) {
    const int* tv = m.tri + 3 * static_cast<int64_t>(t);
    const double* v1 = vertex(m, tv[0]);
    double a[3], b[3], c[3];
    sub3(vertex(m, tv[1]), v1, a);
    sub3(vertex(m, tv[2]), v1, b);
    cross(a, b, c);
    const double area = 0.5 * sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
    double p[3][kMaxGauss];
    gauss_points(m, t, p);
    double sum = 0.0;
    for (int i = 0; i < m.nq; ++i) {
        const double am[3] = {p[0][i] - om[0], p[1][i] - om[1], p[2][i] - om[2]};
        const double an[3] = {p[0][i] - on[0], p[1][i] - on[1], p[2][i] - on[2]};
        sum += 0.5 * dot(am, an) / (2. * area) * m.gauss[4 * i + 3];
    }
    return sum;
}

// The CFIE entries alpha Z_EFIE + (1 - alpha) eta0 Z_MFIE from the sums of
// their 4 triangle pairs: one thread per entry, rows fastest.
__global__ void __launch_bounds__(256) cfie_entries_kernel(const EfieItem* items, Mesh m) {
    const EfieItem item = items[blockIdx.y];
    const int64_t e = static_cast<int64_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (e >= static_cast<int64_t>(item.nrow) * item.ncol) return;
    const int r = static_cast<int>(e % item.nrow);
    const int c = static_cast<int>(e / item.nrow);
    const double k2 = m.k * m.k;
    const int* Er = m.edge + 6 * item.row_edges[r];
    const int* Ec = m.edge + 6 * item.col_edges[c];
    const int* row_tri = item.row_tri + 2 * static_cast<int64_t>(r);
    const int* col_tri = item.col_tri + 2 * static_cast<int64_t>(c);
    dcomplex efie(0.0, 0.0), mfie(0.0, 0.0);
    for (int a = 0; a < 2; ++a) {
        if (row_tri[a] < 0) continue;
        const int t = item.tr[row_tri[a]];
        const int* tv = m.tri + 3 * static_cast<int64_t>(t);
        const double* x3 = vertex(m, tv[2]);
        double F[3][3];  // f1, f2, -u
        sub3(vertex(m, tv[0]), x3, F[0]);
        sub3(vertex(m, tv[1]), x3, F[1]);
        sub3(x3, vertex(m, Er[4 + a]), F[2]);
        const double* nr = m.normals + 3 * static_cast<int64_t>(t);
        for (int b = 0; b < 2; ++b) {
            if (col_tri[b] < 0) continue;
            const int s = item.tc[col_tri[b]];
            const int* sv = m.tri + 3 * static_cast<int64_t>(s);
            const double* y3 = vertex(m, sv[2]);
            double E[3][3];  // e1, e2, -v
            sub3(vertex(m, sv[0]), y3, E[0]);
            sub3(vertex(m, sv[1]), y3, E[1]);
            sub3(y3, vertex(m, Ec[4 + b]), E[2]);
            const double sign = (a == 0 ? 1.0 : -1.0) * (b == 0 ? 1.0 : -1.0);
            const int64_t pair = static_cast<int64_t>(row_tri[a]) * item.ntc + col_tri[b];
            // EFIE (efie_entries_kernel; u = -F[2], v = -E[2])
            const dcomplex* M = item.sums + pair * kPairSums;
            dcomplex vB(0.0, 0.0), uC(0.0, 0.0);
            for (int d = 0; d < 3; ++d) {
                vB += -E[2][d] * (F[0][d] * M[1] + F[1][d] * M[2]);
                uC += -F[2][d] * (E[0][d] * M[3] + E[1][d] * M[4]);
            }
            const dcomplex D = dot(F[0], E[0]) * M[5] + dot(F[0], E[1]) * M[6] + dot(F[1], E[0]) * M[7] +
                               dot(F[1], E[1]) * M[8];
            efie += sign * (k2 * (D - vB - uC + dot(F[2], E[2]) * M[0]) - 4. * M[0]);
            // MFIE
            if (t == s) {
                mfie.re += sign * mfie_identity(m, t, vertex(m, Er[4 + a]), vertex(m, Ec[4 + b]));
                continue;
            }
            const dcomplex* G = item.msums + pair * kMfieSums;
            dcomplex sum(0.0, 0.0);
            for (int p = 0; p < 3; ++p) {
                for (int q = 0; q < 3; ++q) {
                    const dcomplex* g = G + 9 * p + 3 * q;
                    const dcomplex ng = nr[0] * g[0] + nr[1] * g[1] + nr[2] * g[2];
                    const dcomplex fg = F[p][0] * g[0] + F[p][1] * g[1] + F[p][2] * g[2];
                    sum += dot(F[p], E[q]) * ng - dot(nr, E[q]) * fg;
                }
            }
            mfie -= sign * sum;
        }
    }
    const double lm = dist(vertex(m, Er[0]), vertex(m, Er[1]));
    const double ln = dist(vertex(m, Ec[0]), vertex(m, Ec[1]));
    const dcomplex ze = efie_scale(m, ln * lm, efie);
    const dcomplex zm = (mfie * lm) * ln;
    item.out[r + static_cast<int64_t>(c) * item.ld] = m.alpha * ze + ((1. - m.alpha) * m.eta0) * zm;
}

// ---- the application's state: the mesh tables on the device, and the
// edge table on the host for the blocks' triangles

struct Emsurf {
    Mesh mesh{};
    double* reals = nullptr;
    int* ints = nullptr;
    std::vector<int> edges;  // host copy: 6 per edge
};
Emsurf* state = nullptr;

void fail(const char* what, const char* why) {
    std::fprintf(stderr, "emsurf_gpu: %s: %s\n", what, why);
    std::abort();
}
void check(cudaError_t status, const char* what) {
    if (status != cudaSuccess) fail(what, cudaGetErrorString(status));
}

size_t align(size_t bytes) { return (bytes + 255) / 256 * 256; }

// The distinct triangles of n edges (sorted), and each edge's two as
// indices into them (-1: none)
void triangles_of(const std::vector<int>& edge_table, const int64_t* edges, int n, std::vector<int>& tri_of,
                  std::vector<int>& tris) {
    tris.clear();
    for (int i = 0; i < n; ++i) {
        for (int a = 0; a < 2; ++a) {
            const int t = edge_table[6 * static_cast<size_t>(edges[i]) + 2 + a];
            if (t >= 0) tris.push_back(t);
        }
    }
    std::sort(tris.begin(), tris.end());
    tris.erase(std::unique(tris.begin(), tris.end()), tris.end());
    tri_of.resize(2 * static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        for (int a = 0; a < 2; ++a) {
            const int t = edge_table[6 * static_cast<size_t>(edges[i]) + 2 + a];
            tri_of[2 * static_cast<size_t>(i) + a] =
                t < 0 ? -1 : static_cast<int>(std::lower_bound(tris.begin(), tris.end(), t) - tris.begin());
        }
    }
}

// Blocks [b0, b1) of `blocks`: their triangle lists and items go to the
// device in one copy, with the scratch of their pair sums behind them (the
// CFIE's MFIE sums after the EFIE's).  A group that does not fit in the
// device memory the library can give is split in two, and a single block by
// its columns (an entry reads the sums of its own pairs only).
void pair_group(const bpack_gpu_block_batch* batch, const Emsurf& em, const std::vector<bpack_gpu_block>& blocks,
                int b0, int b1, bool cfie) {
    struct Lists {
        std::vector<int> row_tri, col_tri, tr, tc;
    };
    std::vector<Lists> lists(static_cast<size_t>(b1 - b0));
    size_t ints = 0, sums = 0;
    int64_t max_pairs = 0, max_entries = 0;
    for (int b = b0; b < b1; ++b) {
        const bpack_gpu_block& blk = blocks[static_cast<size_t>(b)];
        Lists& l = lists[static_cast<size_t>(b - b0)];
        triangles_of(em.edges, blk.row_ids_host, blk.m, l.row_tri, l.tr);
        triangles_of(em.edges, blk.col_ids_host, blk.n, l.col_tri, l.tc);
        ints += l.row_tri.size() + l.col_tri.size() + l.tr.size() + l.tc.size();
        const int64_t pairs = static_cast<int64_t>(l.tr.size()) * static_cast<int64_t>(l.tc.size());
        sums += align(static_cast<size_t>(pairs) * kPairSums * sizeof(dcomplex));
        if (cfie) sums += align(static_cast<size_t>(pairs) * kMfieSums * sizeof(dcomplex));
        max_pairs = std::max(max_pairs, pairs);
        max_entries = std::max(max_entries, static_cast<int64_t>(blk.m) * blk.n);
    }
    const int count = b1 - b0;
    const size_t items_at = align(ints * sizeof(int));
    const size_t image_bytes = items_at + align(static_cast<size_t>(count) * sizeof(EfieItem));
    char* d = static_cast<char*>(batch->allocate(batch->allocator, image_bytes + sums));
    if (d == nullptr) {
        if (count == 1) {
            const bpack_gpu_block& blk = blocks[static_cast<size_t>(b0)];
            if (blk.n < 2) fail("EMSURF block", "out of device memory for one column's triangle pairs");
            std::vector<bpack_gpu_block> halves(2, blk);
            const int half = blk.n / 2;
            halves[0].n = half;
            halves[1].n = blk.n - half;
            halves[1].out = static_cast<dcomplex*>(blk.out) + static_cast<int64_t>(half) * blk.ld;
            halves[1].col_ids = blk.col_ids + half;
            halves[1].col_ids_host = blk.col_ids_host + half;
            if (blk.col_coords != nullptr) halves[1].col_coords = blk.col_coords + static_cast<int64_t>(half) * batch->dim;
            pair_group(batch, em, halves, 0, 1, cfie);
            pair_group(batch, em, halves, 1, 2, cfie);
            return;
        }
        const int mid = b0 + count / 2;
        pair_group(batch, em, blocks, b0, mid, cfie);
        pair_group(batch, em, blocks, mid, b1, cfie);
        return;
    }
    // the image: the int lists, then the items (device addresses)
    std::vector<char> image(image_bytes);
    int* h_ints = reinterpret_cast<int*>(image.data());
    EfieItem* h_items = reinterpret_cast<EfieItem*>(image.data() + items_at);
    size_t at = 0, sums_at = image_bytes;
    auto put = [&](const std::vector<int>& v) {
        std::copy(v.begin(), v.end(), h_ints + at);
        const int* p = reinterpret_cast<const int*>(d) + at;
        at += v.size();
        return p;
    };
    for (int b = b0; b < b1; ++b) {
        const bpack_gpu_block& blk = blocks[static_cast<size_t>(b)];
        const Lists& l = lists[static_cast<size_t>(b - b0)];
        EfieItem& it = h_items[b - b0];
        it.nrow = blk.m;
        it.ncol = blk.n;
        it.ntr = static_cast<int>(l.tr.size());
        it.ntc = static_cast<int>(l.tc.size());
        it.row_edges = blk.row_ids;
        it.col_edges = blk.col_ids;
        it.row_tri = put(l.row_tri);
        it.col_tri = put(l.col_tri);
        it.tr = put(l.tr);
        it.tc = put(l.tc);
        it.sums = reinterpret_cast<dcomplex*>(d + sums_at);
        sums_at += align(static_cast<size_t>(it.ntr) * it.ntc * kPairSums * sizeof(dcomplex));
        it.msums = nullptr;
        if (cfie) {
            it.msums = reinterpret_cast<dcomplex*>(d + sums_at);
            sums_at += align(static_cast<size_t>(it.ntr) * it.ntc * kMfieSums * sizeof(dcomplex));
        }
        it.out = static_cast<dcomplex*>(blk.out);
        it.ld = blk.ld;
    }
    cudaStream_t stream = static_cast<cudaStream_t>(batch->stream);
    // (pageable: the image is staged before the call returns)
    check(cudaMemcpyAsync(d, image.data(), image_bytes, cudaMemcpyHostToDevice, stream), "cudaMemcpyAsync");
    const EfieItem* d_items = reinterpret_cast<const EfieItem*>(d + items_at);
    for (int i = 0; i < count; i += 65535) {  // (grid.y)
        const unsigned n = static_cast<unsigned>(std::min(count - i, 65535));
        if (max_pairs > 0) {
            efie_pairs_kernel<<<dim3(static_cast<unsigned>((max_pairs + 127) / 128), n), 128, 0, stream>>>(
                d_items + i, em.mesh);
            if (cfie) {
                mfie_pairs_kernel<<<dim3(static_cast<unsigned>((max_pairs + 127) / 128), n), 128, 0, stream>>>(
                    d_items + i, em.mesh);
            }
        }
        if (cfie) {
            cfie_entries_kernel<<<dim3(static_cast<unsigned>((max_entries + 255) / 256), n), 256, 0, stream>>>(
                d_items + i, em.mesh);
        } else {
            efie_entries_kernel<<<dim3(static_cast<unsigned>((max_entries + 255) / 256), n), 256, 0, stream>>>(
                d_items + i, em.mesh);
        }
    }
    check(cudaGetLastError(), cfie ? "CFIE kernels" : "EFIE kernels");
    batch->release(batch->allocator, d);  // (later work on the stream is ordered after the kernels)
}

void efie_blocks(const bpack_gpu_block_batch* batch, void* user) {
    if (batch->count <= 0) return;
    const std::vector<bpack_gpu_block> blocks(batch->blocks, batch->blocks + batch->count);
    pair_group(batch, *static_cast<const Emsurf*>(user), blocks, 0, batch->count, false);
}
void cfie_blocks(const bpack_gpu_block_batch* batch, void* user) {
    if (batch->count <= 0) return;
    const std::vector<bpack_gpu_block> blocks(batch->blocks, batch->blocks + batch->count);
    pair_group(batch, *static_cast<const Emsurf*>(user), blocks, 0, batch->count, true);
}

}  // namespace

extern "C" void emsurf_gpu_register(void** bmat, const double* params, int nparams, const double* reals,
                                    int64_t nreals, const int* ints, int64_t nints, double cfie_alpha) {
    if (nparams < 8) fail("emsurf_gpu_register", "8 parameters expected");
    emsurf_gpu_release();
    state = new Emsurf;
    Emsurf& em = *state;
    check(cudaMalloc(reinterpret_cast<void**>(&em.reals), static_cast<size_t>(nreals) * sizeof(double)), "cudaMalloc");
    check(cudaMalloc(reinterpret_cast<void**>(&em.ints), static_cast<size_t>(nints) * sizeof(int)), "cudaMalloc");
    check(cudaMemcpy(em.reals, reals, static_cast<size_t>(nreals) * sizeof(double), cudaMemcpyHostToDevice),
          "cudaMemcpy");
    check(cudaMemcpy(em.ints, ints, static_cast<size_t>(nints) * sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy");
    const int nq = static_cast<int>(params[4]);
    const int64_t vertices = static_cast<int64_t>(params[5]);
    const int64_t edges = static_cast<int64_t>(params[6]);
    if (nq > kMaxGauss) fail("emsurf_gpu_register", "more than 7 Gauss points per triangle");
    Mesh& m = em.mesh;
    m.xyz = em.reals;
    m.gauss = em.reals + 3 * vertices;
    m.edge = em.ints;
    m.tri = em.ints + 6 * edges;
    m.nq = nq;
    m.k = params[0];
    m.freq = params[1];
    m.eps0 = params[2];
    m.pi = params[3];
    m.alpha = cfie_alpha;
    em.edges.assign(ints, ints + 6 * edges);
    if (cfie_alpha == 1.0) {
        const int flags = BPACK_GPU_SYMMETRIC | BPACK_GPU_HOST_IDS;
        z_c_bpack_set_gpu_block_evaluator(bmat, &efie_blocks, &em, &flags);
    } else {
        const int64_t eta0_at = 3 * vertices + 4 * static_cast<int64_t>(nq);
        m.eta0 = reals[eta0_at];
        m.normals = em.reals + eta0_at + 1;
        const int flags = BPACK_GPU_HOST_IDS;  // (not symmetric)
        z_c_bpack_set_gpu_block_evaluator(bmat, &cfie_blocks, &em, &flags);
    }
}

extern "C" void emsurf_gpu_release(void) {
    if (state == nullptr) return;
    cudaFree(state->reals);
    cudaFree(state->ints);
    delete state;
    state = nullptr;
}
