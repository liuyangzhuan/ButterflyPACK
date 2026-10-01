// GPU evaluators of EMSURF_Driver.cpp (gpu_evaluators.h): the entry Z(m, n)
// of the RWG basis functions of the edges m and n of a triangle mesh, as
// Zelem_EMSURF (EMSURF_Module.f90, e^{-jkr} convention), in the host's order
// of operations (Gauss points, singular self terms by ianalytic and
// ianalytic2).
//
//  - EFIE (cfie_alpha == 1): a block evaluator.  An RWG edge lives on two
//    triangles, so an entry is a sum over 4 triangle pairs, and a triangle
//    pair serves the up to 9 edge pairs on it.  For a block of row edges by
//    column edges, the quadrature of every pair of their triangles runs
//    once, into 9 sums per pair that do not depend on the edges; the entries
//    then combine the sums of their 4 pairs with the edges' geometry.  The
//    block's distinct triangles are found on the host (BPACK_GPU_HOST_IDS).
//
//  - CFIE: alpha Z_EFIE + (1 - alpha) eta0 Z_MFIE, an entry evaluator (not
//    symmetric: HODLR only).
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

// The EFIE entry Z(em, en) (Zelem_EMSURF, value_e), one entry at a time
__device__ __noinline__ dcomplex efie_entry(const Mesh& m, int64_t em, int64_t en) {
    const double k = m.k;
    const int* Em = m.edge + 6 * em;
    const int* En = m.edge + 6 * en;
    const double lm = dist(vertex(m, Em[0]), vertex(m, Em[1]));
    const double ln = dist(vertex(m, En[0]), vertex(m, En[1]));
    double pn[2][3][kMaxGauss];  // source points of both triangles of en
    for (int jj = 0; jj < 2; ++jj) {
        if (En[2 + jj] >= 0) gauss_points(m, En[2 + jj], pn[jj]);
    }
    double pm[3][kMaxGauss];
    dcomplex c1(0.0, 0.0), c2(0.0, 0.0);
    for (int ii = 0; ii < 2; ++ii) {
        const int tm = Em[2 + ii];
        if (tm < 0) continue;
        const double sm = ii == 0 ? 1.0 : -1.0;
        gauss_points(m, tm, pm);
        const double* om = vertex(m, Em[4 + ii]);
        for (int i = 0; i < m.nq; ++i) {
            const double xm[3] = {pm[0][i], pm[1][i], pm[2][i]};
            const double am[3] = {xm[0] - om[0], xm[1] - om[1], xm[2] - om[2]};
            dcomplex bb(0.0, 0.0);
            dcomplex aa[3] = {dcomplex(0.0, 0.0), dcomplex(0.0, 0.0), dcomplex(0.0, 0.0)};
            for (int jj = 0; jj < 2; ++jj) {
                const int tn = En[2 + jj];
                if (tn < 0) continue;
                const double sn = jj == 0 ? 1.0 : -1.0;
                dcomplex imp(0.0, 0.0), imp1(0.0, 0.0), imp2(0.0, 0.0);
                for (int j = 0; j < m.nq; ++j) {
                    const double xn[3] = {pn[jj][0][j], pn[jj][1][j], pn[jj][2][j]};
                    const double d = dist(xm, xn);
                    const double wj = m.gauss[4 * j + 3];
                    const double w1 = m.gauss[4 * j] * wj, w2 = m.gauss[4 * j + 1] * wj;
                    if (d == 0.0) {
                        // wn (-j k) plus the analytic singular parts (the same
                        // triangle, i == j)
                        double a0, a1, a2;
                        analytic(m, tn, xn, a0, a1, a2);
                        imp.im += wj * (-k);
                        imp1.im += w1 * (-k);
                        imp2.im += w2 * (-k);
                        imp.re += a0;
                        imp1.re += a1;
                        imp2.re += a2;
                    } else {
                        double sn_kd, cs_kd;
                        sincos(k * d, &sn_kd, &cs_kd);  // e^{-jkd} = (cos kd, -sin kd)
                        imp += dcomplex(wj * cs_kd / d, wj * -sn_kd / d);
                        imp1 += dcomplex(w1 * cs_kd / d, w1 * -sn_kd / d);
                        imp2 += dcomplex(w2 * cs_kd / d, w2 * -sn_kd / d);
                    }
                }
                const dcomplex imp3 = imp - imp1 - imp2;
                const int* tv = m.tri + 3 * static_cast<int64_t>(tn);
                const double* x1 = vertex(m, tv[0]);
                const double* x2 = vertex(m, tv[1]);
                const double* x3 = vertex(m, tv[2]);
                const double* xo = vertex(m, En[4 + jj]);
                const double coef = sn * (k * k);
                for (int d3 = 0; d3 < 3; ++d3) {
                    const dcomplex s = x1[d3] * imp1 + x2[d3] * imp2 + x3[d3] * imp3 - xo[d3] * imp;
                    aa[d3] += coef * s;
                }
                bb += sn * imp;
            }
            const dcomplex ctemp = aa[0] * am[0] + aa[1] * am[1] + aa[2] * am[2];
            c1 += (sm * ctemp) * m.gauss[4 * i + 3];
            c2 += ((4. * sm) * bb) * m.gauss[4 * i + 3];
        }
    }
    return efie_scale(m, ln * lm, c1 - c2);
}

// The MFIE part value_m of Zelem_EMSURF (times lm ln): the identity term
// 0.5 (am . an) / (2 area) on the triangles the two edges share, and
// n_m x (an x grad G) . am on the others, grad G = (xm - xn) (1 + j k d)
// e^{-jkd} / (4 pi d^3)
__device__ __noinline__ dcomplex mfie_entry(const Mesh& m, int64_t em, int64_t en) {
    const double k = m.k;
    const int* Em = m.edge + 6 * em;
    const int* En = m.edge + 6 * en;
    const double lm = dist(vertex(m, Em[0]), vertex(m, Em[1]));
    const double ln = dist(vertex(m, En[0]), vertex(m, En[1]));
    double pn[2][3][kMaxGauss];
    for (int jj = 0; jj < 2; ++jj) {
        if (En[2 + jj] >= 0) gauss_points(m, En[2 + jj], pn[jj]);
    }
    double pm[3][kMaxGauss];
    dcomplex value(0.0, 0.0);
    for (int ii = 0; ii < 2; ++ii) {
        const int tm = Em[2 + ii];
        if (tm < 0) continue;
        const int sm = ii == 0 ? 1 : -1;
        gauss_points(m, tm, pm);
        const double* om = vertex(m, Em[4 + ii]);
        const double* nr = m.normals + 3 * static_cast<int64_t>(tm);
        for (int i = 0; i < m.nq; ++i) {
            const double xm[3] = {pm[0][i], pm[1][i], pm[2][i]};
            const double am[3] = {xm[0] - om[0], xm[1] - om[1], xm[2] - om[2]};
            const double wi = m.gauss[4 * i + 3];
            for (int jj = 0; jj < 2; ++jj) {
                const int tn = En[2 + jj];
                if (tn < 0) continue;
                const int s = sm * (jj == 0 ? 1 : -1);
                const double* on = vertex(m, En[4 + jj]);
                if (tm == tn) {
                    const int* tv = m.tri + 3 * static_cast<int64_t>(tn);
                    const double* v1 = vertex(m, tv[0]);
                    const double* v2 = vertex(m, tv[1]);
                    const double* v3 = vertex(m, tv[2]);
                    double a[3], b[3], c[3];
                    for (int d = 0; d < 3; ++d) {
                        a[d] = v2[d] - v1[d];
                        b[d] = v3[d] - v1[d];
                    }
                    cross(a, b, c);
                    const double area = 0.5 * sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
                    const double an[3] = {xm[0] - on[0], xm[1] - on[1], xm[2] - on[2]};
                    const double temp = dot(am, an);
                    value.re += s * 0.5 * temp / (2. * area) * wi;
                } else {
                    for (int j = 0; j < m.nq; ++j) {
                        const double xn[3] = {pn[jj][0][j], pn[jj][1][j], pn[jj][2][j]};
                        const double d = dist(xm, xn);
                        const double an[3] = {xn[0] - on[0], xn[1] - on[1], xn[2] - on[2]};
                        double sn_kd, cs_kd;
                        sincos(k * d, &sn_kd, &cs_kd);
                        const double kd = k * d;
                        const double den = 4 * m.pi * (d * d * d);
                        dcomplex dg[3];
                        for (int q = 0; q < 3; ++q) {
                            // (x (1 + j k d)) e^{-jkd} / (4 pi d^3), x = xm - xn
                            const dcomplex f(xm[q] - xn[q], (xm[q] - xn[q]) * kd);
                            const dcomplex g = f * dcomplex(cs_kd, -sn_kd);
                            dg[q] = dcomplex(g.re / den, g.im / den);
                        }
                        // dg1 = an x dg, dg2 = nr x dg1 (rccurl), ctemp = dg2 . am (cscalar)
                        const dcomplex dg1[3] = {an[1] * dg[2] - dg[1] * an[2], -an[0] * dg[2] + dg[0] * an[2],
                                                 an[0] * dg[1] - dg[0] * an[1]};
                        const dcomplex dg2[3] = {nr[1] * dg1[2] - dg1[1] * nr[2], -nr[0] * dg1[2] + dg1[0] * nr[2],
                                                 nr[0] * dg1[1] - dg1[0] * nr[1]};
                        const dcomplex ctemp = dg2[0] * am[0] + dg2[1] * am[1] + dg2[2] * am[2];
                        value -= ((static_cast<double>(s) * ctemp) * wi) * m.gauss[4 * j + 3];
                    }
                }
            }
        }
    }
    return (value * lm) * ln;
}

// ---- CFIE: an entry evaluator

struct CfieEntry {
    using value_type = dcomplex;
    Mesh mesh;

    __device__ dcomplex operator()(const double*, int64_t i, const double*, int64_t j) const {
        const dcomplex ze = efie_entry(mesh, i, j);
        const dcomplex zm = mfie_entry(mesh, i, j);
        return mesh.alpha * ze + ((1. - mesh.alpha) * mesh.eta0) * zm;
    }
};

// ---- EFIE: a block evaluator by triangle pairs
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

constexpr int kPairSums = 9;

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
    dcomplex* sums;  // ntr * ntc * kPairSums
    dcomplex* out;   // out[r + c * ld]
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

// Blocks [b0, b1) of the batch: their triangle lists and items go to the
// device in one copy, with the scratch of their pair sums behind them; a
// group that does not fit in the device memory the library can give is
// split in two.
void efie_group(const bpack_gpu_block_batch* batch, const Emsurf& em, int b0, int b1) {
    struct Lists {
        std::vector<int> row_tri, col_tri, tr, tc;
    };
    std::vector<Lists> lists(static_cast<size_t>(b1 - b0));
    size_t ints = 0, sums = 0;
    int64_t max_pairs = 0, max_entries = 0;
    for (int b = b0; b < b1; ++b) {
        const bpack_gpu_block& blk = batch->blocks[b];
        Lists& l = lists[static_cast<size_t>(b - b0)];
        triangles_of(em.edges, blk.row_ids_host, blk.m, l.row_tri, l.tr);
        triangles_of(em.edges, blk.col_ids_host, blk.n, l.col_tri, l.tc);
        ints += l.row_tri.size() + l.col_tri.size() + l.tr.size() + l.tc.size();
        const int64_t pairs = static_cast<int64_t>(l.tr.size()) * static_cast<int64_t>(l.tc.size());
        sums += align(static_cast<size_t>(pairs) * kPairSums * sizeof(dcomplex));
        max_pairs = std::max(max_pairs, pairs);
        max_entries = std::max(max_entries, static_cast<int64_t>(blk.m) * blk.n);
    }
    const int count = b1 - b0;
    const size_t items_at = align(ints * sizeof(int));
    const size_t image_bytes = items_at + align(static_cast<size_t>(count) * sizeof(EfieItem));
    char* d = static_cast<char*>(batch->allocate(batch->allocator, image_bytes + sums));
    if (d == nullptr) {
        if (count == 1) fail("EFIE block", "out of device memory for one block's triangle pairs");
        const int mid = b0 + count / 2;
        efie_group(batch, em, b0, mid);
        efie_group(batch, em, mid, b1);
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
        const bpack_gpu_block& blk = batch->blocks[b];
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
        }
        efie_entries_kernel<<<dim3(static_cast<unsigned>((max_entries + 255) / 256), n), 256, 0, stream>>>(
            d_items + i, em.mesh);
    }
    check(cudaGetLastError(), "EFIE kernels");
    batch->release(batch->allocator, d);  // (later work on the stream is ordered after the kernels)
}

void efie_blocks(const bpack_gpu_block_batch* batch, void* user) {
    if (batch->count > 0) efie_group(batch, *static_cast<const Emsurf*>(user), 0, batch->count);
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
    if (cfie_alpha == 1.0) {
        em.edges.assign(ints, ints + 6 * edges);
        const int flags = BPACK_GPU_SYMMETRIC | BPACK_GPU_HOST_IDS;
        z_c_bpack_set_gpu_block_evaluator(bmat, &efie_blocks, &em, &flags);
    } else {
        const int64_t eta0_at = 3 * vertices + 4 * static_cast<int64_t>(nq);
        m.eta0 = reals[eta0_at];
        m.normals = em.reals + eta0_at + 1;
        bpack::gpu::set_entry_evaluator(bmat, CfieEntry{m}, 0);  // reads only the ids; not symmetric
    }
}

extern "C" void emsurf_gpu_release(void) {
    if (state == nullptr) return;
    cudaFree(state->reals);
    cudaFree(state->ints);
    delete state;
    state = nullptr;
}
