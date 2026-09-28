#pragma once
// Kind 3 of KernelSpec: the EFIE entry Z(m, n) of the RWG basis functions of
// the edges m and n of a triangle mesh -- EMSURF's Zelem_EMSURF
// (EXAMPLE/EMSURF_Module.f90) with CFIE_alpha = 1, e^{-jkr} convention --
// in the host's order of operations (Gauss points, singular self terms by
// ianalytic / ianalytic2).
//   p[0] wavenumber k, p[1] frequency, p[2] eps0, p[3] pi, p[4] Gauss points
//   per triangle (at most kEmsurfMaxGauss), p[5] vertices, p[6] edges
//   treal: vertex xyz (3 per vertex), then the Gauss rule (ng1, ng2, ng3, w
//          per point)
//   tint:  per edge (6): its two vertices, its two triangles, and their
//          vertices opposite the edge; then per triangle (3): its vertices
//          (all 0-based; a missing triangle is -1)

#include "device_kernels.hpp"

namespace fmm {
namespace gpu {

constexpr int kEmsurfMaxGauss = 7;

namespace emsurf {

struct Mesh {
    const double* xyz;
    const double* gauss;  // ng1, ng2, ng3, w per point
    const int* edge;
    const int* tri;
    int nq;
};

__device__ __forceinline__ Mesh mesh_of(const KernelSpec& spec) {
    Mesh m;
    m.xyz = spec.treal;
    m.gauss = spec.treal + 3 * static_cast<int64_t>(spec.p[5]);
    m.edge = spec.tint;
    m.tri = spec.tint + 6 * static_cast<int64_t>(spec.p[6]);
    m.nq = static_cast<int>(spec.p[4]);
    return m;
}

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
__device__ __forceinline__ void gauss_points(const Mesh& m, int t, double (&p)[3][kEmsurfMaxGauss]) {
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
static __device__ __noinline__ void analytic(const Mesh& m, int t, const double* x, double& i0, double& i1,
                                             double& i2) {
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

// Z(em, en) (Zelem_EMSURF, value_e)
static __device__ __noinline__ dcomplex efie_entry(const KernelSpec& spec, int64_t em, int64_t en) {
    const Mesh m = mesh_of(spec);
    const double k = spec.p[0];
    const int* Em = m.edge + 6 * em;
    const int* En = m.edge + 6 * en;
    const double lm = dist(vertex(m, Em[0]), vertex(m, Em[1]));
    const double ln = dist(vertex(m, En[0]), vertex(m, En[1]));
    double pn[2][3][kEmsurfMaxGauss];  // source points of both triangles of en
    for (int jj = 0; jj < 2; ++jj) {
        if (En[2 + jj] >= 0) gauss_points(m, En[2 + jj], pn[jj]);
    }
    double pm[3][kEmsurfMaxGauss];
    dcomplex c1(0.0, 0.0), c2(0.0, 0.0);
    for (int ii = 0; ii < 2; ++ii) {
        const int tm = Em[2 + ii];
        if (tm < 0) continue;
        const double sm = ii == 0 ? 1.0 : -1.0;  // (-1)**(ii+1), ii = 3, 4
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
    // ln lm j (c1 - c2) / 8 / pi^2 / freq / eps0
    const double t = ln * lm;
    const dcomplex diff = c1 - c2;
    const double pi2 = spec.p[3] * spec.p[3];
    double re = -(t * diff.im), im = t * diff.re;
    re = re / 8. / pi2 / spec.p[1] / spec.p[2];
    im = im / 8. / pi2 / spec.p[1] / spec.p[2];
    return dcomplex(re, im);
}

}  // namespace emsurf
}  // namespace gpu
}  // namespace fmm
