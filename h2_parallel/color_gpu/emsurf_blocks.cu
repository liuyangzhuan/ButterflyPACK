// Kind-3 (EMSURF EFIE) blocks by triangle pairs: see emsurf_blocks.hpp.
//
// Zelem_EMSURF's EFIE entry of test edge m and source edge n is
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

#include "emsurf_blocks.hpp"
#include "emsurf_kernel.cuh"

#include <stdexcept>
#include <string>

namespace fmm {
namespace gpu {

namespace {

void check_launch(const char* what) {
    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA launch failed in ") + what + ": " + cudaGetErrorString(status));
    }
}

__device__ __forceinline__ const int* int_list(const char* meta, int64_t offset) {
    return reinterpret_cast<const int*>(meta + offset);
}

// The 9 sums of each triangle pair (test triangle of the rows, source
// triangle of the columns): one thread per pair.
__global__ void __launch_bounds__(128) efie_pairs_kernel(const EfieBlockItem* items, const char* meta,
                                                         KernelSpec spec) {
    const EfieBlockItem item = items[blockIdx.y];
    const int64_t pair = static_cast<int64_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (pair >= static_cast<int64_t>(item.ntr) * item.ntc) return;
    const int t = int_list(meta, item.tr)[pair / item.ntc];
    const int s = int_list(meta, item.tc)[pair % item.ntc];
    const emsurf::Mesh m = emsurf::mesh_of(spec);
    const double k = spec.p[0];
    const int nq = m.nq;
    double y[3][kEmsurfMaxGauss];
    emsurf::gauss_points(m, s, y);
    const int* tv = m.tri + 3 * static_cast<int64_t>(t);
    const double* x1 = emsurf::vertex(m, tv[0]);
    const double* x2 = emsurf::vertex(m, tv[1]);
    const double* x3 = emsurf::vertex(m, tv[2]);
    dcomplex sum[kEfieSums];
#pragma unroll
    for (int q = 0; q < kEfieSums; ++q) sum[q] = dcomplex(0.0, 0.0);
    for (int i = 0; i < nq; ++i) {
        const double a1 = m.gauss[4 * i], a2 = m.gauss[4 * i + 1], a3 = m.gauss[4 * i + 2], wi = m.gauss[4 * i + 3];
        double xi[3];
        for (int d = 0; d < 3; ++d) {
            xi[d] = __dadd_rn(__dadd_rn(__dmul_rn(a1, x1[d]), __dmul_rn(a2, x2[d])), __dmul_rn(a3, x3[d]));
        }
        dcomplex i0(0.0, 0.0), i1(0.0, 0.0), i2(0.0, 0.0);
#pragma unroll
        for (int j = 0; j < kEmsurfMaxGauss; ++j) {
            if (j >= nq) break;
            const double yj[3] = {y[0][j], y[1][j], y[2][j]};
            const double d = emsurf::dist(xi, yj);
            const double wj = m.gauss[4 * j + 3];
            const double w1 = m.gauss[4 * j] * wj, w2 = m.gauss[4 * j + 1] * wj;
            if (t == s && d == 0.0) {
                // the same triangle, i == j: the host's analytic singular part
                double s0, s1, s2;
                emsurf::analytic(m, s, yj, s0, s1, s2);
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
    dcomplex* out = item.M + pair * kEfieSums;
#pragma unroll
    for (int q = 0; q < kEfieSums; ++q) out[q] = sum[q];
}

__device__ __forceinline__ void sub3(const double* a, const double* b, double* c) {
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
}

// The entries: row edge r tests, column edge c is the source; one thread
// per entry, the index of the smaller output stride fastest.
__global__ void __launch_bounds__(256) efie_entries_kernel(const EfieBlockItem* items, const char* meta,
                                                           KernelSpec spec) {
    const EfieBlockItem item = items[blockIdx.y];
    const int64_t e = static_cast<int64_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (e >= static_cast<int64_t>(item.nrow) * item.ncol) return;
    const bool rows_fast = item.rs <= item.cs;
    const int r = static_cast<int>(rows_fast ? e % item.nrow : e / item.ncol);
    const int c = static_cast<int>(rows_fast ? e / item.nrow : e % item.ncol);
    const emsurf::Mesh m = emsurf::mesh_of(spec);
    const double k2 = spec.p[0] * spec.p[0];
    const int* Er = m.edge + 6 * static_cast<int64_t>(int_list(meta, item.row_edges)[r]);
    const int* Ec = m.edge + 6 * static_cast<int64_t>(int_list(meta, item.col_edges)[c]);
    const int* row_tri = int_list(meta, item.row_tri) + 2 * static_cast<int64_t>(r);
    const int* col_tri = int_list(meta, item.col_tri) + 2 * static_cast<int64_t>(c);
    const int* tr = int_list(meta, item.tr);
    const int* tc = int_list(meta, item.tc);
    dcomplex total(0.0, 0.0);
    for (int a = 0; a < 2; ++a) {
        if (row_tri[a] < 0) continue;
        const int* tv = m.tri + 3 * static_cast<int64_t>(tr[row_tri[a]]);
        const double* x3 = emsurf::vertex(m, tv[2]);
        double f1[3], f2[3], u[3];
        sub3(emsurf::vertex(m, tv[0]), x3, f1);
        sub3(emsurf::vertex(m, tv[1]), x3, f2);
        sub3(emsurf::vertex(m, Er[4 + a]), x3, u);
        for (int b = 0; b < 2; ++b) {
            if (col_tri[b] < 0) continue;
            const int* sv = m.tri + 3 * static_cast<int64_t>(tc[col_tri[b]]);
            const double* y3 = emsurf::vertex(m, sv[2]);
            double e1[3], e2[3], v[3];
            sub3(emsurf::vertex(m, sv[0]), y3, e1);
            sub3(emsurf::vertex(m, sv[1]), y3, e2);
            sub3(emsurf::vertex(m, Ec[4 + b]), y3, v);
            const dcomplex* M = item.M + (static_cast<int64_t>(row_tri[a]) * item.ntc + col_tri[b]) * kEfieSums;
            dcomplex vB(0.0, 0.0), uC(0.0, 0.0);
            for (int d = 0; d < 3; ++d) {
                vB += v[d] * (f1[d] * M[1] + f2[d] * M[2]);
                uC += u[d] * (e1[d] * M[3] + e2[d] * M[4]);
            }
            const dcomplex D = emsurf::dot(f1, e1) * M[5] + emsurf::dot(f1, e2) * M[6] + emsurf::dot(f2, e1) * M[7] +
                               emsurf::dot(f2, e2) * M[8];
            const dcomplex F = k2 * (D - vB - uC + emsurf::dot(u, v) * M[0]) - 4. * M[0];
            total += ((a == 0 ? 1.0 : -1.0) * (b == 0 ? 1.0 : -1.0)) * F;
        }
    }
    // lm ln j total / 8 / pi^2 / freq / eps0, as efie_entry
    const double lengths = emsurf::dist(emsurf::vertex(m, Ec[0]), emsurf::vertex(m, Ec[1])) *
                           emsurf::dist(emsurf::vertex(m, Er[0]), emsurf::vertex(m, Er[1]));
    const double pi2 = spec.p[3] * spec.p[3];
    double re = -(lengths * total.im), im = lengths * total.re;
    re = re / 8. / pi2 / spec.p[1] / spec.p[2];
    im = im / 8. / pi2 / spec.p[1] / spec.p[2];
    item.out[r * item.rs + c * item.cs] = dcomplex(re, im);
}

}  // namespace

void launch_efie_blocks(const EfieBlockItem* items, int count, int64_t max_pairs, int64_t max_entries,
                        const char* meta, KernelSpec spec, cudaStream_t stream) {
    if (count <= 0) return;
    if (spec.kind != 3) throw std::runtime_error("launch_efie_blocks: kind 3 kernel expected");
    if (max_pairs > 0) {
        const dim3 grid(static_cast<unsigned>((max_pairs + 127) / 128), static_cast<unsigned>(count));
        efie_pairs_kernel<<<grid, 128, 0, stream>>>(items, meta, spec);
        check_launch("efie_pairs_kernel");
    }
    if (max_entries > 0) {
        const dim3 grid(static_cast<unsigned>((max_entries + 255) / 256), static_cast<unsigned>(count));
        efie_entries_kernel<<<grid, 256, 0, stream>>>(items, meta, spec);
        check_launch("efie_entries_kernel");
    }
}

}  // namespace gpu
}  // namespace fmm
