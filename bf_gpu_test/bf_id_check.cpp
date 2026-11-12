// Check the batched interpolative decomposition of bf_gpu against a host
// reference and against an invariant that needs no reference at all.
//
// Three checks per panel:
//   (a) the rank matches a host truncated QRCP run with the same rule;
//   (b) A ~ A(:, J) X to the requested tolerance;
//   (c) X(:, J) = I - the interpolation reproduces the columns it chose
//       exactly.  That one holds for any correct ID whatever the pivots are,
//       so it catches a wrong transpose, a wrong permutation or a wrong
//       triangular solve without sharing anything with the reference.  It does
//       NOT apply to a numerically zero panel: there the CPU returns a
//       placeholder - rank 1, a zero factor, the first pivot
//       (Bplus_compress.f90:1918-1922) - which is checked on its own terms.
//
// Pivot order is deliberately NOT compared: the device reduces in a different
// order, so a near-tie can pick the other column (see bf_id.hpp).

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

#include <cuda_runtime.h>

#include "bf_gpu/bf_id.hpp"

namespace {

// A host truncated QRCP with LAPACK's pivot rule and geqp3mod's early exit,
// written straight from the Fortran, then the same two follow-up steps.
// Returns the rank; fills jpvt and X (rank x n, lda = m).
int host_id(const std::vector<double>& panel, int m, int n, bool transposed, double rtol, double atol,
            std::vector<int>& jpvt, std::vector<double>& X) {
  const int pr = transposed ? n : m;
  std::vector<double> A(size_t(m) * n), B(size_t(m) * n);
  for (int j = 0; j < n; j++)
    for (int i = 0; i < m; i++) {
      const double v = transposed ? panel[size_t(j) + size_t(i) * pr] : panel[size_t(i) + size_t(j) * pr];
      A[size_t(i) + size_t(j) * m] = v;
      B[size_t(i) + size_t(j) * m] = v;
    }
  jpvt.resize(n);
  for (int j = 0; j < n; j++) jpvt[j] = j;
  std::vector<double> tau(std::min(m, n), 0.0), nrm(n, 0.0);
  for (int j = 0; j < n; j++) {
    double s = 0;
    for (int i = 0; i < m; i++) s += A[size_t(i) + size_t(j) * m] * A[size_t(i) + size_t(j) * m];
    nrm[j] = std::sqrt(s);
  }
  int rank = 0;
  double r11 = 0;
  const int mn = std::min(m, n);
  for (int c = 0; c < mn; c++) {
    int piv = c;
    for (int j = c; j < n; j++) if (nrm[j] > nrm[piv]) piv = j;
    if (piv != c) {
      for (int i = 0; i < m; i++) std::swap(A[size_t(i) + size_t(c) * m], A[size_t(i) + size_t(piv) * m]);
      std::swap(jpvt[c], jpvt[piv]);
      std::swap(nrm[c], nrm[piv]);
    }
    double tail = 0;
    for (int i = c + 1; i < m; i++) tail += A[size_t(i) + size_t(c) * m] * A[size_t(i) + size_t(c) * m];
    const double alpha = A[size_t(c) + size_t(c) * m];
    double beta = alpha;
    if (tail != 0.0) {
      const double q = std::sqrt(alpha * alpha + tail);
      beta = alpha >= 0 ? -q : q;
      tau[c] = (beta - alpha) / beta;
      const double scale = 1.0 / (alpha - beta);
      for (int i = c + 1; i < m; i++) A[size_t(i) + size_t(c) * m] *= scale;
    }
    A[size_t(c) + size_t(c) * m] = beta;
    if (c == 0) r11 = std::fabs(beta);
    if (r11 < atol) break;
    if (std::fabs(beta) / r11 < rtol) break;
    rank++;
    if (c + 1 >= n) break;
    if (tau[c] != 0.0) {
      for (int j = c + 1; j < n; j++) {
        double* col = &A[size_t(j) * m];
        double w = col[c];
        for (int i = c + 1; i < m; i++) w += A[size_t(i) + size_t(c) * m] * col[i];
        w *= tau[c];
        col[c] -= w;
        for (int i = c + 1; i < m; i++) col[i] -= w * A[size_t(i) + size_t(c) * m];
      }
    }
    for (int j = c + 1; j < n; j++) {
      double s = 0;
      for (int i = c + 1; i < m; i++) s += A[size_t(i) + size_t(j) * m] * A[size_t(i) + size_t(j) * m];
      nrm[j] = std::sqrt(s);
    }
  }
  X.assign(size_t(m) * n, 0.0);
  if (rank == 0) { jpvt[0] = 0; return 1; }
  for (int c = 0; c < rank; c++) {
    if (tau[c] == 0.0) continue;
    for (int j = 0; j < n; j++) {
      double* col = &B[size_t(j) * m];
      double w = col[c];
      for (int i = c + 1; i < m; i++) w += A[size_t(i) + size_t(c) * m] * col[i];
      w *= tau[c];
      col[c] -= w;
      for (int i = c + 1; i < m; i++) col[i] -= w * A[size_t(i) + size_t(c) * m];
    }
  }
  for (int j = 0; j < n; j++) {
    double* col = &B[size_t(j) * m];
    for (int i = rank - 1; i >= 0; i--) {
      double s = col[i];
      for (int k = i + 1; k < rank; k++) s -= A[size_t(i) + size_t(k) * m] * col[k];
      col[i] = s / A[size_t(i) + size_t(i) * m];
    }
  }
  X = B;
  return rank;
}

struct Case {
  int rows, cols, transposed, true_rank;
  const char* what;
};

}  // namespace

int main(int argc, char** argv) {
  const double rtol = argc > 1 ? std::atof(argv[1]) : 1e-6;
  const double atol = 1e-30;  // BPACK_SafeUnderflow is tiny; the rtol branch decides

  // Shapes chosen to cover the branches: wide and tall, both orientations, a
  // numerically deficient panel, a full-rank one, a single column, and a zero
  // panel (the CPU's rank-0 case, which must come back as rank 1).
  const std::vector<Case> cases = {
      {120,  80, 0, 17, "tall, decaying"},
      { 80, 120, 0, 23, "wide, decaying"},
      {120,  80, 1, 11, "tall, transposed"},
      { 80, 120, 1, 29, "wide, transposed"},
      {200, 150, 0, 150, "full rank"},
      { 64,  64, 0, 1,  "rank one"},
      { 40,  30, 0, 0,  "zero panel"},
      {  1,  50, 0, 1,  "single row"},
      { 50,   1, 0, 1,  "single column"},
  };

  // One contiguous entry buffer for the whole batch, as the Fortran hands a
  // level over.
  std::vector<double> entries;
  bf::gpu::IdBatch batch;
  std::vector<int64_t> entry_at;
  std::mt19937_64 rng(20261025);
  std::normal_distribution<double> gauss(0.0, 1.0);
  for (const Case& c : cases) {
    const int pr = c.rows, pc = c.cols;
    std::vector<double> P(size_t(pr) * pc, 0.0);
    if (c.true_rank > 0) {
      // P = sum_k s_k u_k v_k^T with a sharp drop after true_rank, so the
      // numerical rank at rtol is unambiguous
      const int r = std::min(c.true_rank, std::min(pr, pc));
      std::vector<double> u(pr), v(pc);
      for (int k = 0; k < r; k++) {
        for (int i = 0; i < pr; i++) u[i] = gauss(rng);
        for (int j = 0; j < pc; j++) v[j] = gauss(rng);
        const double s = std::pow(10.0, -3.0 * double(k) / std::max(r, 1));
        for (int j = 0; j < pc; j++)
          for (int i = 0; i < pr; i++) P[size_t(i) + size_t(j) * pr] += s * u[i] * v[j];
      }
    }
    bf::gpu::IdPanel panel;
    panel.rows = pr;
    panel.cols = pc;
    panel.transposed = c.transposed;
    panel.entry_offset = int64_t(entries.size());
    entry_at.push_back(panel.entry_offset);
    entries.insert(entries.end(), P.begin(), P.end());
    batch.panels.push_back(panel);
  }

  fmm::gpu::DeviceHeap& heap = fmm::gpu::DeviceHeap::instance();
  heap.ensure_initialized();
  double* d_entries = heap.alloc<double>(entries.size() * sizeof(double), fmm::gpu::ResidentKind::Adaptive);
  cudaMemcpy(d_entries, entries.data(), entries.size() * sizeof(double), cudaMemcpyHostToDevice);
  if (!bf::gpu::plan_id_batch(batch)) {
    std::printf("FAIL: no arena room for the batch\n");
    return 1;
  }
  bf::gpu::run_id_batch(d_entries, batch, rtol, atol, nullptr);
  cudaDeviceSynchronize();
  const cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    std::printf("FAIL: %s\n", cudaGetErrorString(err));
    return 1;
  }

  std::vector<int> rank(batch.count()), jpvt(batch.jpvt_elements);
  std::vector<double> X(batch.mat_elements);
  cudaMemcpy(rank.data(), batch.d_rank, rank.size() * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(jpvt.data(), batch.d_jpvt, jpvt.size() * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(X.data(), batch.d_x, X.size() * sizeof(double), cudaMemcpyDeviceToHost);

  int fails = 0;
  std::printf("rtol=%.1e, %d panels\n", rtol, batch.count());
  for (int p = 0; p < batch.count(); p++) {
    const bf::gpu::IdPanel& pan = batch.panels[p];
    const int m = pan.m(), n = pan.n(), pr = pan.rows;
    const double* P = entries.data() + entry_at[size_t(p)];
    const int r = rank[size_t(p)];
    const int* J = jpvt.data() + pan.jpvt_offset;
    const double* Xp = X.data() + pan.x_offset;

    std::vector<int> hjpvt;
    std::vector<double> hX;
    const int href = host_id(std::vector<double>(P, P + size_t(pr) * pan.cols), m, n, pan.transposed != 0,
                             rtol, atol, hjpvt, hX);

    // (b) A ~ A(:, J) X, and (c) X(:, J) = I
    auto Aij = [&](int i, int j) {
      return pan.transposed ? P[size_t(j) + size_t(i) * pr] : P[size_t(i) + size_t(j) * pr];
    };
    double num = 0, den = 0, eye = 0;
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < m; i++) {
        double s = 0;
        for (int k = 0; k < r; k++) s += Aij(i, J[k]) * Xp[size_t(k) + size_t(j) * m];
        const double a = Aij(i, j);
        num += (s - a) * (s - a);
        den += a * a;
      }
    }
    for (int k = 0; k < r; k++)
      for (int k2 = 0; k2 < r; k2++) {
        const double want = k == k2 ? 1.0 : 0.0;
        eye = std::max(eye, std::fabs(Xp[size_t(k) + size_t(J[k2]) * m] - want));
      }
    const double rel = den > 0 ? std::sqrt(num / den) : std::sqrt(num);

    const bool ok_rank = (r == href);
    const bool ok_rec = den == 0.0 ? (rel == 0.0) : (rel < 200.0 * rtol);
    bool ok_eye;
    if (den == 0.0) {
      // the placeholder the CPU returns for an empty panel: rank 1, X zero,
      // first pivot.  X(:, J) = I does not hold and must not.
      ok_eye = (r == 1) && (J[0] == 0);
      for (int j = 0; j < n && ok_eye; j++) ok_eye = (Xp[size_t(j) * m] == 0.0);
    } else {
      ok_eye = eye < 1e-8;
    }
    const char* verdict = (ok_rank && ok_rec && ok_eye) ? "PASS" : "FAIL";
    if (den == 0.0) {
      std::printf("  %-18s %3dx%-3d%s rank=%-4d host=%-4d  ||A-A(:,J)X||/||A||=%.2e  zero-panel placeholder  %s\n",
                  cases[size_t(p)].what, pan.rows, pan.cols, pan.transposed ? " T" : "  ", r, href, rel,
                  verdict);
    } else {
      std::printf("  %-18s %3dx%-3d%s rank=%-4d host=%-4d  ||A-A(:,J)X||/||A||=%.2e  max|X(:,J)-I|=%.1e  %s\n",
                  cases[size_t(p)].what, pan.rows, pan.cols, pan.transposed ? " T" : "  ", r, href, rel, eye,
                  verdict);
    }
    if (!(ok_rank && ok_rec && ok_eye)) fails++;
  }

  bf::gpu::release_id_batch(batch);
  heap.free(d_entries);
  std::printf("%s\n", fails == 0 ? "BATCHED ID OK" : "BATCHED ID FAILED");
  return fails == 0 ? 0 : 1;
}
