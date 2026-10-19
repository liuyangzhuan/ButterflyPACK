#pragma once
// The drivers' check of c_bpack_inv_mult and the transposed c_bpack_mult
// (--inv_mult_check <ncol>, after the solve; needs a factorization, so not
// with precon 2).  On `ncol` test vectors U (ncol >= 2):
//   - inv_mult 'T' against 'N' (for a transpose-symmetric matrix F^{-T} =
//     F^{-1}; a nonsymmetric one prints their difference), and 'C' against
//     F^{-H} U = conj(F^{-T} conj(U));
//   - the transposes against 'N' for any matrix: u0^T F^{-T} u1 against
//     u1^T F^{-1} u0, and the same with F (bilinear, no conjugate);
//   - the residuals mult(op, inv_mult(op, U)) - U for op = 'N', 'T', 'C';
//   - inv_mult 'N' on the solve's right-hand sides against its solution;
//   - the symmetry of the factorized inverse and of the matrix as applied:
//     u0^T F^{-1} u1 against u1^T F^{-1} u0 (bilinear, no conjugate), and
//     the same with F - what makes 'T' equal 'N' beyond rounding;
//   - the time of `repeats` calls of inv_mult 'N' on the ncol columns (the
//     Schur-sampling loop: many calls on blocks of vectors).
// With precon 2 (no factorization) check_mult_transpose checks the
// transposed mult alone.
// Scalar is double or std::complex<double>; the callables wrap the C API:
//   inv_mult(char trans, const Scalar* in, Scalar* out, int ncol)
//   mult(char trans, const Scalar* in, Scalar* out, int ncol)

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <type_traits>
#include <vector>

namespace bpack_example {

template<typename Scalar>
Scalar inv_mult_check_conj(Scalar value) {
  if constexpr (std::is_same_v<Scalar, double>) {
    return value;
  } else {
    return std::conj(value);
  }
}

// test vectors: smooth in the row, different per rank and column
template<typename Scalar>
std::vector<Scalar> inv_mult_check_vectors(int rank, int local_rows, int ncol) {
  std::vector<Scalar> u(static_cast<size_t>(local_rows) * static_cast<size_t>(ncol));
  for (int column = 0; column < ncol; ++column) {
    for (int row = 0; row < local_rows; ++row) {
      const double phase = 0.37 * (row + 1009.0 * rank) + 1.3 * column;
      Scalar value{};
      if constexpr (!std::is_same_v<Scalar, double>) {
        value = Scalar(std::sin(phase) + 0.1 * (column + 1), std::cos(1.7 * phase));
      } else {
        value = std::sin(phase) + 0.1 * (column + 1);
      }
      u[static_cast<size_t>(row) + static_cast<size_t>(column) * local_rows] = value;
    }
  }
  return u;
}

// column j of a against column k of b over the ranks, bilinear (no conjugate)
template<typename Scalar>
Scalar inv_mult_check_dot(MPI_Comm comm, int local_rows, const std::vector<Scalar>& a, int j,
                          const std::vector<Scalar>& b, int k) {
  Scalar s{};
  for (int row = 0; row < local_rows; ++row) {
    s += a[static_cast<size_t>(row) + static_cast<size_t>(j) * local_rows] *
         b[static_cast<size_t>(row) + static_cast<size_t>(k) * local_rows];
  }
  if constexpr (!std::is_same_v<Scalar, double>) {
    double parts[2] = {s.real(), s.imag()};
    MPI_Allreduce(MPI_IN_PLACE, parts, 2, MPI_DOUBLE, MPI_SUM, comm);
    return Scalar(parts[0], parts[1]);
  } else {
    MPI_Allreduce(MPI_IN_PLACE, &s, 1, MPI_DOUBLE, MPI_SUM, comm);
    return s;
  }
}

template<typename Scalar, typename InvMult, typename Mult>
void check_inv_mult(MPI_Comm comm, int local_rows, int ncol, int repeats,
                    const std::vector<Scalar>& solve_rhs,
                    const std::vector<Scalar>& solve_solution, int solve_nrhs,
                    InvMult inv_mult, Mult mult, const char* label) {
  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  ncol = std::max(ncol, 2);
  const size_t count = static_cast<size_t>(local_rows) * static_cast<size_t>(ncol);

  const std::vector<Scalar> u = inv_mult_check_vectors<Scalar>(rank, local_rows, ncol);
  auto sum_all = [&](double value) {
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_DOUBLE, MPI_SUM, comm);
    return value;
  };
  auto max_all = [&](double value) {
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_DOUBLE, MPI_MAX, comm);
    return value;
  };
  auto norm = [&](const std::vector<Scalar>& v) {
    double s = 0.0;
    for (const Scalar& x : v) s += std::norm(x);
    return std::sqrt(sum_all(s));
  };
  auto max_abs_difference = [&](const std::vector<Scalar>& a, const std::vector<Scalar>& b) {
    double m = 0.0;
    for (size_t i = 0; i < a.size(); ++i) m = std::max(m, std::abs(a[i] - b[i]));
    return max_all(m);
  };
  auto max_abs = [&](const std::vector<Scalar>& a) {
    double m = 0.0;
    for (const Scalar& x : a) m = std::max(m, std::abs(x));
    return max_all(m);
  };
  auto dot = [&](const std::vector<Scalar>& a, int j, const std::vector<Scalar>& b, int k) {
    return inv_mult_check_dot(comm, local_rows, a, j, b, k);
  };
  auto symmetry_defect = [&](const std::vector<Scalar>& applied) {
    const Scalar a01 = dot(u, 0, applied, 1), a10 = dot(u, 1, applied, 0);
    return std::abs(a01 - a10) / std::max(std::abs(a01), std::abs(a10));
  };

  std::vector<Scalar> x_n(count), x_t(count), x_c(count);
  inv_mult('N', u.data(), x_n.data(), ncol);
  inv_mult('T', u.data(), x_t.data(), ncol);
  inv_mult('C', u.data(), x_c.data(), ncol);
  // F^{-H} U against conj(F^{-T} conj(U))
  std::vector<Scalar> u_conj(u), x_c_ref(count);
  for (Scalar& x : u_conj) x = inv_mult_check_conj(x);
  inv_mult('T', u_conj.data(), x_c_ref.data(), ncol);
  for (Scalar& x : x_c_ref) x = inv_mult_check_conj(x);

  const double u_norm = norm(u);
  auto residual = [&](char op, const std::vector<Scalar>& x) {
    std::vector<Scalar> product(count);
    mult(op, x.data(), product.data(), ncol);
    for (size_t i = 0; i < count; ++i) product[i] -= u[i];
    return norm(product) / u_norm;
  };
  const double res_n = residual('N', x_n), res_t = residual('T', x_t), res_c = residual('C', x_c);
  std::vector<Scalar> a_u(count), a_t(count);
  mult('N', u.data(), a_u.data(), ncol);
  mult('T', u.data(), a_t.data(), ncol);

  double solve_difference = -1.0;
  if (solve_nrhs > 0 && !solve_rhs.empty()) {
    std::vector<Scalar> x_s(solve_rhs.size());
    inv_mult('N', solve_rhs.data(), x_s.data(), solve_nrhs);
    solve_difference = max_abs_difference(x_s, solve_solution) / max_abs(solve_solution);
  }

  const double x_scale = max_abs(x_n);
  const double t_vs_n = max_abs_difference(x_t, x_n) / x_scale;
  const double c_vs_ref = max_abs_difference(x_c, x_c_ref) / x_scale;
  const double inverse_symmetry = symmetry_defect(x_n), matrix_symmetry = symmetry_defect(a_u);
  // u0^T (applied_t)_1 against u1^T (applied_n)_0
  auto transpose_defect = [&](const std::vector<Scalar>& applied_t, const std::vector<Scalar>& applied_n) {
    const Scalar a = dot(u, 0, applied_t, 1), b = dot(u, 1, applied_n, 0);
    return std::abs(a - b) / std::max(std::abs(a), std::abs(b));
  };
  const double inverse_transpose = transpose_defect(x_t, x_n), matrix_transpose = transpose_defect(a_t, a_u);

  // the repeated apply
  std::vector<Scalar> scratch(count);
  MPI_Barrier(comm);
  const double t0 = MPI_Wtime();
  for (int r = 0; r < repeats; ++r) inv_mult('N', u.data(), scratch.data(), ncol);
  const double per_call = max_all((MPI_Wtime() - t0) / std::max(repeats, 1));
  if (rank == 0) {
    std::printf("%s inv_mult check (%d columns): 'T' vs 'N' %.3e, 'C' vs conj(T(conj)) %.3e; "
                "residual N %.3e T %.3e C %.3e; vs the solve %.3e; symmetry of F^{-1} %.3e, of F %.3e; "
                "transpose of F^{-1} %.3e, of F %.3e; %d calls %.4f s each (%.5f s a column)\n",
                label, ncol, t_vs_n, c_vs_ref, res_n, res_t, res_c, solve_difference, inverse_symmetry,
                matrix_symmetry, inverse_transpose, matrix_transpose, repeats, per_call, per_call / ncol);
    std::fflush(stdout);
  }
}

// The transposed mult of a compression-only H2 (--inv_mult_check with
// precon 2: there is no inverse): u0^T F^T u1 against u1^T F u0, and F^H U
// against conj(F^T conj(U)); also 'T' against 'N' and the symmetry of F.
template<typename Scalar, typename Mult>
void check_mult_transpose(MPI_Comm comm, int local_rows, int ncol, Mult mult, const char* label) {
  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  ncol = std::max(ncol, 2);
  const size_t count = static_cast<size_t>(local_rows) * static_cast<size_t>(ncol);
  const std::vector<Scalar> u = inv_mult_check_vectors<Scalar>(rank, local_rows, ncol);
  std::vector<Scalar> a_n(count), a_t(count), a_c(count), u_conj(u), a_c_ref(count);
  mult('N', u.data(), a_n.data(), ncol);
  mult('T', u.data(), a_t.data(), ncol);
  mult('C', u.data(), a_c.data(), ncol);
  for (Scalar& x : u_conj) x = inv_mult_check_conj(x);
  mult('T', u_conj.data(), a_c_ref.data(), ncol);
  for (Scalar& x : a_c_ref) x = inv_mult_check_conj(x);
  auto dot = [&](const std::vector<Scalar>& a, int j, const std::vector<Scalar>& b, int k) {
    return inv_mult_check_dot(comm, local_rows, a, j, b, k);
  };
  auto defect = [](Scalar a, Scalar b) { return std::abs(a - b) / std::max(std::abs(a), std::abs(b)); };
  const double transpose = defect(dot(u, 0, a_t, 1), dot(u, 1, a_n, 0));
  const double symmetry = defect(dot(u, 0, a_n, 1), dot(u, 1, a_n, 0));
  double m[3] = {0.0, 0.0, 0.0};  // max |a_c - a_c_ref|, |a_t - a_n|, |a_n|
  for (size_t i = 0; i < count; ++i) {
    m[0] = std::max(m[0], std::abs(a_c[i] - a_c_ref[i]));
    m[1] = std::max(m[1], std::abs(a_t[i] - a_n[i]));
    m[2] = std::max(m[2], std::abs(a_n[i]));
  }
  MPI_Allreduce(MPI_IN_PLACE, m, 3, MPI_DOUBLE, MPI_MAX, comm);
  if (rank == 0) {
    std::printf("%s mult check (%d columns): transpose of F %.3e, 'C' vs conj(T(conj)) %.3e; 'T' vs 'N' %.3e, "
                "symmetry of F %.3e\n",
                label, ncol, transpose, m[0] / m[2], m[1] / m[2], symmetry);
    std::fflush(stdout);
  }
}

}  // namespace bpack_example
