#pragma once
// The drivers' check of c_bpack_inv_mult and the transposed c_bpack_mult
// (--inv_mult_check <ncol>, after the solve; needs a factorization, so not
// with precon 2).  On `ncol` test vectors U (ncol >= 2):
//   - inv_mult 'T' and 'C' against 'N': for a transpose-symmetric matrix
//     F^{-T} = F^{-1}, and F^{-H} U = conj(F^{-1} conj(U));
//   - the residuals mult(op, inv_mult(op, U)) - U for op = 'N', 'T', 'C';
//   - inv_mult 'N' on the solve's right-hand sides against its solution;
//   - the symmetry of the factorized inverse and of the matrix as applied:
//     u0^T F^{-1} u1 against u1^T F^{-1} u0 (bilinear, no conjugate), and
//     the same with F - what makes 'T' equal 'N' beyond rounding;
//   - the time of `repeats` calls of inv_mult 'N' on the ncol columns (the
//     Schur-sampling loop: many calls on blocks of vectors).
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

template<typename Scalar, typename InvMult, typename Mult>
void check_inv_mult(MPI_Comm comm, int local_rows, int ncol, int repeats,
                    const std::vector<Scalar>& solve_rhs,
                    const std::vector<Scalar>& solve_solution, int solve_nrhs,
                    InvMult inv_mult, Mult mult, const char* label) {
  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  ncol = std::max(ncol, 2);
  const size_t count = static_cast<size_t>(local_rows) * static_cast<size_t>(ncol);
  constexpr bool complex_data = !std::is_same_v<Scalar, double>;

  // test vectors: smooth in the row, different per rank and column
  std::vector<Scalar> u(count);
  for (int column = 0; column < ncol; ++column) {
    for (int row = 0; row < local_rows; ++row) {
      const double phase = 0.37 * (row + 1009.0 * rank) + 1.3 * column;
      Scalar value{};
      if constexpr (complex_data) {
        value = Scalar(std::sin(phase) + 0.1 * (column + 1), std::cos(1.7 * phase));
      } else {
        value = std::sin(phase) + 0.1 * (column + 1);
      }
      u[static_cast<size_t>(row) + static_cast<size_t>(column) * local_rows] = value;
    }
  }
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
  // column j of a against column k of b, bilinear (no conjugate)
  auto dot = [&](const std::vector<Scalar>& a, int j, const std::vector<Scalar>& b, int k) {
    Scalar s{};
    for (int row = 0; row < local_rows; ++row) {
      s += a[static_cast<size_t>(row) + static_cast<size_t>(j) * local_rows] *
           b[static_cast<size_t>(row) + static_cast<size_t>(k) * local_rows];
    }
    if constexpr (complex_data) {
      double parts[2] = {s.real(), s.imag()};
      MPI_Allreduce(MPI_IN_PLACE, parts, 2, MPI_DOUBLE, MPI_SUM, comm);
      return Scalar(parts[0], parts[1]);
    } else {
      return sum_all(s);
    }
  };
  auto symmetry_defect = [&](const std::vector<Scalar>& applied) {
    const Scalar a01 = dot(u, 0, applied, 1), a10 = dot(u, 1, applied, 0);
    return std::abs(a01 - a10) / std::max(std::abs(a01), std::abs(a10));
  };

  std::vector<Scalar> x_n(count), x_t(count), x_c(count);
  inv_mult('N', u.data(), x_n.data(), ncol);
  inv_mult('T', u.data(), x_t.data(), ncol);
  inv_mult('C', u.data(), x_c.data(), ncol);
  // F^{-H} U against conj(F^{-1} conj(U))
  std::vector<Scalar> u_conj(u), x_c_ref(count);
  for (Scalar& x : u_conj) x = inv_mult_check_conj(x);
  inv_mult('N', u_conj.data(), x_c_ref.data(), ncol);
  for (Scalar& x : x_c_ref) x = inv_mult_check_conj(x);

  const double u_norm = norm(u);
  auto residual = [&](char op, const std::vector<Scalar>& x) {
    std::vector<Scalar> product(count);
    mult(op, x.data(), product.data(), ncol);
    for (size_t i = 0; i < count; ++i) product[i] -= u[i];
    return norm(product) / u_norm;
  };
  const double res_n = residual('N', x_n), res_t = residual('T', x_t), res_c = residual('C', x_c);
  std::vector<Scalar> a_u(count);
  mult('N', u.data(), a_u.data(), ncol);

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

  // the repeated apply
  std::vector<Scalar> scratch(count);
  MPI_Barrier(comm);
  const double t0 = MPI_Wtime();
  for (int r = 0; r < repeats; ++r) inv_mult('N', u.data(), scratch.data(), ncol);
  const double per_call = max_all((MPI_Wtime() - t0) / std::max(repeats, 1));
  if (rank == 0) {
    std::printf("%s inv_mult check (%d columns): 'T' vs 'N' %.3e, 'C' vs conj(N(conj)) %.3e; "
                "residual N %.3e T %.3e C %.3e; vs the solve %.3e; symmetry of F^{-1} %.3e, of F %.3e; "
                "%d calls %.4f s each (%.5f s a column)\n",
                label, ncol, t_vs_n, c_vs_ref, res_n, res_t, res_c, solve_difference, inverse_symmetry,
                matrix_symmetry, repeats, per_call, per_call / ncol);
    std::fflush(stdout);
  }
}

}  // namespace bpack_example
