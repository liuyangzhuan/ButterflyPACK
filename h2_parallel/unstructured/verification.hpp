#pragma once

#include "core/butterfly_matvec.hpp"

#include "solver_impl.hpp"
#include "partial.hpp"

namespace butterfly {
namespace color_unstructured {

template<typename CoordType, typename DataType, typename KernelType>
double h2_quick_verification_unstructured(
    ParallelTree<CoordType, DataType>* tree,
    KernelType* kernel,               // H2Kernel with evaluate_block_by_index
    int num_src = 20,                 // sparse support size (Fortran uses 20)
    unsigned seed = 12345,            // fixed → deterministic regression check
    bool verbose = true) {


    int rank;
    MPI_Comm_rank(tree->comm, &rank);

    const int leaf_level = tree->num_levels - 1;
    auto& leaf = tree->levels[leaf_level];


    using RealType = std::conditional_t<
        std::is_same_v<DataType, std::complex<double>>, double,
        std::conditional_t<std::is_same_v<DataType, std::complex<float>>, float, DataType>>;

    auto mag2 = [](const DataType& z) -> RealType {
        RealType a = std::abs(z);   // works for real and complex
        return a * a;
    };

    // step 0: set up sparse x_true and b = A_exact*x_true
    const auto verification = make_sparse_mvp_verification_data(
        tree, kernel, num_src, seed);
    const auto& x_true = verification.input;
    const auto& b_vec = verification.exact_output;
    const int64_t Nloc = static_cast<int64_t>(x_true.size());

    // step 1: calculate H2_approx*x_true
    std::vector<std::vector<fmm::SolveDataRequest<CoordType, DataType>>>
        mul_data(tree->num_levels);

    const double t0 = MPI_Wtime();
    hierarchical_mul_unstructured(tree, x_true, mul_data, 1, false);
    double t_mvp = MPI_Wtime() - t0;

    // step 2: compare against the sparse exact product prepared above
    RealType approx_Ax = 0, exact_Ax = 0, diff_mvp = 0, x_true_norm = 0;
    int64_t pos = 0;
    for (int64_t b = 0; b < leaf.num_boxes_local; ++b) {
        const auto& solve_box = mul_data[leaf_level][b];
        for (size_t i = 0; i < solve_box.left_side.size(); ++i, ++pos) {
            const DataType y_ref = b_vec[static_cast<size_t>(pos)];
            const DataType y_h2 = solve_box.left_side[i];     // compressed A x
            approx_Ax += mag2(y_h2);
            exact_Ax += mag2(y_ref);
            diff_mvp += mag2(y_h2 - y_ref);
            x_true_norm += mag2(x_true[static_cast<size_t>(pos)]);
        }
    }

    // step 3: solve H2_approx*x_h2 = b
    std::vector<std::vector<fmm::SolveDataRequest<CoordType, DataType>>>
        solve_data(tree->num_levels);

    const double t1 = MPI_Wtime();
    hierarchical_solve_unstructured(tree, b_vec, solve_data, 1, -1);
    double t_solve = MPI_Wtime() - t1;

    std::vector<DataType> x_h2(static_cast<size_t>(Nloc), DataType{0.0});
    RealType diff_forward = 0;
    pos = 0;
    for (int64_t b = 0; b < leaf.num_boxes_local; ++b) {
        const auto& solve_box = solve_data[leaf_level][b];
        for (size_t i = 0; i < solve_box.left_side.size(); ++i, ++pos) {
            const DataType x_approx = solve_box.left_side[i];
            x_h2[static_cast<size_t>(pos)] = x_approx;
            diff_forward += mag2(x_true[static_cast<size_t>(pos)] - x_approx);
        }
    }

    // step 4: calculate H2_approx*x_h2 for the backward residual
    std::vector<std::vector<fmm::SolveDataRequest<CoordType, DataType>>>
        backward_mul_data(tree->num_levels);

    const double t2 = MPI_Wtime();
    hierarchical_mul_unstructured(
        tree, x_h2, backward_mul_data, 1, false);
    double t_backward_mvp = MPI_Wtime() - t2;

    RealType diff_backward = 0;
    pos = 0;
    for (int64_t b = 0; b < leaf.num_boxes_local; ++b) {
        const auto& solve_box = backward_mul_data[leaf_level][b];
        for (size_t i = 0; i < solve_box.left_side.size(); ++i, ++pos) {
            const DataType h2_x_h2 = solve_box.left_side[i];
            diff_backward += mag2(h2_x_h2 - b_vec[static_cast<size_t>(pos)]);
        }
    }

    // step 5: reduce norms and timings across ranks

    RealType localv[6]  = {
        approx_Ax,
        exact_Ax,
        diff_mvp,
        x_true_norm,
        diff_forward,
        diff_backward
    };
    RealType globalv[6] = { 0, 0, 0, 0, 0, 0 };
    const MPI_Datatype real_mpi =
        std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Allreduce(localv, globalv, 6, real_mpi, MPI_SUM, tree->comm);

    double timings[3] = { t_mvp, t_solve, t_backward_mvp };
    MPI_Allreduce(MPI_IN_PLACE, timings, 3, MPI_DOUBLE, MPI_MAX, tree->comm);

    const RealType nrm_h2  = std::sqrt(globalv[0]);
    const RealType nrm_b = std::sqrt(globalv[1]);
    const RealType acc_mvp = (globalv[1] > RealType(0))
                             ? std::sqrt(globalv[2] / globalv[1]) : RealType(0);
    const RealType acc_forward = (globalv[3] > RealType(0))
                                 ? std::sqrt(globalv[4] / globalv[3]) : RealType(0);
    const RealType acc_backward = (globalv[1] > RealType(0))
                                  ? std::sqrt(globalv[5] / globalv[1]) : RealType(0);

    if (rank == 0 && verbose) {
        std::cout << "H2_CheckError(quick): fnorm: "
                  << std::scientific << std::setprecision(7)
                  << nrm_h2 << "  " << nrm_b
                  << "  acc_mvp: " << acc_mvp
                  << "  acc_forward: " << acc_forward
                  << "  acc_backward: " << acc_backward
                  << "  time_mvp: " << timings[0]
                  << "  time_solve: " << timings[1]
                  << "  time_backward_mvp: " << timings[2] << std::endl;
    }
    return static_cast<double>(acc_mvp);
}

// H2_partial (partial.hpp): a sparse x_true and b = A x_true (exact).  The
// forward partial solve of b gives b2 - F21 F11^-1 b1 = CB x2_true, compared
// with the Schur complement's own product (acc_schur: the two paths agree);
// the backward one with x2_true gives x1_true back (acc_forward, over the
// F11 points: how well F11 is factored).
template<typename CoordType, typename DataType, typename KernelType>
double h2_quick_verification_partial(
    ParallelTree<CoordType, DataType>* tree,
    KernelType* kernel,
    int num_src = 40,
    unsigned seed = 12345,
    bool verbose = true) {
    int rank;
    MPI_Comm_rank(tree->comm, &rank);
    const int leaf_level = tree->num_levels - 1;
    auto& leaf = tree->levels[leaf_level];
    using RealType = std::conditional_t<std::is_same_v<DataType, std::complex<double>>, double, DataType>;
    auto mag2 = [](const DataType& z) -> RealType { RealType a = std::abs(z); return a * a; };
    const auto verification = make_sparse_mvp_verification_data(tree, kernel, num_src, seed);
    const auto& x_true = verification.input;
    const auto& b_vec = verification.exact_output;
    const int64_t pc = tree->partial_count;
    const int64_t n22 = tree->num_points - pc;
    // x2_true, replicated in F22 index order
    std::vector<DataType> x2(static_cast<size_t>(n22), DataType(0.0));
    {
        int64_t pos = 0;
        for (int64_t b = 0; b < leaf.num_boxes_local; ++b) {
            const auto& box = leaf.local_boxes[static_cast<size_t>(b)];
            for (int64_t i = 0; i < box.num_points; ++i, ++pos) {
                const int64_t id = box.point_indices[static_cast<size_t>(i)];
                if (id >= pc) x2[static_cast<size_t>(id - pc)] = x_true[static_cast<size_t>(pos)];
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, x2.data(), static_cast<int>(n22), partial_mpi_type<DataType>(), MPI_SUM, tree->comm);
    }
    std::vector<std::vector<fmm::SolveDataRequest<CoordType, DataType>>> solve_data(tree->num_levels);
    std::vector<DataType> out2;
    const double t0 = MPI_Wtime();
    hierarchical_solve_unstructured(tree, b_vec, solve_data, 1, -1, false, 1,
                                    static_cast<const std::vector<DataType>*>(nullptr), &out2);
    const double t_forward = MPI_Wtime() - t0;
    std::vector<DataType> y2(static_cast<size_t>(n22), DataType(0.0));
    partial_schur_multiply(*tree, x2.data(), y2.data(), 1);
    RealType diff_schur = 0, norm_schur = 0, norm_b2 = 0;
    for (int64_t j = 0; j < n22; ++j) {
        diff_schur += mag2(out2[static_cast<size_t>(j)] - y2[static_cast<size_t>(j)]);
        norm_schur += mag2(out2[static_cast<size_t>(j)]);
    }
    {   // the F22 rows of b = A x_true: the Schur product's scale before the cancellation
        int64_t pos = 0;
        for (int64_t bi = 0; bi < leaf.num_boxes_local; ++bi) {
            const auto& box = leaf.local_boxes[static_cast<size_t>(bi)];
            for (int64_t i = 0; i < box.num_points; ++i, ++pos) {
                if (box.point_indices[static_cast<size_t>(i)] >= pc) norm_b2 += mag2(b_vec[static_cast<size_t>(pos)]);
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, &norm_b2, 1, std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT, MPI_SUM, tree->comm);
    }
    const double t1 = MPI_Wtime();
    hierarchical_solve_unstructured(tree, b_vec, solve_data, 1, -1, false, 2, &x2,
                                    static_cast<std::vector<DataType>*>(nullptr));
    const double t_backward = MPI_Wtime() - t1;
    RealType diff_forward = 0, norm_x1 = 0;
    {
        int64_t pos = 0;
        for (int64_t b = 0; b < leaf.num_boxes_local; ++b) {
            const auto& box = leaf.local_boxes[static_cast<size_t>(b)];
            const auto& solve_box = solve_data[leaf_level][static_cast<size_t>(b)];
            for (int64_t i = 0; i < box.num_points; ++i, ++pos) {
                if (box.point_indices[static_cast<size_t>(i)] >= pc) continue;
                diff_forward += mag2(solve_box.left_side[static_cast<size_t>(i)] - x_true[static_cast<size_t>(pos)]);
                norm_x1 += mag2(x_true[static_cast<size_t>(pos)]);
            }
        }
    }
    RealType localv[2] = {diff_forward, norm_x1}, globalv[2] = {0, 0};
    MPI_Allreduce(localv, globalv, 2, std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT, MPI_SUM, tree->comm);
    double timings[2] = {t_forward, t_backward};
    MPI_Allreduce(MPI_IN_PLACE, timings, 2, MPI_DOUBLE, MPI_MAX, tree->comm);
    const RealType acc_forward = globalv[1] > RealType(0) ? std::sqrt(globalv[0] / globalv[1]) : RealType(0);
    const RealType acc_schur = norm_schur > RealType(0) ? std::sqrt(diff_schur / norm_schur) : RealType(0);
    const RealType acc_schur_b2 = norm_b2 > RealType(0) ? std::sqrt(diff_schur / norm_b2) : RealType(0);
    if (rank == 0 && verbose) {
        std::cout << "H2_CheckError(partial): acc_forward: " << std::scientific << std::setprecision(7) << acc_forward
                  << " (over " << std::sqrt(globalv[1]) << " of x1)  acc_schur: " << acc_schur
                  << " (vs b2: " << acc_schur_b2 << ")"
                  << "  time_forward: " << timings[0] << "  time_backward: " << timings[1] << std::endl;
    }
    return static_cast<double>(acc_forward);
}

}  // namespace color_unstructured
}  // namespace butterfly
