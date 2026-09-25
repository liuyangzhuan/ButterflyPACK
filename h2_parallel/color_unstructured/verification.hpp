#pragma once

#include "solver_impl.hpp"

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

}  // namespace color_unstructured
}  // namespace butterfly
