#pragma once

#include "color_unstructured/factorization.hpp"
#include "color_unstructured/compression.hpp"
#include "color_unstructured/compression_apply.hpp"
#include "color_unstructured/solver.hpp"

namespace butterfly {

template<typename CoordType, typename DataType>
void dispatch_h2_factorization(H2<CoordType, DataType>* solver,
                               double* factorization_time,
                               double* entryeval_time) {
    if (solver->options.unstructured == 1) {
        color_unstructured::factorize(
            solver, factorization_time, entryeval_time);
    } else {
        butterfly::butterfly_factorization_parallel(
            solver, factorization_time, entryeval_time);
    }
}

template<typename CoordType, typename DataType>
void dispatch_h2_compression(H2<CoordType, DataType>* solver,
                             double* compression_time,
                             double* entryeval_time) {
    if (solver->options.unstructured == 1) {
        color_unstructured::compress(
            solver, compression_time, entryeval_time);
    } else {
        butterfly::butterfly_compression_parallel(
            solver, compression_time, entryeval_time);
    }
}

template<typename CoordType, typename DataType>
double dispatch_h2_compression_verification(
    H2<CoordType, DataType>* solver) {
    if (solver->options.unstructured == 1) {
        return color_unstructured::compression_quick_verification(
            solver->tree.get(), &solver->kernel);
    }
    return butterfly::h2_compression_quick_verification(
        solver->tree.get(), &solver->kernel);
}

template<typename CoordType, typename DataType>
void dispatch_h2_iterative_solve(
    H2<CoordType, DataType>* solver,
    const std::vector<DataType>& rhs,
    std::vector<DataType>& solution,
    int nrhs,
    double tolerance,
    int max_iterations,
    int* completed_iterations,
    double* final_relative_residual,
    bool verbose) {
    if (solver->options.unstructured == 1) {
        color_unstructured::compressed_bicgstab(
            solver->tree.get(), rhs, solution, nrhs, tolerance,
            max_iterations, completed_iterations,
            final_relative_residual, verbose);
        return;
    }
    butterfly::hierarchical_h2_bicgstab_parallel(
        solver->tree.get(), rhs, solution, nrhs, tolerance,
        max_iterations, completed_iterations,
        final_relative_residual, verbose);
}

template<typename CoordType, typename DataType>
void dispatch_h2_compressed_multiply(
    H2<CoordType, DataType>* solver,
    const std::vector<DataType>& input,
    std::vector<DataType>& output,
    int nrhs,
    bool verbose) {
    if (solver->options.unstructured == 1) {
        color_unstructured::compressed_multiply(
            solver->tree.get(), input, output, nrhs, verbose);
        return;
    }
    butterfly::hierarchical_h2_mul_parallel(
        solver->tree.get(), input, output, nrhs, verbose);
}

template<typename CoordType, typename DataType>
void dispatch_h2_solve(
    H2<CoordType, DataType>* solver,
    const std::vector<DataType>& rhs,
    std::vector<std::vector<fmm::SolveDataRequest<CoordType, DataType>>>& solve_data,
    int nrhs,
    int verbosity) {
    if (solver->options.unstructured == 1) {
        color_unstructured::solve(
            solver, rhs, solve_data, nrhs, verbosity);
    } else {
        butterfly::hierarchical_solve_parallel(
            solver->tree.get(), rhs, solve_data, nrhs, verbosity);
    }
}

template<typename CoordType, typename DataType>
void dispatch_h2_multiply(
    H2<CoordType, DataType>* solver,
    const std::vector<DataType>& input,
    std::vector<std::vector<fmm::SolveDataRequest<CoordType, DataType>>>& multiply_data,
    int nrhs,
    bool verbose) {
    if (solver->options.unstructured == 1) {
        color_unstructured::multiply(
            solver, input, multiply_data, nrhs, verbose);
    } else {
        butterfly::hierarchical_mul_parallel(
            solver->tree.get(), input, multiply_data, nrhs, verbose);
    }
}

}  // namespace butterfly
