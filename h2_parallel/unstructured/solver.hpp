#pragma once

#include "options.hpp"
#include "solver_impl.hpp"

namespace butterfly {
namespace color_unstructured {

template<typename CoordType, typename DataType>
void solve(
    H2<CoordType, DataType>* solver,
    const std::vector<DataType>& rhs,
    std::vector<std::vector<fmm::SolveDataRequest<CoordType, DataType>>>& solve_data,
    int nrhs,
    int verbosity) {
    validate_tree(solver);
    hierarchical_solve_unstructured(
        solver->tree.get(), rhs, solve_data, nrhs, verbosity);
}

template<typename CoordType, typename DataType>
void multiply(
    H2<CoordType, DataType>* solver,
    const std::vector<DataType>& input,
    std::vector<std::vector<fmm::SolveDataRequest<CoordType, DataType>>>& multiply_data,
    int nrhs,
    bool verbose) {
    validate_tree(solver);
    hierarchical_mul_unstructured(
        solver->tree.get(), input, multiply_data, nrhs, verbose);
}

}  // namespace color_unstructured
}  // namespace butterfly
