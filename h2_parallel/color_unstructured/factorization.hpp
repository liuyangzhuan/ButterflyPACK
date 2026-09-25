#pragma once

#include "options.hpp"
#include "occupied_topology.hpp"
#include "factorization_impl.hpp"

namespace butterfly {
namespace color_unstructured {

template<typename CoordType, typename DataType>
void factorize(H2<CoordType, DataType>* solver,
               double* factorization_time,
               double* entryeval_time) {
    validate_tree(solver);

    if (solver->build_state == H2BuildState::RS_FACTORIZED) {
        if (factorization_time) *factorization_time = 0.0;
        if (entryeval_time) *entryeval_time = 0.0;
        return;
    }
    if (solver->build_state == H2BuildState::H2_COMPRESSED) {
        throw std::runtime_error(
            "color_unstructured: solver already contains a compression-only H2 representation");
    }

    const auto occupied = prepare_occupied_topology(
        solver->tree.get(), solver->options.verbosity);

    if constexpr (std::is_same_v<DataType, double> ||
                  std::is_same_v<DataType, std::complex<double>>) {
    const auto factorization_method = fmm::FactorizationMethod::BUNCH_KAUFMAN;
    fmm::HierarchicalFactorization<
        CoordType, DataType, butterfly::H2Kernel<CoordType, DataType>> factorizer(
            solver->options.N,
            fmm::MatrixProperty::SYMMETRIC,
            &solver->kernel,
            solver->options.dimension,
            factorization_method,
            solver->options.num_proxy);

    solver->kernel.entryeval_time_per_thread.assign(
        omp_get_max_threads(), 0.0);
    auto& entry_times = solver->kernel.entryeval_time_per_thread;

    const double start = MPI_Wtime();
    hierarchical_factorization_unstructured(
        solver->tree.get(),
        &solver->kernel,
        solver->options.tolerance,
        solver->options.use_sketch,
        true,
        false,
        factorization_method,
        factorizer.get_unit_proxy_points(),
        factorizer.get_num_proxy_points(),
        CoordType(0),
        &solver->last_factor_rankmax,
        &solver->factorization_memory,
        solver->options.lazy_schur,
        solver->options.gemm_split,
        0,
        0,
        0,
        solver->options.verbosity,
        occupied);

    double elapsed = MPI_Wtime() - start;
    MPI_Allreduce(
        MPI_IN_PLACE, &elapsed, 1, MPI_DOUBLE, MPI_MAX, solver->comm);
    if (factorization_time) *factorization_time = elapsed;

    double entry_elapsed = 0.0;
    if (!entry_times.empty()) {
        entry_elapsed = *std::max_element(
            entry_times.begin(), entry_times.end());
    }
    MPI_Allreduce(
        MPI_IN_PLACE, &entry_elapsed, 1, MPI_DOUBLE, MPI_MAX, solver->comm);
    if (entryeval_time) *entryeval_time = entry_elapsed;

    MPI_Allreduce(
        MPI_IN_PLACE, &solver->last_factor_rankmax, 1,
        MPI_INT64_T, MPI_MAX, solver->comm);
    solver->factorized = true;
    solver->build_state = H2BuildState::RS_FACTORIZED;
    } else {
        throw std::runtime_error(
            "color_unstructured only supports double precision");
    }
}

}  // namespace color_unstructured
}  // namespace butterfly
