#pragma once

#include "options.hpp"
#include "occupied_topology.hpp"
#include "compression_impl.hpp"

namespace butterfly {
namespace color_unstructured {

template<typename CoordType, typename DataType>
void compress(H2<CoordType, DataType>* solver,
              double* compression_time,
              double* entryeval_time) {
    validate_tree(solver);

    if constexpr (!(std::is_same_v<DataType, double> ||
                    std::is_same_v<DataType, std::complex<double>>)) {
        throw std::runtime_error(
            "color_unstructured compression-only supports double and "
            "complex<double> data");
    } else {
        if (solver->build_state == H2BuildState::H2_COMPRESSED) {
            if (compression_time) *compression_time = 0.0;
            if (entryeval_time) *entryeval_time = 0.0;
            return;
        }
        if (solver->build_state == H2BuildState::RS_FACTORIZED) {
            throw std::runtime_error(
                "color_unstructured: solver already contains an RS-S factorization");
        }

        const auto occupied = prepare_occupied_topology(
            solver->tree.get(), solver->options.verbosity);
        solver->kernel.entryeval_time_per_thread.assign(
            omp_get_max_threads(), 0.0);

        const double start = MPI_Wtime();
        hierarchical_compression_unstructured(
            solver->tree.get(),
            &solver->kernel,
            solver->options.tolerance,
            &solver->last_factor_rankmax,
            &solver->factorization_memory,
            solver->options.use_sketch,
            solver->options.verbosity >= 1,
            occupied);

        double elapsed = MPI_Wtime() - start;
        MPI_Allreduce(
            MPI_IN_PLACE, &elapsed, 1, MPI_DOUBLE, MPI_MAX, solver->comm);
        if (compression_time) *compression_time = elapsed;

        double entry_time = 0.0;
        if (!solver->kernel.entryeval_time_per_thread.empty()) {
            entry_time = *std::max_element(
                solver->kernel.entryeval_time_per_thread.begin(),
                solver->kernel.entryeval_time_per_thread.end());
        }
        MPI_Allreduce(
            MPI_IN_PLACE, &entry_time, 1, MPI_DOUBLE, MPI_MAX, solver->comm);
        if (entryeval_time) *entryeval_time = entry_time;

        MPI_Allreduce(
            MPI_IN_PLACE, &solver->last_factor_rankmax,
            1, MPI_INT64_T, MPI_MAX, solver->comm);
        solver->build_state = H2BuildState::H2_COMPRESSED;
        solver->factorized = false;

        if (solver->options.verbosity >= 0 &&
            solver->tree->mpi_rank == 0) {
            std::cout << "\n========================================\n"
                      << "Unstructured H2 Compression Complete\n"
                      << "  total time: "
                      << std::llround(elapsed * 1000.0) << " ms\n"
                      << "========================================\n"
                      << std::endl;
        }
    }
}

}  // namespace color_unstructured
}  // namespace butterfly
