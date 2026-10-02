#pragma once

// H2 (format-7) integration: per-level MPI communicator sets shared by the
// structured and unstructured factorization and solve drivers.

#include <vector>

#include "tree.hpp"

namespace butterfly {
using namespace fmm;

struct SolveCommunicatorSet {
    std::vector<MPI_Comm> level;
};

template<typename CoordType, typename DataType>
SolveCommunicatorSet make_solve_communicators(
    ParallelTree<CoordType, DataType>* tree) {
    SolveCommunicatorSet communicators;
    communicators.level.assign(
        static_cast<size_t>(tree->num_levels), MPI_COMM_NULL);
    for (int level = 0; level < tree->num_levels; ++level) {
        MPI_Comm_split(
            tree->comm,
            tree->levels[level].is_process_active ? 0 : MPI_UNDEFINED,
            tree->mpi_rank,
            &communicators.level[static_cast<size_t>(level)]);
    }
    return communicators;
}

inline void destroy_solve_communicators(
    SolveCommunicatorSet& communicators) {
    for (MPI_Comm& comm : communicators.level) {
        if (comm != MPI_COMM_NULL) MPI_Comm_free(&comm);
    }
}

struct FactorizationCommunicatorSet {
    std::vector<MPI_Comm> level;
    std::vector<MPI_Comm> transition;
};

template<typename CoordType, typename DataType>
FactorizationCommunicatorSet make_factorization_communicators(
    fmm::ParallelTree<CoordType, DataType>* tree) {
    FactorizationCommunicatorSet communicators;
    communicators.level.assign(
        static_cast<size_t>(tree->num_levels), MPI_COMM_NULL);
    communicators.transition.assign(
        static_cast<size_t>(tree->num_levels), MPI_COMM_NULL);

    for (int level = 0; level < tree->num_levels; ++level) {
        const int color = tree->levels[level].is_process_active
            ? 0
            : MPI_UNDEFINED;
        MPI_Comm_split(
            tree->comm,
            color,
            tree->mpi_rank,
            &communicators.level[static_cast<size_t>(level)]);
    }

    for (int child_level = 1;
         child_level < tree->num_levels;
         ++child_level) {
        const bool participates =
            tree->levels[child_level].is_process_active ||
            tree->levels[child_level - 1].is_process_active;
        MPI_Comm_split(
            tree->comm,
            participates ? 0 : MPI_UNDEFINED,
            tree->mpi_rank,
            &communicators.transition[static_cast<size_t>(child_level)]);
    }

    return communicators;
}

inline void destroy_factorization_communicators(
    FactorizationCommunicatorSet& communicators) {
    for (MPI_Comm& comm : communicators.transition) {
        if (comm != MPI_COMM_NULL) {
            MPI_Comm_free(&comm);
        }
    }
    for (MPI_Comm& comm : communicators.level) {
        if (comm != MPI_COMM_NULL) {
            MPI_Comm_free(&comm);
        }
    }
}

} // namespace butterfly
