#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <mpi.h>

#include "runtime_thread_support.hpp"

namespace butterfly {

/**
 * Distributed ordering metadata shared by the format-neutral 64-bit C API.
 *
 * The input ordering is the caller's original rank-local ordering.  The
 * internal ordering is ButterflyPACK's rank-contiguous, tree-ordered layout.
 * Counts are kept in points; vector redistribution expands them by nrhs.
 */
struct DistributedLayout64 {
    MPI_Comm comm = MPI_COMM_NULL;
    int rank = 0;
    int size = 1;
    int dimension = 0;

    int64_t global_size = 0;
    int64_t input_local_size = 0;
    int64_t internal_local_size = 0;
    int64_t internal_global_start = 1;  // 1-based

    std::vector<int64_t> input_global_ids;       // 1-based original IDs
    std::vector<int64_t> internal_global_ids;    // 1-based original IDs
    std::vector<int> input_internal_owner;
    std::vector<int64_t> input_internal_global_index;  // 1-based

    // Point counts/displacements for input owner -> internal owner.
    std::vector<int> forward_send_counts;
    std::vector<int> forward_send_displacements;
    std::vector<int> forward_recv_counts;
    std::vector<int> forward_recv_displacements;

    // Packing/scattering slots in the same packet order as the arrays above.
    std::vector<int64_t> forward_send_input_slots;
    std::vector<int64_t> forward_recv_internal_slots;

    bool initialized() const noexcept {
        return comm != MPI_COMM_NULL;
    }

    template<typename DataType>
    void input_to_internal(const DataType* input, DataType* internal,
                           int nrhs) const {
        redistribute(input, internal, nrhs, true);
    }

    template<typename DataType>
    void internal_to_input(const DataType* internal, DataType* input,
                           int nrhs) const {
        redistribute(internal, input, nrhs, false);
    }

private:
    static int checked_scaled_count(int count, int nrhs,
                                    const char* caller) {
        if (count < 0 || nrhs <= 0 ||
            static_cast<int64_t>(count) * nrhs >
                std::numeric_limits<int>::max()) {
            throw std::overflow_error(
                std::string(caller) +
                ": MPI count exceeds the current INT_MAX implementation limit");
        }
        return count * nrhs;
    }

    template<typename DataType>
    void redistribute(const DataType* source, DataType* destination,
                      int nrhs, bool forward) const {
        if (!initialized()) {
            throw std::runtime_error(
                "DistributedLayout64::redistribute: layout is not initialized");
        }
        if (nrhs <= 0) {
            throw std::invalid_argument(
                "DistributedLayout64::redistribute: nrhs must be positive");
        }

        const int64_t source_size =
            forward ? input_local_size : internal_local_size;
        const int64_t destination_size =
            forward ? internal_local_size : input_local_size;
        const auto& send_slots = forward
            ? forward_send_input_slots : forward_recv_internal_slots;
        const auto& recv_slots = forward
            ? forward_recv_internal_slots : forward_send_input_slots;
        const auto& point_send_counts = forward
            ? forward_send_counts : forward_recv_counts;
        const auto& point_send_displacements = forward
            ? forward_send_displacements : forward_recv_displacements;
        const auto& point_recv_counts = forward
            ? forward_recv_counts : forward_send_counts;
        const auto& point_recv_displacements = forward
            ? forward_recv_displacements : forward_send_displacements;

        if ((source_size > 0 && source == nullptr) ||
            (destination_size > 0 && destination == nullptr)) {
            throw std::invalid_argument(
                "DistributedLayout64::redistribute: null vector buffer");
        }

        std::vector<int> send_counts(static_cast<size_t>(size));
        std::vector<int> send_displacements(static_cast<size_t>(size));
        std::vector<int> recv_counts(static_cast<size_t>(size));
        std::vector<int> recv_displacements(static_cast<size_t>(size));
        for (int process = 0; process < size; ++process) {
            send_counts[static_cast<size_t>(process)] = checked_scaled_count(
                point_send_counts[static_cast<size_t>(process)], nrhs,
                "DistributedLayout64::redistribute");
            recv_counts[static_cast<size_t>(process)] = checked_scaled_count(
                point_recv_counts[static_cast<size_t>(process)], nrhs,
                "DistributedLayout64::redistribute");
            const int64_t send_displacement =
                static_cast<int64_t>(point_send_displacements[
                    static_cast<size_t>(process)]) * nrhs;
            const int64_t recv_displacement =
                static_cast<int64_t>(point_recv_displacements[
                    static_cast<size_t>(process)]) * nrhs;
            if (send_displacement > std::numeric_limits<int>::max() ||
                recv_displacement > std::numeric_limits<int>::max()) {
                throw std::overflow_error(
                    "DistributedLayout64::redistribute: MPI displacement "
                    "exceeds the current INT_MAX implementation limit");
            }
            send_displacements[static_cast<size_t>(process)] =
                static_cast<int>(send_displacement);
            recv_displacements[static_cast<size_t>(process)] =
                static_cast<int>(recv_displacement);
        }

        const size_t send_points = send_slots.size();
        const size_t recv_points = recv_slots.size();
        if (send_points > std::numeric_limits<size_t>::max() /
                              static_cast<size_t>(nrhs) ||
            recv_points > std::numeric_limits<size_t>::max() /
                              static_cast<size_t>(nrhs)) {
            throw std::overflow_error(
                "DistributedLayout64::redistribute: vector size overflow");
        }
        std::vector<DataType> send_buffer(send_points * nrhs);
        std::vector<DataType> recv_buffer(recv_points * nrhs);

        // Each point packet stores all RHS contiguously.  This uses one
        // collective for the full batch rather than one collective per RHS.
        for (size_t packet = 0; packet < send_points; ++packet) {
            const int64_t slot = send_slots[packet];
            if (slot < 0 || slot >= source_size) {
                throw std::runtime_error(
                    "DistributedLayout64::redistribute: invalid send slot");
            }
            for (int rhs = 0; rhs < nrhs; ++rhs) {
                send_buffer[packet * nrhs + static_cast<size_t>(rhs)] =
                    source[static_cast<size_t>(slot) +
                           static_cast<size_t>(rhs) * source_size];
            }
        }

        MPI_Alltoallv(
            send_buffer.data(), send_counts.data(), send_displacements.data(),
            fmm::mpi_datatype_for<DataType>(),
            recv_buffer.data(), recv_counts.data(), recv_displacements.data(),
            fmm::mpi_datatype_for<DataType>(), comm);

        for (size_t packet = 0; packet < recv_points; ++packet) {
            const int64_t slot = recv_slots[packet];
            if (slot < 0 || slot >= destination_size) {
                throw std::runtime_error(
                    "DistributedLayout64::redistribute: invalid receive slot");
            }
            for (int rhs = 0; rhs < nrhs; ++rhs) {
                destination[static_cast<size_t>(slot) +
                            static_cast<size_t>(rhs) * destination_size] =
                    recv_buffer[packet * nrhs + static_cast<size_t>(rhs)];
            }
        }
    }
};

}  // namespace butterfly
