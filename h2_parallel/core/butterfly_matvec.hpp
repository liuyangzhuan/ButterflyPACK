#pragma once

// H2 (format-7) integration: shared nested-basis machinery - box
// skeletonization, parent-box assembly, one-hop vector exchange, the H2
// matvec and column reductions, and the sparse-MVP verification data the
// compression is checked against.

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "butterfly_types.hpp"
#include "id_decomposition.hpp"
#include "morton.hpp"
#include "runtime_thread_support.hpp"
#include "tree.hpp"

namespace butterfly {
using namespace fmm;

template<typename CoordType, typename DataType>
void ensure_h2_assisting_slots(
    TreeLevel<CoordType, DataType>& level,
    const std::vector<int64_t>& mortons) {

    for (int64_t morton_index : mortons) {
        if (level.find_local_box(morton_index) != nullptr) continue;
        if (level.assisting_box_points_for_kernel_evaluation.count(morton_index)) continue;

        const int64_t slot = static_cast<int64_t>(level.assisting_boxes.size());
        PointDataRequest<CoordType> request;
        request.morton_index = morton_index;
        level.assisting_boxes.push_back(std::move(request));
        level.assisting_box_points_for_kernel_evaluation[morton_index] = slot;
    }
}

template<typename CoordType, typename DataType>
std::vector<int64_t> h2_box_global_indices(
    TreeLevel<CoordType, DataType>& level,
    int64_t morton_index,
    bool skeleton_only) {

    if (auto* box = level.find_local_box(morton_index)) {
        if (!skeleton_only) return box->point_indices;

        std::vector<int64_t> result;
        result.reserve(box->skeleton_indices.size());
        for (int64_t local_index : box->skeleton_indices) {
            if (local_index < 0 || local_index >= static_cast<int64_t>(box->point_indices.size())) {
                throw std::runtime_error("h2_box_global_indices: invalid local skeleton index");
            }
            result.push_back(box->point_indices[static_cast<size_t>(local_index)]);
        }
        return result;
    }

    auto assist_it = level.assisting_box_points_for_kernel_evaluation.find(morton_index);
    if (assist_it == level.assisting_box_points_for_kernel_evaluation.end()) {
        throw std::runtime_error(
            "h2_box_global_indices: missing remote box " + std::to_string(morton_index));
    }
    const auto& request = level.assisting_boxes.at(static_cast<size_t>(assist_it->second));
    if (!skeleton_only) {
        if (request.indices.empty()) {
            throw std::runtime_error(
                "h2_box_global_indices: remote box has no point indices " +
                std::to_string(morton_index));
        }
        return request.indices;
    }

    std::vector<int64_t> result;
    result.reserve(request.skel_indices.size());
    for (int64_t local_index : request.skel_indices) {
        if (local_index < 0 || local_index >= static_cast<int64_t>(request.indices.size())) {
            throw std::runtime_error("h2_box_global_indices: invalid remote skeleton index");
        }
        result.push_back(request.indices[static_cast<size_t>(local_index)]);
    }
    return result;
}

template<typename CoordType, typename DataType, typename KernelType>
void h2_skeletonize_box(
    const ParallelTree<CoordType, DataType>* tree,
    BoxData<CoordType, DataType>* box,
    TreeLevel<CoordType, DataType>& level,
    KernelType* kernel,
    double tolerance,
    int use_sketch) {

    // Adaptive ID rows (H2_ID_proxy 2) with H2_use_sketch 2: selected against
    // the sketch of the ring rows and appended to it unsketched, as the
    // factorization's streamed sketch does (and the device); otherwise
    // against the materialized target, which is then sketched or not.
    const bool adaptive_on_sketch = use_sketch == 2 && tree->id_proxy_mode == 2;
    FactorizationThreadScratch<CoordType, DataType> scratch;
    gather_id_workspace(
        tree,
        box, level, kernel, tolerance,
        static_cast<const CoordType*>(nullptr), 0, CoordType{0}, true,
        scratch.workspace, scratch.workspace_rows, scratch.workspace_cols,
        0, box->on_boundary, false, &scratch.id_adaptive, !adaptive_on_sketch);

    if (scratch.workspace_cols == 0) {
        box->skeleton_indices.clear();
        box->redundant_indices.clear();
        box->interpolation_matrix.allocate(0, 0, MatrixStorage<DataType>::NONE);
        return;
    }
    if (scratch.workspace_rows == 0) {
        if (h2_id_trace_enabled()) {
            h2_id_trace_write("L" + std::to_string(box->level) + " m=" + std::to_string(box->morton_index) +
                              " local=1 ob=" + std::to_string(box->on_boundary) +
                              " n=" + std::to_string(box->num_points) + " k=" + std::to_string(box->num_points) +
                              " wave=-1 | rows=0 d=0 sketch_norm=0");
        }
        box->skeleton_indices.resize(static_cast<size_t>(box->num_points));
        std::iota(box->skeleton_indices.begin(), box->skeleton_indices.end(), int64_t{0});
        box->redundant_indices.clear();
        box->interpolation_matrix.allocate(box->num_points, 0);
        return;
    }

    constexpr double sketch_factor = 1.0;
    constexpr int sketch_nonzeros = 4;
    std::string trace_tail;
    if (h2_id_trace_enabled()) {
        // the norm of the ID input (the sketch, as compute_id_sparse_sketch
        // forms it), for box-by-box comparisons with the device path
        const int64_t m = scratch.workspace_rows, n = scratch.workspace_cols;
        int64_t d = m;
        double ssq = 0.0;
        if (use_sketch) {
            d = std::max<int64_t>(std::min<int64_t>(static_cast<int64_t>(std::ceil(sketch_factor * n)), m), n);
            std::vector<DataType> sketch(static_cast<size_t>(d * n));
            sketch_sparse_random(scratch.workspace.data(), m, n, m, sketch.data(), d, d,
                                 std::min<int>(sketch_nonzeros, static_cast<int>(d)),
                                 static_cast<uint64_t>(box->morton_index + 1));
            for (const DataType& v : sketch) ssq += std::norm(v);
        } else {
            for (const DataType& v : scratch.workspace) ssq += std::norm(v);
        }
        std::ostringstream tail;
        tail << std::setprecision(10) << " | rows=" << m << " d=" << d << " sketch_norm=" << std::sqrt(ssq);
        trace_tail = tail.str();
    }
    IDResult<DataType> id;
    if (adaptive_on_sketch) {
        const int64_t m = scratch.workspace_rows, n = scratch.workspace_cols;
        const int64_t d = std::max<int64_t>(std::min<int64_t>(static_cast<int64_t>(std::ceil(sketch_factor * n)), m), n);
        scratch.sketch_storage.assign(static_cast<size_t>(d * n), DataType{});
        sketch_sparse_random(scratch.workspace.data(), m, n, m, scratch.sketch_storage.data(), d, d,
                             std::min<int>(sketch_nonzeros, static_cast<int>(d)),
                             static_cast<uint64_t>(box->morton_index + 1));
        int64_t rows = d;
        append_adaptive_id_training_rows(tree, box, kernel, tolerance, true, scratch.sketch_storage, rows, n,
                                         &scratch.id_adaptive);
        id = compute_id_complex(scratch.sketch_storage.data(), rows, n, rows, tolerance, 0);
    } else if (use_sketch) {
        id = compute_id_sparse_sketch(
            scratch.workspace.data(), scratch.sketch_storage,
            scratch.workspace_rows, scratch.workspace_cols, scratch.workspace_rows,
            tolerance, sketch_factor, sketch_nonzeros,
            static_cast<uint64_t>(box->morton_index + 1));
    } else {
        id = compute_id_complex(
            scratch.workspace.data(),
            scratch.workspace_rows, scratch.workspace_cols, scratch.workspace_rows,
            tolerance, 0);
    }
    ensure_nonempty_box_id_rank(id, scratch.workspace_cols);
    if (h2_id_trace_enabled()) {
        h2_id_trace_write("L" + std::to_string(box->level) + " m=" + std::to_string(box->morton_index) +
                          " local=1 ob=" + std::to_string(box->on_boundary) +
                          " n=" + std::to_string(box->num_points) +
                          " k=" + std::to_string(id.skeleton_indices.size()) + " wave=-1" + trace_tail +
                          h2_id_adaptive_trace(scratch.id_adaptive));
    }

    box->skeleton_indices = std::move(id.skeleton_indices);
    box->redundant_indices = std::move(id.redundant_indices);
    box->interpolation_matrix = std::move(id.interpolation);
}

template<typename CoordType, typename DataType>
void install_h2_parent_boxes(
    ParallelTree<CoordType, DataType>* tree,
    int child_level_number,
    std::vector<BoxData<CoordType, DataType>> local_parents) {

    auto& child_level = tree->levels[child_level_number];
    auto& parent_level = tree->levels[child_level_number - 1];
    const int rank = tree->mpi_rank;
    const bool reduction =
        parent_level.num_active_processes != child_level.num_active_processes;

    if (!reduction) {
        if (parent_level.is_process_active) {
            parent_level.local_boxes = std::move(local_parents);
        }
        return;
    }

    const bool keep_local =
        child_level.is_process_active && child_level.parent_level_owner == rank;
    const bool send_parents =
        child_level.is_process_active && child_level.parent_level_owner != rank;
    const int size_tag = 7200 + 2 * child_level_number;
    const int data_tag = size_tag + 1;

    if (parent_level.is_process_active) {
        std::vector<BoxData<CoordType, DataType>> gathered;
        for (int child_rank : parent_level.children_senders) {
            if (keep_local && child_rank == rank) {
                gathered.insert(
                    gathered.end(),
                    std::make_move_iterator(local_parents.begin()),
                    std::make_move_iterator(local_parents.end()));
                continue;
            }

            int64_t buffer_size = 0;
            MPI_Status status;
            MPI_Recv(
                &buffer_size, 1, MPI_INT64_T, child_rank, size_tag, tree->comm, &status);
            std::vector<char> buffer(static_cast<size_t>(buffer_size));
            MPI_Recv_large(
                buffer.data(), static_cast<size_t>(buffer_size), MPI_CHAR,
                child_rank, data_tag, tree->comm, &status);
            auto received = deserialize_boxes<CoordType, DataType>(buffer);
            gathered.insert(
                gathered.end(),
                std::make_move_iterator(received.begin()),
                std::make_move_iterator(received.end()));
        }
        std::sort(
            gathered.begin(), gathered.end(),
            [](const auto& lhs, const auto& rhs) {
                return lhs.morton_index < rhs.morton_index;
            });
        parent_level.local_boxes = std::move(gathered);
    }

    if (send_parents) {
        std::vector<char> buffer = serialize_boxes(local_parents);
        const int64_t buffer_size = static_cast<int64_t>(buffer.size());
        MPI_Send(
            &buffer_size, 1, MPI_INT64_T,
            child_level.parent_level_owner, size_tag, tree->comm);
        MPI_Send_large(
            buffer.data(), buffer.size(), MPI_CHAR,
            child_level.parent_level_owner, data_tag, tree->comm);
    }

    if (parent_level.is_process_active &&
        parent_level.local_boxes.size() != static_cast<size_t>(parent_level.num_boxes_local)) {
        throw std::runtime_error(
            "install_h2_parent_boxes: received the wrong number of parent boxes");
    }
}

template<typename CoordType, typename DataType>
std::vector<DataType> h2_local_vector_for_morton(
    TreeLevel<CoordType, DataType>& level,
    const std::vector<SolveDataRequest<CoordType, DataType>>& level_data,
    int64_t morton_index,
    bool skeleton_only) {

    auto* box = level.find_local_box(morton_index);
    if (box == nullptr) {
        throw std::runtime_error(
            "h2_local_vector_for_morton: requested box is not local");
    }
    const int64_t local_index = morton_index - level.local_morton_start;
    if (local_index < 0 || local_index >= static_cast<int64_t>(level_data.size())) {
        throw std::runtime_error(
            "h2_local_vector_for_morton: missing local vector data");
    }
    const auto& data = level_data[static_cast<size_t>(local_index)];

    if (!skeleton_only) return data.right_side;

    const int64_t skeleton_count =
        static_cast<int64_t>(box->skeleton_indices.size());
    std::vector<DataType> result(
        static_cast<size_t>(skeleton_count * data.nrhs));
    for (int64_t column = 0; column < data.nrhs; ++column) {
        for (int64_t i = 0; i < skeleton_count; ++i) {
            result[static_cast<size_t>(i + column * skeleton_count)] =
                data.left_side.at(static_cast<size_t>(
                    box->skeleton_indices[static_cast<size_t>(i)] +
                    column * data.num_points));
        }
    }
    return result;
}

template<typename CoordType, typename DataType>
std::unordered_map<int64_t, std::vector<DataType>> exchange_h2_vectors_onehop(
    ParallelTree<CoordType, DataType>* tree,
    int level_number,
    const std::vector<SolveDataRequest<CoordType, DataType>>& level_data,
    std::vector<int64_t> needed_remote_mortons,
    bool skeleton_only,
    int tag_base) {

    auto& level = tree->levels[level_number];
    std::unordered_map<int64_t, std::vector<DataType>> received_vectors;
    if (!level.is_process_active) return received_vectors;

    const int rank = tree->mpi_rank;
    const auto neighbor_ranks =
        compute_one_hop_neighbor_ranks(tree, level, level_number);
    const std::unordered_set<int> neighbor_set(
        neighbor_ranks.begin(), neighbor_ranks.end());

    const uint32_t grid_size = 1u << level_number;
    auto owner_of_morton = [&](int64_t morton_index) {
        std::vector<uint64_t> one{static_cast<uint64_t>(morton_index)};
        const auto regions = morton::assign_to_processes_nd(
            tree->dimension, one, level.num_active_processes, grid_size);
        return level.morton_to_rank.at(static_cast<int>(regions.front()));
    };

    std::sort(needed_remote_mortons.begin(), needed_remote_mortons.end());
    needed_remote_mortons.erase(
        std::unique(needed_remote_mortons.begin(), needed_remote_mortons.end()),
        needed_remote_mortons.end());

    std::unordered_map<int, std::vector<int64_t>> requests_to_send;
    for (int64_t morton_index : needed_remote_mortons) {
        const int owner = owner_of_morton(morton_index);
        if (owner == rank) continue;
        if (!neighbor_set.count(owner)) {
            throw std::runtime_error(
                "exchange_h2_vectors_onehop: source owner is not a neighboring rank");
        }
        requests_to_send[owner].push_back(morton_index);
    }

    std::vector<int> send_counts(neighbor_ranks.size(), 0);
    std::vector<int> recv_counts(neighbor_ranks.size(), 0);
    std::vector<MPI_Request> mpi_requests;
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        auto it = requests_to_send.find(neighbor_ranks[i]);
        if (it != requests_to_send.end()) {
            send_counts[i] = static_cast<int>(it->second.size());
        }
        MPI_Request request;
        MPI_Irecv(
            &recv_counts[i], 1, MPI_INT, neighbor_ranks[i],
            tag_base, tree->comm, &request);
        mpi_requests.push_back(request);
    }
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        MPI_Send(
            &send_counts[i], 1, MPI_INT, neighbor_ranks[i],
            tag_base, tree->comm);
    }
    if (!mpi_requests.empty()) {
        MPI_Waitall(
            static_cast<int>(mpi_requests.size()), mpi_requests.data(),
            MPI_STATUSES_IGNORE);
        mpi_requests.clear();
    }

    std::vector<std::vector<int64_t>> requests_received(neighbor_ranks.size());
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (recv_counts[i] <= 0) continue;
        requests_received[i].resize(static_cast<size_t>(recv_counts[i]));
        MPI_Request request;
        MPI_Irecv(
            requests_received[i].data(), recv_counts[i], MPI_INT64_T,
            neighbor_ranks[i], tag_base + 1, tree->comm, &request);
        mpi_requests.push_back(request);
    }
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (send_counts[i] <= 0) continue;
        MPI_Send(
            requests_to_send[neighbor_ranks[i]].data(), send_counts[i], MPI_INT64_T,
            neighbor_ranks[i], tag_base + 1, tree->comm);
    }
    if (!mpi_requests.empty()) {
        MPI_Waitall(
            static_cast<int>(mpi_requests.size()), mpi_requests.data(),
            MPI_STATUSES_IGNORE);
        mpi_requests.clear();
    }

    std::vector<std::vector<char>> send_buffers(neighbor_ranks.size());
    std::vector<uint64_t> send_sizes(neighbor_ranks.size(), 0);
    std::vector<uint64_t> recv_sizes(neighbor_ranks.size(), 0);
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        size_t total_size = 0;
        std::vector<std::vector<DataType>> values;
        values.reserve(requests_received[i].size());
        for (int64_t morton_index : requests_received[i]) {
            values.push_back(h2_local_vector_for_morton(
                level, level_data, morton_index, skeleton_only));
            total_size += 2 * sizeof(int64_t) +
                values.back().size() * sizeof(DataType);
        }

        send_buffers[i].resize(total_size);
        char* pointer = send_buffers[i].data();
        for (size_t record = 0; record < requests_received[i].size(); ++record) {
            const int64_t morton_index = requests_received[i][record];
            const int64_t value_count =
                static_cast<int64_t>(values[record].size());
            std::memcpy(pointer, &morton_index, sizeof(int64_t));
            pointer += sizeof(int64_t);
            std::memcpy(pointer, &value_count, sizeof(int64_t));
            pointer += sizeof(int64_t);
            if (value_count > 0) {
                std::memcpy(
                    pointer, values[record].data(),
                    static_cast<size_t>(value_count) * sizeof(DataType));
                pointer += static_cast<size_t>(value_count) * sizeof(DataType);
            }
        }
        send_sizes[i] = static_cast<uint64_t>(total_size);
    }

    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        MPI_Request request;
        MPI_Irecv(
            &recv_sizes[i], 1, MPI_UINT64_T, neighbor_ranks[i],
            tag_base + 2, tree->comm, &request);
        mpi_requests.push_back(request);
    }
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        MPI_Send(
            &send_sizes[i], 1, MPI_UINT64_T, neighbor_ranks[i],
            tag_base + 2, tree->comm);
    }
    if (!mpi_requests.empty()) {
        MPI_Waitall(
            static_cast<int>(mpi_requests.size()), mpi_requests.data(),
            MPI_STATUSES_IGNORE);
        mpi_requests.clear();
    }

    std::vector<std::vector<char>> recv_buffers(neighbor_ranks.size());
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        recv_buffers[i].resize(static_cast<size_t>(recv_sizes[i]));
        if (recv_sizes[i] == 0) continue;
        MPI_Irecv_large(
            recv_buffers[i].data(), recv_buffers[i].size(), MPI_CHAR,
            neighbor_ranks[i], tag_base + 3, tree->comm, mpi_requests);
    }
    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        if (send_sizes[i] == 0) continue;
        MPI_Send_large(
            send_buffers[i].data(), send_buffers[i].size(), MPI_CHAR,
            neighbor_ranks[i], tag_base + 3, tree->comm);
    }
    if (!mpi_requests.empty()) {
        MPI_Waitall(
            static_cast<int>(mpi_requests.size()), mpi_requests.data(),
            MPI_STATUSES_IGNORE);
    }

    for (const auto& buffer : recv_buffers) {
        const char* pointer = buffer.data();
        const char* end = pointer + buffer.size();
        while (pointer < end) {
            if (end - pointer < static_cast<ptrdiff_t>(2 * sizeof(int64_t))) {
                throw std::runtime_error(
                    "exchange_h2_vectors_onehop: truncated record header");
            }
            int64_t morton_index = -1;
            int64_t value_count = -1;
            std::memcpy(&morton_index, pointer, sizeof(int64_t));
            pointer += sizeof(int64_t);
            std::memcpy(&value_count, pointer, sizeof(int64_t));
            pointer += sizeof(int64_t);
            if (value_count < 0 ||
                end - pointer < static_cast<ptrdiff_t>(
                    static_cast<size_t>(value_count) * sizeof(DataType))) {
                throw std::runtime_error(
                    "exchange_h2_vectors_onehop: invalid record size");
            }
            auto& values = received_vectors[morton_index];
            values.resize(static_cast<size_t>(value_count));
            if (value_count > 0) {
                std::memcpy(
                    values.data(), pointer,
                    static_cast<size_t>(value_count) * sizeof(DataType));
                pointer += static_cast<size_t>(value_count) * sizeof(DataType);
            }
        }
    }
    return received_vectors;
}

template<typename DataType>
void h2_matrix_vector_product(
    const MatrixStorage<DataType>& matrix,
    const std::vector<DataType>& input,
    std::vector<DataType>& output,
    int64_t nrhs,
    char transpose = 'N') {

    const bool transposed = transpose != 'N' && transpose != 'n';
    const int64_t input_size = transposed ? matrix.rows : matrix.cols;
    const int64_t output_size = transposed ? matrix.cols : matrix.rows;
    if (nrhs <= 0 || static_cast<int64_t>(input.size()) != input_size * nrhs) {
        throw std::runtime_error("h2_matrix_vector_product: input dimension mismatch");
    }
    output.assign(static_cast<size_t>(output_size * nrhs), DataType{0});
    if (input_size == 0 || output_size == 0) return;

    const int m = static_cast<int>(output_size);
    const int n = static_cast<int>(nrhs);
    const int k = static_cast<int>(input_size);
    const int lda = static_cast<int>(matrix.lda);
    const int ldb = static_cast<int>(input_size);
    const int ldc = static_cast<int>(output_size);
    const DataType alpha = DataType{1};
    const DataType beta = DataType{0};
    const char trans_b = 'N';
    gemm_(&transpose, &trans_b, &m, &n, &k, &alpha,
          matrix.data.data(), &lda, input.data(), &ldb,
          &beta, output.data(), &ldc);
}

template<typename CoordType, typename DataType>
void apply_h2_upward_projection(
    const BoxData<CoordType, DataType>& box,
    SolveDataRequest<CoordType, DataType>& data) {

    if (box.redundant_indices.empty()) return;
    const int64_t redundant_count =
        static_cast<int64_t>(box.redundant_indices.size());
    const int64_t skeleton_count =
        static_cast<int64_t>(box.skeleton_indices.size());
    std::vector<DataType> redundant(
        static_cast<size_t>(redundant_count * data.nrhs));
    for (int64_t column = 0; column < data.nrhs; ++column) {
        for (int64_t i = 0; i < redundant_count; ++i) {
            redundant[static_cast<size_t>(i + column * redundant_count)] =
                data.left_side.at(static_cast<size_t>(
                    box.redundant_indices[static_cast<size_t>(i)] +
                    column * data.num_points));
        }
    }
    std::vector<DataType> projected;
    h2_matrix_vector_product(
        box.interpolation_matrix, redundant, projected, data.nrhs, 'N');
    for (int64_t column = 0; column < data.nrhs; ++column) {
        for (int64_t i = 0; i < skeleton_count; ++i) {
            data.left_side.at(static_cast<size_t>(
                box.skeleton_indices[static_cast<size_t>(i)] +
                column * data.num_points)) +=
                projected[static_cast<size_t>(i + column * skeleton_count)];
        }
    }
}

template<typename CoordType, typename DataType>
void apply_h2_downward_interpolation(
    const BoxData<CoordType, DataType>& box,
    SolveDataRequest<CoordType, DataType>& data) {

    if (box.redundant_indices.empty()) return;
    const int64_t skeleton_count =
        static_cast<int64_t>(box.skeleton_indices.size());
    const int64_t redundant_count =
        static_cast<int64_t>(box.redundant_indices.size());
    std::vector<DataType> skeleton(
        static_cast<size_t>(skeleton_count * data.nrhs));
    for (int64_t column = 0; column < data.nrhs; ++column) {
        for (int64_t i = 0; i < skeleton_count; ++i) {
            skeleton[static_cast<size_t>(i + column * skeleton_count)] =
                data.left_side.at(static_cast<size_t>(
                    box.skeleton_indices[static_cast<size_t>(i)] +
                    column * data.num_points));
        }
    }
    std::vector<DataType> interpolated;
    h2_matrix_vector_product(
        box.interpolation_matrix, skeleton, interpolated, data.nrhs, 'T');
    for (int64_t column = 0; column < data.nrhs; ++column) {
        for (int64_t i = 0; i < redundant_count; ++i) {
            data.left_side.at(static_cast<size_t>(
                box.redundant_indices[static_cast<size_t>(i)] +
                column * data.num_points)) +=
                interpolated[static_cast<size_t>(i + column * redundant_count)];
        }
    }
}

template<typename DataType>
std::vector<DataType> h2_global_column_dots(
    const std::vector<DataType>& lhs,
    const std::vector<DataType>& rhs,
    int nrhs,
    MPI_Comm comm) {

    if (nrhs <= 0 || lhs.size() != rhs.size() ||
        lhs.size() % static_cast<size_t>(nrhs) != 0) {
        throw std::runtime_error(
            "h2_global_column_dots: invalid batched vector dimensions");
    }
    const size_t local_rows = lhs.size() / static_cast<size_t>(nrhs);
    std::vector<DataType> local(static_cast<size_t>(nrhs), DataType{0});
    for (int column = 0; column < nrhs; ++column) {
        const size_t offset = static_cast<size_t>(column) * local_rows;
        for (size_t row = 0; row < local_rows; ++row) {
            if constexpr (is_complex_v<DataType>) {
                local[static_cast<size_t>(column)] +=
                    std::conj(lhs[offset + row]) * rhs[offset + row];
            } else {
                local[static_cast<size_t>(column)] +=
                    lhs[offset + row] * rhs[offset + row];
            }
        }
    }
    std::vector<DataType> global(static_cast<size_t>(nrhs), DataType{0});
    MPI_Allreduce(
        local.data(), global.data(), nrhs, mpi_datatype_for<DataType>(),
        MPI_SUM, comm);
    return global;
}

template<typename DataType>
std::vector<double> h2_global_column_norms(
    const std::vector<DataType>& values,
    int nrhs,
    MPI_Comm comm) {

    if (nrhs <= 0 || values.size() % static_cast<size_t>(nrhs) != 0) {
        throw std::runtime_error(
            "h2_global_column_norms: invalid batched vector dimensions");
    }
    const size_t local_rows = values.size() / static_cast<size_t>(nrhs);
    std::vector<double> local(static_cast<size_t>(nrhs), 0.0);
    for (int column = 0; column < nrhs; ++column) {
        const size_t offset = static_cast<size_t>(column) * local_rows;
        for (size_t row = 0; row < local_rows; ++row) {
            if constexpr (is_complex_v<DataType>) {
                local[static_cast<size_t>(column)] +=
                    std::norm(values[offset + row]);
            } else {
                const double value = static_cast<double>(values[offset + row]);
                local[static_cast<size_t>(column)] += value * value;
            }
        }
    }
    std::vector<double> global(static_cast<size_t>(nrhs), 0.0);
    MPI_Allreduce(
        local.data(), global.data(), nrhs, MPI_DOUBLE, MPI_SUM, comm);
    for (double& value : global) value = std::sqrt(value);
    return global;
}

inline uint64_t splitmix64(uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

template<typename RealType>
RealType centered_uniform_from_hash(uint64_t key) {
    constexpr long double scale =
        1.0L / static_cast<long double>(std::numeric_limits<uint64_t>::max());
    const long double u =
        static_cast<long double>(splitmix64(key)) * scale;
    return static_cast<RealType>(u - 0.5L);
}

template<typename DataType>
DataType make_verification_entry(int64_t global_index, uint64_t seed) {
    const uint64_t base =
        splitmix64(seed ^ static_cast<uint64_t>(global_index));
    if constexpr (std::is_same_v<DataType, std::complex<double>>) {
        return DataType(
            centered_uniform_from_hash<double>(base),
            centered_uniform_from_hash<double>(base ^ 0x6a09e667f3bcc909ULL));
    } else if constexpr (std::is_same_v<DataType, std::complex<float>>) {
        return DataType(
            centered_uniform_from_hash<float>(base),
            centered_uniform_from_hash<float>(base ^ 0x6a09e667f3bcc909ULL));
    } else {
        return centered_uniform_from_hash<DataType>(base);
    }
}

// Step 1: pick `num_src` distinct random columns in [0, N) and give each a random weight.
// Fully deterministic in (N, num_src, seed): every rank calls this and gets the identical
// result, so no MPI_Bcast is needed to agree on x.
template<typename DataType>
SparseTestVector<DataType> make_sparse_test_vector(int64_t N, int64_t Npt_src, uint64_t seed) {
    SparseTestVector<DataType> stv;

    stv.idx.reserve(static_cast<size_t>(Npt_src));
    stv.weight.reserve(static_cast<size_t>(Npt_src));

    std::unordered_set<int64_t> seen;
    seen.reserve(static_cast<size_t>(Npt_src) * 2);

    // draw distinct candidate columns via the hash; dedup until we have Npt_src of them
    for (uint64_t t = 0; static_cast<int64_t>(stv.idx.size()) < Npt_src; ++t) {
        const int64_t cand = static_cast<int64_t>(
            splitmix64(seed ^ (0xD1B54A32D192ED03ULL * (t + 1))) % static_cast<uint64_t>(N));
        if (seen.insert(cand).second) {
            stv.idx.push_back(cand);
            stv.weight.push_back(make_verification_entry<DataType>(cand, seed));
        }
    }
    return stv;
}

template<typename DataType>
struct SparseMvpVerificationData {
    std::vector<DataType> input;
    std::vector<DataType> exact_output;
};

template<typename CoordType, typename DataType, typename KernelType>
SparseMvpVerificationData<DataType> make_sparse_mvp_verification_data(
    ParallelTree<CoordType, DataType>* tree,
    KernelType* kernel,
    int num_src,
    unsigned seed) {

    const int64_t N = tree->num_points;
    const int64_t Npt_src = std::min<int64_t>(num_src, N);
    const SparseTestVector<DataType> stv =
        make_sparse_test_vector<DataType>(N, Npt_src, seed);
    kernel->ensure_coordinates_collective(stv.idx, tree->comm);

    std::unordered_map<int64_t, DataType> sparse_values;
    sparse_values.reserve(static_cast<size_t>(Npt_src) * 2);
    for (int64_t k = 0; k < Npt_src; ++k) {
        sparse_values.emplace(stv.idx[static_cast<size_t>(k)],
                              stv.weight[static_cast<size_t>(k)]);
    }

    const int leaf_level = tree->num_levels - 1;
    auto& leaf = tree->levels[leaf_level];
    SparseMvpVerificationData<DataType> verification;
    std::vector<int64_t> local_rows;
    for (int64_t b = 0; b < leaf.num_boxes_local; ++b) {
        const auto& box = leaf.local_boxes[b];
        for (int64_t i = 0; i < box.num_points; ++i) {
            const int64_t global_row = box.point_indices[i];
            local_rows.push_back(global_row);
            const auto value = sparse_values.find(global_row);
            verification.input.push_back(
                value == sparse_values.end() ? DataType{0.0} : value->second);
        }
    }

    const int64_t Nloc = static_cast<int64_t>(local_rows.size());
    std::vector<DataType> block(
        static_cast<size_t>(Nloc) * static_cast<size_t>(Npt_src));
    if (Nloc > 0 && Npt_src > 0) {
        kernel->evaluate_block_by_index(local_rows.data(), Nloc,
                                        stv.idx.data(), Npt_src,
                                        block.data(), Nloc);
    }

    verification.exact_output.assign(static_cast<size_t>(Nloc), DataType{0.0});
    for (int64_t row = 0; row < Nloc; ++row) {
        for (int64_t k = 0; k < Npt_src; ++k) {
            verification.exact_output[static_cast<size_t>(row)] +=
                stv.weight[static_cast<size_t>(k)] *
                block[static_cast<size_t>(row + k * Nloc)];
        }
    }
    return verification;
}

} // namespace butterfly
