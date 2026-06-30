#pragma once

#include <vector>
#include <complex>
#include <stdexcept>
#include "tree.hpp"
#include "kernel.hpp"
#include "morton.hpp"
#include "serialization.hpp"
#include <mpi.h>
#include <vector>
#include <cstdint>
#include <cmath>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <complex>
#include "blas_declare.hpp"

/*used blas functions*/
// dpotrs_, dgemv_

namespace fmm {

// ============================================================================
// SOLVE PHASE FUNCTIONS
// ============================================================================

/**
 * @brief Apply V^{-1} operator (forward sweep)
 * 
 * Applies V^{-1} = L^{-1} U_T^{-1} P^* to the left-hand side vector.
 * Assumes left_side is already initialized as a copy of right_side.
 * 
 * Mathematical operations:
 * 1. L^{-1}: x[S] += T * x[R]                      (interpolation)
 * 2. U_T^{-1}: x[S] += X_SR * x[R]                 (X_SR already has -X_SR*X_RR^{-1})
 *              x[N] += X_NR * x[R]                 (X_NR already has -X_NR*X_RR^{-1})
 * 
 * @param level Tree level containing factorization data
 * @param solve_data Solve data for current box (contains left_side vector)
 * @param level_solve_data All solve data for this level (for neighbor access)
 * @param matrix_property Matrix symmetry property
 * @param is_ghost True if box is a ghost box (use solve_data's matrices)
 */
// template<typename CoordType, typename DataType>
// void apply_forward_elimination(
//     TreeLevel<CoordType, DataType>& level,
//     SolveDataRequest<CoordType, DataType>& solve_data,
//     std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
//     MatrixProperty matrix_property,
//     bool is_ghost) {
    
//     // Get factorization matrices (from solve_data if ghost, from level if local)
//     const MatrixStorage<DataType>* T = nullptr;
//     const MatrixStorage<DataType>* X_SR = nullptr;
//     const MatrixStorage<DataType>* X_NR = nullptr;
//     const std::vector<int64_t>* skeleton_indices = nullptr;
//     const std::vector<int64_t>* redundant_indices = nullptr;
//     BoxData<CoordType, DataType>* box = nullptr;
    
//     if (is_ghost) {
//         // Ghost box: use factorization data from solve_data
//         T = &solve_data.interpolation_matrix;
//         X_SR = &solve_data.X_SR;
//         X_NR = &solve_data.X_NR;
//         skeleton_indices = &solve_data.skeleton_indices;
//         redundant_indices = &solve_data.redundant_indices;
//     } else {
//         // Local box: find in level's local_boxes
//         box = level.find_local_box(solve_data.morton_index);
//         if (box == nullptr) {
//             throw std::runtime_error(
//                 "apply_forward_elimination: Local box " + 
//                 std::to_string(solve_data.morton_index) + " not found");
//         }
        
//         T = &box->interpolation_matrix;
//         X_SR = &box->X_SR;
//         X_NR = &box->X_NR;
//         skeleton_indices = &box->skeleton_indices;
//         redundant_indices = &box->redundant_indices;
//     }
    
//     // Check if box has been eliminated
//     if (skeleton_indices->empty() || redundant_indices->empty()) {
//         return;  // Box not yet factored, skip
//     }
    
//     int64_t k = skeleton_indices->size();  // Skeleton DOFs
//     int64_t r = redundant_indices->size(); // Redundant DOFs
    
//     // Extract x[R]
//     std::vector<DataType> x_R(r);
//     for (int64_t i = 0; i < r; ++i) {
//         x_R[i] = solve_data.left_side[(*redundant_indices)[i]];
//     }
    
//     // ===== Step 1: Apply L^{-1}: x[S] += T * x[R] =====
    
//     if (T->is_allocated()) {
//         // T is (k × r): maps redundant DOFs to skeleton DOFs
//         // x[S] (size k) += T (k × r) * x[R] (size r)
        
//         std::vector<DataType> result(k, DataType{0.0});
        
//         // DGEMV: y = alpha * A * x + beta * y
//         char trans = 'N';  // ← FIX: Use 'N' not 'T'
//         int m = k, n = r;  // T is (k × r)
//         DataType alpha = 1.0, beta = 0.0;
//         int lda = k, incx = 1, incy = 1;
        
//         if constexpr (std::is_same_v<DataType, double>) {
//             dgemv_(&trans, &m, &n, &alpha,
//                    T->data.data(), &lda,
//                    x_R.data(), &incx,
//                    &beta, result.data(), &incy);
//         } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
//             zgemv_(&trans, &m, &n, &alpha,
//                    T->data.data(), &lda,
//                    x_R.data(), &incx,
//                    &beta, result.data(), &incy);
//         }
        
//         // Add result to x[S]
//         for (int64_t i = 0; i < k; ++i) {
//             solve_data.left_side[(*skeleton_indices)[i]] -= result[i];
//         }
//     }
    
//     // ===== Step 2a: Apply U_T^{-1}: x[S] += X_SR * x[R] =====
    
//     if (X_SR->is_allocated()) {
//         // X_SR is (k × r), already stores -X_SR * X_RR^{-1}
//         std::vector<DataType> result(k, DataType{0.0});
        
//         char trans = 'N';
//         int m = k, n = r;
//         DataType alpha = 1.0, beta = 0.0;
//         int lda = k, incx = 1, incy = 1;
        
//         if constexpr (std::is_same_v<DataType, double>) {
//             dgemv_(&trans, &m, &n, &alpha,
//                    X_SR->data.data(), &lda,
//                    x_R.data(), &incx,
//                    &beta, result.data(), &incy);
//         } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
//             zgemv_(&trans, &m, &n, &alpha,
//                    X_SR->data.data(), &lda,
//                    x_R.data(), &incx,
//                    &beta, result.data(), &incy);
//         }
        
//         // Add result to x[S]
//         for (int64_t i = 0; i < k; ++i) {
//             solve_data.left_side[(*skeleton_indices)[i]] += result[i];
//         }
//     }
    
//     // ===== Step 2b: Apply U_T^{-1}: x[N] += X_NR * x[R] (batched) =====
    
//     if (X_NR->is_allocated() && box != nullptr) {
//         // X_NR stores all neighbor blocks concatenated vertically
//         // X_NR is (total_neighbor_points × r)
        
//         // Compute X_NR * x[R] in one DGEMV call
//         int64_t total_neighbor_points = X_NR->rows;
//         std::vector<DataType> neighbor_updates(total_neighbor_points, DataType{0.0});
        
//         char trans = 'N';
//         int m = total_neighbor_points, n = r;
//         DataType alpha = 1.0, beta = 0.0;
//         int lda = total_neighbor_points, incx = 1, incy = 1;
        
//         if constexpr (std::is_same_v<DataType, double>) {
//             dgemv_(&trans, &m, &n, &alpha,
//                    X_NR->data.data(), &lda,
//                    x_R.data(), &incx,
//                    &beta, neighbor_updates.data(), &incy);
//         } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
//             zgemv_(&trans, &m, &n, &alpha,
//                    X_NR->data.data(), &lda,
//                    x_R.data(), &incx,
//                    &beta, neighbor_updates.data(), &incy);
//         }
        
//         // Distribute updates to neighbors
//         int64_t row_offset = 0;
//         for (const auto& modified_block : box->near_field_modified_interactions) {
//             int64_t neighbor_morton = modified_block.neighbor_morton;
//             int64_t n_neighbor = modified_block.A_NS.rows;  // neighbor × self
            
//             // Find neighbor's solve data
//             SolveDataRequest<CoordType, DataType>* neighbor_data = nullptr;
            
//             // Check local boxes using level_solve_data
//             for (auto& solve_box : level_solve_data) {
//                 if (solve_box.morton_index == neighbor_morton) {
//                     neighbor_data = &solve_box;
//                     break;
//                 }
//             }
            
//             // If not found locally, check ghost/assisting boxes (for distributed case)
//             if (neighbor_data == nullptr) {
//                 auto it = level.ghost_and_assisting_box_points_for_solve_map.find(neighbor_morton);
//                 if (it != level.ghost_and_assisting_box_points_for_solve_map.end()) {
//                     neighbor_data = &level.ghost_and_assisting_boxes_for_solve[it->second];
//                 }
//             }
            
//             if (neighbor_data == nullptr) {
//                 row_offset += n_neighbor;
//                 continue;  // Neighbor not found
//             }
            
//             // Add updates to neighbor's skeleton DOFs
//             const auto& neighbor_skel = neighbor_data->skeleton_indices;
//             for (int64_t i = 0; i < n_neighbor; ++i) {
//                 if (i < neighbor_skel.size()) {
//                     neighbor_data->left_side[neighbor_skel[i]] += neighbor_updates[row_offset + i];
//                 }
//             }
            
//             row_offset += n_neighbor;
//         }
//     }
// }



// /**
//  * @brief Gather skeleton DOFs from children to parent (forward sweep transition)
//  * 
//  * After eliminating boxes at level n, aggregate skeleton DOFs to parent at level n-1.
//  * Each parent's DOFs = concatenation of children's skeleton DOFs.
//  */
// template<typename CoordType, typename DataType>
// void gather_skeleton_to_parent(
//     TreeLevel<CoordType, DataType>& child_level,
//     TreeLevel<CoordType, DataType>& parent_level,
//     std::vector<SolveDataRequest<CoordType, DataType>>& child_solve_data,
//     std::vector<SolveDataRequest<CoordType, DataType>>& parent_solve_data) {
    
//     // For each parent box, gather from its children
//     for (int64_t parent_idx = 0; parent_idx < parent_level.num_boxes_local; ++parent_idx) {
//         auto& parent_box = parent_level.local_boxes[parent_idx];
//         auto& parent_solve = parent_solve_data[parent_idx];
        
//         int64_t parent_offset = 0;
        
//         // Iterate through children
//         for (int32_t child_i = 0; child_i < parent_box.num_children; ++child_i) {
//             int64_t child_morton = parent_box.children_morton[child_i];
            
//             // Find child in child_level
//             BoxData<CoordType, DataType>* child_box = child_level.find_local_box(child_morton);
//             if (child_box == nullptr) {
//                 throw std::runtime_error("gather_skeleton_to_parent: Child box " + 
//                                        std::to_string(child_morton) + " not found");
//             }
            
//             // Find child's solve data
//             SolveDataRequest<CoordType, DataType>* child_solve = nullptr;
//             // for (auto& solve_box : child_solve_data) {
//             //     if (solve_box.morton_index == child_morton) {
//             //         child_solve = &solve_box;
//             //         break;
//             //     }
//             // }
//             if (child_morton >= child_level.local_morton_start && child_morton <= child_level.local_morton_end) {
//                 child_solve = &child_solve_data[child_morton - child_level.local_morton_start];
//             }
            
//             if (child_solve == nullptr) {
//                 throw std::runtime_error("gather_skeleton_to_parent: Child solve data not found");
//             }
            
//             // Copy child's skeleton DOFs to parent
//             const auto& child_skel = child_box->skeleton_indices;
//             for (int64_t i = 0; i < child_skel.size(); ++i) {
//                 parent_solve.left_side[parent_offset + i] = child_solve->left_side[child_skel[i]];
//             }
            
//             parent_offset += child_skel.size();
//         }
        
//         // Sanity check: parent's num_points should equal sum of children's skeleton sizes
//         if (parent_offset != parent_box.num_points) {
//             throw std::runtime_error(
//                 "gather_skeleton_to_parent: Size mismatch for parent " + 
//                 std::to_string(parent_box.morton_index) + 
//                 ": expected " + std::to_string(parent_box.num_points) + 
//                 ", got " + std::to_string(parent_offset));
//         }
//     }
// }



template<typename CoordType, typename DataType>
size_t get_parent_gather_transfer_size(
    const SolveDataRequest<CoordType, DataType>& solve_box) {
    return sizeof(int64_t) +
           sizeof(size_t) + solve_box.left_side.size() * sizeof(DataType) +
           sizeof(size_t) + solve_box.right_side.size() * sizeof(DataType);
}

template<typename CoordType, typename DataType>
char* serialize_parent_gather_transfer(
    const SolveDataRequest<CoordType, DataType>& solve_box,
    char* ptr) {
    std::memcpy(ptr, &solve_box.morton_index, sizeof(int64_t));
    ptr += sizeof(int64_t);

    auto serialize_vector = [&ptr](const auto& vec) {
        size_t vec_size = vec.size();
        std::memcpy(ptr, &vec_size, sizeof(size_t));
        ptr += sizeof(size_t);
        if (vec_size > 0) {
            using ValueType = typename std::decay_t<decltype(vec)>::value_type;
            std::memcpy(ptr, vec.data(), vec_size * sizeof(ValueType));
            ptr += vec_size * sizeof(ValueType);
        }
    };

    serialize_vector(solve_box.left_side);
    serialize_vector(solve_box.right_side);
    return ptr;
}

template<typename CoordType, typename DataType>
const char* deserialize_parent_gather_transfer(
    SolveDataRequest<CoordType, DataType>& solve_box,
    const char* ptr) {
    std::memcpy(&solve_box.morton_index, ptr, sizeof(int64_t));
    ptr += sizeof(int64_t);

    auto deserialize_vector = [&ptr](auto& vec) {
        size_t vec_size;
        std::memcpy(&vec_size, ptr, sizeof(size_t));
        ptr += sizeof(size_t);
        if (vec_size > 0) {
            using ValueType = typename std::decay_t<decltype(vec)>::value_type;
            vec.resize(vec_size);
            std::memcpy(vec.data(), ptr, vec_size * sizeof(ValueType));
            ptr += vec_size * sizeof(ValueType);
        } else {
            vec.clear();
        }
    };

    deserialize_vector(solve_box.left_side);
    deserialize_vector(solve_box.right_side);
    return ptr;
}

template<typename CoordType, typename DataType>
size_t get_parent_scatter_transfer_size(
    const SolveDataRequest<CoordType, DataType>& solve_box) {
    return sizeof(int64_t) +
           sizeof(size_t) + solve_box.left_side.size() * sizeof(DataType);
}

template<typename CoordType, typename DataType>
char* serialize_parent_scatter_transfer(
    const SolveDataRequest<CoordType, DataType>& solve_box,
    char* ptr) {
    std::memcpy(ptr, &solve_box.morton_index, sizeof(int64_t));
    ptr += sizeof(int64_t);

    size_t vec_size = solve_box.left_side.size();
    std::memcpy(ptr, &vec_size, sizeof(size_t));
    ptr += sizeof(size_t);
    if (vec_size > 0) {
        std::memcpy(ptr, solve_box.left_side.data(), vec_size * sizeof(DataType));
        ptr += vec_size * sizeof(DataType);
    }
    return ptr;
}

template<typename CoordType, typename DataType>
const char* deserialize_parent_scatter_transfer(
    SolveDataRequest<CoordType, DataType>& solve_box,
    const char* ptr) {
    std::memcpy(&solve_box.morton_index, ptr, sizeof(int64_t));
    ptr += sizeof(int64_t);

    size_t vec_size;
    std::memcpy(&vec_size, ptr, sizeof(size_t));
    ptr += sizeof(size_t);
    if (vec_size > 0) {
        solve_box.left_side.resize(vec_size);
        std::memcpy(solve_box.left_side.data(), ptr, vec_size * sizeof(DataType));
        ptr += vec_size * sizeof(DataType);
    } else {
        solve_box.left_side.clear();
    }
    solve_box.right_side.clear();
    return ptr;
}

template<typename CoordType, typename DataType>
void gather_skeleton_to_parent(
    TreeLevel<CoordType, DataType>& child_level,
    TreeLevel<CoordType, DataType>& parent_level,
    std::vector<SolveDataRequest<CoordType, DataType>>& child_solve_data,
    std::vector<SolveDataRequest<CoordType, DataType>>& parent_solve_data,
    int dimension) {
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const bool child_active = child_level.is_process_active;
    const bool parent_active = parent_level.is_process_active;
    
    // ===== Step 1: Build parent boxes from LOCAL children =====
    
    int children_per_parent = (dimension == 2) ? 4 : 8;
    
    int64_t num_local_parents = child_level.num_boxes_local / children_per_parent;
    
    std::vector<SolveDataRequest<CoordType, DataType>> local_parent_solve;
    local_parent_solve.resize(num_local_parents);
    
    for (int64_t parent_i = 0; parent_i < num_local_parents; ++parent_i) {
        int64_t child_start_idx = parent_i * children_per_parent;
        
        auto& first_child = child_level.local_boxes[child_start_idx];
        int64_t parent_morton = first_child.parent_morton;
        
        auto& parent_solve = local_parent_solve[parent_i];
        parent_solve.morton_index = parent_morton;
        parent_solve.source_rank = child_level.parent_level_owner;
        assert(child_level.parent_level_owner != -1);
        
        // Calculate total skeleton size
        int64_t total_skel_size = 0;
        for (int child_i = 0; child_i < children_per_parent; ++child_i) {
            auto& child_box = child_level.local_boxes[child_start_idx + child_i];
            total_skel_size += child_box.skeleton_indices.size();
        }
        
        parent_solve.left_side.resize(total_skel_size, DataType{0.0});
        parent_solve.right_side.resize(total_skel_size, DataType{0.0});
        
        // Copy metadata from parent box if we own it
        BoxData<CoordType, DataType>* parent_box = parent_level.find_local_box(parent_morton);
        if (parent_box != nullptr) {
            parent_solve.skeleton_indices = parent_box->skeleton_indices;
            parent_solve.redundant_indices = parent_box->redundant_indices;
            parent_solve.one_hop = parent_box->one_hop;
            parent_solve.use_full_set = parent_box->use_full_set;
        }
        
        // Gather skeleton DOFs from children
        int64_t parent_offset = 0;
        
        for (int child_i = 0; child_i < children_per_parent; ++child_i) {
            int64_t child_idx = child_start_idx + child_i;
            auto& child_box = child_level.local_boxes[child_idx];
            auto& child_solve = child_solve_data[child_idx];
            
            const auto& child_skel = child_box.skeleton_indices;
            for (int64_t i = 0; i < child_skel.size(); ++i) {
                parent_solve.left_side[parent_offset + i] = child_solve.left_side[child_skel[i]];
                parent_solve.right_side[parent_offset + i] = child_solve.right_side[child_skel[i]];
            }
            
            parent_offset += child_skel.size();
        }
    }

    const bool reduction_occurred =
        (parent_level.num_active_processes != child_level.num_active_processes);

    if (!reduction_occurred) {
        if (child_active && parent_active) {
            parent_solve_data = std::move(local_parent_solve);
        } else {
            parent_solve_data.clear();
        }
    } else if (parent_active) {
        std::vector<SolveDataRequest<CoordType, DataType>> all_parent_solve;

        for (int child_rank : parent_level.children_senders) {
            if (child_active && child_rank == rank) {
                all_parent_solve.insert(
                    all_parent_solve.end(),
                    std::make_move_iterator(local_parent_solve.begin()),
                    std::make_move_iterator(local_parent_solve.end())
                );
                continue;
            }

            int64_t buffer_size = 0;
            MPI_Status status;
            MPI_Recv(&buffer_size, 1, MPI_INT64_T, child_rank, 400, MPI_COMM_WORLD, &status);

            std::vector<char> recv_buffer(static_cast<size_t>(buffer_size));
            MPI_Recv(recv_buffer.data(), buffer_size, MPI_CHAR, child_rank, 401, MPI_COMM_WORLD, &status);

            const char* ptr = recv_buffer.data();
            const char* buffer_end = ptr + buffer_size;

            int64_t num_boxes;
            std::memcpy(&num_boxes, ptr, sizeof(int64_t));
            ptr += sizeof(int64_t);

            for (int64_t i = 0; i < num_boxes; ++i) {
                SolveDataRequest<CoordType, DataType> parent_solve_box;
                parent_solve_box.source_rank = rank;
                ptr = deserialize_parent_gather_transfer(parent_solve_box, ptr);
                all_parent_solve.push_back(std::move(parent_solve_box));
            }

            if (ptr != buffer_end) {
                throw std::runtime_error(
                    "gather_skeleton_to_parent: Deserialization size mismatch from rank " +
                    std::to_string(child_rank));
            }
        }

        parent_solve_data = std::move(all_parent_solve);

        for (auto& parent_solve : parent_solve_data) {
            int64_t parent_morton = parent_solve.morton_index;
            BoxData<CoordType, DataType>* parent_box = parent_level.find_local_box(parent_morton);

            if (parent_box == nullptr) {
                throw std::runtime_error(
                    "gather_skeleton_to_parent: Parent box " +
                    std::to_string(parent_morton) +
                    " not found in parent_level on owner rank " + std::to_string(rank));
            }

            parent_solve.skeleton_indices = parent_box->skeleton_indices;
            parent_solve.redundant_indices = parent_box->redundant_indices;
            parent_solve.one_hop = parent_box->one_hop;
            parent_solve.use_full_set = parent_box->use_full_set;
            parent_solve.source_rank = rank;
        }
    } else if (child_active) {
        int64_t num_boxes = local_parent_solve.size();
        size_t total_size = sizeof(int64_t);
        
        for (const auto& solve_box : local_parent_solve) {
            total_size += get_parent_gather_transfer_size(solve_box);
        }
        
        std::vector<char> send_buffer(total_size);
        char* ptr = send_buffer.data();
        
        std::memcpy(ptr, &num_boxes, sizeof(int64_t));
        ptr += sizeof(int64_t);
        
        for (const auto& solve_box : local_parent_solve) {
            ptr = serialize_parent_gather_transfer(solve_box, ptr);
        }
        
        int64_t buffer_size = send_buffer.size();
        MPI_Send(&buffer_size, 1, MPI_INT64_T, child_level.parent_level_owner, 400, MPI_COMM_WORLD);
        MPI_Send(send_buffer.data(), buffer_size, MPI_CHAR, child_level.parent_level_owner, 401, MPI_COMM_WORLD);
        
        parent_solve_data.clear();
    } else {
        parent_solve_data.clear();
    }
}

// /**
//  * @brief Scatter solution from parent to children's skeleton (backward sweep transition)
//  * 
//  * After solving at level n-1, distribute parent's solution to children's skeleton DOFs at level n.
//  */
// template<typename CoordType, typename DataType>
// void scatter_solution_to_children(
//     TreeLevel<CoordType, DataType>& child_level,
//     TreeLevel<CoordType, DataType>& parent_level,
//     std::vector<SolveDataRequest<CoordType, DataType>>& child_solve_data,
//     std::vector<SolveDataRequest<CoordType, DataType>>& parent_solve_data) {
    
//     // For each parent box, scatter to its children
//     for (int64_t parent_idx = 0; parent_idx < parent_level.num_boxes_local; ++parent_idx) {
//         auto& parent_box = parent_level.local_boxes[parent_idx];
//         auto& parent_solve = parent_solve_data[parent_idx];
        
//         int64_t parent_offset = 0;
        
//         // Iterate through children
//         for (int32_t child_i = 0; child_i < parent_box.num_children; ++child_i) {
//             int64_t child_morton = parent_box.children_morton[child_i];
//             // printf("child morton: %d\n", child_morton);
//             // Find child in child_level
//             BoxData<CoordType, DataType>* child_box = child_level.find_local_box(child_morton);
//             if (child_box == nullptr) {
//                 throw std::runtime_error("scatter_solution_to_children: Child box " + 
//                                        std::to_string(child_morton) + " not found");
//             }
            
//             // Find child's solve data
//             SolveDataRequest<CoordType, DataType>* child_solve = nullptr;
//             // for (auto& solve_box : child_solve_data) {
//             //     if (solve_box.morton_index == child_morton) {
//             //         child_solve = &solve_box;
//             //         break;
//             //     }
//             // }
//             if (child_morton >= child_level.local_morton_start && child_morton <= child_level.local_morton_end) {
//                 child_solve = &child_solve_data[child_morton - child_level.local_morton_start];
//             }
            
//             if (child_solve == nullptr) {
//                 throw std::runtime_error("scatter_solution_to_children: Child solve data not found");
//             }
            
//             // Copy parent's solution to child's skeleton DOFs
//             const auto& child_skel = child_box->skeleton_indices;
//             for (int64_t i = 0; i < child_skel.size(); ++i) {
//                 child_solve->left_side[child_skel[i]] = parent_solve.left_side[parent_offset + i];
//             }
            
//             parent_offset += child_skel.size();
//         }
//     }
// }



/**
 * @brief Scatter solution from parent to children's skeleton with MPI process expansion
 * 
 * After solving at level n-1, distribute parent's solution to children's skeleton DOFs at level n.
 * Handles two cases:
 * 
 * 1. No process reduction (same # of processes at both levels):
 *    - Direct local scatter from parent boxes to child boxes
 *    - No MPI communication needed
 * 
 * 2. Process expansion occurred (parent owned by fewer processes):
 *    - Parent owner distributes parent boxes to child processes
 *    - Each child process receives its portion and scatters locally
 *    - Example: [0] (parent) → [0,1,2,3] (children)
 * 
 * Data flow (when reduction occurred):
 *    Step 1: Parent owner splits parent_solve_data among children_senders
 *    Step 2: Parent owner sends/keeps parent solve data (with left_side values)
 *    Step 3: All processes scatter received parent.left_side → child.left_side[skeleton]
 * 
 * @tparam CoordType Coordinate type (float or double)
 * @tparam DataType Data type for solution values
 * 
 * @param child_level Child tree level (level n) - where solution is scattered TO
 * @param parent_level Parent tree level (level n-1) - where solution comes FROM
 * @param child_solve_data Solve data at child level (will be updated with parent values)
 * @param parent_solve_data Solve data at parent level (complete, owned by parent_level_owner)
 * @param dimension Spatial dimension (2 or 3) - determines children per parent (4 or 8)
 * 
 * @note Unlike gather_skeleton_to_parent, no metadata filling needed because parent owner
 *       sends complete data (left_side values) directly to children
 * 
 * @throws std::runtime_error if child boxes not found or deserialization fails
 */
template<typename CoordType, typename DataType>
void scatter_solution_to_children(
    TreeLevel<CoordType, DataType>& child_level,
    TreeLevel<CoordType, DataType>& parent_level,
    std::vector<SolveDataRequest<CoordType, DataType>>& child_solve_data,
    std::vector<SolveDataRequest<CoordType, DataType>>& parent_solve_data, int dimension) {
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const bool child_active = child_level.is_process_active;
    const bool parent_active = parent_level.is_process_active;
    
    bool reduction_occurred =
        (parent_level.num_active_processes != child_level.num_active_processes);
    
    if (!reduction_occurred) {
        if (!child_active || !parent_active) {
            return;
        }
        // No process reduction - simple local scatter
        for (int64_t parent_idx = 0; parent_idx < parent_level.num_boxes_local; ++parent_idx) {
            auto& parent_box = parent_level.local_boxes[parent_idx];
            auto& parent_solve = parent_solve_data[parent_idx];
            
            int64_t parent_offset = 0;
            
            for (int32_t child_i = 0; child_i < parent_box.num_children; ++child_i) {
                int64_t child_morton = parent_box.children_morton[child_i];
                
                BoxData<CoordType, DataType>* child_box = child_level.find_local_box(child_morton);
                if (child_box == nullptr) {
                    throw std::runtime_error(
                        "scatter_solution_to_children: Child box " + 
                        std::to_string(child_morton) + " not found");
                }
                
                int64_t child_idx = child_morton - child_level.local_morton_start;
                auto& child_solve = child_solve_data[child_idx];
                
                const auto& child_skel = child_box->skeleton_indices;
                for (int64_t i = 0; i < child_skel.size(); ++i) {
                    child_solve.left_side[child_skel[i]] = parent_solve.left_side[parent_offset + i];
                }
                
                parent_offset += child_skel.size();
            }
        }
        return;
    }

    
    // ===== Step 2: Reduction occurred - send parent boxes, then scatter locally =====
    
    std::vector<SolveDataRequest<CoordType, DataType>> local_parent_data;
    
    if (parent_active) {
        
        // This process owns the parent - send parent boxes to child processes
        
        int num_children_processes = parent_level.children_senders.size();
        int64_t num_parent_boxes = parent_solve_data.size();
        int64_t boxes_per_child = num_parent_boxes / num_children_processes;

        
        
        for (size_t proc_idx = 0; proc_idx < num_children_processes; ++proc_idx) {
            int child_rank = parent_level.children_senders[proc_idx];
            
            
            // Calculate range of parent boxes for this child process
            int64_t start_idx = proc_idx * boxes_per_child;
            int64_t end_idx = (proc_idx == num_children_processes - 1) ? 
                              num_parent_boxes : (proc_idx + 1) * boxes_per_child;
            
            
            if (child_active && child_rank == rank) {
                // Move local parent data (use move iterator)
                local_parent_data.reserve(end_idx - start_idx);
                for (int64_t parent_idx = start_idx; parent_idx < end_idx; ++parent_idx) {
                    local_parent_data.push_back(std::move(parent_solve_data[parent_idx]));
                }
            } else {
                // Send parent boxes to child process
                int64_t num_boxes = end_idx - start_idx;
                
                // Calculate buffer size
                size_t total_size = sizeof(int64_t);  // num_boxes
                for (int64_t parent_idx = start_idx; parent_idx < end_idx; ++parent_idx) {
                    total_size += get_parent_scatter_transfer_size(parent_solve_data[parent_idx]);
                }
                
                std::vector<char> send_buffer(total_size);
                char* ptr = send_buffer.data();
                
                std::memcpy(ptr, &num_boxes, sizeof(int64_t));
                ptr += sizeof(int64_t);
                
                for (int64_t parent_idx = start_idx; parent_idx < end_idx; ++parent_idx) {
                    ptr = serialize_parent_scatter_transfer(parent_solve_data[parent_idx], ptr);
                }
                
                int64_t buffer_size = send_buffer.size();
                MPI_Send(&buffer_size, 1, MPI_INT64_T, child_rank, 500, MPI_COMM_WORLD);
                MPI_Send(send_buffer.data(), buffer_size, MPI_CHAR, child_rank, 501, MPI_COMM_WORLD);
            }
        }
        
    } else if (child_active) {
        // This process does NOT own the parent - receive parent boxes
        
        int64_t buffer_size = 0;
        MPI_Status status;
        MPI_Recv(&buffer_size, 1, MPI_INT64_T, child_level.parent_level_owner, 500, MPI_COMM_WORLD, &status);
        
        std::vector<char> recv_buffer(buffer_size);
        MPI_Recv(recv_buffer.data(), buffer_size, MPI_CHAR, child_level.parent_level_owner, 501, MPI_COMM_WORLD, &status);
        
        const char* ptr = recv_buffer.data();
        const char* buffer_end = ptr + buffer_size;
        
        int64_t num_boxes;
        std::memcpy(&num_boxes, ptr, sizeof(int64_t));
        ptr += sizeof(int64_t);
        
        local_parent_data.resize(num_boxes);
        for (int64_t i = 0; i < num_boxes; ++i) {
            ptr = deserialize_parent_scatter_transfer(local_parent_data[i], ptr);
        }
        
        if (ptr != buffer_end) {
            throw std::runtime_error(
                "scatter_solution_to_children: Deserialization size mismatch");
        }
    } else {
        return;
    }
    
    // ===== Step 3: All processes now scatter locally =====
    if (!child_active) {
        return;
    }

    int children_per_parent = (dimension == 2) ? 4 : 8;

    for (auto& parent_solve : local_parent_data) {
        int64_t parent_morton = parent_solve.morton_index;
        // printf("parent morton %d from rank: %d\n", parent_morton, rank);
        
        // Compute children morton indices directly
        int64_t first_child_morton = parent_morton * children_per_parent;
        
        int64_t parent_offset = 0;
        
        for (int child_i = 0; child_i < children_per_parent; ++child_i) {
            int64_t child_morton = first_child_morton + child_i;
            
            BoxData<CoordType, DataType>* child_box = child_level.find_local_box(child_morton);
            if (child_box == nullptr) {
                throw std::runtime_error(
                    "scatter_solution_to_children: Child box " + 
                    std::to_string(child_morton) + " not found locally");
            }
            
            int64_t child_idx = child_morton - child_level.local_morton_start;
            auto& child_solve = child_solve_data[child_idx];
            
            const auto& child_skel = child_box->skeleton_indices;
            for (int64_t i = 0; i < child_skel.size(); ++i) {
                child_solve.left_side[child_skel[i]] = parent_solve.left_side[parent_offset + i];
            }
            
            parent_offset += child_skel.size();
        }
    }
}




// Helper: first insertion is "replace", subsequent insertions are "+=" accumulate.
template <typename DataType>
static inline void accumulate_replace_then_add(
    std::unordered_map<int64_t, std::vector<DataType>>& mp,
    int64_t morton,
    const DataType* values,
    int64_t n)
{
    if (n < 0) throw std::runtime_error("accumulate_replace_then_add: n < 0");
    if (n == 0) return;

    auto it = mp.find(morton);
    if (it == mp.end()) {
        mp.emplace(morton, std::vector<DataType>(values, values + n)); // REPLACE on first insert
        return;
    }

    auto& dst = it->second;
    if (static_cast<int64_t>(dst.size()) != n) {
        throw std::runtime_error(
            "accumulate_replace_then_add: size mismatch for morton=" + std::to_string(morton) +
            " existing=" + std::to_string(dst.size()) + " incoming=" + std::to_string(n));
    }

    for (int64_t i = 0; i < n; ++i) dst[static_cast<size_t>(i)] += values[static_cast<size_t>(i)];
}

template <typename DataType>
static inline void merge_pending_solve(
    PendingSolveUpdates<DataType>& dst,
    const PendingSolveUpdates<DataType>& src)
{
    for (const auto& [morton, upd] : src.full_updates) {
        accumulate_replace_then_add(
            dst.full_updates,
            morton,
            upd.data(),
            static_cast<int64_t>(upd.size()));
    }

    for (const auto& [morton, upd] : src.skel_updates) {
        accumulate_replace_then_add(
            dst.skel_updates,
            morton,
            upd.data(),
            static_cast<int64_t>(upd.size()));
    }
}

/**
 * @brief Apply V^{-1} operator (forward sweep)
 * 
 * Applies V^{-1} = L U_T P^* to the vector.
 * 
 * Mathematical operations (applied in order):
 * 1. U_T:  x[R] -= T^T * x[S]              (update redundant using skeleton)
 * 2. L:    x[S] -= X_SR*X_RR^{-1} * x[R]   (update skeleton using redundant)
 *          x[N] -= X_NR*X_RR^{-1} * x[R]   (update neighbors using redundant)
 * 
 * Note: Stored X_SR = -X_SR*X_RR^{-1}, so we ADD instead of subtract
 *       use_full_set[i] indicates whether one_hop[i] should use full RHS (1) or skeleton (0)
 */
template<typename CoordType, typename DataType>
void apply_forward_elimination(
    TreeLevel<CoordType, DataType>& level,
    SolveDataRequest<CoordType, DataType>& solve_data,
    std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
    MatrixProperty matrix_property,
    PendingSolveUpdates<DataType>& pending_updates,
    bool is_ghost,
    bool defer_local_updates = false)
{
    (void)matrix_property;

    // Factorization matrices / metadata
    const MatrixStorage<DataType>* T = nullptr;
    const MatrixStorage<DataType>* X_SR = nullptr;
    const MatrixStorage<DataType>* X_NR = nullptr;
    const std::vector<int64_t>* skeleton_indices = nullptr;
    const std::vector<int64_t>* redundant_indices = nullptr;
    const std::vector<int64_t>* one_hop = nullptr;
    const std::vector<int64_t>* use_full_set = nullptr;

    if (is_ghost) {
        T = &solve_data.interpolation_matrix;
        X_SR = &solve_data.X_SR;
        X_NR = &solve_data.X_NR;
        skeleton_indices = &solve_data.skeleton_indices;
        redundant_indices = &solve_data.redundant_indices;
        one_hop = &solve_data.one_hop;
        use_full_set = &solve_data.use_full_set;
    } else {
        auto* box = level.find_local_box(solve_data.morton_index);
        if (!box) {
            throw std::runtime_error("apply_forward_elimination: Local box " +
                                     std::to_string(solve_data.morton_index) + " not found");
        }
        T = &box->interpolation_matrix;
        X_SR = &box->X_SR;
        X_NR = &box->X_NR;
        skeleton_indices = &box->skeleton_indices;
        redundant_indices = &box->redundant_indices;
        one_hop = &box->one_hop;
        use_full_set = &box->use_full_set;
    }

    if (!skeleton_indices || !redundant_indices || skeleton_indices->empty() || redundant_indices->empty())
        return;

    const int64_t k = static_cast<int64_t>(skeleton_indices->size());
    const int64_t r = static_cast<int64_t>(redundant_indices->size());

    // Extract x[S], x[R]
    std::vector<DataType> x_S(static_cast<size_t>(k));
    for (int64_t i = 0; i < k; ++i)
        x_S[static_cast<size_t>(i)] = solve_data.left_side[static_cast<size_t>((*skeleton_indices)[static_cast<size_t>(i)])];

    std::vector<DataType> x_R(static_cast<size_t>(r));
    for (int64_t i = 0; i < r; ++i)
        x_R[static_cast<size_t>(i)] = solve_data.left_side[static_cast<size_t>((*redundant_indices)[static_cast<size_t>(i)])];

    // ===== Step 1: x[R] -= T^T * x[S] =====
    if (T && T->is_allocated()) {
        std::vector<DataType> tmp(static_cast<size_t>(r), DataType{0});

        char trans = 'T';
        int m = static_cast<int>(k), n = static_cast<int>(r);
        DataType alpha = DataType{1}, beta = DataType{0};
        int lda = static_cast<int>(T->lda), incx = 1, incy = 1;

        if constexpr (std::is_same_v<DataType, double>) {
            dgemv_(&trans, &m, &n, &alpha, T->data.data(), &lda,
                   x_S.data(), &incx, &beta, tmp.data(), &incy);
        } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
            zgemv_(&trans, &m, &n, &alpha, T->data.data(), &lda,
                   x_S.data(), &incx, &beta, tmp.data(), &incy);
        }

        for (int64_t i = 0; i < r; ++i)
            solve_data.left_side[static_cast<size_t>((*redundant_indices)[static_cast<size_t>(i)])] -= tmp[static_cast<size_t>(i)];
    }

    // Re-extract x[R]
    for (int64_t i = 0; i < r; ++i)
        x_R[static_cast<size_t>(i)] =
            solve_data.left_side[static_cast<size_t>((*redundant_indices)[static_cast<size_t>(i)])];

    // ===== Step 2: x[S] += X_SR * x[R] =====
    if (X_SR && X_SR->is_allocated()) {
        std::vector<DataType> tmp(static_cast<size_t>(k), DataType{0});

        char trans = 'N';
        int m = static_cast<int>(k), n = static_cast<int>(r);
        DataType alpha = DataType{1}, beta = DataType{0};
        int lda = static_cast<int>(X_SR->lda), incx = 1, incy = 1;

        if constexpr (std::is_same_v<DataType, double>) {
            dgemv_(&trans, &m, &n, &alpha, X_SR->data.data(), &lda,
                   x_R.data(), &incx, &beta, tmp.data(), &incy);
        } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
            zgemv_(&trans, &m, &n, &alpha, X_SR->data.data(), &lda,
                   x_R.data(), &incx, &beta, tmp.data(), &incy);
        }

        for (int64_t i = 0; i < k; ++i)
            solve_data.left_side[static_cast<size_t>((*skeleton_indices)[static_cast<size_t>(i)])] += tmp[static_cast<size_t>(i)];
    }

    // ===== Step 3: neighbor updates =====
    if (!(X_NR && X_NR->is_allocated() && one_hop && use_full_set && !one_hop->empty()))
        return;

    // We now use level.solve_neighbor_size to split X_NR rows into per-neighbor segments.
    // This is assumed to be stored for LOCAL boxes, in one_hop order, already choosing full/skel length.
    if (is_ghost) {
        throw std::runtime_error("apply_forward_elimination: is_ghost=true but solve_neighbor_size is only defined for local boxes");
    }

    const int64_t local_offset = solve_data.morton_index - level.local_morton_start;
    if (local_offset < 0 || local_offset >= static_cast<int64_t>(level.solve_neighbor_size.size())) {
        throw std::runtime_error("apply_forward_elimination: local_offset out of range for solve_neighbor_size");
    }

    const auto& neighbor_sizes = level.solve_neighbor_size[static_cast<size_t>(local_offset)];
    if (neighbor_sizes.size() != one_hop->size()) {
        throw std::runtime_error(
            "apply_forward_elimination: solve_neighbor_size mismatch: sizes=" +
            std::to_string(neighbor_sizes.size()) + " one_hop=" + std::to_string(one_hop->size()));
    }

    const int64_t total_neighbor_points = X_NR->rows;
    std::vector<DataType> neighbor_updates(static_cast<size_t>(total_neighbor_points), DataType{0});

    {
        char trans = 'N';
        int m = static_cast<int>(total_neighbor_points), n = static_cast<int>(r);
        DataType alpha = DataType{1}, beta = DataType{0};
        int lda = static_cast<int>(X_NR->lda), incx = 1, incy = 1;

        if constexpr (std::is_same_v<DataType, double>) {
            dgemv_(&trans, &m, &n, &alpha, X_NR->data.data(), &lda,
                   x_R.data(), &incx, &beta, neighbor_updates.data(), &incy);
        } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
            zgemv_(&trans, &m, &n, &alpha, X_NR->data.data(), &lda,
                   x_R.data(), &incx, &beta, neighbor_updates.data(), &incy);
        }
    }

    auto is_local_morton = [&](int64_t morton) -> bool {
        return (morton >= level.local_morton_start && morton <= level.local_morton_end);
    };

    int64_t row_offset = 0;
    for (size_t neighbor_idx = 0; neighbor_idx < one_hop->size(); ++neighbor_idx) {
        const int64_t neighbor_morton = (*one_hop)[neighbor_idx];
        const bool use_full = ((*use_full_set)[neighbor_idx] == 1);

        const int64_t n_neighbor = neighbor_sizes[neighbor_idx];
        if (n_neighbor <= 0) {
            throw std::runtime_error("apply_forward_elimination: n_neighbor <= 0 for neighbor=" +
                                     std::to_string(neighbor_morton));
        }
        if (row_offset + n_neighbor > total_neighbor_points) {
            throw std::runtime_error("apply_forward_elimination: X_NR row partition overflow");
        }

        const DataType* seg = neighbor_updates.data() + row_offset;

        if (is_local_morton(neighbor_morton)) {
            auto& neighbor_data =
                level_solve_data[static_cast<size_t>(neighbor_morton - level.local_morton_start)];

            if (use_full) {
                if (static_cast<int64_t>(neighbor_data.left_side.size()) != n_neighbor) {
                    throw std::runtime_error("apply_forward_elimination: local full size mismatch for neighbor=" +
                                             std::to_string(neighbor_morton));
                }
                if (defer_local_updates) {
                    accumulate_replace_then_add(
                        pending_updates.full_updates,
                        neighbor_morton,
                        seg,
                        n_neighbor);
                } else {
                    for (int64_t i = 0; i < n_neighbor; ++i)
                        neighbor_data.left_side[static_cast<size_t>(i)] += seg[static_cast<size_t>(i)];
                }
            } else {
                if (static_cast<int64_t>(neighbor_data.skeleton_indices.size()) != n_neighbor) {
                    throw std::runtime_error("apply_forward_elimination: local skel size mismatch for neighbor=" +
                                             std::to_string(neighbor_morton));
                }
                if (defer_local_updates) {
                    accumulate_replace_then_add(
                        pending_updates.skel_updates,
                        neighbor_morton,
                        seg,
                        n_neighbor);
                } else {
                    const auto& skel = neighbor_data.skeleton_indices;
                    for (int64_t i = 0; i < n_neighbor; ++i)
                        neighbor_data.left_side[static_cast<size_t>(skel[static_cast<size_t>(i)])] += seg[static_cast<size_t>(i)];
                }
            }
        } else {
            // Nonlocal: just accumulate the segment; receiver will apply indexing if needed.
            if (use_full) {
                accumulate_replace_then_add(pending_updates.full_updates, neighbor_morton, seg, n_neighbor);
            } else {
                accumulate_replace_then_add(pending_updates.skel_updates, neighbor_morton, seg, n_neighbor);
            }
        }

        row_offset += n_neighbor;
    }

    if (row_offset != total_neighbor_points) {
        throw std::runtime_error(
            "apply_forward_elimination: row_offset " + std::to_string(row_offset) +
            " != total_neighbor_points " + std::to_string(total_neighbor_points));
    }
}


/**
 * @brief Apply W^{-1} operator (backward sweep)
 * 
 * Applies W^{-1} = P L_T U to the vector.
 * 
 * Mathematical operations (applied in order):
 * 1. U:    x[R] -= X_RR^{-1}X_RS * x[S]    (update redundant using skeleton)
 *          x[R] -= X_RR^{-1}X_RN * x[N]    (update redundant using neighbors)
 * 2. L_T:  x[S] -= T * x[R]                (update skeleton using redundant)
 * 
 * Note: Stored X_SR = -X_SR*X_RR^{-1} (symmetric) or X_RS = -X_RS*X_RR^{-1} (nonsymmetric)
 *       Stored X_NR = -X_NR*X_RR^{-1} (symmetric) or X_RN = -X_RN*X_RR^{-1} (nonsymmetric)
 *       use_full_set[i] indicates whether one_hop[i] should use full left_side (1) or skeleton (0)
 */
template<typename CoordType, typename DataType>
void apply_backward_substitution(
    TreeLevel<CoordType, DataType>& level,
    SolveDataRequest<CoordType, DataType>& solve_data,
    std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
    MatrixProperty matrix_property,
    bool is_ghost) {
    
    // Get factorization matrices
    const MatrixStorage<DataType>* T = nullptr;
    const MatrixStorage<DataType>* X_RR = nullptr;
    const MatrixStorage<DataType>* X_RS = nullptr;  // Only for NONSYMMETRIC
    const MatrixStorage<DataType>* X_RN = nullptr;  // Only for NONSYMMETRIC
    const MatrixStorage<DataType>* X_SR = nullptr;  // For SYMMETRIC (transpose)
    const MatrixStorage<DataType>* X_NR = nullptr;  // For SYMMETRIC (transpose)
    const std::vector<int64_t>* skeleton_indices = nullptr;
    const std::vector<int64_t>* redundant_indices = nullptr;
    const std::vector<int64_t>* one_hop = nullptr;
    const std::vector<int64_t>* use_full_set = nullptr;
    BoxData<CoordType, DataType>* box = nullptr;
    
    if (is_ghost) {
        T = &solve_data.interpolation_matrix;
        X_RR = &solve_data.X_RR;
        X_RS = &solve_data.X_RS;
        X_RN = &solve_data.X_RN;
        X_SR = &solve_data.X_SR;
        X_NR = &solve_data.X_NR;
        skeleton_indices = &solve_data.skeleton_indices;
        redundant_indices = &solve_data.redundant_indices;
        one_hop = &solve_data.one_hop;
        use_full_set = &solve_data.use_full_set;
    } else {
        box = level.find_local_box(solve_data.morton_index);
        if (box == nullptr) {
            throw std::runtime_error(
                "apply_backward_substitution: Local box " + 
                std::to_string(solve_data.morton_index) + " not found");
        }
        
        T = &box->interpolation_matrix;
        X_RR = &box->X_RR;
        X_RS = &box->X_RS;
        X_RN = &box->X_RN;
        X_SR = &box->X_SR;
        X_NR = &box->X_NR;
        skeleton_indices = &box->skeleton_indices;
        redundant_indices = &box->redundant_indices;
        one_hop = &box->one_hop;
        use_full_set = &box->use_full_set;
    }
    
    if (skeleton_indices->empty() || redundant_indices->empty()) {
        return;
    }
    
    int64_t k = skeleton_indices->size();
    int64_t r = redundant_indices->size();
    
    // ===== Step 1: Apply U: x[R] -= X_RR^{-1}X_RS * x[S] - X_RR^{-1}X_RN * x[N] =====
    
    std::vector<DataType> x_R(r);
    for (int64_t i = 0; i < r; ++i) {
        x_R[i] = solve_data.left_side[(*redundant_indices)[i]];
    }
    
    std::vector<DataType> x_S(k);
    for (int64_t i = 0; i < k; ++i) {
        x_S[i] = solve_data.left_side[(*skeleton_indices)[i]];
    }
    
    // U: x[R] -= X_RR^{-1}X_RS * x[S]
    // Stored X_SR = -X_SR*X_RR^{-1}, so: x[R] += stored_X_SR^T * x[S]
    assert(X_SR->is_allocated());
    if (matrix_property == MatrixProperty::SYMMETRIC) {
        if (X_SR->is_allocated()) {
            char trans = 'T';
            int m = k, n = r;
            DataType alpha = 1.0;  // ADD (double negative cancels)
            DataType beta = 1.0;
            int lda = k, incx = 1, incy = 1;
            
            if constexpr (std::is_same_v<DataType, double>) {
                dgemv_(&trans, &m, &n, &alpha,
                       X_SR->data.data(), &lda,
                       x_S.data(), &incx,
                       &beta, x_R.data(), &incy);
            } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                zgemv_(&trans, &m, &n, &alpha,
                       X_SR->data.data(), &lda,
                       x_S.data(), &incx,
                       &beta, x_R.data(), &incy);
            }
        }
    } else {
        // NONSYMMETRIC: Use stored X_RS directly
        if (X_RS->is_allocated()) {
            char trans = 'N';
            int m = r, n = k;
            DataType alpha = 1.0;  // ADD
            DataType beta = 1.0;
            int lda = r, incx = 1, incy = 1;
            
            if constexpr (std::is_same_v<DataType, double>) {
                dgemv_(&trans, &m, &n, &alpha,
                       X_RS->data.data(), &lda,
                       x_S.data(), &incx,
                       &beta, x_R.data(), &incy);
            } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                zgemv_(&trans, &m, &n, &alpha,
                       X_RS->data.data(), &lda,
                       x_S.data(), &incx,
                       &beta, x_R.data(), &incy);
            }
        }
    }
    
    // U: x[R] -= X_RR^{-1}X_RN * x[N]
    // Collect neighbor values in one_hop order, respecting use_full_set
    assert(!one_hop->empty());
    if (one_hop != nullptr && !one_hop->empty()) {
        
        if (matrix_property == MatrixProperty::SYMMETRIC && X_NR->is_allocated()) {
            int64_t total_neighbor_points = X_NR->rows;
            std::vector<DataType> neighbor_values(total_neighbor_points, DataType{0.0});
            
            // Collect neighbor values in one_hop order
            int64_t row_offset = 0;
            for (size_t neighbor_idx = 0; neighbor_idx < one_hop->size(); ++neighbor_idx) {
                int64_t neighbor_morton = (*one_hop)[neighbor_idx];
                bool use_full = ((*use_full_set)[neighbor_idx] == 1);
                
                // Find neighbor's solve data
                SolveDataRequest<CoordType, DataType>* neighbor_data = nullptr;
                
                // for (auto& solve_box : level_solve_data) {
                //     if (solve_box.morton_index == neighbor_morton) {
                //         neighbor_data = &solve_box;
                //         break;
                //     }
                // }
                if (neighbor_morton >= level.local_morton_start && neighbor_morton <= level.local_morton_end) {
                    neighbor_data = &level_solve_data[neighbor_morton - level.local_morton_start];
                }
                
                if (neighbor_data == nullptr) {
                    auto it = level.ghost_and_assisting_box_points_for_solve_map.find(neighbor_morton);
                    if (it != level.ghost_and_assisting_box_points_for_solve_map.end()) {
                        neighbor_data = &level.ghost_and_assisting_boxes_for_solve[it->second];
                    }
                }
                
                if (neighbor_data == nullptr) {
                    throw std::runtime_error(
                        "apply_backward_substitution: Neighbor " + 
                        std::to_string(neighbor_morton) + " not found in solve data");
                }
                
                int64_t n_neighbor;
                
                if (use_full) {
                    // Neighbor was eliminated AFTER current box
                    // X_NR was computed against neighbor's full point set
                    // Read from entire left_side vector
                    n_neighbor = neighbor_data->left_side.size();
                    
                    for (int64_t i = 0; i < n_neighbor; ++i) {
                        neighbor_values[row_offset + i] = neighbor_data->left_side[i];
                    }
                } else {
                    // Neighbor was eliminated BEFORE current box
                    // X_NR was computed against neighbor's skeleton DOFs only
                    // Read from skeleton indices
                    n_neighbor = neighbor_data->skeleton_indices.size();
                    const auto& neighbor_skel = neighbor_data->skeleton_indices;
                    
                    for (int64_t i = 0; i < n_neighbor; ++i) {
                        neighbor_values[row_offset + i] = neighbor_data->left_side[neighbor_skel[i]];
                    }
                }
                
                row_offset += n_neighbor;
            }
            if(row_offset != total_neighbor_points)
            {
                throw std::runtime_error(
                    "row_offset " + std::to_string(row_offset) + " does not match total_neighbor_points " +
                    std::to_string(total_neighbor_points));
            }
            
            // U: x[R] -= X_RR^{-1}X_RN*x[N] = x[R] += stored_X_NR^T*x[N]
            char trans = 'T';
            int m = total_neighbor_points, n = r;
            DataType alpha = 1.0;  // ADD
            DataType beta = 1.0;
            int lda = total_neighbor_points, incx = 1, incy = 1;
            
            if constexpr (std::is_same_v<DataType, double>) {
                dgemv_(&trans, &m, &n, &alpha,
                       X_NR->data.data(), &lda,
                       neighbor_values.data(), &incx,
                       &beta, x_R.data(), &incy);
            } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                zgemv_(&trans, &m, &n, &alpha,
                       X_NR->data.data(), &lda,
                       neighbor_values.data(), &incx,
                       &beta, x_R.data(), &incy);
            }
            
        } else if (matrix_property == MatrixProperty::NONSYMMETRIC && X_RN->is_allocated()) {
            // NONSYMMETRIC case
            int64_t total_neighbor_points = X_RN->cols;
            std::vector<DataType> neighbor_values(total_neighbor_points, DataType{0.0});
            
            // Collect neighbor values in one_hop order
            int64_t col_offset = 0;
            for (size_t neighbor_idx = 0; neighbor_idx < one_hop->size(); ++neighbor_idx) {
                int64_t neighbor_morton = (*one_hop)[neighbor_idx];
                bool use_full = ((*use_full_set)[neighbor_idx] == 1);
                
                // Find neighbor's solve data
                SolveDataRequest<CoordType, DataType>* neighbor_data = nullptr;
                
                // for (auto& solve_box : level_solve_data) {
                //     if (solve_box.morton_index == neighbor_morton) {
                //         neighbor_data = &solve_box;
                //         break;
                //     }
                // }
                if (neighbor_morton >= level.local_morton_start && neighbor_morton <= level.local_morton_end) {
                    neighbor_data = &level_solve_data[neighbor_morton - level.local_morton_start];
                }
                
                if (neighbor_data == nullptr) {
                    auto it = level.ghost_and_assisting_box_points_for_solve_map.find(neighbor_morton);
                    if (it != level.ghost_and_assisting_box_points_for_solve_map.end()) {
                        neighbor_data = &level.ghost_and_assisting_boxes_for_solve[it->second];
                    }
                }
                
                if (neighbor_data == nullptr) {
                    throw std::runtime_error(
                        "apply_backward_substitution: Neighbor " + 
                        std::to_string(neighbor_morton) + " not found in solve data");
                }
                
                int64_t n_neighbor;
                
                if (use_full) {
                    // Neighbor was eliminated AFTER current box
                    // X_RN was computed against neighbor's full point set
                    // Read from entire left_side vector
                    n_neighbor = neighbor_data->left_side.size();
                    
                    for (int64_t i = 0; i < n_neighbor; ++i) {
                        neighbor_values[col_offset + i] = neighbor_data->left_side[i];
                    }
                } else {
                    // Neighbor was eliminated BEFORE current box
                    // X_RN was computed against neighbor's skeleton DOFs only
                    // Read from skeleton indices
                    n_neighbor = neighbor_data->skeleton_indices.size();
                    const auto& neighbor_skel = neighbor_data->skeleton_indices;
                    
                    for (int64_t i = 0; i < n_neighbor; ++i) {
                        neighbor_values[col_offset + i] = neighbor_data->left_side[neighbor_skel[i]];
                    }
                }
                
                col_offset += n_neighbor;
            }

            if(col_offset != total_neighbor_points)
            {
                throw std::runtime_error(
                    "col_offset " + std::to_string(col_offset) + " does not match total_neighbor_points " +
                    std::to_string(total_neighbor_points));
            }
            
            // U: x[R] -= X_RR^{-1}X_RN*x[N] = x[R] += stored_X_RN*x[N]
            char trans = 'N';
            int m = r, n = total_neighbor_points;
            DataType alpha = 1.0;  // ADD
            DataType beta = 1.0;
            int lda = r, incx = 1, incy = 1;
            
            if constexpr (std::is_same_v<DataType, double>) {
                dgemv_(&trans, &m, &n, &alpha,
                       X_RN->data.data(), &lda,
                       neighbor_values.data(), &incx,
                       &beta, x_R.data(), &incy);
            } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
                zgemv_(&trans, &m, &n, &alpha,
                       X_RN->data.data(), &lda,
                       neighbor_values.data(), &incx,
                       &beta, x_R.data(), &incy);
            }
        }
    }
    
    // Store updated x[R]
    for (int64_t i = 0; i < r; ++i) {
        solve_data.left_side[(*redundant_indices)[i]] = x_R[i];
    }
    
    // ===== Step 2: Apply L_T: x[S] -= T * x[R] =====
    assert(T->is_allocated());
    if (T->is_allocated()) {
        // Re-extract x[R] after Step 1 updates
        for (int64_t i = 0; i < r; ++i) {
            x_R[i] = solve_data.left_side[(*redundant_indices)[i]];
        }
        
        std::vector<DataType> result(k, DataType{0.0});
        
        char trans = 'N';  // No transpose: T is (k×r)
        int m = k, n = r;
        DataType alpha = 1.0, beta = 0.0;
        int lda = k, incx = 1, incy = 1;
        
        if constexpr (std::is_same_v<DataType, double>) {
            dgemv_(&trans, &m, &n, &alpha,
                   T->data.data(), &lda,
                   x_R.data(), &incx,
                   &beta, result.data(), &incy);
        } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
            zgemv_(&trans, &m, &n, &alpha,
                   T->data.data(), &lda,
                   x_R.data(), &incx,
                   &beta, result.data(), &incy);
        }
        
        // L_T: x[S] -= T * x[R]
        for (int64_t i = 0; i < k; ++i) {
            solve_data.left_side[(*skeleton_indices)[i]] -= result[i];
        }
    }
}

// template<typename CoordType, typename DataType>
// void apply_backward_substitution(
//     TreeLevel<CoordType, DataType>& level,
//     SolveDataRequest<CoordType, DataType>& solve_data,
//     std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
//     MatrixProperty matrix_property,
//     bool is_ghost) {
    
//     // Get factorization matrices
//     const MatrixStorage<DataType>* T = nullptr;
//     const MatrixStorage<DataType>* X_RR = nullptr;
//     const MatrixStorage<DataType>* X_RS = nullptr;  // Only for NONSYMMETRIC
//     const MatrixStorage<DataType>* X_RN = nullptr;  // Only for NONSYMMETRIC
//     const MatrixStorage<DataType>* X_SR = nullptr;  // For SYMMETRIC (transpose)
//     const MatrixStorage<DataType>* X_NR = nullptr;  // For SYMMETRIC (transpose)
//     const std::vector<int64_t>* skeleton_indices = nullptr;
//     const std::vector<int64_t>* redundant_indices = nullptr;
//     BoxData<CoordType, DataType>* box = nullptr;
    
//     if (is_ghost) {
//         T = &solve_data.interpolation_matrix;
//         X_RR = &solve_data.X_RR;
//         X_RS = &solve_data.X_RS;
//         X_RN = &solve_data.X_RN;
//         X_SR = &solve_data.X_SR;
//         X_NR = &solve_data.X_NR;
//         skeleton_indices = &solve_data.skeleton_indices;
//         redundant_indices = &solve_data.redundant_indices;
//     } else {
//         box = level.find_local_box(solve_data.morton_index);
//         if (box == nullptr) {
//             throw std::runtime_error(
//                 "apply_backward_substitution: Local box " + 
//                 std::to_string(solve_data.morton_index) + " not found");
//         }
        
//         T = &box->interpolation_matrix;
//         X_RR = &box->X_RR;
//         X_RS = &box->X_RS;
//         X_RN = &box->X_RN;
//         X_SR = &box->X_SR;
//         X_NR = &box->X_NR;
//         skeleton_indices = &box->skeleton_indices;
//         redundant_indices = &box->redundant_indices;
//     }
    
//     if (skeleton_indices->empty() || redundant_indices->empty()) {
//         return;
//     }
    
//     int64_t k = skeleton_indices->size();
//     int64_t r = redundant_indices->size();
    
//     // ===== Step 1: Apply U (NOT U^{-1}) =====
//     // U * [x[R]; x[S]; x[N]; x[F]] gives:
//     // x[R]_new = x[R] - X_RR^{-1}X_RS * x[S] - X_RR^{-1}X_RN * x[N]
//     // But stored X_SR = -X_RR^{-1} * A_SR, so we ADD!
    
//     std::vector<DataType> x_R(r);
//     for (int64_t i = 0; i < r; ++i) {
//         x_R[i] = solve_data.left_side[(*redundant_indices)[i]];
//     }
    
//     std::vector<DataType> x_S(k);
//     for (int64_t i = 0; i < k; ++i) {
//         x_S[i] = solve_data.left_side[(*skeleton_indices)[i]];
//     }
    
//     // Add X_SR^T * x[S] (or X_RS * x[S] for nonsymmetric)
    
//     if (matrix_property == MatrixProperty::SYMMETRIC) {
//         if (X_SR->is_allocated()) {
//             // Stored X_SR = -X_RR^{-1} * A_SR
//             // We want: x[R] - X_RR^{-1} * A_RS = x[R] - X_RR^{-1} * A_SR^T
//             //        = x[R] - (-(stored_X_SR))^T = x[R] + stored_X_SR^T
            
//             char trans = 'T';
//             int m = k, n = r;
//             DataType alpha = 1.0;  // ← ADD (double negative cancels)
//             DataType beta = 1.0;
//             int lda = k, incx = 1, incy = 1;
            
//             if constexpr (std::is_same_v<DataType, double>) {
//                 dgemv_(&trans, &m, &n, &alpha,
//                        X_SR->data.data(), &lda,
//                        x_S.data(), &incx,
//                        &beta, x_R.data(), &incy);
//             } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
//                 zgemv_(&trans, &m, &n, &alpha,
//                        X_SR->data.data(), &lda,
//                        x_S.data(), &incx,
//                        &beta, x_R.data(), &incy);
//             }
//         }
//     } else {
//         // NONSYMMETRIC: Use stored X_RS directly
//         if (X_RS->is_allocated()) {
//             // Stored X_RS = -X_RR^{-1} * A_RS
//             // We want: x[R] - X_RR^{-1} * A_RS = x[R] + stored_X_RS
            
//             char trans = 'N';
//             int m = r, n = k;
//             DataType alpha = 1.0;  // ← ADD
//             DataType beta = 1.0;
//             int lda = r, incx = 1, incy = 1;
            
//             if constexpr (std::is_same_v<DataType, double>) {
//                 dgemv_(&trans, &m, &n, &alpha,
//                        X_RS->data.data(), &lda,
//                        x_S.data(), &incx,
//                        &beta, x_R.data(), &incy);
//             } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
//                 zgemv_(&trans, &m, &n, &alpha,
//                        X_RS->data.data(), &lda,
//                        x_S.data(), &incx,
//                        &beta, x_R.data(), &incy);
//             }
//         }
//     }
    
//     // Add X_NR^T * x[N] (or X_RN * x[N] for nonsymmetric) - batched
    
//     if (box != nullptr && !box->near_field_modified_interactions.empty()) {
        
//         if (matrix_property == MatrixProperty::SYMMETRIC && X_NR->is_allocated()) {
//             int64_t total_neighbor_points = X_NR->rows;
//             std::vector<DataType> neighbor_values(total_neighbor_points, DataType{0.0});
            
//             // Collect neighbor skeleton values
//             int64_t row_offset = 0;
//             for (const auto& modified_block : box->near_field_modified_interactions) {
//                 int64_t neighbor_morton = modified_block.neighbor_morton;
//                 int64_t n_neighbor = modified_block.A_NS.rows;
                
//                 SolveDataRequest<CoordType, DataType>* neighbor_data = nullptr;
                
//                 for (auto& solve_box : level_solve_data) {
//                     if (solve_box.morton_index == neighbor_morton) {
//                         neighbor_data = &solve_box;
//                         break;
//                     }
//                 }
                
//                 if (neighbor_data == nullptr) {
//                     auto it = level.ghost_and_assisting_box_points_for_solve_map.find(neighbor_morton);
//                     if (it != level.ghost_and_assisting_box_points_for_solve_map.end()) {
//                         neighbor_data = &level.ghost_and_assisting_boxes_for_solve[it->second];
//                     }
//                 }
                
//                 if (neighbor_data == nullptr) {
//                     row_offset += n_neighbor;
//                     continue;
//                 }
                
//                 const auto& neighbor_skel = neighbor_data->skeleton_indices;
//                 for (int64_t i = 0; i < n_neighbor && i < neighbor_skel.size(); ++i) {
//                     neighbor_values[row_offset + i] = neighbor_data->left_side[neighbor_skel[i]];
//                 }
                
//                 row_offset += n_neighbor;
//             }
            
//             // X_NR^T * neighbor_values
//             char trans = 'T';
//             int m = total_neighbor_points, n = r;
//             DataType alpha = 1.0;  // ← ADD
//             DataType beta = 1.0;
//             int lda = total_neighbor_points, incx = 1, incy = 1;
            
//             if constexpr (std::is_same_v<DataType, double>) {
//                 dgemv_(&trans, &m, &n, &alpha,
//                        X_NR->data.data(), &lda,
//                        neighbor_values.data(), &incx,
//                        &beta, x_R.data(), &incy);
//             } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
//                 zgemv_(&trans, &m, &n, &alpha,
//                        X_NR->data.data(), &lda,
//                        neighbor_values.data(), &incx,
//                        &beta, x_R.data(), &incy);
//             }
            
//         } else if (matrix_property == MatrixProperty::NONSYMMETRIC && X_RN->is_allocated()) {
//             // NONSYMMETRIC case
//             int64_t total_neighbor_points = X_RN->cols;
//             std::vector<DataType> neighbor_values(total_neighbor_points, DataType{0.0});
            
//             int64_t col_offset = 0;
//             for (const auto& modified_block : box->near_field_modified_interactions) {
//                 int64_t neighbor_morton = modified_block.neighbor_morton;
//                 int64_t n_neighbor = modified_block.A_NS.rows;
                
//                 SolveDataRequest<CoordType, DataType>* neighbor_data = nullptr;
                
//                 for (auto& solve_box : level_solve_data) {
//                     if (solve_box.morton_index == neighbor_morton) {
//                         neighbor_data = &solve_box;
//                         break;
//                     }
//                 }
                
//                 if (neighbor_data == nullptr) {
//                     auto it = level.ghost_and_assisting_box_points_for_solve_map.find(neighbor_morton);
//                     if (it != level.ghost_and_assisting_box_points_for_solve_map.end()) {
//                         neighbor_data = &level.ghost_and_assisting_boxes_for_solve[it->second];
//                     }
//                 }
                
//                 if (neighbor_data == nullptr) {
//                     col_offset += n_neighbor;
//                     continue;
//                 }
                
//                 const auto& neighbor_skel = neighbor_data->skeleton_indices;
//                 for (int64_t i = 0; i < n_neighbor && i < neighbor_skel.size(); ++i) {
//                     neighbor_values[col_offset + i] = neighbor_data->left_side[neighbor_skel[i]];
//                 }
                
//                 col_offset += n_neighbor;
//             }
            
//             // X_RN * neighbor_values
//             char trans = 'N';
//             int m = r, n = total_neighbor_points;
//             DataType alpha = 1.0;  // ← ADD
//             DataType beta = 1.0;
//             int lda = r, incx = 1, incy = 1;
            
//             if constexpr (std::is_same_v<DataType, double>) {
//                 dgemv_(&trans, &m, &n, &alpha,
//                        X_RN->data.data(), &lda,
//                        neighbor_values.data(), &incx,
//                        &beta, x_R.data(), &incy);
//             } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
//                 zgemv_(&trans, &m, &n, &alpha,
//                        X_RN->data.data(), &lda,
//                        neighbor_values.data(), &incx,
//                        &beta, x_R.data(), &incy);
//             }
//         }
//     }
    
//     // ===== NO X_RR solve! Already done in diagonal phase =====
    
//     // Store updated x[R]
//     for (int64_t i = 0; i < r; ++i) {
//         solve_data.left_side[(*redundant_indices)[i]] = x_R[i];
//     }
    
   

//     if (T->is_allocated()) {
//         // Re-extract x[R] after Step 1 updates
//         for (int64_t i = 0; i < r; ++i) {
//             x_R[i] = solve_data.left_side[(*redundant_indices)[i]];
//         }
        
//         std::vector<DataType> result(k, DataType{0.0});  // ← size k (skeleton), not r!
        
//         char trans = 'N';  // ← No transpose! T is (k×r)
//         int m = k, n = r;
//         DataType alpha = 1.0, beta = 0.0;
//         int lda = k, incx = 1, incy = 1;
        
//         if constexpr (std::is_same_v<DataType, double>) {
//             dgemv_(&trans, &m, &n, &alpha,
//                 T->data.data(), &lda,
//                 x_R.data(), &incx,  // ← Use x_R, not x_S!
//                 &beta, result.data(), &incy);
//         } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
//             zgemv_(&trans, &m, &n, &alpha,
//                 T->data.data(), &lda,
//                 x_R.data(), &incx,
//                 &beta, result.data(), &incy);
//         }
        
//         // L_T: x[S] -= T * x[R]
//         for (int64_t i = 0; i < k; ++i) {
//             solve_data.left_side[(*skeleton_indices)[i]] -= result[i];  // ← Update skeleton!
//         }
//     }
// }



 // // ===== Step 2: Apply L_T: x[R] += T^T * x[S] =====
    
    // if (T->is_allocated()) {
    //     std::vector<DataType> result(r, DataType{0.0});
        
    //     char trans = 'T';
    //     int m = k, n = r;
    //     DataType alpha = 1.0, beta = 0.0;
    //     int lda = k, incx = 1, incy = 1;
        
    //     if constexpr (std::is_same_v<DataType, double>) {
    //         dgemv_(&trans, &m, &n, &alpha,
    //                T->data.data(), &lda,
    //                x_S.data(), &incx,
    //                &beta, result.data(), &incy);
    //     } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
    //         zgemv_(&trans, &m, &n, &alpha,
    //                T->data.data(), &lda,
    //                x_S.data(), &incx,
    //                &beta, result.data(), &incy);
    //     }
        
    //     for (int64_t i = 0; i < r; ++i) {
    //         solve_data.left_side[(*redundant_indices)[i]] -= result[i];
    //     }
    // }
    
    // ===== Step 2: Apply L_T: x[S] -= T * x[R] =====
    // L_T has -T in (S,R) block
    // This updates skeleton DOFs using redundant DOFs


// template<typename CoordType, typename DataType>
// void apply_backward_substitution(
//     TreeLevel<CoordType, DataType>& level,
//     SolveDataRequest<CoordType, DataType>& solve_data,
//     std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
//     MatrixProperty matrix_property,
//     bool is_ghost) {
    
//     // Get factorization matrices
//     const MatrixStorage<DataType>* T = nullptr;
//     const MatrixStorage<DataType>* X_RR = nullptr;
//     const MatrixStorage<DataType>* X_RS = nullptr;  // Only for NONSYMMETRIC
//     const MatrixStorage<DataType>* X_RN = nullptr;  // Only for NONSYMMETRIC
//     const MatrixStorage<DataType>* X_SR = nullptr;  // For SYMMETRIC (transpose)
//     const MatrixStorage<DataType>* X_NR = nullptr;  // For SYMMETRIC (transpose)
//     const std::vector<int64_t>* skeleton_indices = nullptr;
//     const std::vector<int64_t>* redundant_indices = nullptr;
//     BoxData<CoordType, DataType>* box = nullptr;
    
//     if (is_ghost) {
//         T = &solve_data.interpolation_matrix;
//         X_RR = &solve_data.X_RR;
//         X_RS = &solve_data.X_RS;
//         X_RN = &solve_data.X_RN;
//         X_SR = &solve_data.X_SR;
//         X_NR = &solve_data.X_NR;
//         skeleton_indices = &solve_data.skeleton_indices;
//         redundant_indices = &solve_data.redundant_indices;
//     } else {
//         box = level.find_local_box(solve_data.morton_index);
//         if (box == nullptr) {
//             throw std::runtime_error(
//                 "apply_backward_substitution: Local box " + 
//                 std::to_string(solve_data.morton_index) + " not found");
//         }
        
//         T = &box->interpolation_matrix;
//         X_RR = &box->X_RR;
//         X_RS = &box->X_RS;
//         X_RN = &box->X_RN;
//         X_SR = &box->X_SR;
//         X_NR = &box->X_NR;
//         skeleton_indices = &box->skeleton_indices;
//         redundant_indices = &box->redundant_indices;
//     }
    
//     if (skeleton_indices->empty() || redundant_indices->empty()) {
//         return;
//     }
    
//     int64_t k = skeleton_indices->size();
//     int64_t r = redundant_indices->size();
    
//     // ===== Step 1: Apply U^{-1} =====
    
//     // Start with temp = x[R]
//     std::vector<DataType> temp(r);
//     for (int64_t i = 0; i < r; ++i) {
//         temp[i] = solve_data.left_side[(*redundant_indices)[i]];
//     }
    
//     // ===== Add X_RS * x[S] (or X_SR^T * x[S] for symmetric) =====
    
//     std::vector<DataType> x_S(k);
//     for (int64_t i = 0; i < k; ++i) {
//         x_S[i] = solve_data.left_side[(*skeleton_indices)[i]];
//     }
    
//     if (matrix_property == MatrixProperty::SYMMETRIC) {
//         // Use X_SR^T instead of X_RS
//         if (X_SR->is_allocated()) {
//             // X_SR is (k × r), X_SR^T is (r × k)
//             // temp (r) += X_SR^T (r × k) * x_S (k)
            
//             char trans = 'T';  // Transpose
//             int m = k, n = r;  // X_SR dimensions
//             DataType alpha = 1.0, beta = 1.0;  // beta=1 to accumulate
//             int lda = k, incx = 1, incy = 1;
            
//             if constexpr (std::is_same_v<DataType, double>) {
//                 dgemv_(&trans, &m, &n, &alpha,
//                        X_SR->data.data(), &lda,
//                        x_S.data(), &incx,
//                        &beta, temp.data(), &incy);
//             } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
//                 zgemv_(&trans, &m, &n, &alpha,
//                        X_SR->data.data(), &lda,
//                        x_S.data(), &incx,
//                        &beta, temp.data(), &incy);
//             }
//         }
//     } else {
//         // NONSYMMETRIC: Use X_RS directly
//         if (X_RS->is_allocated()) {
//             // X_RS is (r × k)
//             // temp (r) += X_RS (r × k) * x_S (k)
            
//             char trans = 'N';
//             int m = r, n = k;
//             DataType alpha = 1.0, beta = 1.0;  // beta=1 to accumulate
//             int lda = r, incx = 1, incy = 1;
            
//             if constexpr (std::is_same_v<DataType, double>) {
//                 dgemv_(&trans, &m, &n, &alpha,
//                        X_RS->data.data(), &lda,
//                        x_S.data(), &incx,
//                        &beta, temp.data(), &incy);
//             } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
//                 zgemv_(&trans, &m, &n, &alpha,
//                        X_RS->data.data(), &lda,
//                        x_S.data(), &incx,
//                        &beta, temp.data(), &incy);
//             }
//         }
//     }
    
//     // ===== Add X_RN * x[N] (or X_NR^T * x[N] for symmetric) - batched =====
    
//     if (box != nullptr && !box->near_field_modified_interactions.empty()) {
        
//         if (matrix_property == MatrixProperty::SYMMETRIC && X_NR->is_allocated()) {
//             // Use X_NR^T instead of X_RN
//             // X_NR is (total_neighbor_points × r)
//             // X_NR^T is (r × total_neighbor_points)
            
//             int64_t total_neighbor_points = X_NR->rows;
            
//             // Collect all neighbor skeleton values
//             std::vector<DataType> neighbor_values(total_neighbor_points, DataType{0.0});
            
//             int64_t row_offset = 0;
//             for (const auto& modified_block : box->near_field_modified_interactions) {
//                 int64_t neighbor_morton = modified_block.neighbor_morton;
//                 int64_t n_neighbor = modified_block.A_NS.rows;
                
//                 // Find neighbor's solve data
//                 SolveDataRequest<CoordType, DataType>* neighbor_data = nullptr;
                
//                 for (auto& solve_box : level_solve_data) {
//                     if (solve_box.morton_index == neighbor_morton) {
//                         neighbor_data = &solve_box;
//                         break;
//                     }
//                 }
                
//                 if (neighbor_data == nullptr) {
//                     auto it = level.ghost_and_assisting_box_points_for_solve_map.find(neighbor_morton);
//                     if (it != level.ghost_and_assisting_box_points_for_solve_map.end()) {
//                         neighbor_data = &level.ghost_and_assisting_boxes_for_solve[it->second];
//                     }
//                 }
                
//                 if (neighbor_data == nullptr) {
//                     row_offset += n_neighbor;
//                     continue;
//                 }
                
//                 // Copy neighbor's skeleton values
//                 const auto& neighbor_skel = neighbor_data->skeleton_indices;
//                 for (int64_t i = 0; i < n_neighbor && i < neighbor_skel.size(); ++i) {
//                     neighbor_values[row_offset + i] = neighbor_data->left_side[neighbor_skel[i]];
//                 }
                
//                 row_offset += n_neighbor;
//             }
            
//             // Compute X_NR^T * neighbor_values: (r × n) * (n × 1) = (r × 1)
//             char trans = 'T';
//             int m = total_neighbor_points, n = r;  // X_NR dimensions
//             DataType alpha = 1.0, beta = 1.0;  // beta=1 to accumulate
//             int lda = total_neighbor_points, incx = 1, incy = 1;
            
//             if constexpr (std::is_same_v<DataType, double>) {
//                 dgemv_(&trans, &m, &n, &alpha,
//                        X_NR->data.data(), &lda,
//                        neighbor_values.data(), &incx,
//                        &beta, temp.data(), &incy);
//             } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
//                 zgemv_(&trans, &m, &n, &alpha,
//                        X_NR->data.data(), &lda,
//                        neighbor_values.data(), &incx,
//                        &beta, temp.data(), &incy);
//             }
            
//         } else if (matrix_property == MatrixProperty::NONSYMMETRIC && X_RN->is_allocated()) {
//             // Use X_RN directly
//             int64_t total_neighbor_points = X_RN->cols;
            
//             // Collect all neighbor skeleton values
//             std::vector<DataType> neighbor_values(total_neighbor_points, DataType{0.0});
            
//             int64_t col_offset = 0;
//             for (const auto& modified_block : box->near_field_modified_interactions) {
//                 int64_t neighbor_morton = modified_block.neighbor_morton;
//                 int64_t n_neighbor = modified_block.A_NS.rows;
                
//                 // Find neighbor's solve data
//                 SolveDataRequest<CoordType, DataType>* neighbor_data = nullptr;
                
//                 for (auto& solve_box : level_solve_data) {
//                     if (solve_box.morton_index == neighbor_morton) {
//                         neighbor_data = &solve_box;
//                         break;
//                     }
//                 }
                
//                 if (neighbor_data == nullptr) {
//                     auto it = level.ghost_and_assisting_box_points_for_solve_map.find(neighbor_morton);
//                     if (it != level.ghost_and_assisting_box_points_for_solve_map.end()) {
//                         neighbor_data = &level.ghost_and_assisting_boxes_for_solve[it->second];
//                     }
//                 }
                
//                 if (neighbor_data == nullptr) {
//                     col_offset += n_neighbor;
//                     continue;
//                 }
                
//                 // Copy neighbor's skeleton values
//                 const auto& neighbor_skel = neighbor_data->skeleton_indices;
//                 for (int64_t i = 0; i < n_neighbor && i < neighbor_skel.size(); ++i) {
//                     neighbor_values[col_offset + i] = neighbor_data->left_side[neighbor_skel[i]];
//                 }
                
//                 col_offset += n_neighbor;
//             }
            
//             // Compute X_RN * neighbor_values: (r × n) * (n × 1) = (r × 1)
//             char trans = 'N';
//             int m = r, n = total_neighbor_points;
//             DataType alpha = 1.0, beta = 1.0;  // beta=1 to accumulate
//             int lda = r, incx = 1, incy = 1;
            
//             if constexpr (std::is_same_v<DataType, double>) {
//                 dgemv_(&trans, &m, &n, &alpha,
//                        X_RN->data.data(), &lda,
//                        neighbor_values.data(), &incx,
//                        &beta, temp.data(), &incy);
//             } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
//                 zgemv_(&trans, &m, &n, &alpha,
//                        X_RN->data.data(), &lda,
//                        neighbor_values.data(), &incx,
//                        &beta, temp.data(), &incy);
//             }
//         }
//     }
    
//     // Solve X_RR * x[R] = temp using Cholesky factorization
//     std::vector<DataType> x_R = temp;
    
//     if (X_RR->format == MatrixStorage<DataType>::CHOLESKY_L) {
//         char uplo = 'L';
//         int n = r, nrhs = 1;
//         int lda = r, ldb = r, info = 0;
        
//         if constexpr (std::is_same_v<DataType, double>) {
//             dpotrs_(&uplo, &n, &nrhs,
//                     X_RR->data.data(), &lda,
//                     x_R.data(), &ldb, &info);
//         } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
//             zpotrs_(&uplo, &n, &nrhs,
//                     X_RR->data.data(), &lda,
//                     x_R.data(), &ldb, &info);
//         }
        
//         if (info != 0) {
//             throw std::runtime_error("Cholesky solve failed in backward substitution");
//         }
//     }
    
//     // Store result in x[R]
//     for (int64_t i = 0; i < r; ++i) {
//         solve_data.left_side[(*redundant_indices)[i]] = x_R[i];
//     }
    
//     // ===== Step 2: Apply L_T^{-1}: x[R] += T * x[S] =====
    
//     if (T->is_allocated()) {
//         std::vector<DataType> x_S(k);
//         for (int64_t i = 0; i < k; ++i) {
//             x_S[i] = solve_data.left_side[(*skeleton_indices)[i]];
//         }
        
//         // Compute T^T * x_S: (r × k) * (k × 1) = (r × 1)
//         std::vector<DataType> result(r, DataType{0.0});
        
//         char trans = 'T';  // T is (k × r), we need T^T
//         int m = k, n = r;
//         DataType alpha = 1.0, beta = 0.0;
//         int lda = k, incx = 1, incy = 1;
        
//         if constexpr (std::is_same_v<DataType, double>) {
//             dgemv_(&trans, &m, &n, &alpha,
//                    T->data.data(), &lda,
//                    x_S.data(), &incx,
//                    &beta, result.data(), &incy);
//         } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
//             zgemv_(&trans, &m, &n, &alpha,
//                    T->data.data(), &lda,
//                    x_S.data(), &incx,
//                    &beta, result.data(), &incy);
//         }
        
//         // Add result to x[R]
//         for (int64_t i = 0; i < r; ++i) {
//             solve_data.left_side[(*redundant_indices)[i]] += result[i];
//         }
//     }
// }

/**
 * @brief Solve diagonal block using X_RR factorization
 * 
 * Solves X_RR * x = b where X_RR is in Cholesky form (L * L^T).
 * Uses dpotrs which performs two triangular solves: L^{-1} then L^{-T}.
 * 
 * At level 0 (root): All DOFs are "skeleton" but we solve with X_RR
 * At other levels: Would be used for redundant DOFs (if applicable)
 * 
 * @param level Tree level
 * @param solve_data Solve data for the box
 * @param is_ghost True if box is a ghost box
 */
template<typename CoordType, typename DataType>
void apply_diagonal_solve(
    TreeLevel<CoordType, DataType>& level,
    SolveDataRequest<CoordType, DataType>& solve_data,
    bool is_ghost) {
    
    const MatrixStorage<DataType>* X_RR = nullptr;
    const std::vector<int>* X_RR_pivots = nullptr;
    const std::vector<int64_t>* skeleton_indices = nullptr;
    
    if (is_ghost) {
        X_RR = &solve_data.X_RR;
        X_RR_pivots = &solve_data.X_RR_pivots;
        skeleton_indices = &solve_data.skeleton_indices;
    } else {
        BoxData<CoordType, DataType>* box = level.find_local_box(solve_data.morton_index);
        if (box == nullptr) {
            throw std::runtime_error(
                "apply_diagonal_solve: Local box " + 
                std::to_string(solve_data.morton_index) + " not found");
        }
        
        X_RR = &box->X_RR;
        X_RR_pivots = &box->X_RR_pivots;
        skeleton_indices = &box->skeleton_indices;
    }
    
    if (!X_RR->is_allocated()) {
        throw std::runtime_error("X_RR not allocated for diagonal solve");
    }
    
    int64_t k = skeleton_indices->size();
    
    
    // Extract x[S] from left_side
    // (At level 0, "skeleton" indices represent all DOFs)
    std::vector<DataType> x_S(k);
    for (int64_t i = 0; i < k; ++i) {
        x_S[i] = solve_data.left_side[(*skeleton_indices)[i]];
    }
    
    int n = static_cast<int>(k);
    int nrhs = 1;
    int lda = static_cast<int>(k);
    int ldb = static_cast<int>(k);
    int info = 0;

    if (X_RR->format == MatrixStorage<DataType>::CHOLESKY_L) {
        char uplo = 'L';

        if constexpr (std::is_same_v<DataType, double>) {
            dpotrs_(&uplo, &n, &nrhs,
                    X_RR->data.data(), &lda,
                    x_S.data(), &ldb, &info);
        } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
            zsychol_solve_(&uplo, &n, &nrhs,
                        X_RR->data.data(), &lda,
                        x_S.data(), &ldb, &info);
        }

        if (info != 0) {
            throw std::runtime_error(
                "Diagonal solve (dpotrs) failed with info = " + std::to_string(info));
        }
    } else if (X_RR->format == MatrixStorage<DataType>::LU_FACTORED) {
        if (X_RR_pivots == nullptr || X_RR_pivots->size() < static_cast<size_t>(k)) {
            throw std::runtime_error("X_RR LU solve is missing pivot data");
        }

        char trans = 'N';
        if constexpr (std::is_same_v<DataType, double>) {
            dgetrs_(&trans, &n, &nrhs,
                    const_cast<double*>(X_RR->data.data()), &lda,
                    const_cast<int*>(X_RR_pivots->data()),
                    x_S.data(), &ldb, &info);
        } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
            zgetrs_(&trans, &n, &nrhs,
                    const_cast<std::complex<double>*>(X_RR->data.data()), &lda,
                    const_cast<int*>(X_RR_pivots->data()),
                    x_S.data(), &ldb, &info);
        }

        if (info != 0) {
            throw std::runtime_error(
                "Diagonal solve (dgetrs) failed with info = " + std::to_string(info));
        }
    } else if (X_RR->format == MatrixStorage<DataType>::BUNCH_KAUFMAN) {
        if (X_RR_pivots == nullptr || X_RR_pivots->size() < static_cast<size_t>(k)) {
            throw std::runtime_error("X_RR Bunch-Kaufman solve is missing pivot data");
        }

        if constexpr (std::is_same_v<DataType, std::complex<double>>) {
            char uplo = 'L';
            zsytrs_(&uplo, &n, &nrhs,
                    X_RR->data.data(), &lda,
                    const_cast<int*>(X_RR_pivots->data()),
                    x_S.data(), &ldb, &info);
        } else {
            throw std::runtime_error("BUNCH_KAUFMAN format only supported for complex<double>");
        }

        if (info != 0) {
            throw std::runtime_error(
                "Diagonal solve (zsytrs) failed with info = " + std::to_string(info));
        }
    } else {
        throw std::runtime_error("Unsupported X_RR format for diagonal solve");
    }
    
    // Store result back to left_side[S]
    for (int64_t i = 0; i < k; ++i) {
        solve_data.left_side[(*skeleton_indices)[i]] = x_S[i];
    }
}

} // namespace fmm
