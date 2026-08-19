#pragma once

#include <vector>
#include <complex>
#include <stdexcept>
#include <cassert>
#include <unordered_map>
#include "tree.hpp"
#include "blas_declare.hpp"

/**
 * @file apply_mul.hpp
 * @brief Hierarchical factorization MULTIPLY operators (F * x).
 *
 * The factorization gives:
 *   F = (prod_i V_i) * P_l * D * P_l^* * (prod_i W_i)   ~  K
 *
 * So F * x is applied as:
 *   Phase 1 (forward, leaf->root): apply W_n, W_{n-1}, ..., W_1
 *   Phase 2 (diagonal):            apply D  (multiply by X_RR blocks)
 *   Phase 3 (backward, root->leaf): apply V_1, V_2, ..., V_n
 *
 * Each operator is the sign-flipped, step-order-reversed version of
 * the corresponding solve operator:
 *   W  = U^{-1} L_T^{-1} P^*    (solve uses W^{-1} = P L_T U)
 *   V  = P U_T^{-1} L^{-1}      (solve uses V^{-1} = L U_T P^*)
 *   D                             (solve uses D^{-1})
 */

namespace fmm {

// Color has no ghost/assisting boxes, so both "any" and "local-or-ghost"
// resolvers reduce to a simple range check on the local slab.
template<typename CoordType, typename DataType>
inline SolveDataRequest<CoordType, DataType>* resolve_any_solve_data_for_morton(
    TreeLevel<CoordType, DataType>& level,
    std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
    int64_t morton) {
    if (morton >= level.local_morton_start && morton <= level.local_morton_end) {
        return &level_solve_data[static_cast<size_t>(morton - level.local_morton_start)];
    }
    return nullptr;
}

template<typename CoordType, typename DataType>
inline SolveDataRequest<CoordType, DataType>* resolve_local_or_ghost_solve_data_for_morton(
    TreeLevel<CoordType, DataType>& level,
    std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
    int64_t morton) {
    return resolve_any_solve_data_for_morton(level, level_solve_data, morton);
}

// ============================================================================
// Phase 1: Forward sweep — apply W operator (leaf -> root)
// ============================================================================
//
// W = U^{-1} L_T^{-1} P^*
//   Step 1 (L_T^{-1}):  x[S] += T * x[R]             (self-box only)
//   Step 2 (U^{-1}):    x[R] -= stored_X_SR^T * x[S]  (self-box, uses updated x[S])
//   Step 3 (U^{-1}):    x[R] -= stored_X_NR^T * x[N]  (reads from neighbors)
//
// Reads from neighbors but does NOT write to them -> no pending updates needed.
// Same communication pattern as the solve's backward sweep (W^{-1}).

template<typename CoordType, typename DataType>
void apply_mul_forward_W(
    TreeLevel<CoordType, DataType>& level,
    SolveDataRequest<CoordType, DataType>& solve_data,
    std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
    MatrixProperty matrix_property,
    bool is_ghost) {

    const MatrixStorage<DataType>* T = nullptr;
    const MatrixStorage<DataType>* X_SR = nullptr;
    const MatrixStorage<DataType>* X_NR = nullptr;
    const MatrixStorage<DataType>* X_RS = nullptr;
    const MatrixStorage<DataType>* X_RN = nullptr;
    const std::vector<int64_t>* skeleton_indices = nullptr;
    const std::vector<int64_t>* redundant_indices = nullptr;
    const std::vector<int64_t>* one_hop = nullptr;
    const std::vector<int64_t>* use_full_set = nullptr;
    BoxData<CoordType, DataType>* box = nullptr;

    if (is_ghost) {
        box = level.find_ghost_box(solve_data.morton_index);
        if (box == nullptr) {
            T = &solve_data.interpolation_matrix;
            X_SR = &solve_data.X_SR;
            X_NR = &solve_data.X_NR;
            X_RS = &solve_data.X_RS;
            X_RN = &solve_data.X_RN;
            skeleton_indices = &solve_data.skeleton_indices;
            redundant_indices = &solve_data.redundant_indices;
            one_hop = &solve_data.one_hop;
            use_full_set = &solve_data.use_full_set;
        } else {
            T = &box->interpolation_matrix;
            X_SR = &box->X_SR;
            X_NR = &box->X_NR;
            X_RS = &box->X_RS;
            X_RN = &box->X_RN;
            skeleton_indices = &box->skeleton_indices;
            redundant_indices = &box->redundant_indices;
            one_hop = &box->one_hop;
            use_full_set = &box->use_full_set;
        }
    } else {
        box = level.find_local_box(solve_data.morton_index);
        if (box == nullptr) {
            throw std::runtime_error(
                "apply_mul_forward_W: Local box " +
                std::to_string(solve_data.morton_index) + " not found");
        }
        T = &box->interpolation_matrix;
        X_SR = &box->X_SR;
        X_NR = &box->X_NR;
        X_RS = &box->X_RS;
        X_RN = &box->X_RN;
        skeleton_indices = &box->skeleton_indices;
        redundant_indices = &box->redundant_indices;
        one_hop = &box->one_hop;
        use_full_set = &box->use_full_set;
    }

    if (skeleton_indices->empty() || redundant_indices->empty()) {
        return;
    }

    int64_t k = static_cast<int64_t>(skeleton_indices->size());
    int64_t r = static_cast<int64_t>(redundant_indices->size());
    const int64_t nrhs = solve_data.nrhs;

    // Extract x[R]
    std::vector<DataType> x_R(static_cast<size_t>(r * nrhs));
    for (int64_t column = 0; column < nrhs; ++column) {
        for (int64_t i = 0; i < r; ++i) {
            x_R[static_cast<size_t>(i + column * r)] =
                solve_data.left_side[static_cast<size_t>(
                    (*redundant_indices)[static_cast<size_t>(i)] +
                    column * solve_data.num_points)];
        }
    }

    // ===== Step 1: L_T^{-1}: x[S] += T * x[R] =====
    // (Compare solve backward: L_T: x[S] -= T * x[R])
    assert(T->is_allocated());
    if (T->is_allocated()) {
        std::vector<DataType> result(
            static_cast<size_t>(k * nrhs), DataType{0.0});

        char trans = 'N';
        int m = static_cast<int>(k), n = static_cast<int>(nrhs);
        int inner = static_cast<int>(r);
        DataType alpha = 1.0, beta = 0.0;
        int lda = static_cast<int>(T->lda), ldb = static_cast<int>(r);
        int ldc = static_cast<int>(k);
        gemm_(&trans, &trans, &m, &n, &inner, &alpha,
              T->data.data(), &lda, x_R.data(), &ldb,
              &beta, result.data(), &ldc);

        for (int64_t column = 0; column < nrhs; ++column) {
            for (int64_t i = 0; i < k; ++i) {
                solve_data.left_side[static_cast<size_t>(
                    (*skeleton_indices)[static_cast<size_t>(i)] +
                    column * solve_data.num_points)] +=
                    result[static_cast<size_t>(i + column * k)];
            }
        }
    }

    // Re-extract x[S] after update (needed for Step 2)
    std::vector<DataType> x_S(static_cast<size_t>(k * nrhs));
    for (int64_t column = 0; column < nrhs; ++column) {
        for (int64_t i = 0; i < k; ++i) {
            x_S[static_cast<size_t>(i + column * k)] =
                solve_data.left_side[static_cast<size_t>(
                    (*skeleton_indices)[static_cast<size_t>(i)] +
                    column * solve_data.num_points)];
        }
    }

    // ===== Step 2: U^{-1}: x[R] -= stored_X_SR^T * x[S] =====
    // (Compare solve backward: U: x[R] += stored_X_SR^T * x[S])
    // stored_X_SR = -X_SR*X_RR^{-1}
    // U^{-1} wants +X_RR^{-1}*X_RS = -(stored_X_SR)^T for symmetric
    // So x[R] += -(stored_X_SR)^T * x[S] = x[R] -= stored_X_SR^T * x[S]
    if (matrix_property == MatrixProperty::SYMMETRIC) {
        assert(X_SR->is_allocated());
        if (X_SR->is_allocated()) {
            char trans = 'T';
            char trans_b = 'N';
            int m = static_cast<int>(r), n = static_cast<int>(nrhs);
            int inner = static_cast<int>(k);
            DataType alpha = -1.0;
            DataType beta = 1.0;
            int lda = static_cast<int>(X_SR->lda), ldb = static_cast<int>(k);
            int ldc = static_cast<int>(r);
            gemm_(&trans, &trans_b, &m, &n, &inner, &alpha,
                  X_SR->data.data(), &lda, x_S.data(), &ldb,
                  &beta, x_R.data(), &ldc);
        }
    } else {
        // NONSYMMETRIC: stored_X_RS = -X_RR^{-1}*X_RS
        // U^{-1} wants +X_RR^{-1}*X_RS = -stored_X_RS
        if (X_RS->is_allocated()) {
            char trans = 'N';
            int m = static_cast<int>(r), n = static_cast<int>(nrhs);
            int inner = static_cast<int>(k);
            DataType alpha = -1.0;
            DataType beta = 1.0;
            int lda = static_cast<int>(X_RS->lda), ldb = static_cast<int>(k);
            int ldc = static_cast<int>(r);
            gemm_(&trans, &trans, &m, &n, &inner, &alpha,
                  X_RS->data.data(), &lda, x_S.data(), &ldb,
                  &beta, x_R.data(), &ldc);
        }
    }

    // ===== Step 3: U^{-1}: x[R] -= stored_X_NR^T * x[N] =====
    // (Compare solve backward: U: x[R] += stored_X_NR^T * x[N])
    if (one_hop != nullptr && !one_hop->empty()) {

        if (matrix_property == MatrixProperty::SYMMETRIC && X_NR->is_allocated()) {
            int64_t total_neighbor_points = X_NR->rows;
            std::vector<DataType> neighbor_values(
                static_cast<size_t>(total_neighbor_points * nrhs), DataType{0.0});

            // Collect neighbor values in one_hop order
            int64_t row_offset = 0;
            for (size_t neighbor_idx = 0; neighbor_idx < one_hop->size(); ++neighbor_idx) {
                int64_t neighbor_morton = (*one_hop)[neighbor_idx];
                bool use_full = ((*use_full_set)[neighbor_idx] == 1);

                SolveDataRequest<CoordType, DataType>* neighbor_data = nullptr;

                if (neighbor_morton >= level.local_morton_start && neighbor_morton <= level.local_morton_end) {
                    neighbor_data = &level_solve_data[static_cast<size_t>(neighbor_morton - level.local_morton_start)];
                }

                if (neighbor_data == nullptr) {
                    auto it = level.ghost_and_assisting_box_points_for_solve_map.find(neighbor_morton);
                    if (it != level.ghost_and_assisting_box_points_for_solve_map.end()) {
                        neighbor_data = &level.ghost_and_assisting_boxes_for_solve[it->second];
                    }
                }
                if (neighbor_data == nullptr) {
                    throw std::runtime_error(
                        "apply_mul_forward_W: Neighbor " +
                        std::to_string(neighbor_morton) + " not found in solve data");
                }
                if (neighbor_data->nrhs != nrhs) {
                    throw std::runtime_error(
                        "apply_mul_forward_W: neighbor RHS count mismatch");
                }

                int64_t n_neighbor;

                if (use_full) {
                    n_neighbor = neighbor_data->num_points;
                    for (int64_t column = 0; column < nrhs; ++column) {
                        for (int64_t i = 0; i < n_neighbor; ++i) {
                            neighbor_values[static_cast<size_t>(
                                row_offset + i + column * total_neighbor_points)] =
                                neighbor_data->left_side[static_cast<size_t>(
                                    i + column * neighbor_data->num_points)];
                        }
                    }
                } else {
                    n_neighbor = static_cast<int64_t>(neighbor_data->skeleton_indices.size());
                    const auto& neighbor_skel = neighbor_data->skeleton_indices;
                    for (int64_t column = 0; column < nrhs; ++column) {
                        for (int64_t i = 0; i < n_neighbor; ++i) {
                            neighbor_values[static_cast<size_t>(
                                row_offset + i + column * total_neighbor_points)] =
                                neighbor_data->left_side[static_cast<size_t>(
                                    neighbor_skel[static_cast<size_t>(i)] +
                                    column * neighbor_data->num_points)];
                        }
                    }
                }

                row_offset += n_neighbor;
            }
            if (row_offset != total_neighbor_points) {
                throw std::runtime_error(
                    "apply_mul_forward_W: row_offset " + std::to_string(row_offset) +
                    " does not match total_neighbor_points " + std::to_string(total_neighbor_points));
            }

            // x[R] -= stored_X_NR^T * x[N]  (sign flip from solve's +=)
            char trans_a = 'T', trans_b = 'N';
            int m = static_cast<int>(r), n = static_cast<int>(nrhs);
            int inner = static_cast<int>(total_neighbor_points);
            DataType alpha = -1.0;
            DataType beta = 1.0;
            int lda = static_cast<int>(X_NR->lda);
            int ldb = static_cast<int>(total_neighbor_points);
            int ldc = static_cast<int>(r);
            gemm_(&trans_a, &trans_b, &m, &n, &inner, &alpha,
                  X_NR->data.data(), &lda, neighbor_values.data(), &ldb,
                  &beta, x_R.data(), &ldc);

        } else if (matrix_property == MatrixProperty::NONSYMMETRIC && X_RN->is_allocated()) {
            int64_t total_neighbor_points = X_RN->cols;
            std::vector<DataType> neighbor_values(
                static_cast<size_t>(total_neighbor_points * nrhs), DataType{0.0});

            int64_t col_offset = 0;
            for (size_t neighbor_idx = 0; neighbor_idx < one_hop->size(); ++neighbor_idx) {
                int64_t neighbor_morton = (*one_hop)[neighbor_idx];
                bool use_full = ((*use_full_set)[neighbor_idx] == 1);

                SolveDataRequest<CoordType, DataType>* neighbor_data = nullptr;
                if (neighbor_morton >= level.local_morton_start && neighbor_morton <= level.local_morton_end) {
                    neighbor_data = &level_solve_data[static_cast<size_t>(neighbor_morton - level.local_morton_start)];
                }
                if (neighbor_data == nullptr) {
                    auto it = level.ghost_and_assisting_box_points_for_solve_map.find(neighbor_morton);
                    if (it != level.ghost_and_assisting_box_points_for_solve_map.end()) {
                        neighbor_data = &level.ghost_and_assisting_boxes_for_solve[it->second];
                    }
                }
                if (neighbor_data == nullptr) {
                    throw std::runtime_error(
                        "apply_mul_forward_W: Neighbor " +
                        std::to_string(neighbor_morton) + " not found");
                }
                if (neighbor_data->nrhs != nrhs) {
                    throw std::runtime_error(
                        "apply_mul_forward_W: neighbor RHS count mismatch");
                }

                int64_t n_neighbor;
                if (use_full) {
                    n_neighbor = neighbor_data->num_points;
                    for (int64_t column = 0; column < nrhs; ++column) {
                        for (int64_t i = 0; i < n_neighbor; ++i) {
                            neighbor_values[static_cast<size_t>(
                                col_offset + i + column * total_neighbor_points)] =
                                neighbor_data->left_side[static_cast<size_t>(
                                    i + column * neighbor_data->num_points)];
                        }
                    }
                } else {
                    n_neighbor = static_cast<int64_t>(neighbor_data->skeleton_indices.size());
                    const auto& neighbor_skel = neighbor_data->skeleton_indices;
                    for (int64_t column = 0; column < nrhs; ++column) {
                        for (int64_t i = 0; i < n_neighbor; ++i) {
                            neighbor_values[static_cast<size_t>(
                                col_offset + i + column * total_neighbor_points)] =
                                neighbor_data->left_side[static_cast<size_t>(
                                    neighbor_skel[static_cast<size_t>(i)] +
                                    column * neighbor_data->num_points)];
                        }
                    }
                }
                col_offset += n_neighbor;
            }
            if (col_offset != total_neighbor_points) {
                throw std::runtime_error(
                    "apply_mul_forward_W: col_offset mismatch");
            }

            // x[R] -= stored_X_RN * x[N]  (sign flip from solve's +=)
            char trans = 'N';
            int m = static_cast<int>(r), n = static_cast<int>(nrhs);
            int inner = static_cast<int>(total_neighbor_points);
            DataType alpha = -1.0;
            DataType beta = 1.0;
            int lda = static_cast<int>(X_RN->lda);
            int ldb = static_cast<int>(total_neighbor_points);
            int ldc = static_cast<int>(r);
            gemm_(&trans, &trans, &m, &n, &inner, &alpha,
                  X_RN->data.data(), &lda, neighbor_values.data(), &ldb,
                  &beta, x_R.data(), &ldc);
        }
    }

    // Store updated x[R]
    for (int64_t column = 0; column < nrhs; ++column) {
        for (int64_t i = 0; i < r; ++i) {
            solve_data.left_side[static_cast<size_t>(
                (*redundant_indices)[static_cast<size_t>(i)] +
                column * solve_data.num_points)] =
                x_R[static_cast<size_t>(i + column * r)];
        }
    }
}


// ============================================================================
// Bunch-Kaufman multiply: y = A*x where A = L*D*L^T (with interleaved pivots)
// ============================================================================
//
// xSYTRF with UPLO='L' produces A = L*D*L^T where:
//   L = product of permutation and unit lower triangular elementary matrices
//   D = block diagonal with 1x1 and 2x2 blocks
//   IPIV encodes the pivot structure (positive=1x1, negative=2x2)
//
// To multiply y = A*x, we reverse the xSYTRS solve steps:
//   Step 1: Apply L^T with forward permutations (k = 0..n-1)
//   Step 2: Apply D (block diagonal multiply)
//   Step 3: Apply L with backward permutations (k = n-1..0)

template <typename DataType>
inline void bunch_kaufman_multiply(
    int n,
    const DataType* A, int lda,
    const int* ipiv,
    DataType* x, int ldx, int nrhs)
{
    if (n == 0 || nrhs == 0) {
        return;
    }

    const DataType one = DataType{1};
    const char trans = 'T';
    const char no_trans = 'N';

    // Step 1: L^T multiply with permutations (forward, k=0..n-1)
    // This is the inverse of the SYTRS backward pass
    int k = 0;
    while (k < n) {
        if (ipiv[k] > 0) {
            // 1x1 pivot
            int kp = ipiv[k] - 1;  // convert to 0-based
            if (kp != k) {
                for (int column = 0; column < nrhs; ++column) {
                    std::swap(x[k + column * ldx], x[kp + column * ldx]);
                }
            }
            // X[k,:] += L(k+1:n,k)^T * X[k+1:n,:]
            int tail = n - k - 1;
            if (tail > 0) {
                int one_row = 1;
                gemm_(&trans, &no_trans, &one_row, &nrhs, &tail,
                      &one, A + (k + 1) + k * lda, &lda,
                      x + k + 1, &ldx, &one, x + k, &ldx);
            }
            k += 1;
        } else {
            // 2x2 pivot: ipiv[k] < 0 and ipiv[k+1] < 0
            int kp = -ipiv[k] - 1;  // convert to 0-based
            if (kp != k + 1) {
                for (int column = 0; column < nrhs; ++column) {
                    std::swap(x[k + 1 + column * ldx],
                              x[kp + column * ldx]);
                }
            }
            // X[k:k+1,:] += L(k+2:n,k:k+1)^T * X[k+2:n,:]
            int tail = n - k - 2;
            if (tail > 0) {
                int two_rows = 2;
                gemm_(&trans, &no_trans, &two_rows, &nrhs, &tail,
                      &one, A + (k + 2) + k * lda, &lda,
                      x + k + 2, &ldx, &one, x + k, &ldx);
            }
            k += 2;
        }
    }

    // Step 2: D multiply (block diagonal)
    k = 0;
    while (k < n) {
        if (ipiv[k] > 0) {
            // 1x1 block: D(k,k) = A(k,k)
            for (int column = 0; column < nrhs; ++column) {
                x[k + column * ldx] *= A[k + k * lda];
            }
            k += 1;
        } else {
            // 2x2 block: D = [A(k,k), A(k+1,k); A(k+1,k), A(k+1,k+1)]
            DataType d11 = A[k + k * lda];
            DataType d21 = A[(k + 1) + k * lda];
            DataType d22 = A[(k + 1) + (k + 1) * lda];
            for (int column = 0; column < nrhs; ++column) {
                DataType* x_column = x + column * ldx;
                DataType t0 = x_column[k], t1 = x_column[k + 1];
                x_column[k] = d11 * t0 + d21 * t1;
                x_column[k + 1] = d21 * t0 + d22 * t1;
            }
            k += 2;
        }
    }

    // Step 3: L multiply with permutations (backward, k=n-1..0)
    // This is the inverse of the SYTRS forward pass
    k = n - 1;
    while (k >= 0) {
        if (ipiv[k] > 0) {
            // 1x1 pivot
            // X[k+1:n,:] += L(k+1:n,k) * X[k,:]
            int tail = n - k - 1;
            if (tail > 0) {
                int one_column = 1;
                gemm_(&no_trans, &no_trans, &tail, &nrhs, &one_column,
                      &one, A + (k + 1) + k * lda, &lda,
                      x + k, &ldx, &one, x + k + 1, &ldx);
            }
            int kp = ipiv[k] - 1;
            if (kp != k) {
                for (int column = 0; column < nrhs; ++column) {
                    std::swap(x[k + column * ldx], x[kp + column * ldx]);
                }
            }
            k -= 1;
        } else {
            // 2x2 pivot: ipiv[k] < 0 and ipiv[k-1] < 0
            // x[k+1:n] += L(k+1:n, k-1) * x[k-1] + L(k+1:n, k) * x[k]
            // Note: the 2x2 block is at (k-1, k), so we process k-1 and k together
            int tail = n - k - 1;
            if (tail > 0) {
                int two_columns = 2;
                gemm_(&no_trans, &no_trans, &tail, &nrhs, &two_columns,
                      &one, A + (k + 1) + (k - 1) * lda, &lda,
                      x + k - 1, &ldx, &one, x + k + 1, &ldx);
            }
            int kp = -ipiv[k - 1] - 1;
            if (kp != k) {
                for (int column = 0; column < nrhs; ++column) {
                    std::swap(x[k + column * ldx], x[kp + column * ldx]);
                }
            }
            k -= 2;
        }
    }
}

// ============================================================================
// Phase 2: Diagonal multiply — multiply by D (X_RR blocks)
// ============================================================================
//
// Instead of solving X_RR^{-1} * b, we compute X_RR * x.
// X_RR is stored in factored form, so we reconstruct the multiply:
//   CHOLESKY_L:  X_RR = L * L^T  ->  y = L^T * x, then z = L * y  (dtrmv)
//   LU_FACTORED: X_RR = P * L * U ->  y = U*x, z = L*y, w = P*z   (dtrmv + laswp)
//   BUNCH_KAUFMAN: A = L*D*L^T with interleaved pivots (xSYTRF)

template<typename CoordType, typename DataType>
void apply_diagonal_multiply(
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
                "apply_diagonal_multiply: Local box " +
                std::to_string(solve_data.morton_index) + " not found");
        }
        X_RR = &box->X_RR;
        X_RR_pivots = &box->X_RR_pivots;
        skeleton_indices = &box->skeleton_indices;
    }

    if (!X_RR->is_allocated()) {
        throw std::runtime_error("apply_diagonal_multiply: X_RR not allocated");
    }

    int64_t k = static_cast<int64_t>(skeleton_indices->size());
    const int64_t rhs_count = solve_data.nrhs;

    // Extract x[S]
    std::vector<DataType> x_S(static_cast<size_t>(k * rhs_count));
    for (int64_t column = 0; column < rhs_count; ++column) {
        for (int64_t i = 0; i < k; ++i) {
            x_S[static_cast<size_t>(i + column * k)] =
                solve_data.left_side[static_cast<size_t>(
                    (*skeleton_indices)[static_cast<size_t>(i)] +
                    column * solve_data.num_points)];
        }
    }

    int n = static_cast<int>(k);
    int nrhs = static_cast<int>(rhs_count);
    int ldb = n;
    char side = 'L';
    DataType one = DataType{1};

    if (X_RR->format == MatrixStorage<DataType>::CHOLESKY_L) {
        // X_RR = L * L^T (real symmetric) or L * L^T (complex symmetric via zsychol_)
        // Multiply: first y = L^T * x, then z = L * y
        char uplo = 'L';
        char diag = 'N';
        int lda = static_cast<int>(k);

        char trans_T = 'T';
        trmm_(&side, &uplo, &trans_T, &diag, &n, &nrhs, &one,
              X_RR->data.data(), &lda, x_S.data(), &ldb);
        char trans_N = 'N';
        trmm_(&side, &uplo, &trans_N, &diag, &n, &nrhs, &one,
              X_RR->data.data(), &lda, x_S.data(), &ldb);

    } else if (X_RR->format == MatrixStorage<DataType>::LU_FACTORED) {
        // X_RR = P * L * U
        // Multiply: y = U*x, z = L*y, w = P*z
        if (X_RR_pivots == nullptr || X_RR_pivots->size() < static_cast<size_t>(k)) {
            throw std::runtime_error("apply_diagonal_multiply: missing LU pivots");
        }

        int lda = static_cast<int>(k);

        char uplo_U = 'U', trans_N = 'N', diag_N = 'N';
        trmm_(&side, &uplo_U, &trans_N, &diag_N, &n, &nrhs, &one,
              X_RR->data.data(), &lda, x_S.data(), &ldb);
        char uplo_L = 'L', diag_U = 'U';
        trmm_(&side, &uplo_L, &trans_N, &diag_U, &n, &nrhs, &one,
              X_RR->data.data(), &lda, x_S.data(), &ldb);
        int k1 = 1, k2 = n, inc_rev = -1;
        if constexpr (std::is_same_v<DataType, double>) {
            dlaswp_(&nrhs, x_S.data(), &n,
                    &k1, &k2, X_RR_pivots->data(), &inc_rev);
        } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
            zlaswp_(&nrhs, x_S.data(), &n,
                    &k1, &k2, X_RR_pivots->data(), &inc_rev);
        }

    } else if (X_RR->format == MatrixStorage<DataType>::BUNCH_KAUFMAN) {
        if (X_RR_pivots == nullptr || X_RR_pivots->size() < static_cast<size_t>(k)) {
            throw std::runtime_error("apply_diagonal_multiply: missing Bunch-Kaufman pivots");
        }
        bunch_kaufman_multiply(
            n, X_RR->data.data(), n, X_RR_pivots->data(),
            x_S.data(), ldb, nrhs);
    } else {
        throw std::runtime_error(
            "apply_diagonal_multiply: unsupported X_RR format");
    }

    // Store result
    for (int64_t column = 0; column < rhs_count; ++column) {
        for (int64_t i = 0; i < k; ++i) {
            solve_data.left_side[static_cast<size_t>(
                (*skeleton_indices)[static_cast<size_t>(i)] +
                column * solve_data.num_points)] =
                x_S[static_cast<size_t>(i + column * k)];
        }
    }
}


// ============================================================================
// Phase 3: Backward sweep — apply V operator (root -> leaf)
// ============================================================================
//
// V = P U_T^{-1} L^{-1}
//   Step 1 (L^{-1}):   x[S] -= stored_X_SR * x[R]  (self-box)
//                       x[N] -= stored_X_NR * x[R]  (WRITES to neighbors)
//   Step 2 (U_T^{-1}): x[R] += T^T * x[S]          (self-box, uses updated x[S])
//
// Writes to neighbors -> needs pending updates (like solve's forward sweep).
// stored_X_SR = -X_SR*X_RR^{-1}, so:
//   L^{-1} has +X_SR*X_RR^{-1} on (S,R), meaning x[S] += X_SR*X_RR^{-1}*x[R]
//   = x[S] += (-stored_X_SR)*x[R] = x[S] -= stored_X_SR*x[R]
// Similarly for X_NR.

template<typename CoordType, typename DataType>
void apply_mul_backward_V_with_pending(
    TreeLevel<CoordType, DataType>& level,
    SolveDataRequest<CoordType, DataType>& solve_data,
    std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
    MatrixProperty matrix_property,
    PendingSolveUpdates<DataType>& pending_updates,
    bool is_ghost,
    bool local_or_ghost_targets_only = false) {
    (void)matrix_property;

    const MatrixStorage<DataType>* T = nullptr;
    const MatrixStorage<DataType>* X_SR = nullptr;
    const MatrixStorage<DataType>* X_NR = nullptr;
    const std::vector<int64_t>* skeleton_indices = nullptr;
    const std::vector<int64_t>* redundant_indices = nullptr;
    const std::vector<int64_t>* one_hop = nullptr;
    const std::vector<int64_t>* use_full_set = nullptr;
    BoxData<CoordType, DataType>* box = nullptr;

    if (is_ghost) {
        box = level.find_ghost_box(solve_data.morton_index);
        if (box == nullptr) {
            T = &solve_data.interpolation_matrix;
            X_SR = &solve_data.X_SR;
            X_NR = &solve_data.X_NR;
            skeleton_indices = &solve_data.skeleton_indices;
            redundant_indices = &solve_data.redundant_indices;
            one_hop = &solve_data.one_hop;
            use_full_set = &solve_data.use_full_set;
        } else {
            T = &box->interpolation_matrix;
            X_SR = &box->X_SR;
            X_NR = &box->X_NR;
            skeleton_indices = &box->skeleton_indices;
            redundant_indices = &box->redundant_indices;
            one_hop = &box->one_hop;
            use_full_set = &box->use_full_set;
        }
    } else {
        box = level.find_local_box(solve_data.morton_index);
        if (box == nullptr) {
            throw std::runtime_error(
                "apply_mul_backward_V_with_pending: local box not found");
        }
        T = &box->interpolation_matrix;
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

    const int64_t k = static_cast<int64_t>(skeleton_indices->size());
    const int64_t r = static_cast<int64_t>(redundant_indices->size());
    const int64_t nrhs = solve_data.nrhs;

    // Extract x[R] (read-only for L^{-1})
    std::vector<DataType> x_R(static_cast<size_t>(r * nrhs));
    for (int64_t column = 0; column < nrhs; ++column) {
        for (int64_t i = 0; i < r; ++i) {
            x_R[static_cast<size_t>(i + column * r)] =
                solve_data.left_side[static_cast<size_t>(
                    (*redundant_indices)[static_cast<size_t>(i)] +
                    column * solve_data.num_points)];
        }
    }

    // ===== Step 1a: L^{-1}: x[S] -= stored_X_SR * x[R] =====
    // (Compare solve forward V^{-1}: L: x[S] += stored_X_SR * x[R])
    if (X_SR->is_allocated()) {
        std::vector<DataType> result(
            static_cast<size_t>(k * nrhs), DataType{0.0});
        char trans = 'N';
        int m = static_cast<int>(k), n = static_cast<int>(nrhs);
        int inner = static_cast<int>(r);
        int lda = static_cast<int>(X_SR->lda), ldb = static_cast<int>(r);
        int ldc = static_cast<int>(k);
        DataType alpha = 1.0, beta = 0.0;
        gemm_(&trans, &trans, &m, &n, &inner, &alpha,
              X_SR->data.data(), &lda, x_R.data(), &ldb,
              &beta, result.data(), &ldc);
        for (int64_t column = 0; column < nrhs; ++column) {
            for (int64_t i = 0; i < k; ++i) {
                solve_data.left_side[static_cast<size_t>(
                    (*skeleton_indices)[static_cast<size_t>(i)] +
                    column * solve_data.num_points)] -=
                    result[static_cast<size_t>(i + column * k)];
            }
        }
    }

    // ===== Step 1b: L^{-1}: x[N] -= stored_X_NR * x[R] (writes to neighbors via pending) =====
    // (Compare solve forward V^{-1}: L: x[N] += stored_X_NR * x[R])
    if (!(X_NR->is_allocated() && one_hop != nullptr && use_full_set != nullptr && !one_hop->empty())) {
        // No neighbor interaction; still do Step 2 below.
    } else {
        const int64_t total_neighbor_points = X_NR->rows;
        std::vector<DataType> neighbor_updates(
            static_cast<size_t>(total_neighbor_points * nrhs), DataType{0.0});

        char trans = 'N';
        int m = static_cast<int>(total_neighbor_points);
        int n = static_cast<int>(nrhs), inner = static_cast<int>(r);
        int lda = static_cast<int>(X_NR->lda), ldb = static_cast<int>(r);
        int ldc = static_cast<int>(total_neighbor_points);
        DataType alpha = 1.0, beta = 0.0;
        gemm_(&trans, &trans, &m, &n, &inner, &alpha,
              X_NR->data.data(), &lda, x_R.data(), &ldb,
              &beta, neighbor_updates.data(), &ldc);

        const auto neighbor_sizes = solve_neighbor_sizes_for_box(
            level, level_solve_data, solve_data, *one_hop, *use_full_set,
            is_ghost);

        // Distribute updates to neighbors (sign flipped: subtract instead of add).
        // Remote neighbors are queued to pending_updates and dispatched by the transport step.
        int64_t row_offset = 0;
        for (size_t neighbor_idx = 0; neighbor_idx < one_hop->size(); ++neighbor_idx) {
            const int64_t neighbor_morton = (*one_hop)[neighbor_idx];
            const bool use_full = ((*use_full_set)[neighbor_idx] == 1);

            const int64_t n_neighbor = neighbor_sizes[neighbor_idx];
            if (n_neighbor <= 0) {
                throw std::runtime_error(
                    "apply_mul_backward_V_with_pending: n_neighbor <= 0 for neighbor=" +
                    std::to_string(neighbor_morton));
            }
            if (row_offset + n_neighbor > total_neighbor_points) {
                throw std::runtime_error(
                    "apply_mul_backward_V_with_pending: row partition overflow");
            }

            // Negate the updates before accumulating (sign flip from solve's forward)
            std::vector<DataType> neg_seg(
                static_cast<size_t>(n_neighbor * nrhs));
            for (int64_t column = 0; column < nrhs; ++column) {
                for (int64_t i = 0; i < n_neighbor; ++i) {
                    neg_seg[static_cast<size_t>(i + column * n_neighbor)] =
                        -neighbor_updates[static_cast<size_t>(
                            row_offset + i + column * total_neighbor_points)];
                }
            }
            const int64_t segment_size = n_neighbor * nrhs;

            auto* neighbor_target = resolve_solve_data_for_morton(
                level, level_solve_data, neighbor_morton,
                /*local_or_ghost_only=*/true);
            if (!local_or_ghost_targets_only || neighbor_target != nullptr) {
                if (use_full) {
                    accumulate_replace_then_add(
                        pending_updates.full_updates,
                        neighbor_morton,
                        neg_seg.data(),
                        segment_size);
                } else {
                    accumulate_replace_then_add(
                        pending_updates.skel_updates,
                        neighbor_morton,
                        neg_seg.data(),
                        segment_size);
                }
            }

            row_offset += n_neighbor;
        }

        if (row_offset != total_neighbor_points) {
            throw std::runtime_error(
                "apply_mul_backward_V_with_pending: row_offset " +
                std::to_string(row_offset) + " mismatch");
        }
    }

    // ===== Step 2: U_T^{-1}: x[R] += T^T * x[S] =====
    // (Compare solve forward V^{-1}: U_T: x[R] -= T^T * x[S])
    // Uses the UPDATED x[S] from Step 1a.
    assert(T->is_allocated());
    if (T->is_allocated()) {
        // Re-extract x[S] after Step 1a update
        std::vector<DataType> x_S(static_cast<size_t>(k * nrhs));
        for (int64_t column = 0; column < nrhs; ++column) {
            for (int64_t i = 0; i < k; ++i) {
                x_S[static_cast<size_t>(i + column * k)] =
                    solve_data.left_side[static_cast<size_t>(
                        (*skeleton_indices)[static_cast<size_t>(i)] +
                        column * solve_data.num_points)];
            }
        }

        std::vector<DataType> result(
            static_cast<size_t>(r * nrhs), DataType{0.0});
        char trans_a = 'T', trans_b = 'N';
        int m = static_cast<int>(r), n = static_cast<int>(nrhs);
        int inner = static_cast<int>(k);
        int lda = static_cast<int>(T->lda), ldb = static_cast<int>(k);
        int ldc = static_cast<int>(r);
        DataType alpha = 1.0, beta = 0.0;
        gemm_(&trans_a, &trans_b, &m, &n, &inner, &alpha,
              T->data.data(), &lda, x_S.data(), &ldb,
              &beta, result.data(), &ldc);

        for (int64_t column = 0; column < nrhs; ++column) {
            for (int64_t i = 0; i < r; ++i) {
                solve_data.left_side[static_cast<size_t>(
                    (*redundant_indices)[static_cast<size_t>(i)] +
                    column * solve_data.num_points)] +=
                    result[static_cast<size_t>(i + column * r)];
            }
        }
    }
}

} // namespace fmm
