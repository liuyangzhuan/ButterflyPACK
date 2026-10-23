#pragma once
// The partial factorization (H2_partial, PARTIAL_FACTORIZATION_PLAN.md):
// the F22 phase of a level, the root, the partial solves' F22 steps and the
// Schur complement CB = F22 - F21 F11^-1 F12 by entries.
//
// Variables.  A box b's F22 ID gives T22 with A(:, R22) ~ A(:, S22) T22 (and
// the rows likewise, the joint ID).  With P_b = [I T; 0 I] on the box's
// (S', R22) slots and x~ = P x, the matrix in the new variables is
// M' = P^-T A P^-1 over every box: its R22 rows and columns are the
// sparsified blocks kept here (X_RR, X_RS, X_SR over the box; X_RN, X_NR
// over the one-hop neighbours' raw F22 points, raw since every box's phase
// reads its neighbours before any of them transforms), its S' rows and
// columns go on to the parent.  Nothing is eliminated.  After the F11
// eliminations CB = P22^T CB' P22 with CB' the Schur complement in the
// transformed variables: the sparsified blocks at every level plus the
// root's dense block on the root's F22 slots.  An entry CB(i, j) is
// t_i^T CB' t_j with t = P e the point's chain of vectors up the tree, which
// splits per level into the sparsified part p = r[R22] and w = r[S'] + T p that
// goes up (unstructured/extraction.hpp's walk), so
//   CB(i, j) = sum over levels and near pairs of
//       p_i^T X_RN_b(R_b, b' raw) r_j  +  w_i^T X_NR_b'(S_b raw rows, R_b') p_j
//       (b = b': p_i X_RR p_j + p_i X_RS w_j + w_i X_SR p_j)
//     + r_i(root)^T CB_root r_j(root).
// The solves: forward, after a level's waves, x_R22 -= T22^T x_S' (P^-T);
// backward, before a level's reversed waves, x_S' -= T22 x_R22 (P^-1); the
// forward partial solve returns b2 - F21 F11^-1 b1 = P22^T y~2 through a
// down-pass (x_R22 += T22^T x_S'), the backward one takes x2 through an
// up-pass (x_S' += T22 x_R22).
//
// Milestone P2a: one rank (the neighbours and parents local).

#include "core/tree.hpp"
#include "core/blas_declare.hpp"
#include "core/factorization.hpp"
#include "core/id_decomposition.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace butterfly {
namespace color_unstructured {

template<typename DataType>
MPI_Datatype partial_mpi_type() {
    if constexpr (std::is_same_v<DataType, double>) return MPI_DOUBLE;
    else if constexpr (std::is_same_v<DataType, float>) return MPI_FLOAT;
    else if constexpr (std::is_same_v<DataType, std::complex<float>>) return MPI_C_FLOAT_COMPLEX;
    else return MPI_C_DOUBLE_COMPLEX;
}

// C = alpha op(A) op(B) + beta C, column-major, A m x k (lda), B k x n (ldb).
template<typename DataType>
void partial_gemm(char ta, char tb, int64_t m, int64_t n, int64_t k, DataType alpha,
                  const DataType* a, int64_t lda, const DataType* b, int64_t ldb, DataType beta,
                  DataType* c, int64_t ldc) {
    if (m == 0 || n == 0) return;
    int M = static_cast<int>(m), N = static_cast<int>(n), K = static_cast<int>(k);
    int LDA = static_cast<int>(std::max<int64_t>(lda, 1)), LDB = static_cast<int>(std::max<int64_t>(ldb, 1)),
        LDC = static_cast<int>(std::max<int64_t>(ldc, 1));
    if (k == 0) {
        for (int64_t j = 0; j < n; ++j)
            for (int64_t i = 0; i < m; ++i) c[i + j * ldc] *= beta;
        return;
    }
    if constexpr (std::is_same_v<DataType, double>) {
        dgemm_(&ta, &tb, &M, &N, &K, &alpha, const_cast<double*>(a), &LDA, const_cast<double*>(b), &LDB, &beta, c,
               &LDC);
    } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
        zgemm_(&ta, &tb, &M, &N, &K, &alpha, const_cast<std::complex<double>*>(a), &LDA,
               const_cast<std::complex<double>*>(b), &LDB, &beta, c, &LDC);
    } else {
        throw std::logic_error("partial factorization: double precision only");
    }
}

template<typename CoordType, typename DataType>
void partial_require_single_rank(const fmm::ParallelTree<CoordType, DataType>& tree, const char* what) {
    if (tree.mpi_size != 1) {
        throw std::runtime_error(std::string("H2_partial: ") + what + " runs on one rank for now (milestone P2a)");
    }
}

// A box's offset in its parent's point list (the children's skeletons in Morton order).
template<typename CoordType, typename DataType>
int64_t partial_child_offset(fmm::TreeLevel<CoordType, DataType>& child_level, int64_t morton, int dimension) {
    const int64_t parent = morton >> dimension;
    int64_t offset = 0;
    for (int64_t sib = parent << dimension; sib < morton; ++sib) {
        if (const auto* sb = child_level.find_local_box(sib)) offset += static_cast<int64_t>(sb->skeleton_indices.size());
    }
    return offset;
}

// =====================================================================
// The F22 phase of a level (after its F11 waves): every local box's F22
// part is skeletonized against the far ring and the near field's live F11
// points, its R22 sparsified, the box's live set reduced to S' = S11, S22.
// =====================================================================
template<typename CoordType, typename DataType, typename KernelType>
void partial_f22_phase(fmm::ParallelTree<CoordType, DataType>* tree, fmm::TreeLevel<CoordType, DataType>& level,
                       KernelType* kernel, double tolerance, const std::vector<int64_t>& occupied_local,
                       bool print_detail) {
    using Box = fmm::BoxData<CoordType, DataType>;
    using MS = fmm::MatrixStorage<DataType>;
    partial_require_single_rank(*tree, "the F22 phase");
    const double abs_tol = fmm::id_absolute_tolerance(tree, tolerance);
    // the live sets before the phase, for every box (the records read the
    // neighbours' raw layouts while they may already have been transformed)
    for (int64_t li : occupied_local) {
        Box& b = level.local_boxes[static_cast<size_t>(li)];
        b.skeleton_before_f22 = b.skeleton_indices;
        b.f22 = typename Box::F22Part{};
    }
    auto a_ns_block = [&](Box& b, int64_t nb) -> fmm::ModifiedBlock<DataType>* {  // rows nb, cols b
        auto it = b.near_field_interaction_map.find(nb);
        return it == b.near_field_interaction_map.end() ? nullptr : &b.near_field_modified_interactions[it->second];
    };
    auto a_sn_block = [&](Box& b, int64_t nb) -> fmm::ModifiedBlock<DataType>* {  // rows b, cols nb
        auto it = b.near_field_interaction_map_nonsymmetry.find(nb);
        return it == b.near_field_interaction_map_nonsymmetry.end() ? nullptr
                                                                    : &b.near_field_modified_interactions[it->second];
    };
    std::exception_ptr failure;
    int64_t sparsified_boxes = 0, sparsified_points = 0, kept_points = 0;
    // the boxes' kept positions within S_old; the shared state (live sets,
    // the boxes' own blocks) is changed after the loop, the loop only reads
    // the neighbours and rings
    std::vector<std::vector<int64_t>> kept(occupied_local.size());
#pragma omp parallel
    {
        fmm::FactorizationThreadScratch<CoordType, DataType> scratch;
#pragma omp for schedule(dynamic, 1) reduction(+ : sparsified_boxes, sparsified_points, kept_points)
        for (int64_t oi = 0; oi < static_cast<int64_t>(occupied_local.size()); ++oi) {
            try {
                Box& b = level.local_boxes[static_cast<size_t>(occupied_local[static_cast<size_t>(oi)])];
                const std::vector<int64_t> S_old = b.skeleton_before_f22;
                const int64_t k_old = static_cast<int64_t>(S_old.size());
                std::vector<int64_t> pos11, pos22;  // positions within S_old
                for (int64_t q = 0; q < k_old; ++q) {
                    if (tree->partial_f11(b.point_indices[static_cast<size_t>(S_old[static_cast<size_t>(q)])])) pos11.push_back(q);
                    else pos22.push_back(q);
                }
                const int64_t n22 = static_cast<int64_t>(pos22.size()), k11 = static_cast<int64_t>(pos11.size());
                if (n22 == 0) continue;
                if (!b.schur_complement.is_allocated() || b.schur_complement.rows != k_old) {
                    throw std::runtime_error("F22 phase: box " + std::to_string(b.morton_index) +
                                             " has no Schur block over its skeleton");
                }
                const MS& A = b.schur_complement;
                auto Aat = [&](int64_t i, int64_t j) { return A.data[static_cast<size_t>(i + j * A.lda)]; };
                // ---- the ID target: the far ring's sketch (columns = S_old),
                // the near field's live F11 rows and the box's own, both orientations
                // (the gather resets use_full_set for eliminated neighbours; the
                // F11 step's flags describe its X_NR rows and must survive)
                const std::vector<int64_t> use_full_at_f11 = b.use_full_set;
                fmm::gather_id_target_streamed_nonsymmetric(tree, &b, level, kernel, scratch, b.on_boundary, tolerance,
                                                             &S_old);
                b.use_full_set = use_full_at_f11;
                if (!scratch.streamed_sketch_valid || scratch.workspace_cols != k_old) {
                    throw std::runtime_error("F22 phase: no streamed sketch for box " + std::to_string(b.morton_index));
                }
                const int64_t d = scratch.streamed_sketch_rows;
                struct NearF11 { Box* nb; std::vector<int64_t> q11, q22; };  // positions within nb's S_old
                std::vector<NearF11> near(b.one_hop.size());
                int64_t near_rows = 0;
                for (size_t idx = 0; idx < b.one_hop.size(); ++idx) {
                    Box* nb = level.find_local_box(b.one_hop[idx]);
                    if (nb == nullptr) throw std::runtime_error("F22 phase: a neighbour is not local");
                    near[idx].nb = nb;
                    const auto& nS = nb->skeleton_before_f22;
                    for (int64_t q = 0; q < static_cast<int64_t>(nS.size()); ++q) {
                        if (tree->partial_f11(nb->point_indices[static_cast<size_t>(nS[static_cast<size_t>(q)])])) near[idx].q11.push_back(q);
                        else near[idx].q22.push_back(q);
                    }
                    near_rows += static_cast<int64_t>(near[idx].q11.size());
                }
                const int64_t rows = d + 2 * (near_rows + k11);
                std::vector<DataType> target(static_cast<size_t>(rows * n22), DataType(0));
                for (int64_t c = 0; c < n22; ++c) {
                    const int64_t col = pos22[static_cast<size_t>(c)];
                    std::copy_n(scratch.sketch_storage.data() + col * d, d, target.data() + c * rows);
                }
                int64_t at = d;
                for (size_t idx = 0; idx < b.one_hop.size(); ++idx) {  // K(F11_nb, F22_b)
                    auto* blk = a_ns_block(b, b.one_hop[idx]);
                    if (blk == nullptr || !blk->a_ns_is_allocated()) continue;
                    if (blk->a_ns_cols() != k_old) throw std::runtime_error("F22 phase: a neighbour's rows are not over the skeleton");
                    if (blk->a_ns_rows() != static_cast<int64_t>(near[idx].q11.size() + near[idx].q22.size())) {
                        throw std::runtime_error("F22 phase: a neighbour's rows are not over its skeleton (" +
                                                 std::to_string(blk->a_ns_rows()) + " rows, box " +
                                                 std::to_string(b.morton_index) + ")");
                    }
                    for (int64_t c = 0; c < n22; ++c)
                        for (size_t i = 0; i < near[idx].q11.size(); ++i)
                            target[static_cast<size_t>(at + static_cast<int64_t>(i) + c * rows)] =
                                blk->a_ns(near[idx].q11[i], pos22[static_cast<size_t>(c)]);
                    at += static_cast<int64_t>(near[idx].q11.size());
                }
                for (int64_t c = 0; c < n22; ++c)  // the box's own S11 rows
                    for (int64_t i = 0; i < k11; ++i)
                        target[static_cast<size_t>(at + i + c * rows)] = Aat(pos11[static_cast<size_t>(i)], pos22[static_cast<size_t>(c)]);
                at += k11;
                for (size_t idx = 0; idx < b.one_hop.size(); ++idx) {  // K(F22_b, F11_nb)^T
                    auto* blk = a_sn_block(b, b.one_hop[idx]);
                    if (blk == nullptr || !blk->A_SN.is_allocated()) continue;
                    if (blk->A_SN.rows != k_old) throw std::runtime_error("F22 phase: the box's rows are not over the skeleton");
                    if (blk->A_SN.cols != static_cast<int64_t>(near[idx].q11.size() + near[idx].q22.size())) {
                        throw std::runtime_error("F22 phase: a neighbour's columns are not over its skeleton");
                    }
                    for (int64_t c = 0; c < n22; ++c)
                        for (size_t i = 0; i < near[idx].q11.size(); ++i)
                            target[static_cast<size_t>(at + static_cast<int64_t>(i) + c * rows)] =
                                blk->A_SN(pos22[static_cast<size_t>(c)], near[idx].q11[i]);
                    at += static_cast<int64_t>(near[idx].q11.size());
                }
                for (int64_t c = 0; c < n22; ++c)
                    for (int64_t i = 0; i < k11; ++i)
                        target[static_cast<size_t>(at + i + c * rows)] = Aat(pos22[static_cast<size_t>(c)], pos11[static_cast<size_t>(i)]);
                at += k11;
                if (at != rows) throw std::logic_error("F22 phase: target rows");
                fmm::IDResult<DataType> id = fmm::compute_id_complex(target.data(), rows, n22, rows, tolerance, 0, abs_tol);
                fmm::ensure_nonempty_box_id_rank(id, n22);
                const int64_t k22 = static_cast<int64_t>(id.skeleton_indices.size());
                const int64_t r22 = static_cast<int64_t>(id.redundant_indices.size());
                // ---- the new live set S' = S11 (its order), S22 (the ID's) and R22, as positions and slots
                std::vector<int64_t> posS(pos11), posR;
                for (int64_t s : id.skeleton_indices) posS.push_back(pos22[static_cast<size_t>(s)]);
                for (int64_t s : id.redundant_indices) posR.push_back(pos22[static_cast<size_t>(s)]);
                const int64_t kp = static_cast<int64_t>(posS.size());
                auto& f = b.f22;
                f.sparsified = true;
                for (int64_t s : id.skeleton_indices) f.skeleton.push_back(S_old[static_cast<size_t>(pos22[static_cast<size_t>(s)])]);
                for (int64_t s : id.redundant_indices) f.redundant.push_back(S_old[static_cast<size_t>(pos22[static_cast<size_t>(s)])]);
                // T over S' rows (zero on S11)
                f.T.allocate(kp, r22, MS::FULL);
                for (int64_t j = 0; j < r22; ++j)
                    for (int64_t i = 0; i < k22; ++i) f.T.data[static_cast<size_t>(k11 + i + j * kp)] = id.interpolation(i, j);
                const DataType* T = f.T.data.data();
                // ---- the box's own blocks over (S', R): X_RR = A_RR - T^T A_SR - A_RS T + T^T A_SS T,
                // X_RS = A_RS - T^T A_SS, X_SR = A_SR - A_SS T
                std::vector<DataType> A_SS(static_cast<size_t>(kp * kp)), A_SR(static_cast<size_t>(kp * r22)),
                    A_RS(static_cast<size_t>(r22 * kp)), A_RR(static_cast<size_t>(r22 * r22));
                for (int64_t j = 0; j < kp; ++j)
                    for (int64_t i = 0; i < kp; ++i) A_SS[static_cast<size_t>(i + j * kp)] = Aat(posS[static_cast<size_t>(i)], posS[static_cast<size_t>(j)]);
                for (int64_t j = 0; j < r22; ++j)
                    for (int64_t i = 0; i < kp; ++i) A_SR[static_cast<size_t>(i + j * kp)] = Aat(posS[static_cast<size_t>(i)], posR[static_cast<size_t>(j)]);
                for (int64_t j = 0; j < kp; ++j)
                    for (int64_t i = 0; i < r22; ++i) A_RS[static_cast<size_t>(i + j * r22)] = Aat(posR[static_cast<size_t>(i)], posS[static_cast<size_t>(j)]);
                for (int64_t j = 0; j < r22; ++j)
                    for (int64_t i = 0; i < r22; ++i) A_RR[static_cast<size_t>(i + j * r22)] = Aat(posR[static_cast<size_t>(i)], posR[static_cast<size_t>(j)]);
                f.X_SR.allocate(kp, r22, MS::FULL);
                f.X_SR.data = A_SR;
                partial_gemm<DataType>('N', 'N', kp, r22, kp, DataType(-1), A_SS.data(), kp, T, kp, DataType(1), f.X_SR.data.data(), kp);
                f.X_RS.allocate(r22, kp, MS::FULL);
                f.X_RS.data = A_RS;
                partial_gemm<DataType>('T', 'N', r22, kp, kp, DataType(-1), T, kp, A_SS.data(), kp, DataType(1), f.X_RS.data.data(), r22);
                f.X_RR.allocate(r22, r22, MS::FULL);
                f.X_RR.data = A_RR;
                // X_RR = A_RR - T^T A_SR - X_RS T  (X_RS T = A_RS T - T^T A_SS T)
                partial_gemm<DataType>('T', 'N', r22, r22, kp, DataType(-1), T, kp, A_SR.data(), kp, DataType(1), f.X_RR.data.data(), r22);
                partial_gemm<DataType>('N', 'N', r22, r22, kp, DataType(-1), f.X_RS.data.data(), r22, T, kp, DataType(1), f.X_RR.data.data(), r22);
                // ---- the neighbours' raw F22 points: X_NR = A_NR - A_NS T, X_RN = A_RN - T^T A_SN
                f.neighbor_f22_counts.assign(b.one_hop.size(), 0);
                f.neighbor_f22_slots.assign(b.one_hop.size(), {});
                int64_t total = 0;
                for (size_t idx = 0; idx < b.one_hop.size(); ++idx) {
                    const auto& nS = near[idx].nb->skeleton_before_f22;
                    for (int64_t q : near[idx].q22) f.neighbor_f22_slots[idx].push_back(nS[static_cast<size_t>(q)]);
                    f.neighbor_f22_counts[idx] = static_cast<int64_t>(near[idx].q22.size());
                    total += f.neighbor_f22_counts[idx];
                }
                f.X_NR.allocate(total, r22, MS::FULL);
                f.X_RN.allocate(r22, total, MS::FULL);
                int64_t row0 = 0;
                for (size_t idx = 0; idx < b.one_hop.size(); ++idx) {
                    const int64_t m = f.neighbor_f22_counts[idx];
                    if (m == 0) continue;
                    auto* ns = a_ns_block(b, b.one_hop[idx]);
                    auto* sn = a_sn_block(b, b.one_hop[idx]);
                    if (ns == nullptr || sn == nullptr || !ns->a_ns_is_allocated() || !sn->A_SN.is_allocated()) {
                        throw std::runtime_error("F22 phase: a near pair's views are missing (box " +
                                                 std::to_string(b.morton_index) + ")");
                    }
                    std::vector<DataType> A_NR(static_cast<size_t>(m * r22)), A_NS(static_cast<size_t>(m * kp)),
                        A_RN(static_cast<size_t>(r22 * m)), A_SN(static_cast<size_t>(kp * m));
                    for (int64_t i = 0; i < m; ++i) {
                        const int64_t q = near[idx].q22[static_cast<size_t>(i)];
                        for (int64_t j = 0; j < r22; ++j) A_NR[static_cast<size_t>(i + j * m)] = ns->a_ns(q, posR[static_cast<size_t>(j)]);
                        for (int64_t j = 0; j < kp; ++j) A_NS[static_cast<size_t>(i + j * m)] = ns->a_ns(q, posS[static_cast<size_t>(j)]);
                        for (int64_t j = 0; j < r22; ++j) A_RN[static_cast<size_t>(j + i * r22)] = sn->A_SN(posR[static_cast<size_t>(j)], q);
                        for (int64_t j = 0; j < kp; ++j) A_SN[static_cast<size_t>(j + i * kp)] = sn->A_SN(posS[static_cast<size_t>(j)], q);
                    }
                    partial_gemm<DataType>('N', 'N', m, r22, kp, DataType(-1), A_NS.data(), m, T, kp, DataType(1), A_NR.data(), m);
                    partial_gemm<DataType>('T', 'N', r22, m, kp, DataType(-1), T, kp, A_SN.data(), kp, DataType(1), A_RN.data(), r22);
                    for (int64_t j = 0; j < r22; ++j)
                        std::copy_n(A_NR.data() + j * m, m, f.X_NR.data.data() + row0 + j * total);
                    for (int64_t i = 0; i < m; ++i)
                        std::copy_n(A_RN.data() + i * r22, r22, f.X_RN.data.data() + (row0 + i) * r22);
                    row0 += m;
                }
                kept[static_cast<size_t>(oi)] = std::move(posS);
                ++sparsified_boxes;
                sparsified_points += r22;
                kept_points += k22;
            } catch (...) {
#pragma omp critical(partial_f22_phase_failure)
                if (!failure) failure = std::current_exception();
            }
        }
    }
    if (failure) std::rethrow_exception(failure);
    // ---- the boxes' live sets and their own blocks over S'
    for (size_t oi = 0; oi < occupied_local.size(); ++oi) {
        Box& b = level.local_boxes[static_cast<size_t>(occupied_local[oi])];
        if (!b.f22.sparsified) continue;
        const auto& posS = kept[oi];
        const int64_t kp = static_cast<int64_t>(posS.size());
        {
            const MS& A = b.schur_complement;
            MS sliced;
            sliced.allocate(kp, kp, MS::FULL);
            for (int64_t j = 0; j < kp; ++j)
                for (int64_t i = 0; i < kp; ++i)
                    sliced.data[static_cast<size_t>(i + j * kp)] =
                        A.data[static_cast<size_t>(posS[static_cast<size_t>(i)] + posS[static_cast<size_t>(j)] * A.lda)];
            b.schur_complement = std::move(sliced);
        }
        for (size_t idx = 0; idx < b.one_hop.size(); ++idx) {
            if (auto* ns = a_ns_block(b, b.one_hop[idx]); ns != nullptr && ns->a_ns_is_allocated()) {
                const int64_t m = ns->a_ns_rows();
                std::vector<DataType> v(static_cast<size_t>(m * kp));
                for (int64_t j = 0; j < kp; ++j)
                    for (int64_t i = 0; i < m; ++i) v[static_cast<size_t>(i + j * m)] = ns->a_ns(i, posS[static_cast<size_t>(j)]);
                ns->set_a_ns_owned(m, kp, std::move(v), MS::FULL);
            }
            if (auto* sn = a_sn_block(b, b.one_hop[idx]); sn != nullptr && sn->A_SN.is_allocated()) {
                const int64_t m = sn->A_SN.cols;
                std::vector<DataType> v(static_cast<size_t>(kp * m));
                for (int64_t j = 0; j < m; ++j)
                    for (int64_t i = 0; i < kp; ++i) v[static_cast<size_t>(i + j * kp)] = sn->A_SN(posS[static_cast<size_t>(i)], j);
                sn->A_SN.set_owned(kp, m, std::move(v), MS::FULL);
            }
        }
        b.skeleton_indices.clear();
        for (int64_t q : posS) b.skeleton_indices.push_back(b.skeleton_before_f22[static_cast<size_t>(q)]);
    }
    // ---- the neighbours' copies of the views: the sparsified boxes' rows
    // and columns to their new live sets
    for (int64_t li : occupied_local) {
        Box& c = level.local_boxes[static_cast<size_t>(li)];
        for (size_t idx = 0; idx < c.one_hop.size(); ++idx) {
            Box* nb = level.find_local_box(c.one_hop[idx]);
            if (nb == nullptr || !nb->f22.sparsified) continue;
            std::vector<int64_t> pos;  // nb's new skeleton within its old one
            const auto& old = nb->skeleton_before_f22;
            for (int64_t s : nb->skeleton_indices) {
                const auto it = std::find(old.begin(), old.end(), s);
                if (it == old.end()) throw std::logic_error("F22 phase: a skeleton slot outside the old skeleton");
                pos.push_back(static_cast<int64_t>(it - old.begin()));
            }
            const int64_t kn = static_cast<int64_t>(pos.size());
            if (auto* sn = a_sn_block(c, nb->morton_index); sn != nullptr && sn->A_SN.is_allocated()) {  // c rows x nb cols
                if (sn->A_SN.cols != static_cast<int64_t>(old.size())) {
                    throw std::runtime_error("F22 phase: a neighbour's view is not over its old skeleton (columns)");
                }
                const int64_t m = sn->A_SN.rows;
                std::vector<DataType> v(static_cast<size_t>(m * kn));
                for (int64_t j = 0; j < kn; ++j)
                    for (int64_t i = 0; i < m; ++i) v[static_cast<size_t>(i + j * m)] = sn->A_SN(i, pos[static_cast<size_t>(j)]);
                sn->A_SN.set_owned(m, kn, std::move(v), MS::FULL);
            }
            if (auto* ns = a_ns_block(c, nb->morton_index); ns != nullptr && ns->a_ns_is_allocated()) {  // nb rows x c cols
                if (ns->a_ns_rows() != static_cast<int64_t>(old.size())) {
                    throw std::runtime_error("F22 phase: a neighbour's view is not over its old skeleton (rows)");
                }
                const int64_t m = ns->a_ns_cols();
                std::vector<DataType> v(static_cast<size_t>(kn * m));
                for (int64_t j = 0; j < m; ++j)
                    for (int64_t i = 0; i < kn; ++i) v[static_cast<size_t>(i + j * kn)] = ns->a_ns(pos[static_cast<size_t>(i)], j);
                ns->set_a_ns_owned(kn, m, std::move(v), MS::FULL);
            }
        }
    }
    if (print_detail) {
        std::printf("  F22 phase: %lld boxes, %lld points sparsified, %lld kept\n",
                    static_cast<long long>(sparsified_boxes), static_cast<long long>(sparsified_points),
                    static_cast<long long>(kept_points));
    }
}

// =====================================================================
// The root: the root matrix split by point id, A11 = P L U in the root
// box's X_RR (the logdet of F11 reads it), W = A11^-1 A12, A21 and CB_root
// = A22 - A21 W on the root's F22 slots.
// =====================================================================
template<typename DataType>
void partial_root_solve11(int64_t n11, const fmm::MatrixStorage<DataType>& lu, const std::vector<int>& pivots,
                          DataType* b, int64_t ncol, int64_t ldb) {
    if (n11 == 0 || ncol == 0) return;
    char trans = 'N';
    int n = static_cast<int>(n11), nrhs = static_cast<int>(ncol), lda = static_cast<int>(lu.lda),
        LDB = static_cast<int>(ldb), info = 0;
    if constexpr (std::is_same_v<DataType, double>) {
        dgetrs_(&trans, &n, &nrhs, const_cast<double*>(lu.data.data()), &lda, const_cast<int*>(pivots.data()), b,
                &LDB, &info);
    } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
        zgetrs_(&trans, &n, &nrhs, const_cast<std::complex<double>*>(lu.data.data()), &lda,
                const_cast<int*>(pivots.data()), b, &LDB, &info);
    } else {
        throw std::logic_error("partial factorization: double precision only");
    }
    if (info != 0) throw std::runtime_error("partial root: getrs failed, info " + std::to_string(info));
}

template<typename CoordType, typename DataType>
void partial_root_factorize(fmm::ParallelTree<CoordType, DataType>& tree,
                            fmm::BoxData<CoordType, DataType>& root_box) {
    using MS = fmm::MatrixStorage<DataType>;
    auto& pr = tree.partial_root;
    const int64_t n = root_box.num_points;
    const MS& A = root_box.schur_complement;
    if (!A.is_allocated() || A.rows != n || A.cols != n) {
        throw std::runtime_error("partial root: the root's Schur complement does not cover its points");
    }
    pr = fmm::PartialRootFactor<DataType>{};
    pr.active = true;
    for (int64_t s = 0; s < n; ++s) {
        const int64_t id = root_box.point_indices[static_cast<size_t>(s)];
        if (tree.partial_f11(id)) pr.slots11.push_back(s);
        else pr.slots22.push_back(s);
    }
    pr.n11 = static_cast<int64_t>(pr.slots11.size());
    pr.n22 = static_cast<int64_t>(pr.slots22.size());
    pr.ids22.resize(static_cast<size_t>(pr.n22));
    for (int64_t j = 0; j < pr.n22; ++j) {
        pr.ids22[static_cast<size_t>(j)] = root_box.point_indices[static_cast<size_t>(pr.slots22[static_cast<size_t>(j)])];
    }
    pr.n22_global = tree.num_points - tree.partial_count;
    auto gather = [&](const std::vector<int64_t>& rows, const std::vector<int64_t>& cols, MS& out) {
        const int64_t m = static_cast<int64_t>(rows.size()), k = static_cast<int64_t>(cols.size());
        out.allocate(m, k, MS::FULL);
        for (int64_t j = 0; j < k; ++j)
            for (int64_t i = 0; i < m; ++i)
                out.data[static_cast<size_t>(i + j * m)] =
                    A.data[static_cast<size_t>(rows[static_cast<size_t>(i)] + cols[static_cast<size_t>(j)] * A.lda)];
    };
    gather(pr.slots11, pr.slots11, root_box.X_RR);
    root_box.X_RR_pivots.assign(static_cast<size_t>(std::max<int64_t>(pr.n11, 1)), 0);
    if (pr.n11 > 0) {
        int nn = static_cast<int>(pr.n11), info = 0;
        if constexpr (std::is_same_v<DataType, double>) {
            dgetrf_(&nn, &nn, root_box.X_RR.data.data(), &nn, root_box.X_RR_pivots.data(), &info);
        } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
            zgetrf_(&nn, &nn, root_box.X_RR.data.data(), &nn, root_box.X_RR_pivots.data(), &info);
        } else {
            throw std::logic_error("partial factorization: double precision only");
        }
        if (info != 0) {
            throw std::runtime_error("partial root: LU of the F11 block failed at pivot " + std::to_string(info));
        }
    }
    root_box.X_RR.format = MS::LU_FACTORED;
    gather(pr.slots11, pr.slots22, pr.w);
    partial_root_solve11(pr.n11, root_box.X_RR, root_box.X_RR_pivots, pr.w.data.data(), pr.n22, std::max<int64_t>(pr.n11, 1));
    gather(pr.slots22, pr.slots11, pr.a21);
    gather(pr.slots22, pr.slots22, pr.cb);
    partial_gemm<DataType>('N', 'N', pr.n22, pr.n22, pr.n11, DataType(-1.0), pr.a21.data.data(), std::max<int64_t>(pr.n22, 1),
                           pr.w.data.data(), std::max<int64_t>(pr.n11, 1), DataType(1.0), pr.cb.data.data(),
                           std::max<int64_t>(pr.n22, 1));
    pr.cb_released = false;
}

// The forward partial solve on the root: v = [y1; y2] -> [A11^-1 y1; y2 - A21 A11^-1 y1].
template<typename DataType>
void partial_root_forward(const fmm::PartialRootFactor<DataType>& pr, const fmm::MatrixStorage<DataType>& lu,
                          const std::vector<int>& pivots, std::vector<DataType>& v, int64_t n, int64_t ncol) {
    std::vector<DataType> b1(static_cast<size_t>(std::max<int64_t>(pr.n11, 1)) * static_cast<size_t>(ncol));
    std::vector<DataType> b2(static_cast<size_t>(std::max<int64_t>(pr.n22, 1)) * static_cast<size_t>(ncol));
    for (int64_t c = 0; c < ncol; ++c) {
        for (int64_t i = 0; i < pr.n11; ++i)
            b1[static_cast<size_t>(i + c * pr.n11)] = v[static_cast<size_t>(pr.slots11[static_cast<size_t>(i)] + c * n)];
        for (int64_t i = 0; i < pr.n22; ++i)
            b2[static_cast<size_t>(i + c * pr.n22)] = v[static_cast<size_t>(pr.slots22[static_cast<size_t>(i)] + c * n)];
    }
    partial_root_solve11(pr.n11, lu, pivots, b1.data(), ncol, std::max<int64_t>(pr.n11, 1));
    partial_gemm<DataType>('N', 'N', pr.n22, ncol, pr.n11, DataType(-1.0), pr.a21.data.data(), std::max<int64_t>(pr.n22, 1),
                           b1.data(), std::max<int64_t>(pr.n11, 1), DataType(1.0), b2.data(), std::max<int64_t>(pr.n22, 1));
    for (int64_t c = 0; c < ncol; ++c) {
        for (int64_t i = 0; i < pr.n11; ++i)
            v[static_cast<size_t>(pr.slots11[static_cast<size_t>(i)] + c * n)] = b1[static_cast<size_t>(i + c * pr.n11)];
        for (int64_t i = 0; i < pr.n22; ++i)
            v[static_cast<size_t>(pr.slots22[static_cast<size_t>(i)] + c * n)] = b2[static_cast<size_t>(i + c * pr.n22)];
    }
}

// The backward partial solve on the root: v = [y1; x2] -> [y1 - W x2; x2].
template<typename DataType>
void partial_root_backward(const fmm::PartialRootFactor<DataType>& pr, std::vector<DataType>& v, int64_t n,
                           int64_t ncol) {
    std::vector<DataType> x1(static_cast<size_t>(std::max<int64_t>(pr.n11, 1)) * static_cast<size_t>(ncol));
    std::vector<DataType> x2(static_cast<size_t>(std::max<int64_t>(pr.n22, 1)) * static_cast<size_t>(ncol));
    for (int64_t c = 0; c < ncol; ++c) {
        for (int64_t i = 0; i < pr.n11; ++i)
            x1[static_cast<size_t>(i + c * pr.n11)] = v[static_cast<size_t>(pr.slots11[static_cast<size_t>(i)] + c * n)];
        for (int64_t i = 0; i < pr.n22; ++i)
            x2[static_cast<size_t>(i + c * pr.n22)] = v[static_cast<size_t>(pr.slots22[static_cast<size_t>(i)] + c * n)];
    }
    partial_gemm<DataType>('N', 'N', pr.n11, ncol, pr.n22, DataType(-1.0), pr.w.data.data(), std::max<int64_t>(pr.n11, 1),
                           x2.data(), std::max<int64_t>(pr.n22, 1), DataType(1.0), x1.data(), std::max<int64_t>(pr.n11, 1));
    for (int64_t c = 0; c < ncol; ++c)
        for (int64_t i = 0; i < pr.n11; ++i)
            v[static_cast<size_t>(pr.slots11[static_cast<size_t>(i)] + c * n)] = x1[static_cast<size_t>(i + c * pr.n11)];
}

template<typename CoordType, typename DataType>
int partial_root_rank(const fmm::ParallelTree<CoordType, DataType>& tree) {
    int rank = 0;
    MPI_Comm_rank(tree.comm, &rank);
    int mine = tree.partial_root.active ? rank : -1;
    int owner = -1;
    MPI_Allreduce(&mine, &owner, 1, MPI_INT, MPI_MAX, tree.comm);
    if (owner < 0) throw std::runtime_error("partial root: no rank holds the root");
    return owner;
}

// =====================================================================
// The solves' F22 steps (a level's local boxes' vectors, column-major with
// num_points rows a column).
// =====================================================================
// forward, after the level's waves: x_R22 -= T22^T x_S'  (P^-T)
template<typename CoordType, typename DataType>
void partial_f22_forward_level(fmm::TreeLevel<CoordType, DataType>& level,
                               std::vector<fmm::SolveDataRequest<CoordType, DataType>>& data, int nrhs) {
    for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
        const auto& box = level.local_boxes[static_cast<size_t>(bi)];
        if (!box.f22.sparsified || box.f22.redundant.empty()) continue;
        auto& v = data[static_cast<size_t>(bi)].left_side;
        const int64_t n = box.num_points, kp = box.f22.T.rows, r = box.f22.T.cols;
        for (int c = 0; c < nrhs; ++c) {
            for (int64_t j = 0; j < r; ++j) {
                DataType acc(0);
                for (int64_t i = 0; i < kp; ++i)
                    acc += box.f22.T.data[static_cast<size_t>(i + j * kp)] *
                           v[static_cast<size_t>(box.skeleton_indices[static_cast<size_t>(i)] + c * n)];
                v[static_cast<size_t>(box.f22.redundant[static_cast<size_t>(j)] + c * n)] -= acc;
            }
        }
    }
}
// backward, before the level's reversed waves: x_S' -= T22 x_R22  (P^-1)
template<typename CoordType, typename DataType>
void partial_f22_backward_level(fmm::TreeLevel<CoordType, DataType>& level,
                                std::vector<fmm::SolveDataRequest<CoordType, DataType>>& data, int nrhs) {
    for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
        const auto& box = level.local_boxes[static_cast<size_t>(bi)];
        if (!box.f22.sparsified || box.f22.redundant.empty()) continue;
        auto& v = data[static_cast<size_t>(bi)].left_side;
        const int64_t n = box.num_points, kp = box.f22.T.rows, r = box.f22.T.cols;
        for (int c = 0; c < nrhs; ++c) {
            for (int64_t i = 0; i < kp; ++i) {
                DataType acc(0);
                for (int64_t j = 0; j < r; ++j)
                    acc += box.f22.T.data[static_cast<size_t>(i + j * kp)] *
                           v[static_cast<size_t>(box.f22.redundant[static_cast<size_t>(j)] + c * n)];
                v[static_cast<size_t>(box.skeleton_indices[static_cast<size_t>(i)] + c * n)] -= acc;
            }
        }
    }
}

// The forward partial solve's output: from the swept vectors (kept), a
// down-pass applying P22^T (x_R22 += T22^T x_S') to a copy, the leaves' F22
// entries in F22 index order (replicated).
template<typename CoordType, typename DataType>
void partial_forward_output(fmm::ParallelTree<CoordType, DataType>* tree,
                            const std::vector<std::vector<fmm::SolveDataRequest<CoordType, DataType>>>& solve_data,
                            int nrhs, std::vector<DataType>& out2) {
    partial_require_single_rank(*tree, "the forward partial solve's output");
    const int64_t pc = tree->partial_count, n22 = tree->num_points - pc;
    const int dimension = tree->dimension;
    const int leaf = tree->num_levels - 1;
    std::vector<std::vector<std::vector<DataType>>> v(static_cast<size_t>(tree->num_levels));
    for (int l = 0; l <= leaf; ++l) {
        auto& level = tree->levels[static_cast<size_t>(l)];
        v[static_cast<size_t>(l)].resize(static_cast<size_t>(level.num_boxes_local));
        for (int64_t bi = 0; bi < level.num_boxes_local; ++bi)
            v[static_cast<size_t>(l)][static_cast<size_t>(bi)] = solve_data[static_cast<size_t>(l)][static_cast<size_t>(bi)].left_side;
    }
    for (int l = 1; l <= leaf; ++l) {
        auto& level = tree->levels[static_cast<size_t>(l)];
        auto& parent_level = tree->levels[static_cast<size_t>(l) - 1];
        for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
            const auto& box = level.local_boxes[static_cast<size_t>(bi)];
            if (box.skeleton_indices.empty()) continue;
            const int64_t parent = box.morton_index >> dimension;
            const int64_t pbi = parent - parent_level.local_morton_start;
            const auto* pbox = parent_level.find_local_box(parent);
            if (pbox == nullptr) throw std::logic_error("partial_forward_output: a parent is not local");
            const int64_t offset = partial_child_offset(level, box.morton_index, dimension);
            auto& x = v[static_cast<size_t>(l)][static_cast<size_t>(bi)];
            const auto& px = v[static_cast<size_t>(l) - 1][static_cast<size_t>(pbi)];
            const int64_t k = static_cast<int64_t>(box.skeleton_indices.size());
            for (int c = 0; c < nrhs; ++c)
                for (int64_t i = 0; i < k; ++i)
                    x[static_cast<size_t>(box.skeleton_indices[static_cast<size_t>(i)] + c * box.num_points)] =
                        px[static_cast<size_t>(offset + i + c * pbox->num_points)];
            if (box.f22.sparsified && !box.f22.redundant.empty()) {  // x_R22 += T22^T x_S'
                const int64_t kp = box.f22.T.rows, r = box.f22.T.cols;
                for (int c = 0; c < nrhs; ++c)
                    for (int64_t j = 0; j < r; ++j) {
                        DataType acc(0);
                        for (int64_t i = 0; i < kp; ++i)
                            acc += box.f22.T.data[static_cast<size_t>(i + j * kp)] *
                                   x[static_cast<size_t>(box.skeleton_indices[static_cast<size_t>(i)] + c * box.num_points)];
                        x[static_cast<size_t>(box.f22.redundant[static_cast<size_t>(j)] + c * box.num_points)] += acc;
                    }
            }
        }
    }
    out2.assign(static_cast<size_t>(n22) * static_cast<size_t>(nrhs), DataType(0));
    auto& leaf_level = tree->levels[static_cast<size_t>(leaf)];
    for (int64_t bi = 0; bi < leaf_level.num_boxes_local; ++bi) {
        const auto& box = leaf_level.local_boxes[static_cast<size_t>(bi)];
        const auto& x = v[static_cast<size_t>(leaf)][static_cast<size_t>(bi)];
        for (int64_t s = 0; s < box.num_points; ++s) {
            const int64_t id = box.point_indices[static_cast<size_t>(s)];
            if (id < pc) continue;
            for (int c = 0; c < nrhs; ++c)
                out2[static_cast<size_t>((id - pc) + c * n22)] = x[static_cast<size_t>(s + c * box.num_points)];
        }
    }
}

// The backward partial solve's input: x2 (F22 index order) placed at the
// leaves and taken up with P22 (x_S' += T22 x_R22) into every level's F22
// slots, the root's included; the F11 entries are left as they are.
template<typename CoordType, typename DataType>
void partial_backward_input(fmm::ParallelTree<CoordType, DataType>* tree,
                            std::vector<std::vector<fmm::SolveDataRequest<CoordType, DataType>>>& solve_data,
                            int nrhs, const std::vector<DataType>& x2) {
    partial_require_single_rank(*tree, "the backward partial solve's input");
    const int64_t pc = tree->partial_count, n22 = tree->num_points - pc;
    const int dimension = tree->dimension;
    const int leaf = tree->num_levels - 1;
    auto& leaf_level = tree->levels[static_cast<size_t>(leaf)];
    for (int64_t bi = 0; bi < leaf_level.num_boxes_local; ++bi) {
        const auto& box = leaf_level.local_boxes[static_cast<size_t>(bi)];
        auto& x = solve_data[static_cast<size_t>(leaf)][static_cast<size_t>(bi)].left_side;
        for (int64_t s = 0; s < box.num_points; ++s) {
            const int64_t id = box.point_indices[static_cast<size_t>(s)];
            if (id < pc) continue;
            for (int c = 0; c < nrhs; ++c)
                x[static_cast<size_t>(s + c * box.num_points)] = x2[static_cast<size_t>((id - pc) + c * n22)];
        }
    }
    for (int l = leaf; l >= 1; --l) {
        auto& level = tree->levels[static_cast<size_t>(l)];
        auto& parent_level = tree->levels[static_cast<size_t>(l) - 1];
        for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
            const auto& box = level.local_boxes[static_cast<size_t>(bi)];
            if (box.skeleton_indices.empty()) continue;
            auto& x = solve_data[static_cast<size_t>(l)][static_cast<size_t>(bi)].left_side;
            if (box.f22.sparsified && !box.f22.redundant.empty()) {  // x_S' += T22 x_R22
                const int64_t kp = box.f22.T.rows, r = box.f22.T.cols;
                for (int c = 0; c < nrhs; ++c)
                    for (int64_t i = 0; i < kp; ++i) {
                        DataType acc(0);
                        for (int64_t j = 0; j < r; ++j)
                            acc += box.f22.T.data[static_cast<size_t>(i + j * kp)] *
                                   x[static_cast<size_t>(box.f22.redundant[static_cast<size_t>(j)] + c * box.num_points)];
                        x[static_cast<size_t>(box.skeleton_indices[static_cast<size_t>(i)] + c * box.num_points)] += acc;
                    }
            }
            const int64_t parent = box.morton_index >> dimension;
            const int64_t pbi = parent - parent_level.local_morton_start;
            const auto* pbox = parent_level.find_local_box(parent);
            if (pbox == nullptr) throw std::logic_error("partial_backward_input: a parent is not local");
            const int64_t offset = partial_child_offset(level, box.morton_index, dimension);
            auto& px = solve_data[static_cast<size_t>(l) - 1][static_cast<size_t>(pbi)].left_side;
            const int64_t k = static_cast<int64_t>(box.skeleton_indices.size());
            for (int64_t i = 0; i < k; ++i) {
                const int64_t s = box.skeleton_indices[static_cast<size_t>(i)];
                if (box.point_indices[static_cast<size_t>(s)] < pc) continue;
                for (int c = 0; c < nrhs; ++c)
                    px[static_cast<size_t>(offset + i + c * pbox->num_points)] = x[static_cast<size_t>(s + c * box.num_points)];
            }
        }
    }
}

// =====================================================================
// CB(rows, cols), batched: rows and cols F22 indices, out column-major
// nrows x ncols (replicated).  The requested points' chains (header) go
// up the tree as blocks: the points of a box at a level that came up
// through one child share their raw vectors' support (that child's F22
// skeleton positions in the box's point list; at the leaf, their slots),
// so a block keeps them as one dense matrix over that support.  At each
// level a block splits into its p part (rows in R22) and its w part (the
// S22 rows, dense over S22 once T p is added), and every term of the
// entry formula becomes a gathered sub-block of a record times two
// blocks' parts: per pair of a row block and a column block in the same
// or one-hop boxes,
//   same box:  vp_i^T (X_RR vp_j + X_RS wv_j) + wv_i^T X_SR vp_j
//   one hop:   vp_i^T X_RN_bi(:, bj raw) v_j + wv_i^T X_NR_bj(bi raw, :) vp_j
// and at the root v_i^T CB_root v_j over the F22 slots.  The chains are
// zero on F11 slots (T has no S11 rows), so F11 rows never enter a block.
// Row blocks partition the requested rows, so a level's terms run in
// parallel over them.
// =====================================================================
template<typename DataType>
struct PartialChainBlock {
    int64_t box = -1;                ///< local box index at the level
    std::vector<int64_t> members;    ///< indices into the side's distinct points
    std::vector<int64_t> rows;       ///< raw rows: positions in the box's point list
    std::vector<DataType> v;         ///< |rows| x |members|
    std::vector<int64_t> ridx;       ///< rows in R22: indices into f22.redundant
    std::vector<DataType> vp;        ///< |ridx| x |members|
    std::vector<int64_t> wrows;      ///< the w part's rows: S' indices
    std::vector<DataType> wv;        ///< |wrows| x |members|
    std::vector<int32_t> rowpos;     ///< box slot -> position in rows (-1)
    std::vector<int32_t> wpos;       ///< S' index -> position in wrows (-1)
};

template<typename CoordType, typename DataType>
void partial_chain_parts(const fmm::BoxData<CoordType, DataType>& box, const std::vector<int32_t>& s_of,
                         const std::vector<int32_t>& r_of, PartialChainBlock<DataType>& b) {
    const int64_t n = static_cast<int64_t>(b.members.size());
    const int64_t nr = static_cast<int64_t>(b.rows.size());
    const int64_t kp = static_cast<int64_t>(box.skeleton_indices.size());
    const auto& f = box.f22;
    const bool own = f.sparsified && !f.redundant.empty();
    std::vector<int64_t> rp, sp, sidx;  // row positions of the R / S' rows, and the S' indices
    for (int64_t t = 0; t < nr; ++t) {
        const int64_t s = b.rows[static_cast<size_t>(t)];
        if (own && r_of[static_cast<size_t>(s)] >= 0) {
            rp.push_back(t);
            b.ridx.push_back(r_of[static_cast<size_t>(s)]);
        } else if (s_of[static_cast<size_t>(s)] >= 0) {
            sp.push_back(t);
            sidx.push_back(s_of[static_cast<size_t>(s)]);
        }
    }
    const int64_t nrp = static_cast<int64_t>(rp.size());
    b.vp.resize(static_cast<size_t>(nrp * n));
    for (int64_t c = 0; c < n; ++c)
        for (int64_t i = 0; i < nrp; ++i)
            b.vp[static_cast<size_t>(i + c * nrp)] = b.v[static_cast<size_t>(rp[static_cast<size_t>(i)] + c * nr)];
    if (nrp > 0) {
        // dense over S22 (S' = S11 then S22; T has no S11 rows)
        const int64_t k22 = static_cast<int64_t>(f.skeleton.size()), k11 = kp - k22;
        b.wrows.resize(static_cast<size_t>(k22));
        for (int64_t i = 0; i < k22; ++i) b.wrows[static_cast<size_t>(i)] = k11 + i;
        b.wv.assign(static_cast<size_t>(k22 * n), DataType(0));
        for (size_t q = 0; q < sp.size(); ++q) {
            const int64_t i = sidx[q] - k11;
            if (i < 0) continue;  // an S11 row: zero in an F22 chain
            for (int64_t c = 0; c < n; ++c)
                b.wv[static_cast<size_t>(i + c * k22)] = b.v[static_cast<size_t>(sp[q] + c * nr)];
        }
        std::vector<DataType> tsub(static_cast<size_t>(k22 * nrp));
        for (int64_t j = 0; j < nrp; ++j)
            std::copy_n(f.T.data.data() + k11 + b.ridx[static_cast<size_t>(j)] * kp, k22, tsub.data() + j * k22);
        partial_gemm<DataType>('N', 'N', k22, n, nrp, DataType(1), tsub.data(), k22, b.vp.data(), nrp, DataType(1),
                               b.wv.data(), k22);
    } else {
        const int64_t ns = static_cast<int64_t>(sp.size());
        b.wrows = sidx;
        b.wv.resize(static_cast<size_t>(ns * n));
        for (int64_t c = 0; c < n; ++c)
            for (int64_t i = 0; i < ns; ++i)
                b.wv[static_cast<size_t>(i + c * ns)] = b.v[static_cast<size_t>(sp[static_cast<size_t>(i)] + c * nr)];
    }
    b.rowpos.assign(static_cast<size_t>(box.num_points), -1);
    for (int64_t t = 0; t < nr; ++t) b.rowpos[static_cast<size_t>(b.rows[static_cast<size_t>(t)])] = static_cast<int32_t>(t);
    b.wpos.assign(static_cast<size_t>(kp), -1);
    for (size_t t = 0; t < b.wrows.size(); ++t) b.wpos[static_cast<size_t>(b.wrows[t])] = static_cast<int32_t>(t);
}

template<typename CoordType, typename DataType>
void partial_schur_entries(fmm::ParallelTree<CoordType, DataType>& tree, const int64_t* rows, int64_t nrows,
                           const int64_t* cols, int64_t ncols, DataType* out) {
    using Block = PartialChainBlock<DataType>;
    partial_require_single_rank(tree, "the Schur complement's entries");
    auto& pr = tree.partial_root;
    if (pr.cb_released) throw std::runtime_error("partial Schur: released (c_bpack_partial_schur_release)");
    const int64_t pc = tree.partial_count, n22 = tree.num_points - pc;
    const int dimension = tree.dimension;
    const int leaf = tree.num_levels - 1;
    std::fill(out, out + nrows * ncols, DataType(0));
    if (nrows == 0 || ncols == 0) return;
    if (pr.leaf_box.empty() && n22 > 0) {
        pr.leaf_box.assign(static_cast<size_t>(n22), -1);
        pr.leaf_slot.assign(static_cast<size_t>(n22), -1);
        const auto& level = tree.levels[static_cast<size_t>(leaf)];
        for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
            const auto& box = level.local_boxes[static_cast<size_t>(bi)];
            for (int64_t s = 0; s < box.num_points; ++s) {
                const int64_t id = box.point_indices[static_cast<size_t>(s)];
                if (id < pc) continue;
                pr.leaf_box[static_cast<size_t>(id - pc)] = bi;
                pr.leaf_slot[static_cast<size_t>(id - pc)] = s;
            }
        }
    }
    // the two sides: distinct points, their requested positions, the blocks of the current level
    struct Side {
        std::vector<int64_t> points;
        std::vector<std::vector<int64_t>> positions;
        std::vector<Block> blocks;
    };
    auto make_side = [&](const int64_t* idx, int64_t n) {
        Side side;
        std::unordered_map<int64_t, int64_t> of;
        for (int64_t q = 0; q < n; ++q) {
            const int64_t f = idx[q];
            if (f < 0 || f >= n22) throw std::runtime_error("partial Schur: an index out of range");
            auto it = of.find(f);
            if (it == of.end()) {
                it = of.emplace(f, static_cast<int64_t>(side.points.size())).first;
                side.points.push_back(f);
                side.positions.emplace_back();
            }
            side.positions[static_cast<size_t>(it->second)].push_back(q);
        }
        std::unordered_map<int64_t, size_t> block_of;  // leaf box -> block
        for (size_t m = 0; m < side.points.size(); ++m) {
            const int64_t f = side.points[m];
            const int64_t b = pr.leaf_box[static_cast<size_t>(f)];
            if (b < 0) throw std::runtime_error("partial Schur: a point is not in a local leaf");
            auto it = block_of.find(b);
            if (it == block_of.end()) {
                it = block_of.emplace(b, side.blocks.size()).first;
                side.blocks.emplace_back();
                side.blocks.back().box = b;
            }
            auto& blk = side.blocks[it->second];
            blk.members.push_back(static_cast<int64_t>(m));
            blk.rows.push_back(pr.leaf_slot[static_cast<size_t>(f)]);
        }
        for (auto& blk : side.blocks) {
            const int64_t k = static_cast<int64_t>(blk.members.size());
            blk.v.assign(static_cast<size_t>(k * k), DataType(0));
            for (int64_t t = 0; t < k; ++t) blk.v[static_cast<size_t>(t + t * k)] = DataType(1);
        }
        return side;
    };
    Side R = make_side(rows, nrows), C = make_side(cols, ncols);
    // G (|members_i| x |members_j|) into out at the requested positions
    auto scatter = [&](const Block& bi, const Block& bj, const std::vector<DataType>& G) {
        const int64_t ni = static_cast<int64_t>(bi.members.size());
        for (size_t c = 0; c < bj.members.size(); ++c)
            for (int64_t jj : C.positions[static_cast<size_t>(bj.members[c])])
                for (int64_t r = 0; r < ni; ++r) {
                    const DataType g = G[static_cast<size_t>(r + static_cast<int64_t>(c) * ni)];
                    for (int64_t ii : R.positions[static_cast<size_t>(bi.members[static_cast<size_t>(r)])])
                        out[ii + jj * nrows] += g;
                }
    };
    for (int l = leaf; l >= 0; --l) {
        auto& level = tree.levels[static_cast<size_t>(l)];
        if (l == 0) {
            // the root: v_i^T CB_root v_j over the root's F22 slots
            const auto& root = level.local_boxes[0];
            std::vector<int64_t> pos22(static_cast<size_t>(root.num_points), -1);
            for (int64_t a = 0; a < pr.n22; ++a) pos22[static_cast<size_t>(pr.slots22[static_cast<size_t>(a)])] = a;
            auto f22_part = [&](const Block& b, std::vector<int64_t>& a, std::vector<DataType>& v) {
                const int64_t nr = static_cast<int64_t>(b.rows.size()), n = static_cast<int64_t>(b.members.size());
                std::vector<int64_t> t;
                for (int64_t q = 0; q < nr; ++q) {
                    const int64_t p = pos22[static_cast<size_t>(b.rows[static_cast<size_t>(q)])];
                    if (p >= 0) { a.push_back(p); t.push_back(q); }
                }
                const int64_t na = static_cast<int64_t>(a.size());
                v.resize(static_cast<size_t>(na * n));
                for (int64_t c = 0; c < n; ++c)
                    for (int64_t i = 0; i < na; ++i) v[static_cast<size_t>(i + c * na)] = b.v[static_cast<size_t>(t[static_cast<size_t>(i)] + c * nr)];
            };
            std::vector<std::vector<int64_t>> ca(C.blocks.size());
            std::vector<std::vector<DataType>> cv(C.blocks.size());
            for (size_t j = 0; j < C.blocks.size(); ++j) f22_part(C.blocks[j], ca[j], cv[j]);
#pragma omp parallel for schedule(dynamic, 1)
            for (size_t i = 0; i < R.blocks.size(); ++i) {
                std::vector<int64_t> ai;
                std::vector<DataType> vi, M, tmp, G;
                f22_part(R.blocks[i], ai, vi);
                const int64_t na = static_cast<int64_t>(ai.size()), ni = static_cast<int64_t>(R.blocks[i].members.size());
                if (na == 0) continue;
                for (size_t j = 0; j < C.blocks.size(); ++j) {
                    const int64_t nb = static_cast<int64_t>(ca[j].size()), nj = static_cast<int64_t>(C.blocks[j].members.size());
                    if (nb == 0) continue;
                    M.resize(static_cast<size_t>(na * nb));
                    for (int64_t c = 0; c < nb; ++c)
                        for (int64_t r = 0; r < na; ++r)
                            M[static_cast<size_t>(r + c * na)] =
                                pr.cb.data[static_cast<size_t>(ai[static_cast<size_t>(r)] + ca[j][static_cast<size_t>(c)] * pr.cb.lda)];
                    tmp.resize(static_cast<size_t>(na * nj));
                    partial_gemm<DataType>('N', 'N', na, nj, nb, DataType(1), M.data(), na, cv[j].data(), nb, DataType(0), tmp.data(), na);
                    G.assign(static_cast<size_t>(ni * nj), DataType(0));
                    partial_gemm<DataType>('T', 'N', ni, nj, na, DataType(1), vi.data(), na, tmp.data(), na, DataType(0), G.data(), ni);
                    scatter(R.blocks[i], C.blocks[j], G);
                }
            }
            break;
        }
        // the slot maps of the boxes the blocks touch
        std::unordered_map<int64_t, std::pair<std::vector<int32_t>, std::vector<int32_t>>> maps;
        for (const Side* side : {&R, &C})
            for (const auto& b : side->blocks) maps.emplace(b.box, std::pair<std::vector<int32_t>, std::vector<int32_t>>{});
        std::vector<int64_t> touched;
        for (const auto& kv : maps) touched.push_back(kv.first);
#pragma omp parallel for schedule(dynamic, 1)
        for (size_t q = 0; q < touched.size(); ++q) {
            const auto& box = level.local_boxes[static_cast<size_t>(touched[q])];
            auto& m = maps.at(touched[q]);
            m.first.assign(static_cast<size_t>(box.num_points), -1);
            m.second.assign(static_cast<size_t>(box.num_points), -1);
            for (size_t i = 0; i < box.skeleton_indices.size(); ++i) m.first[static_cast<size_t>(box.skeleton_indices[i])] = static_cast<int32_t>(i);
            if (box.f22.sparsified)
                for (size_t i = 0; i < box.f22.redundant.size(); ++i) m.second[static_cast<size_t>(box.f22.redundant[i])] = static_cast<int32_t>(i);
        }
        for (Side* side : {&R, &C}) {
            auto& blocks = side->blocks;
#pragma omp parallel for schedule(dynamic, 1)
            for (size_t i = 0; i < blocks.size(); ++i) {
                const auto& m = maps.at(blocks[i].box);
                partial_chain_parts(level.local_boxes[static_cast<size_t>(blocks[i].box)], m.first, m.second, blocks[i]);
            }
        }
        // the level's terms
        if (l >= 2) {
            std::unordered_map<int64_t, std::vector<size_t>> col_by_box;
            for (size_t j = 0; j < C.blocks.size(); ++j) col_by_box[C.blocks[j].box].push_back(j);
#pragma omp parallel for schedule(dynamic, 1)
            for (size_t i = 0; i < R.blocks.size(); ++i) {
                const Block& bi = R.blocks[i];
                const auto& boxI = level.local_boxes[static_cast<size_t>(bi.box)];
                const auto& fi = boxI.f22;
                const bool own_i = fi.sparsified && !fi.redundant.empty() && !bi.ridx.empty();
                const int64_t ni = static_cast<int64_t>(bi.members.size()), ri = static_cast<int64_t>(bi.ridx.size());
                std::vector<DataType> X, tmp, G, V;
                auto pair_terms = [&](int64_t nb_local, int64_t idx_in_i) {
                    const auto it = col_by_box.find(nb_local);
                    if (it == col_by_box.end()) return;
                    const auto& boxJ = level.local_boxes[static_cast<size_t>(nb_local)];
                    const auto& fj = boxJ.f22;
                    for (size_t jb : it->second) {
                        const Block& bj = C.blocks[jb];
                        const int64_t nj = static_cast<int64_t>(bj.members.size()), rj = static_cast<int64_t>(bj.ridx.size());
                        const bool own_j = fj.sparsified && !fj.redundant.empty() && rj > 0;
                        G.assign(static_cast<size_t>(ni * nj), DataType(0));
                        bool any = false;
                        if (nb_local == bi.box) {
                            const int64_t r = static_cast<int64_t>(fi.redundant.size()), kp = static_cast<int64_t>(boxI.skeleton_indices.size());
                            if (own_i) {
                                // tmp = X_RR(ri, rj) vp_j + X_RS(ri, wj) wv_j;  G += vp_i^T tmp
                                tmp.assign(static_cast<size_t>(ri * nj), DataType(0));
                                if (rj > 0) {
                                    X.resize(static_cast<size_t>(ri * rj));
                                    for (int64_t c = 0; c < rj; ++c)
                                        for (int64_t a = 0; a < ri; ++a)
                                            X[static_cast<size_t>(a + c * ri)] = fi.X_RR.data[static_cast<size_t>(
                                                bi.ridx[static_cast<size_t>(a)] + bj.ridx[static_cast<size_t>(c)] * r)];
                                    partial_gemm<DataType>('N', 'N', ri, nj, rj, DataType(1), X.data(), ri, bj.vp.data(), rj, DataType(1), tmp.data(), ri);
                                }
                                const int64_t wj = static_cast<int64_t>(bj.wrows.size());
                                if (wj > 0) {
                                    X.resize(static_cast<size_t>(ri * wj));
                                    for (int64_t c = 0; c < wj; ++c)
                                        for (int64_t a = 0; a < ri; ++a)
                                            X[static_cast<size_t>(a + c * ri)] = fi.X_RS.data[static_cast<size_t>(
                                                bi.ridx[static_cast<size_t>(a)] + bj.wrows[static_cast<size_t>(c)] * r)];
                                    partial_gemm<DataType>('N', 'N', ri, nj, wj, DataType(1), X.data(), ri, bj.wv.data(), wj, DataType(1), tmp.data(), ri);
                                }
                                partial_gemm<DataType>('T', 'N', ni, nj, ri, DataType(1), bi.vp.data(), ri, tmp.data(), ri, DataType(1), G.data(), ni);
                                any = true;
                            }
                            const int64_t wi = static_cast<int64_t>(bi.wrows.size());
                            if (fi.sparsified && !fi.redundant.empty() && rj > 0 && wi > 0) {
                                // G += wv_i^T X_SR(wi, rj) vp_j
                                X.resize(static_cast<size_t>(wi * rj));
                                for (int64_t c = 0; c < rj; ++c)
                                    for (int64_t a = 0; a < wi; ++a)
                                        X[static_cast<size_t>(a + c * wi)] = fi.X_SR.data[static_cast<size_t>(
                                            bi.wrows[static_cast<size_t>(a)] + bj.ridx[static_cast<size_t>(c)] * kp)];
                                tmp.resize(static_cast<size_t>(wi * nj));
                                partial_gemm<DataType>('N', 'N', wi, nj, rj, DataType(1), X.data(), wi, bj.vp.data(), rj, DataType(0), tmp.data(), wi);
                                partial_gemm<DataType>('T', 'N', ni, nj, wi, DataType(1), bi.wv.data(), wi, tmp.data(), wi, DataType(1), G.data(), ni);
                                any = true;
                            }
                        } else {
                            if (own_i && fi.neighbor_f22_counts[static_cast<size_t>(idx_in_i)] > 0) {
                                // G += vp_i^T X_RN_bi(ri, bj raw slots in bj's rows) v_j
                                const int64_t r = static_cast<int64_t>(fi.redundant.size());
                                const int64_t m = fi.neighbor_f22_counts[static_cast<size_t>(idx_in_i)];
                                int64_t col0 = 0;
                                for (int64_t q = 0; q < idx_in_i; ++q) col0 += fi.neighbor_f22_counts[static_cast<size_t>(q)];
                                const auto& slots = fi.neighbor_f22_slots[static_cast<size_t>(idx_in_i)];
                                std::vector<int64_t> cpos, tpos;
                                for (int64_t c = 0; c < m; ++c) {
                                    const int32_t t = bj.rowpos[static_cast<size_t>(slots[static_cast<size_t>(c)])];
                                    if (t >= 0) { cpos.push_back(col0 + c); tpos.push_back(t); }
                                }
                                const int64_t nc = static_cast<int64_t>(cpos.size());
                                if (nc > 0) {
                                    X.resize(static_cast<size_t>(ri * nc));
                                    for (int64_t c = 0; c < nc; ++c)
                                        for (int64_t a = 0; a < ri; ++a)
                                            X[static_cast<size_t>(a + c * ri)] = fi.X_RN.data[static_cast<size_t>(
                                                bi.ridx[static_cast<size_t>(a)] + cpos[static_cast<size_t>(c)] * r)];
                                    const int64_t nrj = static_cast<int64_t>(bj.rows.size());
                                    V.resize(static_cast<size_t>(nc * nj));
                                    for (int64_t c = 0; c < nj; ++c)
                                        for (int64_t q = 0; q < nc; ++q)
                                            V[static_cast<size_t>(q + c * nc)] = bj.v[static_cast<size_t>(tpos[static_cast<size_t>(q)] + c * nrj)];
                                    tmp.resize(static_cast<size_t>(ri * nj));
                                    partial_gemm<DataType>('N', 'N', ri, nj, nc, DataType(1), X.data(), ri, V.data(), nc, DataType(0), tmp.data(), ri);
                                    partial_gemm<DataType>('T', 'N', ni, nj, ri, DataType(1), bi.vp.data(), ri, tmp.data(), ri, DataType(1), G.data(), ni);
                                    any = true;
                                }
                            }
                            if (own_j && !bi.wrows.empty()) {
                                // G += wv_i^T X_NR_bj(bi raw slots in bi's w rows, rj) vp_j
                                const auto itb = std::find(boxJ.one_hop.begin(), boxJ.one_hop.end(), boxI.morton_index);
                                if (itb != boxJ.one_hop.end()) {
                                    const size_t idi = static_cast<size_t>(itb - boxJ.one_hop.begin());
                                    const int64_t m = fj.neighbor_f22_counts[idi], total = fj.X_NR.rows;
                                    int64_t row0 = 0;
                                    for (size_t q = 0; q < idi; ++q) row0 += fj.neighbor_f22_counts[q];
                                    const auto& slots = fj.neighbor_f22_slots[idi];
                                    const auto& s_of = maps.at(bi.box).first;
                                    std::vector<int64_t> crow, wposl;
                                    for (int64_t c = 0; c < m; ++c) {
                                        const int32_t sp = s_of[static_cast<size_t>(slots[static_cast<size_t>(c)])];
                                        if (sp < 0) continue;  // in R_bi: the X_RN term
                                        const int32_t t = bi.wpos[static_cast<size_t>(sp)];
                                        if (t >= 0) { crow.push_back(row0 + c); wposl.push_back(t); }
                                    }
                                    const int64_t nc = static_cast<int64_t>(crow.size());
                                    if (nc > 0) {
                                        X.resize(static_cast<size_t>(nc * rj));
                                        for (int64_t c = 0; c < rj; ++c)
                                            for (int64_t q = 0; q < nc; ++q)
                                                X[static_cast<size_t>(q + c * nc)] = fj.X_NR.data[static_cast<size_t>(
                                                    crow[static_cast<size_t>(q)] + bj.ridx[static_cast<size_t>(c)] * total)];
                                        const int64_t nwi = static_cast<int64_t>(bi.wrows.size());
                                        V.resize(static_cast<size_t>(nc * ni));
                                        for (int64_t c = 0; c < ni; ++c)
                                            for (int64_t q = 0; q < nc; ++q)
                                                V[static_cast<size_t>(q + c * nc)] = bi.wv[static_cast<size_t>(wposl[static_cast<size_t>(q)] + c * nwi)];
                                        tmp.resize(static_cast<size_t>(nc * nj));
                                        partial_gemm<DataType>('N', 'N', nc, nj, rj, DataType(1), X.data(), nc, bj.vp.data(), rj, DataType(0), tmp.data(), nc);
                                        partial_gemm<DataType>('T', 'N', ni, nj, nc, DataType(1), V.data(), nc, tmp.data(), nc, DataType(1), G.data(), ni);
                                        any = true;
                                    }
                                }
                            }
                        }
                        if (any) scatter(bi, bj, G);
                    }
                };
                pair_terms(bi.box, -1);
                for (size_t q = 0; q < boxI.one_hop.size(); ++q) {
                    const int64_t nb = boxI.one_hop[q];
                    if (!level.is_box_on_process(nb)) throw std::runtime_error("partial Schur: a neighbour is not local");
                    pair_terms(nb - level.local_morton_start, static_cast<int64_t>(q));
                }
            }
        }
        // the parent blocks: each block's w part at its box's offset; the dense
        // parts of one box's blocks share their rows and merge
        auto& parent_level = tree.levels[static_cast<size_t>(l) - 1];
        for (Side* side : {&R, &C}) {
            std::vector<Block> up;
            std::unordered_map<int64_t, size_t> dense_of;  // box at this level -> merged block at the parent
            for (auto& b : side->blocks) {
                const auto& box = level.local_boxes[static_cast<size_t>(b.box)];
                if (b.wrows.empty()) continue;  // nothing left of these chains above
                const int64_t morton = box.morton_index;
                const int64_t parent = morton >> dimension;
                if (!parent_level.is_box_on_process(parent)) throw std::runtime_error("partial Schur: a parent is not local");
                const int64_t off = partial_child_offset(level, morton, dimension);
                const bool dense = !b.ridx.empty();
                const int64_t nw = static_cast<int64_t>(b.wrows.size());
                if (dense) {
                    auto it = dense_of.find(b.box);
                    if (it != dense_of.end()) {
                        auto& u = up[it->second];
                        u.members.insert(u.members.end(), b.members.begin(), b.members.end());
                        u.v.insert(u.v.end(), b.wv.begin(), b.wv.end());
                        continue;
                    }
                    dense_of.emplace(b.box, up.size());
                }
                Block u;
                u.box = parent - parent_level.local_morton_start;
                u.members = std::move(b.members);
                u.rows.resize(static_cast<size_t>(nw));
                for (int64_t t = 0; t < nw; ++t) u.rows[static_cast<size_t>(t)] = off + b.wrows[static_cast<size_t>(t)];
                u.v = std::move(b.wv);
                up.push_back(std::move(u));
            }
            side->blocks = std::move(up);
        }
    }
}

// y2 = CB x2 = P22^T CB' P22 x2 (header): an up-pass with P (w = x_S' +
// T x_R, the raw vectors kept), the sparsified blocks and the root's block
// applied level by level, a down-pass with P^T (a_S' = g_S' + a_parent,
// a_R = g_R + T^T a_S').  x2, y2 column-major n22 x ncol in F22 index order
// (replicated).  Every pass runs over a level's boxes in parallel: a box
// writes only its own vectors and its range of its parent's.
template<typename CoordType, typename DataType>
void partial_schur_multiply(fmm::ParallelTree<CoordType, DataType>& tree, const DataType* x2, DataType* y2,
                            int64_t ncol) {
    partial_require_single_rank(tree, "the Schur complement's product");
    const auto& pr = tree.partial_root;
    if (pr.cb_released) throw std::runtime_error("partial Schur: released (c_bpack_partial_schur_release)");
    const int64_t pc = tree.partial_count, n22 = tree.num_points - pc;
    const int dimension = tree.dimension;
    const int leaf = tree.num_levels - 1;
    if (ncol <= 0) return;
    auto own = [](const fmm::BoxData<CoordType, DataType>& box) {
        return box.f22.sparsified && !box.f22.redundant.empty();
    };
    // per level and box: raw (the box's points x ncol), p = raw[R], w = raw[S'] + T p,
    // and the accumulator h over the box's points
    using Vecs = std::vector<std::vector<DataType>>;
    std::vector<Vecs> raw(static_cast<size_t>(tree.num_levels)), pv(static_cast<size_t>(tree.num_levels)),
        h(static_cast<size_t>(tree.num_levels));
    std::vector<std::vector<int64_t>> offset(static_cast<size_t>(tree.num_levels));
    for (int l = 0; l <= leaf; ++l) {
        auto& level = tree.levels[static_cast<size_t>(l)];
        const int64_t nb = level.num_boxes_local;
        raw[static_cast<size_t>(l)].resize(static_cast<size_t>(nb));
        pv[static_cast<size_t>(l)].resize(static_cast<size_t>(nb));
        h[static_cast<size_t>(l)].resize(static_cast<size_t>(nb));
        offset[static_cast<size_t>(l)].assign(static_cast<size_t>(nb), 0);
#pragma omp parallel for schedule(dynamic, 16)
        for (int64_t bi = 0; bi < nb; ++bi) {
            const auto& box = level.local_boxes[static_cast<size_t>(bi)];
            const size_t n = static_cast<size_t>(box.num_points) * static_cast<size_t>(ncol);
            raw[static_cast<size_t>(l)][static_cast<size_t>(bi)].assign(n, DataType(0));
            h[static_cast<size_t>(l)][static_cast<size_t>(bi)].assign(n, DataType(0));
            if (l > 0 && !box.skeleton_indices.empty())
                offset[static_cast<size_t>(l)][static_cast<size_t>(bi)] = partial_child_offset(level, box.morton_index, dimension);
        }
    }
    {
        auto& level = tree.levels[static_cast<size_t>(leaf)];
#pragma omp parallel for schedule(dynamic, 16)
        for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
            const auto& box = level.local_boxes[static_cast<size_t>(bi)];
            auto& v = raw[static_cast<size_t>(leaf)][static_cast<size_t>(bi)];
            for (int64_t s = 0; s < box.num_points; ++s) {
                const int64_t id = box.point_indices[static_cast<size_t>(s)];
                if (id < pc) continue;
                for (int64_t c = 0; c < ncol; ++c)
                    v[static_cast<size_t>(s + c * box.num_points)] = x2[static_cast<size_t>((id - pc) + c * n22)];
            }
        }
    }
    // up: p = raw[R], w = raw[S'] + T p into the parent's raw vector at the box's offset
    for (int l = leaf; l >= 1; --l) {
        auto& level = tree.levels[static_cast<size_t>(l)];
        auto& parent_level = tree.levels[static_cast<size_t>(l) - 1];
        std::exception_ptr failure;
#pragma omp parallel for schedule(dynamic, 1)
        for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
            const auto& box = level.local_boxes[static_cast<size_t>(bi)];
            if (box.skeleton_indices.empty()) continue;
            const auto& x = raw[static_cast<size_t>(l)][static_cast<size_t>(bi)];
            const int64_t n = box.num_points, kp = static_cast<int64_t>(box.skeleton_indices.size());
            std::vector<DataType> w(static_cast<size_t>(kp * ncol));
            for (int64_t c = 0; c < ncol; ++c)
                for (int64_t i = 0; i < kp; ++i)
                    w[static_cast<size_t>(i + c * kp)] = x[static_cast<size_t>(box.skeleton_indices[static_cast<size_t>(i)] + c * n)];
            if (own(box)) {
                const int64_t r = static_cast<int64_t>(box.f22.redundant.size());
                auto& p = pv[static_cast<size_t>(l)][static_cast<size_t>(bi)];
                p.resize(static_cast<size_t>(r * ncol));
                for (int64_t c = 0; c < ncol; ++c)
                    for (int64_t j = 0; j < r; ++j)
                        p[static_cast<size_t>(j + c * r)] = x[static_cast<size_t>(box.f22.redundant[static_cast<size_t>(j)] + c * n)];
                partial_gemm<DataType>('N', 'N', kp, ncol, r, DataType(1), box.f22.T.data.data(), kp, p.data(), r, DataType(1),
                                       w.data(), kp);
            }
            const int64_t parent = box.morton_index >> dimension;
            const auto* pbox = parent_level.find_local_box(parent);
            if (pbox == nullptr) {
#pragma omp critical(partial_schur_multiply_failure)
                if (!failure) failure = std::make_exception_ptr(std::logic_error("partial_schur_multiply: a parent is not local"));
                continue;
            }
            auto& px = raw[static_cast<size_t>(l) - 1][static_cast<size_t>(parent - parent_level.local_morton_start)];
            const int64_t off = offset[static_cast<size_t>(l)][static_cast<size_t>(bi)];
            for (int64_t c = 0; c < ncol; ++c)
                std::copy_n(w.data() + c * kp, kp, px.data() + off + c * pbox->num_points);
        }
        if (failure) std::rethrow_exception(failure);
    }
    // the root: g = CB_root v on the root's F22 slots, in row chunks
    {
        const auto& box = tree.levels[0].local_boxes[0];
        const auto& x = raw[0][0];
        auto& out = h[0][0];
        const int64_t n = box.num_points, m = pr.n22;
        if (m > 0) {
            std::vector<DataType> v(static_cast<size_t>(m * ncol)), g(static_cast<size_t>(m * ncol));
            for (int64_t c = 0; c < ncol; ++c)
                for (int64_t a = 0; a < m; ++a)
                    v[static_cast<size_t>(a + c * m)] = x[static_cast<size_t>(pr.slots22[static_cast<size_t>(a)] + c * n)];
            constexpr int64_t chunk = 256;
            const int64_t nchunks = (m + chunk - 1) / chunk;
#pragma omp parallel for schedule(static)
            for (int64_t q = 0; q < nchunks; ++q) {
                const int64_t a0 = q * chunk, rows = std::min(chunk, m - a0);
                partial_gemm<DataType>('N', 'N', rows, ncol, m, DataType(1), pr.cb.data.data() + a0, pr.cb.lda, v.data(), m,
                                       DataType(0), g.data() + a0, m);
            }
            for (int64_t c = 0; c < ncol; ++c)
                for (int64_t a = 0; a < m; ++a)
                    out[static_cast<size_t>(pr.slots22[static_cast<size_t>(a)] + c * n)] = g[static_cast<size_t>(a + c * m)];
        }
    }
    // the records, level by level: a box gathers its own terms and its sparsified neighbours' X_NR rows over its S'
    for (int l = leaf; l >= 2; --l) {
        auto& level = tree.levels[static_cast<size_t>(l)];
#pragma omp parallel for schedule(dynamic, 1)
        for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
            const auto& box = level.local_boxes[static_cast<size_t>(bi)];
            const auto& f = box.f22;
            if (box.skeleton_indices.empty()) continue;
            const bool mine = own(box);
            const int64_t n = box.num_points, kp = static_cast<int64_t>(box.skeleton_indices.size());
            const int64_t r = mine ? static_cast<int64_t>(f.redundant.size()) : 0;
            const auto& x = raw[static_cast<size_t>(l)][static_cast<size_t>(bi)];
            std::vector<DataType> gR(static_cast<size_t>(r * ncol), DataType(0)), gS(static_cast<size_t>(kp * ncol), DataType(0));
            if (mine) {
                // the box's own blocks: gR = X_RR p + X_RS w, gS = X_SR p (w = raw[S'] + T p again, kept local)
                const auto& p = pv[static_cast<size_t>(l)][static_cast<size_t>(bi)];
                std::vector<DataType> w(static_cast<size_t>(kp * ncol));
                for (int64_t c = 0; c < ncol; ++c)
                    for (int64_t i = 0; i < kp; ++i)
                        w[static_cast<size_t>(i + c * kp)] = x[static_cast<size_t>(box.skeleton_indices[static_cast<size_t>(i)] + c * n)];
                partial_gemm<DataType>('N', 'N', kp, ncol, r, DataType(1), f.T.data.data(), kp, p.data(), r, DataType(1), w.data(), kp);
                partial_gemm<DataType>('N', 'N', r, ncol, r, DataType(1), f.X_RR.data.data(), r, p.data(), r, DataType(1), gR.data(), r);
                partial_gemm<DataType>('N', 'N', r, ncol, kp, DataType(1), f.X_RS.data.data(), r, w.data(), kp, DataType(1), gR.data(), r);
                partial_gemm<DataType>('N', 'N', kp, ncol, r, DataType(1), f.X_SR.data.data(), kp, p.data(), r, DataType(1), gS.data(), kp);
                // the neighbours' raw vectors against X_RN
                const int64_t total = f.X_RN.cols;
                if (total > 0) {
                    std::vector<DataType> rn(static_cast<size_t>(total * ncol));
                    int64_t col0 = 0;
                    for (size_t idx = 0; idx < box.one_hop.size(); ++idx) {
                        const int64_t m = f.neighbor_f22_counts[idx];
                        if (m == 0) continue;
                        const int64_t nl = box.one_hop[idx] - level.local_morton_start;
                        const int64_t nn = level.local_boxes[static_cast<size_t>(nl)].num_points;
                        const auto& xn = raw[static_cast<size_t>(l)][static_cast<size_t>(nl)];
                        const auto& slots = f.neighbor_f22_slots[idx];
                        for (int64_t c = 0; c < ncol; ++c)
                            for (int64_t i = 0; i < m; ++i)
                                rn[static_cast<size_t>(col0 + i + c * total)] = xn[static_cast<size_t>(slots[static_cast<size_t>(i)] + c * nn)];
                        col0 += m;
                    }
                    partial_gemm<DataType>('N', 'N', r, ncol, total, DataType(1), f.X_RN.data.data(), r, rn.data(), total, DataType(1),
                                           gR.data(), r);
                }
            }
            // the sparsified neighbours' p against their X_NR rows that are this box's S' slots
            std::vector<int32_t> s_of;
            std::vector<DataType> t;
            for (size_t idx = 0; idx < box.one_hop.size(); ++idx) {
                const int64_t nl = box.one_hop[idx] - level.local_morton_start;
                const auto& nb = level.local_boxes[static_cast<size_t>(nl)];
                if (!own(nb)) continue;
                const auto& fn = nb.f22;
                const auto it = std::find(nb.one_hop.begin(), nb.one_hop.end(), box.morton_index);
                if (it == nb.one_hop.end()) continue;
                const size_t idn = static_cast<size_t>(it - nb.one_hop.begin());
                const int64_t m = fn.neighbor_f22_counts[idn], rn_ = static_cast<int64_t>(fn.redundant.size());
                if (m == 0) continue;
                int64_t row0 = 0;
                for (size_t q = 0; q < idn; ++q) row0 += fn.neighbor_f22_counts[q];
                const auto& pn = pv[static_cast<size_t>(l)][static_cast<size_t>(nl)];
                t.resize(static_cast<size_t>(m * ncol));
                partial_gemm<DataType>('N', 'N', m, ncol, rn_, DataType(1), fn.X_NR.data.data() + row0, fn.X_NR.rows, pn.data(), rn_,
                                       DataType(0), t.data(), m);
                if (s_of.empty()) {
                    s_of.assign(static_cast<size_t>(n), -1);
                    for (int64_t i = 0; i < kp; ++i) s_of[static_cast<size_t>(box.skeleton_indices[static_cast<size_t>(i)])] = static_cast<int32_t>(i);
                }
                const auto& slots = fn.neighbor_f22_slots[idn];
                for (int64_t i = 0; i < m; ++i) {
                    const int32_t si = s_of[static_cast<size_t>(slots[static_cast<size_t>(i)])];
                    if (si < 0) continue;  // this box's R rows: in X_RN above
                    for (int64_t c = 0; c < ncol; ++c) gS[static_cast<size_t>(si + c * kp)] += t[static_cast<size_t>(i + c * m)];
                }
            }
            auto& out = h[static_cast<size_t>(l)][static_cast<size_t>(bi)];
            for (int64_t c = 0; c < ncol; ++c) {
                for (int64_t i = 0; i < kp; ++i) out[static_cast<size_t>(box.skeleton_indices[static_cast<size_t>(i)] + c * n)] += gS[static_cast<size_t>(i + c * kp)];
                for (int64_t j = 0; j < r; ++j) out[static_cast<size_t>(f.redundant[static_cast<size_t>(j)] + c * n)] += gR[static_cast<size_t>(j + c * r)];
            }
        }
    }
    // down: a_S' = g_S' + a_parent, a_R = g_R + T^T a_S'
    for (int l = 1; l <= leaf; ++l) {
        auto& level = tree.levels[static_cast<size_t>(l)];
        auto& parent_level = tree.levels[static_cast<size_t>(l) - 1];
#pragma omp parallel for schedule(dynamic, 1)
        for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
            const auto& box = level.local_boxes[static_cast<size_t>(bi)];
            if (box.skeleton_indices.empty()) continue;
            auto& a = h[static_cast<size_t>(l)][static_cast<size_t>(bi)];
            const int64_t n = box.num_points, kp = static_cast<int64_t>(box.skeleton_indices.size());
            const int64_t parent = box.morton_index >> dimension;
            const int64_t pl = parent - parent_level.local_morton_start;
            const int64_t pn = parent_level.local_boxes[static_cast<size_t>(pl)].num_points;
            const auto& pa = h[static_cast<size_t>(l) - 1][static_cast<size_t>(pl)];
            const int64_t off = offset[static_cast<size_t>(l)][static_cast<size_t>(bi)];
            std::vector<DataType> aS(static_cast<size_t>(kp * ncol));
            for (int64_t c = 0; c < ncol; ++c)
                for (int64_t i = 0; i < kp; ++i) {
                    DataType& ai = a[static_cast<size_t>(box.skeleton_indices[static_cast<size_t>(i)] + c * n)];
                    ai += pa[static_cast<size_t>(off + i + c * pn)];
                    aS[static_cast<size_t>(i + c * kp)] = ai;
                }
            if (own(box)) {
                const int64_t r = static_cast<int64_t>(box.f22.redundant.size());
                std::vector<DataType> aR(static_cast<size_t>(r * ncol));
                partial_gemm<DataType>('T', 'N', r, ncol, kp, DataType(1), box.f22.T.data.data(), kp, aS.data(), kp, DataType(0),
                                       aR.data(), r);
                for (int64_t c = 0; c < ncol; ++c)
                    for (int64_t j = 0; j < r; ++j)
                        a[static_cast<size_t>(box.f22.redundant[static_cast<size_t>(j)] + c * n)] += aR[static_cast<size_t>(j + c * r)];
            }
        }
    }
    {
        auto& level = tree.levels[static_cast<size_t>(leaf)];
#pragma omp parallel for schedule(dynamic, 16)
        for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
            const auto& box = level.local_boxes[static_cast<size_t>(bi)];
            const auto& a = h[static_cast<size_t>(leaf)][static_cast<size_t>(bi)];
            for (int64_t s = 0; s < box.num_points; ++s) {
                const int64_t id = box.point_indices[static_cast<size_t>(s)];
                if (id < pc) continue;
                for (int64_t c = 0; c < ncol; ++c)
                    y2[static_cast<size_t>((id - pc) + c * n22)] = a[static_cast<size_t>(s + c * box.num_points)];
            }
        }
    }
}

template<typename CoordType, typename DataType>
void partial_schur_release(fmm::ParallelTree<CoordType, DataType>& tree) {
    auto& pr = tree.partial_root;
    std::vector<DataType>().swap(pr.cb.data);
    pr.cb.rows = pr.cb.cols = pr.cb.lda = 0;
    pr.cb_released = true;
    for (auto& level : tree.levels) {
        for (auto& box : level.local_boxes) {
            auto& f = box.f22;
            if (!f.sparsified) continue;
            for (auto* m : {&f.X_RR, &f.X_RS, &f.X_SR, &f.X_RN, &f.X_NR}) {
                std::vector<DataType>().swap(m->data);
                m->rows = m->cols = m->lda = 0;
            }
            f.released = true;
        }
    }
}

}  // namespace color_unstructured
}  // namespace butterfly
