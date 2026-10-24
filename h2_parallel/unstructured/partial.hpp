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
// Many ranks: a level's neighbours on other ranks are read through their
// copies (ghosts on a CA level, assisting boxes on a Color one), the F22
// phase publishes each sparsified box's kept slots to the ranks holding a
// copy, and the CB entries and products are collective over the level
// communicators, their results replicated.

#include "core/tree.hpp"
#include "core/blas_declare.hpp"
#include "core/factorization.hpp"
#include "core/id_decomposition.hpp"
#include "core/morton.hpp"
#include "core/solver.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <functional>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <set>
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
// The F22 phase's new live sets to the other ranks: each sparsified box's
// kept positions within its skeleton before the phase go to the owners of
// the boxes within three hops (its ghost copies are one-hop neighbours'
// ranks, its assisting copies up to two hops for the IDs and three for
// the parent transition, which reads a parent neighbour's children); a
// receiver brings its ghost and assisting copies to the new live sets
// (the old ones kept, for the lazy fill's rows) and adds the boxes to
// kept_of (morton -> the old live set's size and the kept positions).
// Collective over level_comm (the level's active ranks).
// =====================================================================
template<typename CoordType, typename DataType>
void partial_f22_publish(fmm::ParallelTree<CoordType, DataType>* tree, fmm::TreeLevel<CoordType, DataType>& level,
                         MPI_Comm level_comm, std::unordered_map<int64_t, std::pair<int64_t, std::vector<int64_t>>>& kept_of) {
    int me = 0, lsize = 0;
    MPI_Comm_rank(tree->comm, &me);
    MPI_Comm_size(level_comm, &lsize);
    MPI_Group tree_group, level_group;
    MPI_Comm_group(tree->comm, &tree_group);
    MPI_Comm_group(level_comm, &level_group);
    std::unordered_map<int, int> level_rank_of;
    auto to_level = [&](int r) {
        auto it = level_rank_of.find(r);
        if (it != level_rank_of.end()) return it->second;
        int out = MPI_UNDEFINED;
        MPI_Group_translate_ranks(tree_group, 1, &r, level_group, &out);
        if (out == MPI_UNDEFINED) throw std::runtime_error("F22 phase: a box owner outside the level's communicator");
        level_rank_of.emplace(r, out);
        return out;
    };
    const int dimension = tree->dimension;
    const int64_t side = int64_t(1) << level.level;  // boxes per dimension
    std::vector<std::vector<int64_t>> send(static_cast<size_t>(lsize));
    for (const auto& kv : kept_of) {
        const int64_t morton = kv.first;
        const int64_t n_old = kv.second.first;
        const auto& kept = kv.second.second;
        uint32_t c[3] = {0, 0, 0};
        morton::decode_nd(dimension, static_cast<uint64_t>(morton), c[0], c[1], c[2]);
        std::set<int> dest;
        const int lo[3] = {-3, dimension > 1 ? -3 : 0, dimension > 2 ? -3 : 0};
        const int hi[3] = {3, dimension > 1 ? 3 : 0, dimension > 2 ? 3 : 0};
        for (int dz = lo[2]; dz <= hi[2]; ++dz)
            for (int dy = lo[1]; dy <= hi[1]; ++dy)
                for (int dx = lo[0]; dx <= hi[0]; ++dx) {
                    const int64_t x = int64_t(c[0]) + dx, y = int64_t(c[1]) + dy, z = int64_t(c[2]) + dz;
                    if (x < 0 || y < 0 || z < 0 || x >= side || y >= side || z >= side) continue;
                    const int64_t m = static_cast<int64_t>(morton::encode_nd(
                        dimension, static_cast<uint32_t>(x), static_cast<uint32_t>(y), static_cast<uint32_t>(z)));
                    const int r = level.owner_rank_of(m);
                    if (r != me) dest.insert(r);
                }
        const auto* box = level.find_local_box(morton);
        for (int r : dest) {
            auto& buf = send[static_cast<size_t>(to_level(r))];
            buf.push_back(morton);
            buf.push_back(n_old);
            buf.push_back(static_cast<int64_t>(kept.size()));
            buf.insert(buf.end(), box->skeleton_before_f22.begin(), box->skeleton_before_f22.end());
            buf.insert(buf.end(), kept.begin(), kept.end());
        }
    }
    MPI_Group_free(&tree_group);
    MPI_Group_free(&level_group);
    std::vector<int> scount(static_cast<size_t>(lsize)), rcount(static_cast<size_t>(lsize)), sdispl(static_cast<size_t>(lsize)),
        rdispl(static_cast<size_t>(lsize));
    for (int q = 0; q < lsize; ++q) scount[static_cast<size_t>(q)] = static_cast<int>(send[static_cast<size_t>(q)].size());
    MPI_Alltoall(scount.data(), 1, MPI_INT, rcount.data(), 1, MPI_INT, level_comm);
    std::vector<int64_t> sbuf, rbuf;
    for (int q = 0; q < lsize; ++q) {
        sdispl[static_cast<size_t>(q)] = static_cast<int>(sbuf.size());
        sbuf.insert(sbuf.end(), send[static_cast<size_t>(q)].begin(), send[static_cast<size_t>(q)].end());
    }
    int rtotal = 0;
    for (int q = 0; q < lsize; ++q) {
        rdispl[static_cast<size_t>(q)] = rtotal;
        rtotal += rcount[static_cast<size_t>(q)];
    }
    rbuf.resize(static_cast<size_t>(rtotal));
    MPI_Alltoallv(sbuf.data(), scount.data(), sdispl.data(), MPI_INT64_T, rbuf.data(), rcount.data(), rdispl.data(),
                  MPI_INT64_T, level_comm);
    // a copy's live sets from its owner's: S_old before the phase, S' after
    auto reduce = [&](std::vector<int64_t>& now, std::vector<int64_t>& before, const int64_t* old, int64_t n_old,
                      const int64_t* kept, int64_t nk) {
        before.assign(old, old + n_old);
        now.resize(static_cast<size_t>(nk));
        for (int64_t i = 0; i < nk; ++i) now[static_cast<size_t>(i)] = before[static_cast<size_t>(kept[i])];
    };
    for (size_t at = 0; at < rbuf.size();) {
        const int64_t morton = rbuf[at], n_old = rbuf[at + 1], nk = rbuf[at + 2];
        const int64_t* old = rbuf.data() + at + 3;
        const int64_t* kept = old + n_old;
        at += 3 + static_cast<size_t>(n_old + nk);
        bool held = false;
        if (auto* g = level.find_ghost_box(morton)) {
            reduce(g->skeleton_indices, g->skeleton_before_f22, old, n_old, kept, nk);
            held = true;
        }
        auto ait = level.assisting_box_points_for_kernel_evaluation.find(morton);
        if (ait != level.assisting_box_points_for_kernel_evaluation.end()) {
            auto& a = level.assisting_boxes[static_cast<size_t>(ait->second)];
            reduce(a.skel_indices, a.skel_before, old, n_old, kept, nk);
            held = true;
        }
        if (held) kept_of[morton] = {n_old, std::vector<int64_t>(kept, kept + nk)};
    }
}

// =====================================================================
// The F22 phase of a level (after its F11 waves): every local box's F22
// part is skeletonized against the far ring and the near field's live F11
// points, its R22 sparsified, the box's live set reduced to S' = S11, S22.
// With `device` (H2_use_gpu: LevelEliminator::f22_boxes) the boxes' steps
// run there from the plans made here; `blocks_back` then brings the level's
// blocks to the host sliced to the new live sets (the boxes the step
// sparsified, this rank's and its copies').
// =====================================================================
template<typename CoordType, typename DataType, typename KernelType>
void partial_f22_phase(fmm::ParallelTree<CoordType, DataType>* tree, fmm::TreeLevel<CoordType, DataType>& level,
                       KernelType* kernel, double tolerance, const std::vector<int64_t>& occupied_local,
                       bool print_detail, MPI_Comm level_comm,
                       const std::function<void(std::vector<fmm::PartialF22Work<DataType>>&)>& device = {},
                       const std::function<void(const fmm::PartialF22Kept&)>& blocks_back = {}) {
    using Box = fmm::BoxData<CoordType, DataType>;
    using MS = fmm::MatrixStorage<DataType>;
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
    // BPACK_TRACE=phase: the step's parts, thread-seconds summed over the boxes
    double t_sketch = 0.0, t_target = 0.0, t_id = 0.0, t_records = 0.0;
    const double t_loop0 = MPI_Wtime();
    // the boxes' kept positions within S_old; the shared state (live sets,
    // the boxes' own blocks) is changed after the loop, the loop only reads
    // the neighbours and rings
    std::vector<std::vector<int64_t>> kept(occupied_local.size());
    // morton -> the old live set's size and the kept positions in it: this
    // rank's sparsified boxes, then (several ranks) the other ranks' boxes
    // this rank holds copies of
    fmm::PartialF22Kept kept_of;
    bool sliced_on_device = false;
    if (device) {
        // ---- the device takes the boxes' steps: their plans, then the results
        std::vector<fmm::PartialF22Work<DataType>> works;
        std::vector<size_t> work_box;  // the box (index of occupied_local) of each work
        for (size_t oi = 0; oi < occupied_local.size(); ++oi) {
            Box& b = level.local_boxes[static_cast<size_t>(occupied_local[oi])];
            fmm::PartialF22Work<DataType> w;
            w.morton = b.morton_index;
            w.S_old = b.skeleton_before_f22;
            for (int64_t q = 0; q < static_cast<int64_t>(w.S_old.size()); ++q) {
                (tree->partial_f11(b.point_indices[static_cast<size_t>(w.S_old[static_cast<size_t>(q)])]) ? w.pos11 : w.pos22)
                    .push_back(q);
            }
            if (w.pos22.empty()) continue;
            w.near.resize(b.one_hop.size());
            for (size_t idx = 0; idx < b.one_hop.size(); ++idx) {
                const int64_t nm = b.one_hop[idx];
                const std::vector<int64_t>* nS = nullptr;
                const std::vector<int64_t>* ids = nullptr;
                if (const Box* nb = level.find_local_box(nm)) {
                    nS = &nb->skeleton_before_f22;
                    ids = &nb->point_indices;
                } else if (const Box* g = level.find_ghost_box(nm)) {
                    nS = &g->skeleton_indices;
                    ids = &g->point_indices;
                } else {
                    auto ait = level.assisting_box_points_for_kernel_evaluation.find(nm);
                    if (ait == level.assisting_box_points_for_kernel_evaluation.end()) {
                        throw std::runtime_error("F22 phase: neighbour " + std::to_string(nm) + " of box " +
                                                 std::to_string(b.morton_index) + " has no local or remote copy");
                    }
                    const auto& a = level.assisting_boxes[static_cast<size_t>(ait->second)];
                    nS = &a.skel_indices;
                    ids = &a.indices;
                }
                if (nS->empty()) {
                    throw std::runtime_error("F22 phase: neighbour " + std::to_string(nm) + " of box " +
                                             std::to_string(b.morton_index) + " has no published skeleton");
                }
                auto& nr = w.near[idx];
                nr.morton = nm;
                nr.S_old = nS;
                for (int64_t q = 0; q < static_cast<int64_t>(nS->size()); ++q) {
                    (tree->partial_f11((*ids)[static_cast<size_t>((*nS)[static_cast<size_t>(q)])]) ? nr.q11 : nr.q22)
                        .push_back(q);
                }
            }
            works.push_back(std::move(w));
            work_box.push_back(oi);
        }
        const double t_device0 = MPI_Wtime();
        t_target += t_device0 - t_loop0;  // (the plans)
        device(works);
        t_sketch += MPI_Wtime() - t_device0;  // (the device's whole step)
        const double t_results0 = MPI_Wtime();
        for (size_t wi = 0; wi < works.size(); ++wi) {
            auto& w = works[wi];
            if (!w.done) throw std::runtime_error("F22 phase: the device skipped box " + std::to_string(w.morton));
            for (int64_t i : w.skel22)
                if (i < 0 || i >= static_cast<int64_t>(w.pos22.size())) throw std::logic_error("F22 phase: device skeleton index " + std::to_string(i) + " of " + std::to_string(w.pos22.size()));
            for (int64_t i : w.red22)
                if (i < 0 || i >= static_cast<int64_t>(w.pos22.size())) throw std::logic_error("F22 phase: device redundant index " + std::to_string(i) + " of " + std::to_string(w.pos22.size()));
            Box& b = level.local_boxes[static_cast<size_t>(occupied_local[work_box[wi]])];
            const int64_t k22 = static_cast<int64_t>(w.skel22.size()), r22 = static_cast<int64_t>(w.red22.size());
            if (r22 == 0) {
                kept_points += k22;
                continue;
            }
            std::vector<int64_t> posS(w.pos11), posR;
            for (int64_t i : w.skel22) posS.push_back(w.pos22[static_cast<size_t>(i)]);
            for (int64_t j : w.red22) posR.push_back(w.pos22[static_cast<size_t>(j)]);
            auto& f = b.f22;
            f.sparsified = true;
            const size_t k11 = w.pos11.size();
            for (size_t i = k11; i < posS.size(); ++i) f.skeleton.push_back(w.S_old[static_cast<size_t>(posS[i])]);
            for (int64_t q : posR) f.redundant.push_back(w.S_old[static_cast<size_t>(q)]);
            f.T = std::move(w.T);
            f.X_SR = std::move(w.X_SR);
            f.X_RS = std::move(w.X_RS);
            f.X_RR = std::move(w.X_RR);
            f.X_NR = std::move(w.X_NR);
            f.X_RN = std::move(w.X_RN);
            f.neighbor_f22_counts.assign(b.one_hop.size(), 0);
            f.neighbor_f22_slots.assign(b.one_hop.size(), {});
            for (size_t idx = 0; idx < b.one_hop.size(); ++idx) {
                const auto& nr = w.near[idx];
                for (int64_t q : nr.q22) f.neighbor_f22_slots[idx].push_back((*nr.S_old)[static_cast<size_t>(q)]);
                f.neighbor_f22_counts[idx] = static_cast<int64_t>(nr.q22.size());
            }
            kept[work_box[wi]] = std::move(posS);
            ++sparsified_boxes;
            sparsified_points += r22;
            kept_points += k22;
        }
        t_records += MPI_Wtime() - t_results0;  // (the results into the boxes)
        // the new live sets, then the blocks back from the device sliced to them
        for (size_t oi = 0; oi < occupied_local.size(); ++oi) {
            const Box& b = level.local_boxes[static_cast<size_t>(occupied_local[oi])];
            if (b.f22.sparsified) kept_of.emplace(b.morton_index, std::make_pair(static_cast<int64_t>(b.skeleton_before_f22.size()), kept[oi]));
        }
        if (tree->mpi_size > 1) partial_f22_publish(tree, level, level_comm, kept_of);
        const double t_back0 = MPI_Wtime();
        if (blocks_back) blocks_back(kept_of);
        t_id += MPI_Wtime() - t_back0;  // (the blocks back)
        for (size_t oi = 0; oi < occupied_local.size(); ++oi) {
            Box& b = level.local_boxes[static_cast<size_t>(occupied_local[oi])];
            if (!b.f22.sparsified) continue;
            const auto& posS = kept[oi];
            const int64_t kp = static_cast<int64_t>(posS.size());
            b.skeleton_indices.clear();
            for (int64_t q : posS) b.skeleton_indices.push_back(b.skeleton_before_f22[static_cast<size_t>(q)]);
            if (!b.schur_complement.is_allocated()) {
                // a self block the device never formed (nothing eliminated, no
                // update): the kernel's over S'
                std::vector<int64_t> ids(static_cast<size_t>(kp));
                for (int64_t i = 0; i < kp; ++i)
                    ids[static_cast<size_t>(i)] = b.point_indices[static_cast<size_t>(b.skeleton_indices[static_cast<size_t>(i)])];
                MS block;
                block.allocate(kp, kp, MS::FULL);
                kernel->evaluate_block_by_index(ids.data(), kp, ids.data(), kp, block.data.data(), kp);
                b.schur_complement = std::move(block);
            } else if (b.schur_complement.rows != kp || b.schur_complement.cols != kp) {
                throw std::runtime_error("F22 phase: box " + std::to_string(b.morton_index) + "'s self block came back " +
                                         std::to_string(b.schur_complement.rows) + " x " +
                                         std::to_string(b.schur_complement.cols) + ", not over S'");
            }
        }
        sliced_on_device = true;
    } else {
#pragma omp parallel
    {
        fmm::FactorizationThreadScratch<CoordType, DataType> scratch;
#pragma omp for schedule(dynamic, 1) reduction(+ : sparsified_boxes, sparsified_points, kept_points, t_sketch, t_target, t_id, t_records)
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
                if (!b.schur_complement.is_allocated() && k_old == b.num_points) {
                    // a box nothing eliminated and no update reached keeps no
                    // self block on the device path: it is the kernel's
                    std::vector<int64_t> ids(static_cast<size_t>(k_old));
                    for (int64_t q = 0; q < k_old; ++q)
                        ids[static_cast<size_t>(q)] = b.point_indices[static_cast<size_t>(S_old[static_cast<size_t>(q)])];
                    MS block;
                    block.allocate(k_old, k_old, MS::FULL);
                    kernel->evaluate_block_by_index(ids.data(), k_old, ids.data(), k_old, block.data.data(), k_old);
                    b.schur_complement = std::move(block);
                }
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
                double t0 = MPI_Wtime();
                fmm::gather_id_target_streamed_nonsymmetric(tree, &b, level, kernel, scratch, b.on_boundary, tolerance,
                                                             &S_old);
                b.use_full_set = use_full_at_f11;
                t_sketch += MPI_Wtime() - t0;
                t0 = MPI_Wtime();
                if (!scratch.streamed_sketch_valid || scratch.workspace_cols != k_old) {
                    throw std::runtime_error("F22 phase: no streamed sketch for box " + std::to_string(b.morton_index));
                }
                const int64_t d = scratch.streamed_sketch_rows;
                // a neighbour's skeleton before the phase and its point ids: a
                // local box's, or its remote copy's (a ghost on a CA level, an
                // assisting box on a Color one; the copies change only after
                // the phase)
                // K(nb, b) and K(b, nb) over the two boxes' skeletons before the
                // phase: this box's views, else the neighbour's views of this box
                // (a local neighbour), else the kernel - a pair no elimination
                // touched keeps no view, and the transition reads it the same way
                struct NearF11 {
                    const std::vector<int64_t>* S_old = nullptr;
                    const std::vector<int64_t>* ids = nullptr;
                    std::vector<int64_t> q11, q22;
                    const fmm::ModifiedBlock<DataType>* nb_b = nullptr;  // K(nb, b)
                    bool nb_b_own = true;                                // b's A_NS[nb], else nb's A_SN[b]
                    const fmm::ModifiedBlock<DataType>* b_nb = nullptr;  // K(b, nb)
                    bool b_nb_own = true;                                // b's A_SN[nb], else nb's A_NS[b]
                };
                std::vector<NearF11> near(b.one_hop.size());
                int64_t near_rows = 0;
                for (size_t idx = 0; idx < b.one_hop.size(); ++idx) {
                    const int64_t nm = b.one_hop[idx];
                    const std::vector<int64_t>* nS_ptr = nullptr;
                    const std::vector<int64_t>* ids = nullptr;
                    if (const Box* nb = level.find_local_box(nm)) {
                        nS_ptr = &nb->skeleton_before_f22;
                        ids = &nb->point_indices;
                    } else if (const Box* g = level.find_ghost_box(nm)) {
                        nS_ptr = &g->skeleton_indices;
                        ids = &g->point_indices;
                    } else {
                        auto ait = level.assisting_box_points_for_kernel_evaluation.find(nm);
                        if (ait == level.assisting_box_points_for_kernel_evaluation.end()) {
                            throw std::runtime_error("F22 phase: neighbour " + std::to_string(nm) + " of box " +
                                                     std::to_string(b.morton_index) + " has no local or remote copy");
                        }
                        const auto& a = level.assisting_boxes[static_cast<size_t>(ait->second)];
                        nS_ptr = &a.skel_indices;
                        ids = &a.indices;
                    }
                    if (nS_ptr->empty()) {
                        throw std::runtime_error("F22 phase: neighbour " + std::to_string(nm) + " of box " +
                                                 std::to_string(b.morton_index) + " has no published skeleton");
                    }
                    near[idx].S_old = nS_ptr;
                    near[idx].ids = ids;
                    {
                        const int64_t n_nb = static_cast<int64_t>(nS_ptr->size());
                        auto* ns = a_ns_block(b, nm);
                        auto* sn = a_sn_block(b, nm);
                        Box* nb_local = level.find_local_box(nm);
                        if (ns != nullptr && ns->a_ns_is_allocated()) {
                            if (ns->a_ns_rows() != n_nb || ns->a_ns_cols() != k_old) {
                                throw std::runtime_error("F22 phase: box " + std::to_string(b.morton_index) + "'s view of " +
                                                         std::to_string(nm) + " is not over the two skeletons");
                            }
                            near[idx].nb_b = ns;
                        } else if (nb_local != nullptr) {
                            if (auto* o = a_sn_block(*nb_local, b.morton_index); o != nullptr && o->A_SN.is_allocated()) {
                                if (o->A_SN.rows != n_nb || o->A_SN.cols != k_old) {
                                    throw std::runtime_error("F22 phase: box " + std::to_string(nm) + "'s view of " +
                                                             std::to_string(b.morton_index) + " is not over the two skeletons");
                                }
                                near[idx].nb_b = o;
                                near[idx].nb_b_own = false;
                            }
                        }
                        if (sn != nullptr && sn->A_SN.is_allocated()) {
                            if (sn->A_SN.rows != k_old || sn->A_SN.cols != n_nb) {
                                throw std::runtime_error("F22 phase: box " + std::to_string(b.morton_index) + "'s view of " +
                                                         std::to_string(nm) + " is not over the two skeletons");
                            }
                            near[idx].b_nb = sn;
                        } else if (nb_local != nullptr) {
                            if (auto* o = a_ns_block(*nb_local, b.morton_index); o != nullptr && o->a_ns_is_allocated()) {
                                if (o->a_ns_rows() != k_old || o->a_ns_cols() != n_nb) {
                                    throw std::runtime_error("F22 phase: box " + std::to_string(nm) + "'s view of " +
                                                             std::to_string(b.morton_index) + " is not over the two skeletons");
                                }
                                near[idx].b_nb = o;
                                near[idx].b_nb_own = false;
                            }
                        }
                    }
                    const auto& nS = *nS_ptr;
                    for (int64_t q = 0; q < static_cast<int64_t>(nS.size()); ++q) {
                        if (tree->partial_f11((*ids)[static_cast<size_t>(nS[static_cast<size_t>(q)])])) near[idx].q11.push_back(q);
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
                // the point ids of a neighbour's positions and of this box's
                auto nb_ids_of = [&](size_t idx, const std::vector<int64_t>& qs) {
                    std::vector<int64_t> out(qs.size());
                    for (size_t i = 0; i < qs.size(); ++i)
                        out[i] = (*near[idx].ids)[static_cast<size_t>((*near[idx].S_old)[static_cast<size_t>(qs[i])])];
                    return out;
                };
                auto b_ids_of = [&](const std::vector<int64_t>& ps) {
                    std::vector<int64_t> out(ps.size());
                    for (size_t i = 0; i < ps.size(); ++i) out[i] = b.point_indices[static_cast<size_t>(S_old[static_cast<size_t>(ps[i])])];
                    return out;
                };
                // K(nb, b)(qs, ps) and K(b, nb)(ps, qs), column-major
                auto coupling_nb_b = [&](size_t idx, const std::vector<int64_t>& qs, const std::vector<int64_t>& ps) {
                    const int64_t m = static_cast<int64_t>(qs.size()), n = static_cast<int64_t>(ps.size());
                    std::vector<DataType> out(static_cast<size_t>(m * n));
                    const auto* v = near[idx].nb_b;
                    if (v == nullptr) {
                        if (m > 0 && n > 0) {
                            const auto ri = nb_ids_of(idx, qs), ci = b_ids_of(ps);
                            kernel->evaluate_block_by_index(ri.data(), m, ci.data(), n, out.data(), m);
                        }
                        return out;
                    }
                    for (int64_t c = 0; c < n; ++c)
                        for (int64_t i = 0; i < m; ++i)
                            out[static_cast<size_t>(i + c * m)] = near[idx].nb_b_own ? v->a_ns(qs[static_cast<size_t>(i)], ps[static_cast<size_t>(c)])
                                                                                     : v->A_SN(qs[static_cast<size_t>(i)], ps[static_cast<size_t>(c)]);
                    return out;
                };
                auto coupling_b_nb = [&](size_t idx, const std::vector<int64_t>& ps, const std::vector<int64_t>& qs) {
                    const int64_t m = static_cast<int64_t>(ps.size()), n = static_cast<int64_t>(qs.size());
                    std::vector<DataType> out(static_cast<size_t>(m * n));
                    const auto* v = near[idx].b_nb;
                    if (v == nullptr) {
                        if (m > 0 && n > 0) {
                            const auto ri = b_ids_of(ps), ci = nb_ids_of(idx, qs);
                            kernel->evaluate_block_by_index(ri.data(), m, ci.data(), n, out.data(), m);
                        }
                        return out;
                    }
                    for (int64_t c = 0; c < n; ++c)
                        for (int64_t i = 0; i < m; ++i)
                            out[static_cast<size_t>(i + c * m)] = near[idx].b_nb_own ? v->A_SN(ps[static_cast<size_t>(i)], qs[static_cast<size_t>(c)])
                                                                                     : v->a_ns(ps[static_cast<size_t>(i)], qs[static_cast<size_t>(c)]);
                    return out;
                };
                for (size_t idx = 0; idx < b.one_hop.size(); ++idx) {  // K(F11_nb, F22_b)
                    const int64_t m = static_cast<int64_t>(near[idx].q11.size());
                    if (m == 0) continue;
                    const auto blk = coupling_nb_b(idx, near[idx].q11, pos22);
                    for (int64_t c = 0; c < n22; ++c)
                        for (int64_t i = 0; i < m; ++i)
                            target[static_cast<size_t>(at + i + c * rows)] = blk[static_cast<size_t>(i + c * m)];
                    at += m;
                }
                for (int64_t c = 0; c < n22; ++c)  // the box's own S11 rows
                    for (int64_t i = 0; i < k11; ++i)
                        target[static_cast<size_t>(at + i + c * rows)] = Aat(pos11[static_cast<size_t>(i)], pos22[static_cast<size_t>(c)]);
                at += k11;
                for (size_t idx = 0; idx < b.one_hop.size(); ++idx) {  // K(F22_b, F11_nb)^T
                    const int64_t m = static_cast<int64_t>(near[idx].q11.size());
                    if (m == 0) continue;
                    const auto blk = coupling_b_nb(idx, pos22, near[idx].q11);  // n22 x m
                    for (int64_t c = 0; c < n22; ++c)
                        for (int64_t i = 0; i < m; ++i)
                            target[static_cast<size_t>(at + i + c * rows)] = blk[static_cast<size_t>(c + i * n22)];
                    at += m;
                }
                for (int64_t c = 0; c < n22; ++c)
                    for (int64_t i = 0; i < k11; ++i)
                        target[static_cast<size_t>(at + i + c * rows)] = Aat(pos22[static_cast<size_t>(c)], pos11[static_cast<size_t>(i)]);
                at += k11;
                if (at != rows) throw std::logic_error("F22 phase: target rows");
                t_target += MPI_Wtime() - t0;
                t0 = MPI_Wtime();
                fmm::IDResult<DataType> id = fmm::compute_id_complex(target.data(), rows, n22, rows, tolerance, 0, abs_tol);
                fmm::ensure_nonempty_box_id_rank(id, n22);
                t_id += MPI_Wtime() - t0;
                t0 = MPI_Wtime();
                const int64_t k22 = static_cast<int64_t>(id.skeleton_indices.size());
                const int64_t r22 = static_cast<int64_t>(id.redundant_indices.size());
                // nothing to sparsify: the box stays as it is (its live set
                // unchanged, no records)
                if (r22 == 0) {
                    kept_points += k22;
                    continue;
                }
                // the ID's skeleton and redundant columns in their order in S_old
                // (T's rows and columns permuted alike): S_old is S11 then F22
                // (compute_id_partial), so S' = S11, S22 keeps S_old's order, and
                // with nothing sparsified S' is S_old - a live set the same size
                // as before is then the same list, as the lazy fill's row lookup
                // assumes
                std::vector<int64_t> ps(static_cast<size_t>(k22)), pr_(static_cast<size_t>(r22));
                for (int64_t i = 0; i < k22; ++i) ps[static_cast<size_t>(i)] = i;
                for (int64_t j = 0; j < r22; ++j) pr_[static_cast<size_t>(j)] = j;
                std::sort(ps.begin(), ps.end(), [&](int64_t x, int64_t y) { return id.skeleton_indices[static_cast<size_t>(x)] < id.skeleton_indices[static_cast<size_t>(y)]; });
                std::sort(pr_.begin(), pr_.end(), [&](int64_t x, int64_t y) { return id.redundant_indices[static_cast<size_t>(x)] < id.redundant_indices[static_cast<size_t>(y)]; });
                // ---- the new live set S' = S11, S22 and R22, as positions and slots
                std::vector<int64_t> posS(pos11), posR;
                for (int64_t i = 0; i < k22; ++i) posS.push_back(pos22[static_cast<size_t>(id.skeleton_indices[static_cast<size_t>(ps[static_cast<size_t>(i)])])]);
                for (int64_t j = 0; j < r22; ++j) posR.push_back(pos22[static_cast<size_t>(id.redundant_indices[static_cast<size_t>(pr_[static_cast<size_t>(j)])])]);
                const int64_t kp = static_cast<int64_t>(posS.size());
                auto& f = b.f22;
                f.sparsified = true;
                for (int64_t i = 0; i < k22; ++i) f.skeleton.push_back(S_old[static_cast<size_t>(posS[static_cast<size_t>(k11 + i)])]);
                for (int64_t j = 0; j < r22; ++j) f.redundant.push_back(S_old[static_cast<size_t>(posR[static_cast<size_t>(j)])]);
                // T over S' rows (zero on S11)
                f.T.allocate(kp, r22, MS::FULL);
                for (int64_t j = 0; j < r22; ++j)
                    for (int64_t i = 0; i < k22; ++i)
                        f.T.data[static_cast<size_t>(k11 + i + j * kp)] = id.interpolation(ps[static_cast<size_t>(i)], pr_[static_cast<size_t>(j)]);
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
                    const auto& nS = *near[idx].S_old;
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
                    // K(nb raw F22, b) over (posR, posS), K(b, nb raw F22) over (posR, posS)
                    const auto& q22 = near[idx].q22;
                    const auto KR = coupling_nb_b(idx, q22, posR), KS = coupling_nb_b(idx, q22, posS);    // m x r22, m x kp
                    const auto LR = coupling_b_nb(idx, posR, q22), LS = coupling_b_nb(idx, posS, q22);    // r22 x m, kp x m
                    std::vector<DataType> A_NR = KR, A_NS = KS, A_RN = LR, A_SN = LS;
                    partial_gemm<DataType>('N', 'N', m, r22, kp, DataType(-1), A_NS.data(), m, T, kp, DataType(1), A_NR.data(), m);
                    partial_gemm<DataType>('T', 'N', r22, m, kp, DataType(-1), T, kp, A_SN.data(), kp, DataType(1), A_RN.data(), r22);
                    for (int64_t j = 0; j < r22; ++j)
                        std::copy_n(A_NR.data() + j * m, m, f.X_NR.data.data() + row0 + j * total);
                    for (int64_t i = 0; i < m; ++i)
                        std::copy_n(A_RN.data() + i * r22, r22, f.X_RN.data.data() + (row0 + i) * r22);
                    row0 += m;
                }
                kept[static_cast<size_t>(oi)] = std::move(posS);
                t_records += MPI_Wtime() - t0;
                ++sparsified_boxes;
                sparsified_points += r22;
                kept_points += k22;
            } catch (...) {
#pragma omp critical(partial_f22_phase_failure)
                if (!failure) failure = std::current_exception();
            }
        }
    }
    }
    if (failure) std::rethrow_exception(failure);
    const double t_loop = MPI_Wtime() - t_loop0;
    if (!sliced_on_device) {
    // ---- the boxes' live sets and their own blocks over S'
    for (size_t oi = 0; oi < occupied_local.size(); ++oi) {
        Box& b = level.local_boxes[static_cast<size_t>(occupied_local[oi])];
        if (!b.f22.sparsified) continue;
        const auto& posS = kept[oi];
        const int64_t kp = static_cast<int64_t>(posS.size());
        if (!b.schur_complement.is_allocated()) {
            // a self block the device never formed (nothing eliminated, no
            // update): the kernel's over S'
            std::vector<int64_t> ids(static_cast<size_t>(kp));
            for (int64_t i = 0; i < kp; ++i)
                ids[static_cast<size_t>(i)] =
                    b.point_indices[static_cast<size_t>(b.skeleton_before_f22[static_cast<size_t>(posS[static_cast<size_t>(i)])])];
            MS block;
            block.allocate(kp, kp, MS::FULL);
            kernel->evaluate_block_by_index(ids.data(), kp, ids.data(), kp, block.data.data(), kp);
            b.schur_complement = std::move(block);
        } else {
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
    // ---- the new live sets: this rank's boxes, and (several ranks) the
    // other ranks' boxes this rank holds copies of
    for (size_t oi = 0; oi < occupied_local.size(); ++oi) {
        const Box& b = level.local_boxes[static_cast<size_t>(occupied_local[oi])];
        if (b.f22.sparsified) kept_of.emplace(b.morton_index, std::make_pair(static_cast<int64_t>(b.skeleton_before_f22.size()), kept[oi]));
    }
    if (tree->mpi_size > 1) partial_f22_publish(tree, level, level_comm, kept_of);
    // ---- every box's views of its sparsified neighbours: their rows and
    // columns to the new live sets
    for (int64_t li : occupied_local) {
        Box& c = level.local_boxes[static_cast<size_t>(li)];
        for (size_t idx = 0; idx < c.one_hop.size(); ++idx) {
            const int64_t nm = c.one_hop[idx];
            const auto kit = kept_of.find(nm);
            if (kit == kept_of.end()) continue;
            const int64_t n_old = kit->second.first;
            const auto& pos = kit->second.second;
            const int64_t kn = static_cast<int64_t>(pos.size());
            if (auto* sn = a_sn_block(c, nm); sn != nullptr && sn->A_SN.is_allocated()) {  // c rows x nb cols
                if (sn->A_SN.cols != n_old) {
                    throw std::runtime_error("F22 phase: box " + std::to_string(c.morton_index) + "'s view of " +
                                             std::to_string(nm) + " is not over its old skeleton (columns)");
                }
                const int64_t m = sn->A_SN.rows;
                std::vector<DataType> v(static_cast<size_t>(m * kn));
                for (int64_t j = 0; j < kn; ++j)
                    for (int64_t i = 0; i < m; ++i) v[static_cast<size_t>(i + j * m)] = sn->A_SN(i, pos[static_cast<size_t>(j)]);
                sn->A_SN.set_owned(m, kn, std::move(v), MS::FULL);
            }
            if (auto* ns = a_ns_block(c, nm); ns != nullptr && ns->a_ns_is_allocated()) {  // nb rows x c cols
                if (ns->a_ns_rows() != n_old) {
                    throw std::runtime_error("F22 phase: box " + std::to_string(c.morton_index) + "'s view of " +
                                             std::to_string(nm) + " is not over its old skeleton (rows)");
                }
                const int64_t m = ns->a_ns_cols();
                std::vector<DataType> v(static_cast<size_t>(kn * m));
                for (int64_t j = 0; j < m; ++j)
                    for (int64_t i = 0; i < kn; ++i) v[static_cast<size_t>(i + j * kn)] = ns->a_ns(pos[static_cast<size_t>(i)], j);
                ns->set_a_ns_owned(kn, m, std::move(v), MS::FULL);
            }
        }
    }
    }
    int64_t counts[3] = {sparsified_boxes, sparsified_points, kept_points};
    if (tree->mpi_size > 1) MPI_Allreduce(MPI_IN_PLACE, counts, 3, MPI_INT64_T, MPI_SUM, level_comm);
    if (print_detail && fmm::env::trace("phase")) {
        if (device) {
            std::printf("  F22 phase parts (this rank): plans %.2f s, device %.2f s, results %.2f s, blocks back %.2f s, "
                        "then %.2f s serial\n", t_target, t_sketch, t_records, t_id, MPI_Wtime() - t_loop0 - t_loop);
        } else {
            std::printf("  F22 phase parts (this rank): loop %.2f s, then %.2f s serial; thread-seconds: sketch %.2f, "
                        "target %.2f, ID %.2f, records %.2f\n", t_loop, MPI_Wtime() - t_loop0 - t_loop, t_sketch,
                        t_target, t_id, t_records);
        }
    }
    if (print_detail) {
        std::printf("  F22 phase: %lld boxes, %lld points sparsified, %lld kept\n", static_cast<long long>(counts[0]),
                    static_cast<long long>(counts[1]), static_cast<long long>(counts[2]));
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

// f(c0, c1) over column ranges [c0, c1) of n columns, in parallel: the
// root's products a column block a thread (the BLAS is sequential)
template<typename F>
void partial_columns_parallel(int64_t n, F&& f) {
    const int64_t threads = std::max(1, omp_get_max_threads());
    const int64_t chunk = std::max<int64_t>(1, (n + 4 * threads - 1) / (4 * threads));
    std::exception_ptr failure;
#pragma omp parallel for schedule(dynamic, 1)
    for (int64_t c0 = 0; c0 < n; c0 += chunk) {
        try {
            f(c0, std::min(n, c0 + chunk));
        } catch (...) {
#pragma omp critical(partial_columns_parallel_failure)
            if (!failure) failure = std::current_exception();
        }
    }
    if (failure) std::rethrow_exception(failure);
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
        partial_columns_parallel(k, [&](int64_t j0, int64_t j1) {
            for (int64_t j = j0; j < j1; ++j)
                for (int64_t i = 0; i < m; ++i)
                    out.data[static_cast<size_t>(i + j * m)] =
                        A.data[static_cast<size_t>(rows[static_cast<size_t>(i)] + cols[static_cast<size_t>(j)] * A.lda)];
        });
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
    partial_columns_parallel(pr.n22, [&](int64_t c0, int64_t c1) {
        partial_root_solve11(pr.n11, root_box.X_RR, root_box.X_RR_pivots, pr.w.data.data() + c0 * pr.n11, c1 - c0,
                             std::max<int64_t>(pr.n11, 1));
    });
    gather(pr.slots22, pr.slots11, pr.a21);
    gather(pr.slots22, pr.slots22, pr.cb);
    partial_columns_parallel(pr.n22, [&](int64_t c0, int64_t c1) {  // CB -= A21 W
        partial_gemm<DataType>('N', 'N', pr.n22, c1 - c0, pr.n11, DataType(-1.0), pr.a21.data.data(),
                               std::max<int64_t>(pr.n22, 1), pr.w.data.data() + c0 * pr.n11, std::max<int64_t>(pr.n11, 1),
                               DataType(1.0), pr.cb.data.data() + c0 * pr.n22, std::max<int64_t>(pr.n22, 1));
    });
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

// Scratch solve records over a level's local boxes (the solve's layout:
// the box's points x nrhs), zero or a copy of the solve's vectors.
template<typename CoordType, typename DataType>
std::vector<std::vector<fmm::SolveDataRequest<CoordType, DataType>>> partial_scratch(
    fmm::ParallelTree<CoordType, DataType>* tree,
    const std::vector<std::vector<fmm::SolveDataRequest<CoordType, DataType>>>* copy_of, int nrhs) {
    int rank = 0;
    MPI_Comm_rank(tree->comm, &rank);
    std::vector<std::vector<fmm::SolveDataRequest<CoordType, DataType>>> v(static_cast<size_t>(tree->num_levels));
    for (int l = 0; l < tree->num_levels; ++l) {
        auto& level = tree->levels[static_cast<size_t>(l)];
        if (!level.is_process_active) continue;
        v[static_cast<size_t>(l)].resize(static_cast<size_t>(level.num_boxes_local));
        for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
            const auto& box = level.local_boxes[static_cast<size_t>(bi)];
            auto& r = v[static_cast<size_t>(l)][static_cast<size_t>(bi)];
            r.initialize(box.morton_index, rank, box.num_points, nrhs);
            if (copy_of != nullptr) r.left_side = (*copy_of)[static_cast<size_t>(l)][static_cast<size_t>(bi)].left_side;
        }
    }
    return v;
}

// The forward partial solve's output: from the swept vectors (kept), a
// down-pass applying P22^T (x_R22 += T22^T x_S') to a copy, through the
// solve's parent-to-children scatter (which places a parent's values at
// its children's live sets, on their ranks); the leaves' F22 entries in
// F22 index order, replicated.
template<typename CoordType, typename DataType>
void partial_forward_output(fmm::ParallelTree<CoordType, DataType>* tree,
                            const std::vector<std::vector<fmm::SolveDataRequest<CoordType, DataType>>>& solve_data,
                            int nrhs, std::vector<DataType>& out2) {
    const int64_t pc = tree->partial_count, n22 = tree->num_points - pc;
    const int leaf = tree->num_levels - 1;
    auto v = partial_scratch(tree, &solve_data, nrhs);
    for (int l = 1; l <= leaf; ++l) {
        auto& level = tree->levels[static_cast<size_t>(l)];
        fmm::scatter_solution_to_children(level, tree->levels[static_cast<size_t>(l) - 1], v[static_cast<size_t>(l)],
                                          v[static_cast<size_t>(l) - 1], tree->dimension, tree->comm);
        if (!level.is_process_active) continue;
#pragma omp parallel for schedule(dynamic, 16)
        for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
            const auto& box = level.local_boxes[static_cast<size_t>(bi)];
            if (!box.f22.sparsified || box.f22.redundant.empty()) continue;
            auto& x = v[static_cast<size_t>(l)][static_cast<size_t>(bi)].left_side;
            const int64_t n = box.num_points, kp = box.f22.T.rows, r = box.f22.T.cols;
            std::vector<DataType> xs(static_cast<size_t>(kp * nrhs)), xr(static_cast<size_t>(r * nrhs));
            for (int c = 0; c < nrhs; ++c)
                for (int64_t i = 0; i < kp; ++i)
                    xs[static_cast<size_t>(i + c * kp)] = x[static_cast<size_t>(box.skeleton_indices[static_cast<size_t>(i)] + c * n)];
            partial_gemm<DataType>('T', 'N', r, nrhs, kp, DataType(1), box.f22.T.data.data(), kp, xs.data(), kp, DataType(0),
                                   xr.data(), r);
            for (int c = 0; c < nrhs; ++c)
                for (int64_t j = 0; j < r; ++j)
                    x[static_cast<size_t>(box.f22.redundant[static_cast<size_t>(j)] + c * n)] += xr[static_cast<size_t>(j + c * r)];
        }
    }
    out2.assign(static_cast<size_t>(n22) * static_cast<size_t>(nrhs), DataType(0));
    auto& leaf_level = tree->levels[static_cast<size_t>(leaf)];
    if (leaf_level.is_process_active) {
        for (int64_t bi = 0; bi < leaf_level.num_boxes_local; ++bi) {
            const auto& box = leaf_level.local_boxes[static_cast<size_t>(bi)];
            const auto& x = v[static_cast<size_t>(leaf)][static_cast<size_t>(bi)].left_side;
            for (int64_t s = 0; s < box.num_points; ++s) {
                const int64_t id = box.point_indices[static_cast<size_t>(s)];
                if (id < pc) continue;
                for (int c = 0; c < nrhs; ++c)
                    out2[static_cast<size_t>((id - pc) + c * n22)] = x[static_cast<size_t>(s + c * box.num_points)];
            }
        }
    }
    if (tree->mpi_size > 1) {
        MPI_Allreduce(MPI_IN_PLACE, out2.data(), static_cast<int>(out2.size()), partial_mpi_type<DataType>(), MPI_SUM,
                      tree->comm);
    }
}

// The backward partial solve's input: x2 (F22 index order) placed at the
// leaves and taken up with P22 (x_S' += T22 x_R) through the solve's
// children-to-parent gather on scratch vectors, then every level's F22
// entries (the root's included) into the solve's vectors; the F11 entries
// are left as the forward solve left them.
template<typename CoordType, typename DataType>
void partial_backward_input(fmm::ParallelTree<CoordType, DataType>* tree,
                            std::vector<std::vector<fmm::SolveDataRequest<CoordType, DataType>>>& solve_data,
                            int nrhs, const std::vector<DataType>& x2) {
    const int64_t pc = tree->partial_count, n22 = tree->num_points - pc;
    const int leaf = tree->num_levels - 1;
    auto t = partial_scratch<CoordType, DataType>(tree, nullptr, nrhs);
    auto& leaf_level = tree->levels[static_cast<size_t>(leaf)];
    if (leaf_level.is_process_active) {
        for (int64_t bi = 0; bi < leaf_level.num_boxes_local; ++bi) {
            const auto& box = leaf_level.local_boxes[static_cast<size_t>(bi)];
            auto& x = t[static_cast<size_t>(leaf)][static_cast<size_t>(bi)].left_side;
            for (int64_t s = 0; s < box.num_points; ++s) {
                const int64_t id = box.point_indices[static_cast<size_t>(s)];
                if (id < pc) continue;
                for (int c = 0; c < nrhs; ++c)
                    x[static_cast<size_t>(s + c * box.num_points)] = x2[static_cast<size_t>((id - pc) + c * n22)];
            }
        }
    }
    for (int l = leaf; l >= 1; --l) {
        auto& level = tree->levels[static_cast<size_t>(l)];
        if (level.is_process_active) {
#pragma omp parallel for schedule(dynamic, 16)
            for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
                const auto& box = level.local_boxes[static_cast<size_t>(bi)];
                if (!box.f22.sparsified || box.f22.redundant.empty()) continue;
                auto& x = t[static_cast<size_t>(l)][static_cast<size_t>(bi)].left_side;
                const int64_t n = box.num_points, kp = box.f22.T.rows, r = box.f22.T.cols;
                std::vector<DataType> xr(static_cast<size_t>(r * nrhs)), xs(static_cast<size_t>(kp * nrhs));
                for (int c = 0; c < nrhs; ++c)
                    for (int64_t j = 0; j < r; ++j)
                        xr[static_cast<size_t>(j + c * r)] = x[static_cast<size_t>(box.f22.redundant[static_cast<size_t>(j)] + c * n)];
                partial_gemm<DataType>('N', 'N', kp, nrhs, r, DataType(1), box.f22.T.data.data(), kp, xr.data(), r, DataType(0),
                                       xs.data(), kp);
                for (int c = 0; c < nrhs; ++c)
                    for (int64_t i = 0; i < kp; ++i)
                        x[static_cast<size_t>(box.skeleton_indices[static_cast<size_t>(i)] + c * n)] += xs[static_cast<size_t>(i + c * kp)];
            }
        }
        fmm::gather_skeleton_to_parent(level, tree->levels[static_cast<size_t>(l) - 1], t[static_cast<size_t>(l)],
                                       t[static_cast<size_t>(l) - 1], tree->dimension, tree->comm);
    }
    for (int l = 0; l <= leaf; ++l) {
        auto& level = tree->levels[static_cast<size_t>(l)];
        if (!level.is_process_active) continue;
        for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
            const auto& box = level.local_boxes[static_cast<size_t>(bi)];
            const auto& from = t[static_cast<size_t>(l)][static_cast<size_t>(bi)].left_side;
            auto& to = solve_data[static_cast<size_t>(l)][static_cast<size_t>(bi)].left_side;
            for (int64_t s = 0; s < box.num_points; ++s) {
                if (box.point_indices[static_cast<size_t>(s)] < pc) continue;
                for (int c = 0; c < nrhs; ++c)
                    to[static_cast<size_t>(s + c * box.num_points)] = from[static_cast<size_t>(s + c * box.num_points)];
            }
        }
    }
}

// A byte exchange over comm: send[r] goes to rank r, the result's [r]
// came from rank r (collective; empty buffers allowed).
inline std::vector<std::vector<char>> partial_exchange(MPI_Comm comm, const std::vector<std::vector<char>>& send, int tag) {
    int size = 1;
    MPI_Comm_size(comm, &size);
    std::vector<int64_t> scount(static_cast<size_t>(size)), rcount(static_cast<size_t>(size));
    for (int r = 0; r < size; ++r) {
        scount[static_cast<size_t>(r)] = static_cast<int64_t>(send[static_cast<size_t>(r)].size());
        if (scount[static_cast<size_t>(r)] > INT32_MAX) throw std::runtime_error("partial exchange: a message over 2 GB");
    }
    MPI_Alltoall(scount.data(), 1, MPI_INT64_T, rcount.data(), 1, MPI_INT64_T, comm);
    std::vector<std::vector<char>> recv(static_cast<size_t>(size));
    std::vector<MPI_Request> reqs;
    for (int r = 0; r < size; ++r) {
        if (rcount[static_cast<size_t>(r)] == 0) continue;
        recv[static_cast<size_t>(r)].resize(static_cast<size_t>(rcount[static_cast<size_t>(r)]));
        reqs.emplace_back();
        MPI_Irecv(recv[static_cast<size_t>(r)].data(), static_cast<int>(rcount[static_cast<size_t>(r)]), MPI_BYTE, r, tag, comm,
                  &reqs.back());
    }
    for (int r = 0; r < size; ++r) {
        if (scount[static_cast<size_t>(r)] == 0) continue;
        reqs.emplace_back();
        MPI_Isend(send[static_cast<size_t>(r)].data(), static_cast<int>(scount[static_cast<size_t>(r)]), MPI_BYTE, r, tag, comm,
                  &reqs.back());
    }
    if (!reqs.empty()) MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
    return recv;
}

template<typename T>
void partial_put(std::vector<char>& buf, const T* v, size_t n) {
    const size_t at = buf.size();
    buf.resize(at + n * sizeof(T));
    std::memcpy(buf.data() + at, v, n * sizeof(T));
}
template<typename T>
const char* partial_get(const char* p, T* v, size_t n) {
    std::memcpy(v, p, n * sizeof(T));
    return p + n * sizeof(T);
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
    // Collective over the tree's communicator, every rank passing the same
    // rows and cols, the result replicated.  A rank holds the blocks of its
    // boxes; a level's terms are computed once: a box's own terms and its
    // X_RN terms by its owner (the neighbours' column blocks received), its
    // X_NR terms by the owner of the column box (the neighbours' row blocks
    // received), the root's term by the root's owner; the blocks go to the
    // parents' owner between levels, and a sum assembles the result.
    using Block = PartialChainBlock<DataType>;
    auto& pr = tree.partial_root;
    if (pr.cb_released) throw std::runtime_error("partial Schur: released (c_bpack_partial_schur_release)");
    const int64_t pc = tree.partial_count, n22 = tree.num_points - pc;
    const int dimension = tree.dimension;
    const int leaf = tree.num_levels - 1;
    int me = 0, nranks = 1;
    MPI_Comm_rank(tree.comm, &me);
    MPI_Comm_size(tree.comm, &nranks);
    std::fill(out, out + nrows * ncols, DataType(0));
    if (nrows == 0 || ncols == 0) return;
    if (pr.leaf_box.empty() && n22 > 0) {
        pr.leaf_box.assign(static_cast<size_t>(n22), -1);
        pr.leaf_slot.assign(static_cast<size_t>(n22), -1);
        const auto& level = tree.levels[static_cast<size_t>(leaf)];
        if (level.is_process_active) {
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
    }
    // the two sides: distinct points (the same on every rank), their requested
    // positions, the blocks of this rank's boxes at the current level
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
            if (b < 0) {
                if (nranks == 1) throw std::runtime_error("partial Schur: a point is not in a local leaf");
                continue;  // another rank's
            }
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
    auto scatter = [&](const std::vector<int64_t>& mi, const std::vector<int64_t>& mj, const std::vector<DataType>& G) {
        const int64_t ni = static_cast<int64_t>(mi.size());
        for (size_t c = 0; c < mj.size(); ++c)
            for (int64_t jj : C.positions[static_cast<size_t>(mj[c])])
                for (int64_t r = 0; r < ni; ++r) {
                    const DataType g = G[static_cast<size_t>(r + static_cast<int64_t>(c) * ni)];
                    for (int64_t ii : R.positions[static_cast<size_t>(mi[static_cast<size_t>(r)])])
                        out[ii + jj * nrows] += g;
                }
    };
    // a block of another rank's box: its members, raw rows and values, w part, and the box's live set
    struct Remote {
        int64_t morton = -1, num_points = 0;
        std::vector<int64_t> members, rows, wrows, skeleton;
        std::vector<DataType> v, wv;
        std::vector<int32_t> rowpos, wpos;
    };
    for (int l = leaf; l >= 0; --l) {
        auto& level = tree.levels[static_cast<size_t>(l)];
        const bool active = level.is_process_active;
        if (l == 0) {
            // the root: v_i^T CB_root v_j over the root's F22 slots
            if (pr.active) {
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
                        scatter(R.blocks[i].members, C.blocks[j].members, G);
                    }
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
            // this rank's blocks to the owners of their boxes' remote neighbours
            std::vector<std::vector<char>> send(static_cast<size_t>(nranks));
            if (active) {
                for (int side_id = 0; side_id < 2; ++side_id) {
                    const Side& side = side_id == 0 ? R : C;
                    for (const auto& b : side.blocks) {
                        const auto& box = level.local_boxes[static_cast<size_t>(b.box)];
                        std::set<int> dest;
                        for (int64_t nm : box.one_hop)
                            if (!level.is_box_on_process(nm)) dest.insert(level.owner_rank_of(nm));
                        if (dest.empty()) continue;
                        std::vector<char> msg;
                        const int64_t head[7] = {side_id, box.morton_index, box.num_points, static_cast<int64_t>(b.members.size()),
                                                 static_cast<int64_t>(b.rows.size()), static_cast<int64_t>(b.wrows.size()),
                                                 static_cast<int64_t>(box.skeleton_indices.size())};
                        partial_put(msg, head, 7);
                        partial_put(msg, b.members.data(), b.members.size());
                        partial_put(msg, b.rows.data(), b.rows.size());
                        partial_put(msg, b.v.data(), b.v.size());
                        partial_put(msg, b.wrows.data(), b.wrows.size());
                        partial_put(msg, b.wv.data(), b.wv.size());
                        partial_put(msg, box.skeleton_indices.data(), box.skeleton_indices.size());
                        for (int r : dest) send[static_cast<size_t>(r)].insert(send[static_cast<size_t>(r)].end(), msg.begin(), msg.end());
                    }
                }
            }
            const auto recv = partial_exchange(tree.comm, send, 620);
            std::vector<Remote> remote_rows, remote_cols;
            for (const auto& buf : recv) {
                for (const char* p = buf.data(); p < buf.data() + buf.size();) {
                    int64_t head[7];
                    p = partial_get(p, head, 7);
                    Remote rb;
                    rb.morton = head[1];
                    rb.num_points = head[2];
                    const int64_t nm = head[3], nr = head[4], nw = head[5], kp = head[6];
                    rb.members.resize(static_cast<size_t>(nm));
                    rb.rows.resize(static_cast<size_t>(nr));
                    rb.v.resize(static_cast<size_t>(nr * nm));
                    rb.wrows.resize(static_cast<size_t>(nw));
                    rb.wv.resize(static_cast<size_t>(nw * nm));
                    rb.skeleton.resize(static_cast<size_t>(kp));
                    p = partial_get(p, rb.members.data(), rb.members.size());
                    p = partial_get(p, rb.rows.data(), rb.rows.size());
                    p = partial_get(p, rb.v.data(), rb.v.size());
                    p = partial_get(p, rb.wrows.data(), rb.wrows.size());
                    p = partial_get(p, rb.wv.data(), rb.wv.size());
                    p = partial_get(p, rb.skeleton.data(), rb.skeleton.size());
                    rb.rowpos.assign(static_cast<size_t>(rb.num_points), -1);
                    for (int64_t t = 0; t < nr; ++t) rb.rowpos[static_cast<size_t>(rb.rows[static_cast<size_t>(t)])] = static_cast<int32_t>(t);
                    rb.wpos.assign(static_cast<size_t>(kp), -1);
                    for (int64_t t = 0; t < nw; ++t) rb.wpos[static_cast<size_t>(rb.wrows[static_cast<size_t>(t)])] = static_cast<int32_t>(t);
                    (head[0] == 0 ? remote_rows : remote_cols).push_back(std::move(rb));
                }
            }
            if (active) {
                // the column blocks by box (Morton): this rank's (index) and received (index into remote_cols)
                std::unordered_map<int64_t, std::vector<size_t>> col_by_box, remote_col_by_box;
                for (size_t j = 0; j < C.blocks.size(); ++j)
                    col_by_box[level.local_boxes[static_cast<size_t>(C.blocks[j].box)].morton_index].push_back(j);
                for (size_t j = 0; j < remote_cols.size(); ++j) remote_col_by_box[remote_cols[j].morton].push_back(j);
                // X_RN_bi(ri, bj raw slots in the column block's rows) times its raw values: G += vp_i^T (X V)
                auto rn_term = [&](const Block& bi, const typename fmm::BoxData<CoordType, DataType>::F22Part& fi, int64_t idx_in_i,
                                   const std::vector<int32_t>& rowpos,
                                   const std::vector<DataType>& vj, int64_t nrj, int64_t nj, std::vector<DataType>& G) {
                    const int64_t ni = static_cast<int64_t>(bi.members.size()), ri = static_cast<int64_t>(bi.ridx.size());
                    const int64_t r = static_cast<int64_t>(fi.redundant.size());
                    const int64_t m = fi.neighbor_f22_counts[static_cast<size_t>(idx_in_i)];
                    if (m == 0 || ri == 0) return false;
                    int64_t col0 = 0;
                    for (int64_t q = 0; q < idx_in_i; ++q) col0 += fi.neighbor_f22_counts[static_cast<size_t>(q)];
                    const auto& slots = fi.neighbor_f22_slots[static_cast<size_t>(idx_in_i)];
                    std::vector<int64_t> cpos, tpos;
                    for (int64_t c = 0; c < m; ++c) {
                        const int32_t t = rowpos[static_cast<size_t>(slots[static_cast<size_t>(c)])];
                        if (t >= 0) { cpos.push_back(col0 + c); tpos.push_back(t); }
                    }
                    const int64_t nc = static_cast<int64_t>(cpos.size());
                    if (nc == 0) return false;
                    std::vector<DataType> X(static_cast<size_t>(ri * nc)), V(static_cast<size_t>(nc * nj)), tmp(static_cast<size_t>(ri * nj));
                    for (int64_t c = 0; c < nc; ++c)
                        for (int64_t a = 0; a < ri; ++a)
                            X[static_cast<size_t>(a + c * ri)] = fi.X_RN.data[static_cast<size_t>(bi.ridx[static_cast<size_t>(a)] + cpos[static_cast<size_t>(c)] * r)];
                    for (int64_t c = 0; c < nj; ++c)
                        for (int64_t q = 0; q < nc; ++q) V[static_cast<size_t>(q + c * nc)] = vj[static_cast<size_t>(tpos[static_cast<size_t>(q)] + c * nrj)];
                    partial_gemm<DataType>('N', 'N', ri, nj, nc, DataType(1), X.data(), ri, V.data(), nc, DataType(0), tmp.data(), ri);
                    partial_gemm<DataType>('T', 'N', ni, nj, ri, DataType(1), bi.vp.data(), ri, tmp.data(), ri, DataType(1), G.data(), ni);
                    return true;
                };
                // X_NR_bj(bi raw slots in the row block's w rows, rj) times vp_j: G += wv_i^T (X vp_j)
                auto nr_term = [&](const std::vector<int32_t>& s_of_i, const std::vector<int32_t>& wpos_i, const std::vector<DataType>& wv_i,
                                   int64_t nwi, int64_t ni, int64_t morton_i, const fmm::BoxData<CoordType, DataType>& boxJ, const Block& bj,
                                   std::vector<DataType>& G) {
                    const auto& fj = boxJ.f22;
                    const int64_t nj = static_cast<int64_t>(bj.members.size()), rj = static_cast<int64_t>(bj.ridx.size());
                    const auto itb = std::find(boxJ.one_hop.begin(), boxJ.one_hop.end(), morton_i);
                    if (itb == boxJ.one_hop.end()) return false;
                    const size_t idi = static_cast<size_t>(itb - boxJ.one_hop.begin());
                    const int64_t m = fj.neighbor_f22_counts[idi], total = fj.X_NR.rows;
                    int64_t row0 = 0;
                    for (size_t q = 0; q < idi; ++q) row0 += fj.neighbor_f22_counts[q];
                    const auto& slots = fj.neighbor_f22_slots[idi];
                    std::vector<int64_t> crow, wposl;
                    for (int64_t c = 0; c < m; ++c) {
                        const int32_t sp = s_of_i[static_cast<size_t>(slots[static_cast<size_t>(c)])];
                        if (sp < 0) continue;  // in R_bi: the X_RN term
                        const int32_t t = wpos_i[static_cast<size_t>(sp)];
                        if (t >= 0) { crow.push_back(row0 + c); wposl.push_back(t); }
                    }
                    const int64_t nc = static_cast<int64_t>(crow.size());
                    if (nc == 0) return false;
                    std::vector<DataType> X(static_cast<size_t>(nc * rj)), V(static_cast<size_t>(nc * ni)), tmp(static_cast<size_t>(nc * nj));
                    for (int64_t c = 0; c < rj; ++c)
                        for (int64_t q = 0; q < nc; ++q)
                            X[static_cast<size_t>(q + c * nc)] = fj.X_NR.data[static_cast<size_t>(crow[static_cast<size_t>(q)] + bj.ridx[static_cast<size_t>(c)] * total)];
                    for (int64_t c = 0; c < ni; ++c)
                        for (int64_t q = 0; q < nc; ++q) V[static_cast<size_t>(q + c * nc)] = wv_i[static_cast<size_t>(wposl[static_cast<size_t>(q)] + c * nwi)];
                    partial_gemm<DataType>('N', 'N', nc, nj, rj, DataType(1), X.data(), nc, bj.vp.data(), rj, DataType(0), tmp.data(), nc);
                    partial_gemm<DataType>('T', 'N', ni, nj, nc, DataType(1), V.data(), nc, tmp.data(), nc, DataType(1), G.data(), ni);
                    return true;
                };
                // this rank's row blocks: own terms, X_RN terms with any column block, X_NR terms with this rank's
#pragma omp parallel for schedule(dynamic, 1)
                for (size_t i = 0; i < R.blocks.size(); ++i) {
                    const Block& bi = R.blocks[i];
                    const auto& boxI = level.local_boxes[static_cast<size_t>(bi.box)];
                    const auto& fi = boxI.f22;
                    const bool own_i = fi.sparsified && !fi.redundant.empty() && !bi.ridx.empty();
                    const int64_t ni = static_cast<int64_t>(bi.members.size()), ri = static_cast<int64_t>(bi.ridx.size());
                    const auto& s_of_i = maps.at(bi.box).first;
                    std::vector<DataType> X, tmp, G;
                    // the same box
                    if (const auto it = col_by_box.find(boxI.morton_index); it != col_by_box.end()) {
                        const int64_t r = static_cast<int64_t>(fi.redundant.size()), kp = static_cast<int64_t>(boxI.skeleton_indices.size());
                        for (size_t jb : it->second) {
                            const Block& bj = C.blocks[jb];
                            const int64_t nj = static_cast<int64_t>(bj.members.size()), rj = static_cast<int64_t>(bj.ridx.size());
                            G.assign(static_cast<size_t>(ni * nj), DataType(0));
                            bool any = false;
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
                            if (any) scatter(bi.members, bj.members, G);
                        }
                    }
                    // the one-hop neighbours
                    for (size_t q = 0; q < boxI.one_hop.size(); ++q) {
                        const int64_t nm = boxI.one_hop[q];
                        if (level.is_box_on_process(nm)) {
                            const auto it = col_by_box.find(nm);
                            if (it == col_by_box.end()) continue;
                            const auto& boxJ = level.local_boxes[static_cast<size_t>(nm - level.local_morton_start)];
                            const bool sp_j = boxJ.f22.sparsified && !boxJ.f22.redundant.empty();
                            for (size_t jb : it->second) {
                                const Block& bj = C.blocks[jb];
                                const int64_t nj = static_cast<int64_t>(bj.members.size());
                                G.assign(static_cast<size_t>(ni * nj), DataType(0));
                                bool any = false;
                                if (own_i) any |= rn_term(bi, fi, static_cast<int64_t>(q), bj.rowpos, bj.v, static_cast<int64_t>(bj.rows.size()), nj, G);
                                if (sp_j && !bj.ridx.empty() && !bi.wrows.empty())
                                    any |= nr_term(s_of_i, bi.wpos, bi.wv, static_cast<int64_t>(bi.wrows.size()), ni, boxI.morton_index, boxJ, bj, G);
                                if (any) scatter(bi.members, bj.members, G);
                            }
                        } else if (own_i) {
                            // a remote column box: its X_NR terms are its owner's
                            const auto it = remote_col_by_box.find(nm);
                            if (it == remote_col_by_box.end()) continue;
                            for (size_t jb : it->second) {
                                const Remote& bj = remote_cols[jb];
                                const int64_t nj = static_cast<int64_t>(bj.members.size());
                                G.assign(static_cast<size_t>(ni * nj), DataType(0));
                                if (rn_term(bi, fi, static_cast<int64_t>(q), bj.rowpos, bj.v, static_cast<int64_t>(bj.rows.size()), nj, G))
                                    scatter(bi.members, bj.members, G);
                            }
                        }
                    }
                }
                // received row blocks: their X_NR terms with this rank's column blocks
#pragma omp parallel for schedule(dynamic, 1)
                for (size_t i = 0; i < remote_rows.size(); ++i) {
                    const Remote& bi = remote_rows[i];
                    if (bi.wrows.empty()) continue;
                    std::vector<int32_t> s_of_i(static_cast<size_t>(bi.num_points), -1);
                    for (size_t t = 0; t < bi.skeleton.size(); ++t) s_of_i[static_cast<size_t>(bi.skeleton[t])] = static_cast<int32_t>(t);
                    const int64_t ni = static_cast<int64_t>(bi.members.size());
                    std::vector<DataType> G;
                    // the column boxes of this rank among the row box's neighbours (they list it)
                    for (const auto& kv : col_by_box) {
                        const auto& boxJ = level.local_boxes[static_cast<size_t>(kv.first - level.local_morton_start)];
                        if (!(boxJ.f22.sparsified && !boxJ.f22.redundant.empty())) continue;
                        if (std::find(boxJ.one_hop.begin(), boxJ.one_hop.end(), bi.morton) == boxJ.one_hop.end()) continue;
                        for (size_t jb : kv.second) {
                            const Block& bj = C.blocks[jb];
                            if (bj.ridx.empty()) continue;
                            G.assign(static_cast<size_t>(ni * static_cast<int64_t>(bj.members.size())), DataType(0));
                            if (nr_term(s_of_i, bi.wpos, bi.wv, static_cast<int64_t>(bi.wrows.size()), ni, bi.morton, boxJ, bj, G))
                                scatter(bi.members, bj.members, G);
                        }
                    }
                }
            }
        }
        // the parent blocks: each block's w part at its box's offset; the dense
        // parts of one box's blocks share their rows and merge; blocks whose
        // parent is another rank's go to it
        auto& parent_level = tree.levels[static_cast<size_t>(l) - 1];
        std::vector<std::vector<char>> send_up(static_cast<size_t>(nranks));
        for (int side_id = 0; side_id < 2; ++side_id) {
            Side* side = side_id == 0 ? &R : &C;
            std::vector<Block> up;
            std::vector<int64_t> up_parent;  // the parents' Morton indices
            std::unordered_map<int64_t, size_t> dense_of;  // box at this level -> merged block at the parent
            for (auto& b : side->blocks) {
                const auto& box = level.local_boxes[static_cast<size_t>(b.box)];
                if (b.wrows.empty()) continue;  // nothing left of these chains above
                const int64_t morton = box.morton_index;
                const int64_t parent = morton >> dimension;
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
                u.members = std::move(b.members);
                u.rows.resize(static_cast<size_t>(nw));
                for (int64_t t = 0; t < nw; ++t) u.rows[static_cast<size_t>(t)] = off + b.wrows[static_cast<size_t>(t)];
                u.v = std::move(b.wv);
                up.push_back(std::move(u));
                up_parent.push_back(parent);
            }
            side->blocks.clear();
            for (size_t k = 0; k < up.size(); ++k) {
                const int64_t parent = up_parent[k];
                if (parent_level.is_process_active && parent_level.is_box_on_process(parent)) {
                    up[k].box = parent - parent_level.local_morton_start;
                    side->blocks.push_back(std::move(up[k]));
                    continue;
                }
                const int dest = level.parent_level_owner;
                if (dest < 0) throw std::runtime_error("partial Schur: no owner of a parent box");
                auto& buf = send_up[static_cast<size_t>(dest)];
                const int64_t head[4] = {side_id, parent, static_cast<int64_t>(up[k].members.size()), static_cast<int64_t>(up[k].rows.size())};
                partial_put(buf, head, 4);
                partial_put(buf, up[k].members.data(), up[k].members.size());
                partial_put(buf, up[k].rows.data(), up[k].rows.size());
                partial_put(buf, up[k].v.data(), up[k].v.size());
            }
        }
        const auto recv_up = partial_exchange(tree.comm, send_up, 621);
        for (const auto& buf : recv_up) {
            for (const char* p = buf.data(); p < buf.data() + buf.size();) {
                int64_t head[4];
                p = partial_get(p, head, 4);
                Block u;
                u.box = head[1] - parent_level.local_morton_start;
                u.members.resize(static_cast<size_t>(head[2]));
                u.rows.resize(static_cast<size_t>(head[3]));
                u.v.resize(static_cast<size_t>(head[2] * head[3]));
                p = partial_get(p, u.members.data(), u.members.size());
                p = partial_get(p, u.rows.data(), u.rows.size());
                p = partial_get(p, u.v.data(), u.v.size());
                (head[0] == 0 ? R : C).blocks.push_back(std::move(u));
            }
        }
    }
    if (nranks > 1) MPI_Allreduce(MPI_IN_PLACE, out, static_cast<int>(nrows * ncols), partial_mpi_type<DataType>(), MPI_SUM, tree.comm);
}

// y2 = CB x2 = P22^T CB' P22 x2 (header): an up-pass with P (w = x_S' +
// T x_R, the raw vectors kept), the sparsified blocks and the root's block
// applied level by level, a down-pass with P^T (a_S' = g_S' + a_parent,
// a_R = g_R + T^T a_S').  x2, y2 column-major n22 x ncol in F22 index order
// (replicated); collective over the tree's communicator.  The passes go
// through the solve's children-to-parent gather and parent-to-children
// scatter; at a level, a box's raw vector goes to its remote neighbours'
// owners (the X_RN terms), and a sparsified box's X_NR rows of a remote
// neighbour, applied to its p, go to that neighbour's owner (the records
// stay where they are).  Every pass runs over a rank's boxes in parallel.
template<typename CoordType, typename DataType>
void partial_schur_multiply(fmm::ParallelTree<CoordType, DataType>& tree, const DataType* x2, DataType* y2,
                            int64_t ncol) {
    const auto& pr = tree.partial_root;
    if (pr.cb_released) throw std::runtime_error("partial Schur: released (c_bpack_partial_schur_release)");
    const int64_t pc = tree.partial_count, n22 = tree.num_points - pc;
    const int dimension = tree.dimension;
    const int leaf = tree.num_levels - 1;
    if (ncol <= 0) return;
    const int nrhs = static_cast<int>(ncol);
    int me = 0, nranks = 1;
    MPI_Comm_rank(tree.comm, &me);
    MPI_Comm_size(tree.comm, &nranks);
    auto own = [](const fmm::BoxData<CoordType, DataType>& box) {
        return box.f22.sparsified && !box.f22.redundant.empty();
    };
    using SD = fmm::SolveDataRequest<CoordType, DataType>;
    // raw: per level and local box, its points x ncol; pv: p = raw[R]
    auto raw = partial_scratch<CoordType, DataType>(&tree, nullptr, nrhs);
    std::vector<std::vector<std::vector<DataType>>> pv(static_cast<size_t>(tree.num_levels));
    {
        auto& level = tree.levels[static_cast<size_t>(leaf)];
        if (level.is_process_active) {
#pragma omp parallel for schedule(dynamic, 16)
            for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
                const auto& box = level.local_boxes[static_cast<size_t>(bi)];
                auto& v = raw[static_cast<size_t>(leaf)][static_cast<size_t>(bi)].left_side;
                for (int64_t s = 0; s < box.num_points; ++s) {
                    const int64_t id = box.point_indices[static_cast<size_t>(s)];
                    if (id < pc) continue;
                    for (int64_t c = 0; c < ncol; ++c)
                        v[static_cast<size_t>(s + c * box.num_points)] = x2[static_cast<size_t>((id - pc) + c * n22)];
                }
            }
        }
    }
    // w = raw[S'] + T p for every local box of a level (column-major kp x ncol)
    auto w_of = [&](const fmm::BoxData<CoordType, DataType>& box, const std::vector<DataType>& x, const std::vector<DataType>* p) {
        const int64_t n = box.num_points, kp = static_cast<int64_t>(box.skeleton_indices.size());
        std::vector<DataType> w(static_cast<size_t>(kp * ncol));
        for (int64_t c = 0; c < ncol; ++c)
            for (int64_t i = 0; i < kp; ++i)
                w[static_cast<size_t>(i + c * kp)] = x[static_cast<size_t>(box.skeleton_indices[static_cast<size_t>(i)] + c * n)];
        if (p != nullptr) {
            const int64_t r = static_cast<int64_t>(box.f22.redundant.size());
            partial_gemm<DataType>('N', 'N', kp, ncol, r, DataType(1), box.f22.T.data.data(), kp, p->data(), r, DataType(1), w.data(), kp);
        }
        return w;
    };
    // up: p = raw[R], and the records with S' = w gathered into the parents' raw vectors
    for (int l = leaf; l >= 1; --l) {
        auto& level = tree.levels[static_cast<size_t>(l)];
        std::vector<SD> up;
        if (level.is_process_active) {
            pv[static_cast<size_t>(l)].assign(static_cast<size_t>(level.num_boxes_local), {});
            up.resize(static_cast<size_t>(level.num_boxes_local));
#pragma omp parallel for schedule(dynamic, 1)
            for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
                const auto& box = level.local_boxes[static_cast<size_t>(bi)];
                const auto& x = raw[static_cast<size_t>(l)][static_cast<size_t>(bi)].left_side;
                up[static_cast<size_t>(bi)].initialize(box.morton_index, me, box.num_points, nrhs);
                auto& u = up[static_cast<size_t>(bi)].left_side;
                u = x;
                if (box.skeleton_indices.empty()) continue;
                const std::vector<DataType>* p = nullptr;
                if (own(box)) {
                    const int64_t r = static_cast<int64_t>(box.f22.redundant.size()), n = box.num_points;
                    auto& pb = pv[static_cast<size_t>(l)][static_cast<size_t>(bi)];
                    pb.resize(static_cast<size_t>(r * ncol));
                    for (int64_t c = 0; c < ncol; ++c)
                        for (int64_t j = 0; j < r; ++j)
                            pb[static_cast<size_t>(j + c * r)] = x[static_cast<size_t>(box.f22.redundant[static_cast<size_t>(j)] + c * n)];
                    p = &pb;
                }
                const auto w = w_of(box, x, p);
                const int64_t n = box.num_points, kp = static_cast<int64_t>(box.skeleton_indices.size());
                for (int64_t c = 0; c < ncol; ++c)
                    for (int64_t i = 0; i < kp; ++i)
                        u[static_cast<size_t>(box.skeleton_indices[static_cast<size_t>(i)] + c * n)] = w[static_cast<size_t>(i + c * kp)];
            }
        }
        fmm::gather_skeleton_to_parent(level, tree.levels[static_cast<size_t>(l) - 1], up, raw[static_cast<size_t>(l) - 1],
                                       dimension, tree.comm);
    }
    // h: the accumulators; the root's block on the root's owner
    auto h = partial_scratch<CoordType, DataType>(&tree, nullptr, nrhs);
    if (pr.active && pr.n22 > 0) {
        const auto& box = tree.levels[0].local_boxes[0];
        const auto& x = raw[0][0].left_side;
        auto& out = h[0][0].left_side;
        const int64_t n = box.num_points, m = pr.n22;
        std::vector<DataType> v(static_cast<size_t>(m * ncol)), g(static_cast<size_t>(m * ncol));
        for (int64_t c = 0; c < ncol; ++c)
            for (int64_t a = 0; a < m; ++a) v[static_cast<size_t>(a + c * m)] = x[static_cast<size_t>(pr.slots22[static_cast<size_t>(a)] + c * n)];
        constexpr int64_t chunk = 256;
        const int64_t nchunks = (m + chunk - 1) / chunk;
#pragma omp parallel for schedule(static)
        for (int64_t q = 0; q < nchunks; ++q) {
            const int64_t a0 = q * chunk, rows = std::min(chunk, m - a0);
            partial_gemm<DataType>('N', 'N', rows, ncol, m, DataType(1), pr.cb.data.data() + a0, pr.cb.lda, v.data(), m, DataType(0),
                                   g.data() + a0, m);
        }
        for (int64_t c = 0; c < ncol; ++c)
            for (int64_t a = 0; a < m; ++a) out[static_cast<size_t>(pr.slots22[static_cast<size_t>(a)] + c * n)] = g[static_cast<size_t>(a + c * m)];
    }
    // the records, level by level
    for (int l = leaf; l >= 2; --l) {
        auto& level = tree.levels[static_cast<size_t>(l)];
        const bool active = level.is_process_active;
        auto is_local = [&](int64_t m) { return active && level.is_box_on_process(m); };
        // (1) raw vectors to the owners of remote neighbours
        std::vector<std::vector<char>> send_raw(static_cast<size_t>(nranks));
        if (active) {
            for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
                const auto& box = level.local_boxes[static_cast<size_t>(bi)];
                if (box.skeleton_indices.empty()) continue;
                std::set<int> dest;
                for (int64_t nm : box.one_hop)
                    if (!level.is_box_on_process(nm)) dest.insert(level.owner_rank_of(nm));
                const auto& x = raw[static_cast<size_t>(l)][static_cast<size_t>(bi)].left_side;
                const int64_t head[2] = {box.morton_index, box.num_points};
                for (int r : dest) {
                    partial_put(send_raw[static_cast<size_t>(r)], head, 2);
                    partial_put(send_raw[static_cast<size_t>(r)], x.data(), x.size());
                }
            }
        }
        const auto recv_raw = partial_exchange(tree.comm, send_raw, 610);
        std::unordered_map<int64_t, std::pair<int64_t, std::vector<DataType>>> remote_raw;
        for (const auto& buf : recv_raw) {
            for (const char* p = buf.data(); p < buf.data() + buf.size();) {
                int64_t head[2];
                p = partial_get(p, head, 2);
                std::vector<DataType> x(static_cast<size_t>(head[1] * ncol));
                p = partial_get(p, x.data(), x.size());
                remote_raw[head[0]] = {head[1], std::move(x)};
            }
        }
        // (2) a sparsified box's X_NR rows of a remote neighbour applied to its p, to that neighbour's owner
        std::vector<std::vector<char>> send_nr(static_cast<size_t>(nranks));
        if (active) {
            std::vector<std::vector<std::pair<int, std::vector<char>>>> parts(static_cast<size_t>(level.num_boxes_local));
#pragma omp parallel for schedule(dynamic, 1)
            for (int64_t ni = 0; ni < level.num_boxes_local; ++ni) {
                const auto& nb = level.local_boxes[static_cast<size_t>(ni)];
                if (!own(nb)) continue;
                const auto& fn = nb.f22;
                const auto& pn = pv[static_cast<size_t>(l)][static_cast<size_t>(ni)];
                const int64_t rn = static_cast<int64_t>(fn.redundant.size());
                int64_t row0 = 0;
                for (size_t idx = 0; idx < nb.one_hop.size(); ++idx) {
                    const int64_t m = fn.neighbor_f22_counts[idx];
                    const int64_t target = nb.one_hop[idx];
                    if (m > 0 && !level.is_box_on_process(target)) {
                        std::vector<DataType> t(static_cast<size_t>(m * ncol));
                        partial_gemm<DataType>('N', 'N', m, ncol, rn, DataType(1), fn.X_NR.data.data() + row0, fn.X_NR.rows, pn.data(), rn,
                                               DataType(0), t.data(), m);
                        std::vector<char> msg;
                        const int64_t head[2] = {target, m};
                        partial_put(msg, head, 2);
                        partial_put(msg, fn.neighbor_f22_slots[idx].data(), static_cast<size_t>(m));
                        partial_put(msg, t.data(), t.size());
                        parts[static_cast<size_t>(ni)].emplace_back(level.owner_rank_of(target), std::move(msg));
                    }
                    row0 += m;
                }
            }
            for (auto& list : parts)
                for (auto& [r, msg] : list) send_nr[static_cast<size_t>(r)].insert(send_nr[static_cast<size_t>(r)].end(), msg.begin(), msg.end());
        }
        const auto recv_nr = partial_exchange(tree.comm, send_nr, 611);
        std::unordered_map<int64_t, std::vector<std::pair<std::vector<int64_t>, std::vector<DataType>>>> pushed;
        for (const auto& buf : recv_nr) {
            for (const char* p = buf.data(); p < buf.data() + buf.size();) {
                int64_t head[2];
                p = partial_get(p, head, 2);
                std::vector<int64_t> slots(static_cast<size_t>(head[1]));
                p = partial_get(p, slots.data(), slots.size());
                std::vector<DataType> t(static_cast<size_t>(head[1] * ncol));
                p = partial_get(p, t.data(), t.size());
                pushed[head[0]].emplace_back(std::move(slots), std::move(t));
            }
        }
        if (!active) continue;
        // (3) every local box's terms
        std::exception_ptr failure;
#pragma omp parallel for schedule(dynamic, 1)
        for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
            try {
                const auto& box = level.local_boxes[static_cast<size_t>(bi)];
                const auto& f = box.f22;
                if (box.skeleton_indices.empty()) continue;
                const bool mine = own(box);
                const int64_t n = box.num_points, kp = static_cast<int64_t>(box.skeleton_indices.size());
                const int64_t r = mine ? static_cast<int64_t>(f.redundant.size()) : 0;
                const auto& x = raw[static_cast<size_t>(l)][static_cast<size_t>(bi)].left_side;
                std::vector<DataType> gR(static_cast<size_t>(r * ncol), DataType(0)), gS(static_cast<size_t>(kp * ncol), DataType(0));
                if (mine) {
                    const auto& p = pv[static_cast<size_t>(l)][static_cast<size_t>(bi)];
                    const auto w = w_of(box, x, &p);
                    partial_gemm<DataType>('N', 'N', r, ncol, r, DataType(1), f.X_RR.data.data(), r, p.data(), r, DataType(1), gR.data(), r);
                    partial_gemm<DataType>('N', 'N', r, ncol, kp, DataType(1), f.X_RS.data.data(), r, w.data(), kp, DataType(1), gR.data(), r);
                    partial_gemm<DataType>('N', 'N', kp, ncol, r, DataType(1), f.X_SR.data.data(), kp, p.data(), r, DataType(1), gS.data(), kp);
                    // the neighbours' raw vectors (local, or received) against X_RN
                    const int64_t total = f.X_RN.cols;
                    if (total > 0) {
                        std::vector<DataType> rn(static_cast<size_t>(total * ncol));
                        int64_t col0 = 0;
                        for (size_t idx = 0; idx < box.one_hop.size(); ++idx) {
                            const int64_t m = f.neighbor_f22_counts[idx];
                            if (m == 0) continue;
                            const int64_t nm = box.one_hop[idx];
                            const std::vector<DataType>* xn = nullptr;
                            int64_t nn = 0;
                            if (is_local(nm)) {
                                const int64_t nl = nm - level.local_morton_start;
                                xn = &raw[static_cast<size_t>(l)][static_cast<size_t>(nl)].left_side;
                                nn = level.local_boxes[static_cast<size_t>(nl)].num_points;
                            } else {
                                const auto it = remote_raw.find(nm);
                                if (it == remote_raw.end()) throw std::runtime_error("partial Schur product: no vector of neighbour " + std::to_string(nm));
                                nn = it->second.first;
                                xn = &it->second.second;
                            }
                            const auto& slots = f.neighbor_f22_slots[idx];
                            for (int64_t c = 0; c < ncol; ++c)
                                for (int64_t i = 0; i < m; ++i)
                                    rn[static_cast<size_t>(col0 + i + c * total)] = (*xn)[static_cast<size_t>(slots[static_cast<size_t>(i)] + c * nn)];
                            col0 += m;
                        }
                        partial_gemm<DataType>('N', 'N', r, ncol, total, DataType(1), f.X_RN.data.data(), r, rn.data(), total, DataType(1),
                                               gR.data(), r);
                    }
                }
                // the sparsified neighbours' p against their X_NR rows over this box's S' (local: here; remote: received)
                std::vector<int32_t> s_of(static_cast<size_t>(n), -1);
                for (int64_t i = 0; i < kp; ++i) s_of[static_cast<size_t>(box.skeleton_indices[static_cast<size_t>(i)])] = static_cast<int32_t>(i);
                auto add_rows = [&](const std::vector<int64_t>& slots, const DataType* t, int64_t m) {
                    for (int64_t i = 0; i < m; ++i) {
                        const int32_t si = s_of[static_cast<size_t>(slots[static_cast<size_t>(i)])];
                        if (si < 0) continue;  // this box's R rows: in X_RN above
                        for (int64_t c = 0; c < ncol; ++c) gS[static_cast<size_t>(si + c * kp)] += t[static_cast<size_t>(i + c * m)];
                    }
                };
                std::vector<DataType> t;
                for (size_t idx = 0; idx < box.one_hop.size(); ++idx) {
                    const int64_t nm = box.one_hop[idx];
                    if (!is_local(nm)) continue;
                    const int64_t nl = nm - level.local_morton_start;
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
                    add_rows(fn.neighbor_f22_slots[idn], t.data(), m);
                }
                if (const auto it = pushed.find(box.morton_index); it != pushed.end())
                    for (const auto& [slots, tv] : it->second) add_rows(slots, tv.data(), static_cast<int64_t>(slots.size()));
                auto& out = h[static_cast<size_t>(l)][static_cast<size_t>(bi)].left_side;
                for (int64_t c = 0; c < ncol; ++c) {
                    for (int64_t i = 0; i < kp; ++i) out[static_cast<size_t>(box.skeleton_indices[static_cast<size_t>(i)] + c * n)] += gS[static_cast<size_t>(i + c * kp)];
                    for (int64_t j = 0; j < r; ++j) out[static_cast<size_t>(f.redundant[static_cast<size_t>(j)] + c * n)] += gR[static_cast<size_t>(j + c * r)];
                }
            } catch (...) {
#pragma omp critical(partial_schur_multiply_failure)
                if (!failure) failure = std::current_exception();
            }
        }
        if (failure) std::rethrow_exception(failure);
    }
    // down: a_S' = g_S' + a_parent (the scatter puts the parent's values at the children's S'), a_R = g_R + T^T a_S'
    for (int l = 1; l <= leaf; ++l) {
        auto& level = tree.levels[static_cast<size_t>(l)];
        std::vector<SD> from_parent;
        if (level.is_process_active) {
            from_parent.resize(static_cast<size_t>(level.num_boxes_local));
            for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
                const auto& box = level.local_boxes[static_cast<size_t>(bi)];
                from_parent[static_cast<size_t>(bi)].initialize(box.morton_index, me, box.num_points, nrhs);
            }
        }
        fmm::scatter_solution_to_children(level, tree.levels[static_cast<size_t>(l) - 1], from_parent, h[static_cast<size_t>(l) - 1],
                                          dimension, tree.comm);
        if (!level.is_process_active) continue;
#pragma omp parallel for schedule(dynamic, 1)
        for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
            const auto& box = level.local_boxes[static_cast<size_t>(bi)];
            if (box.skeleton_indices.empty()) continue;
            auto& a = h[static_cast<size_t>(l)][static_cast<size_t>(bi)].left_side;
            const auto& pa = from_parent[static_cast<size_t>(bi)].left_side;
            const int64_t n = box.num_points, kp = static_cast<int64_t>(box.skeleton_indices.size());
            std::vector<DataType> aS(static_cast<size_t>(kp * ncol));
            for (int64_t c = 0; c < ncol; ++c)
                for (int64_t i = 0; i < kp; ++i) {
                    const size_t q = static_cast<size_t>(box.skeleton_indices[static_cast<size_t>(i)] + c * n);
                    a[q] += pa[q];
                    aS[static_cast<size_t>(i + c * kp)] = a[q];
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
    std::fill(y2, y2 + n22 * ncol, DataType(0));
    {
        auto& level = tree.levels[static_cast<size_t>(leaf)];
        if (level.is_process_active) {
#pragma omp parallel for schedule(dynamic, 16)
            for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
                const auto& box = level.local_boxes[static_cast<size_t>(bi)];
                const auto& a = h[static_cast<size_t>(leaf)][static_cast<size_t>(bi)].left_side;
                for (int64_t s = 0; s < box.num_points; ++s) {
                    const int64_t id = box.point_indices[static_cast<size_t>(s)];
                    if (id < pc) continue;
                    for (int64_t c = 0; c < ncol; ++c)
                        y2[static_cast<size_t>((id - pc) + c * n22)] = a[static_cast<size_t>(s + c * box.num_points)];
                }
            }
        }
    }
    if (nranks > 1) MPI_Allreduce(MPI_IN_PLACE, y2, static_cast<int>(n22 * ncol), partial_mpi_type<DataType>(), MPI_SUM, tree.comm);
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
