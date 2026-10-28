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
#include "core/memory_diagnostics.hpp"
#include "core/partial_record_codec.hpp"
#include "core/morton.hpp"
#include "core/solver.hpp"
#ifdef H2_HAVE_GPU
#include "h2_gpu/partial_device.hpp"
#endif

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <functional>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
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

// The records kept on the device only (PartialRecord::device_only) given
// host copies, for the paths that read them on the host (the host CB entries
// and product, H2_partial_extract host, the BPACK_CHECK comparisons) and
// before an eviction of the device data (drop_device: their device pointers
// go, the store frees the blocks).  Local; nothing without kept records.
template<typename CoordType, typename DataType>
void partial_records_to_host(fmm::ParallelTree<CoordType, DataType>& tree, bool drop_device) {
#ifdef H2_HAVE_GPU
    if (!tree.partial_records_on_device) return;
    for (auto& level : tree.levels) {
        for (auto& box : level.local_boxes) {
            auto& f = box.f22;
            if (!f.sparsified) continue;
            auto fetch = [&](fmm::PartialRecord<DataType>& m) {
                if (m.device_only()) {
                    if (m.compressed()) {
                        size_t n = 0;
                        for (uint32_t st : m.strides) n += st;
                        m.packed.resize(n);
                        fmm::gpu::partial_download(m.packed.data(), m.device, n);
                    } else {
                        m.data.resize(static_cast<size_t>(m.rows * m.cols));
                        fmm::gpu::partial_download(m.data.data(), m.device, m.data.size() * sizeof(DataType));
                    }
                }
                if (drop_device) m.device = nullptr;
            };
            for (auto* m : {&f.X_SR, &f.X_RS, &f.X_RR, &f.X_RN}) fetch(*m);
            for (auto& m : f.X_NR) fetch(m);
        }
    }
    if (drop_device) tree.partial_records_on_device = false;
#else
    (void)tree;
    (void)drop_device;
#endif
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
// Two passes: every box's ID and own records (T, X_SR, X_RS, X_RR, X_RN),
// then - the kept sets known and published to the ranks holding copies -
// every box's X_NR, which keeps only the rows of its neighbours' F22 points
// that stay live: a sparsified neighbour's R22 rows would never be read,
// those couplings are in that neighbour's X_RN (at 64^3, 17.5 of the 49.5 GB
// of records).  With `ids` and `records` (H2_use_gpu: LevelEliminator::
// f22_ids, f22_records) the passes run there from the plans made here,
// the own records' copies to the host behind the IDs; `blocks_back` then has the
// level's blocks sliced to the new live sets (the boxes the step
// sparsified, this rank's and its copies') where they are: the host takes
// nothing from them.
// =====================================================================
template<typename CoordType, typename DataType, typename KernelType>
void partial_f22_phase(fmm::ParallelTree<CoordType, DataType>* tree, fmm::TreeLevel<CoordType, DataType>& level,
                       KernelType* kernel, double tolerance, const std::vector<int64_t>& occupied_local,
                       bool print_detail, MPI_Comm level_comm,
                       const std::function<void(std::vector<fmm::PartialF22Work<DataType>>&)>& ids = {},
                       const std::function<void(std::vector<fmm::PartialF22Work<DataType>>&)>& records = {},
                       const std::function<void(const fmm::PartialF22Kept&)>& blocks_back = {}) {
    using Box = fmm::BoxData<CoordType, DataType>;
    using MS = fmm::MatrixStorage<DataType>;
    const double abs_tol = fmm::id_absolute_tolerance(tree, tolerance);
    const bool device = static_cast<bool>(ids);
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
    // a neighbour's live set before the phase and its point ids: a local
    // box's, or its remote copy's (a ghost on a CA level, an assisting box on
    // a Color one), that one copied now: the publish between the passes
    // reduces the copies' live sets
    struct Old {
        const std::vector<int64_t>* S_old = nullptr;
        const std::vector<int64_t>* ids = nullptr;
    };
    std::unordered_map<int64_t, std::vector<int64_t>> copies_old;
    std::unordered_map<int64_t, Old> old_of;
    for (int64_t li : occupied_local) {
        const Box& b = level.local_boxes[static_cast<size_t>(li)];
        for (int64_t nm : b.one_hop) {
            if (old_of.count(nm) != 0) continue;
            Old o;
            if (const Box* nb = level.find_local_box(nm)) {
                o = Old{&nb->skeleton_before_f22, &nb->point_indices};
            } else if (const Box* g = level.find_ghost_box(nm)) {
                o = Old{&(copies_old[nm] = g->skeleton_indices), &g->point_indices};
            } else {
                auto ait = level.assisting_box_points_for_kernel_evaluation.find(nm);
                if (ait == level.assisting_box_points_for_kernel_evaluation.end()) {
                    throw std::runtime_error("F22 phase: neighbour " + std::to_string(nm) + " of box " +
                                             std::to_string(b.morton_index) + " has no local or remote copy");
                }
                const auto& a = level.assisting_boxes[static_cast<size_t>(ait->second)];
                o = Old{&(copies_old[nm] = a.skel_indices), &a.indices};
            }
            if (o.S_old->empty()) {
                throw std::runtime_error("F22 phase: neighbour " + std::to_string(nm) + " of box " +
                                         std::to_string(b.morton_index) + " has no published skeleton");
            }
            old_of.emplace(nm, o);
        }
    }
    std::exception_ptr failure;
    int64_t sparsified_boxes = 0, sparsified_points = 0, kept_points = 0;
    // BPACK_TRACE=phase: the step's parts, thread-seconds summed over the boxes
    double t_sketch = 0.0, t_target = 0.0, t_id = 0.0, t_records = 0.0;
    const double t_loop0 = MPI_Wtime();
    // the boxes' kept positions within S_old; the shared state (live sets,
    // the boxes' own blocks) is changed after the passes, they only read
    // the neighbours and rings
    std::vector<std::vector<int64_t>> kept(occupied_local.size());
    // morton -> the old live set's size and the kept positions in it: this
    // rank's sparsified boxes, then (several ranks) the other ranks' boxes
    // this rank holds copies of
    fmm::PartialF22Kept kept_of;
    auto publish_kept = [&] {
        for (size_t oi = 0; oi < occupied_local.size(); ++oi) {
            const Box& b = level.local_boxes[static_cast<size_t>(occupied_local[oi])];
            if (b.f22.sparsified) kept_of.emplace(b.morton_index, std::make_pair(static_cast<int64_t>(b.skeleton_before_f22.size()), kept[oi]));
        }
        if (tree->mpi_size > 1) partial_f22_publish(tree, level, level_comm, kept_of);
    };
    // X_NR's rows of a neighbour: the F22 positions it keeps live (all of
    // them when the step did not sparsify it)
    std::unordered_map<int64_t, std::vector<char>> live_of;
    auto live_rows = [&](int64_t nm, const std::vector<int64_t>& q22) {
        const auto it = live_of.find(nm);
        if (it == live_of.end()) return q22;
        std::vector<int64_t> out;
        for (int64_t q : q22)
            if (it->second[static_cast<size_t>(q)] != 0) out.push_back(q);
        return out;
    };
    auto index_live = [&] {
        for (const auto& [nm, kv] : kept_of) {
            auto& live = live_of[nm];
            live.assign(static_cast<size_t>(kv.first), 0);
            for (int64_t q : kv.second) live[static_cast<size_t>(q)] = 1;
        }
    };
    bool sliced_on_device = false;
    if (device) {
        // ---- the plans; the device takes the two passes
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
                const Old& o = old_of.at(nm);
                auto& nr = w.near[idx];
                nr.morton = nm;
                nr.S_old = o.S_old;
                for (int64_t q = 0; q < static_cast<int64_t>(o.S_old->size()); ++q) {
                    (tree->partial_f11((*o.ids)[static_cast<size_t>((*o.S_old)[static_cast<size_t>(q)])]) ? nr.q11 : nr.q22)
                        .push_back(q);
                }
            }
            works.push_back(std::move(w));
            work_box.push_back(oi);
        }
        const double t_ids0 = MPI_Wtime();
        t_target += t_ids0 - t_loop0;  // (the plans)
        ids(works);
        const double t_kept0 = MPI_Wtime();
        t_sketch += t_kept0 - t_ids0;  // (the device's passes)
        // ---- the kept sets, published, and the records' X_NR rows
        for (size_t wi = 0; wi < works.size(); ++wi) {
            auto& w = works[wi];
            if (!w.done) throw std::runtime_error("F22 phase: the device skipped box " + std::to_string(w.morton));
            for (int64_t i : w.skel22)
                if (i < 0 || i >= static_cast<int64_t>(w.pos22.size())) throw std::logic_error("F22 phase: device skeleton index " + std::to_string(i) + " of " + std::to_string(w.pos22.size()));
            for (int64_t i : w.red22)
                if (i < 0 || i >= static_cast<int64_t>(w.pos22.size())) throw std::logic_error("F22 phase: device redundant index " + std::to_string(i) + " of " + std::to_string(w.pos22.size()));
            Box& b = level.local_boxes[static_cast<size_t>(occupied_local[work_box[wi]])];
            const int64_t k22 = static_cast<int64_t>(w.skel22.size()), r22 = static_cast<int64_t>(w.red22.size());
            kept_points += k22;
            if (r22 == 0) continue;
            std::vector<int64_t> posS(w.pos11);
            for (int64_t i : w.skel22) posS.push_back(w.pos22[static_cast<size_t>(i)]);
            auto& f = b.f22;
            f.sparsified = true;
            const size_t k11 = w.pos11.size();
            for (size_t i = k11; i < posS.size(); ++i) f.skeleton.push_back(w.S_old[static_cast<size_t>(posS[i])]);
            for (int64_t j : w.red22) f.redundant.push_back(w.S_old[static_cast<size_t>(w.pos22[static_cast<size_t>(j)])]);
            kept[work_box[wi]] = std::move(posS);
            ++sparsified_boxes;
            sparsified_points += r22;
        }
        publish_kept();
        index_live();
        for (auto& w : works)
            for (auto& nr : w.near) nr.q22_nr = live_rows(nr.morton, nr.q22);
        const double t_records0 = MPI_Wtime();
        t_records += t_records0 - t_kept0;  // (the kept sets)
        records(works);
        const double t_results0 = MPI_Wtime();
        t_sketch += t_results0 - t_records0;
        // ---- the records into the boxes
        for (size_t wi = 0; wi < works.size(); ++wi) {
            auto& w = works[wi];
            Box& b = level.local_boxes[static_cast<size_t>(occupied_local[work_box[wi]])];
            auto& f = b.f22;
            if (!f.sparsified) continue;
            f.T = std::move(w.T);
            f.X_SR = std::move(w.X_SR);
            f.X_RS = std::move(w.X_RS);
            f.X_RR = std::move(w.X_RR);
            f.X_NR = std::move(w.X_NR);
            f.X_RN = std::move(w.X_RN);
            f.neighbor_f22_counts.assign(b.one_hop.size(), 0);
            f.neighbor_f22_slots.assign(b.one_hop.size(), {});
            f.neighbor_nr_counts.assign(b.one_hop.size(), 0);
            f.neighbor_nr_slots.assign(b.one_hop.size(), {});
            for (size_t idx = 0; idx < b.one_hop.size(); ++idx) {
                const auto& nr = w.near[idx];
                for (int64_t q : nr.q22) f.neighbor_f22_slots[idx].push_back((*nr.S_old)[static_cast<size_t>(q)]);
                for (int64_t q : nr.q22_nr) f.neighbor_nr_slots[idx].push_back((*nr.S_old)[static_cast<size_t>(q)]);
                f.neighbor_f22_counts[idx] = static_cast<int64_t>(nr.q22.size());
                f.neighbor_nr_counts[idx] = static_cast<int64_t>(nr.q22_nr.size());
            }
        }
        for (size_t oi = 0; oi < occupied_local.size(); ++oi) {  // the new live sets
            Box& b = level.local_boxes[static_cast<size_t>(occupied_local[oi])];
            if (!b.f22.sparsified) continue;
            b.skeleton_indices.clear();
            for (int64_t q : kept[oi]) b.skeleton_indices.push_back(b.skeleton_before_f22[static_cast<size_t>(q)]);
        }
        t_records += MPI_Wtime() - t_results0;  // (the results into the boxes)
        const double t_back0 = MPI_Wtime();
        if (blocks_back) blocks_back(kept_of);
        t_id += MPI_Wtime() - t_back0;  // (the blocks sliced)
        sliced_on_device = true;
    } else {
        // a box between the two passes: its ID's split of S_old and its
        // neighbours.  K(nb, b) and K(b, nb) over the two boxes' live sets
        // before the phase: this box's views, else the neighbour's views of
        // this box (a local neighbour), else the kernel - a pair no
        // elimination touched keeps no view, and the transition reads it the
        // same way
        struct NearF11 {
            const std::vector<int64_t>* S_old = nullptr;
            const std::vector<int64_t>* ids = nullptr;
            std::vector<int64_t> q11, q22, q22_nr;
            const fmm::ModifiedBlock<DataType>* nb_b = nullptr;  // K(nb, b)
            bool nb_b_own = true;                                // b's A_NS[nb], else nb's A_SN[b]
            const fmm::ModifiedBlock<DataType>* b_nb = nullptr;  // K(b, nb)
            bool b_nb_own = true;                                // b's A_SN[nb], else nb's A_NS[b]
        };
        struct Step {
            std::vector<int64_t> pos11, posS, posR;
            std::vector<NearF11> near;
        };
        std::vector<Step> steps(occupied_local.size());
        // the point ids of a neighbour's positions and of this box's
        auto nb_ids_of = [](const NearF11& nr, const std::vector<int64_t>& qs) {
            std::vector<int64_t> out(qs.size());
            for (size_t i = 0; i < qs.size(); ++i) out[i] = (*nr.ids)[static_cast<size_t>((*nr.S_old)[static_cast<size_t>(qs[i])])];
            return out;
        };
        auto b_ids_of = [](const Box& b, const std::vector<int64_t>& ps) {
            std::vector<int64_t> out(ps.size());
            for (size_t i = 0; i < ps.size(); ++i)
                out[i] = b.point_indices[static_cast<size_t>(b.skeleton_before_f22[static_cast<size_t>(ps[i])])];
            return out;
        };
        // K(nb, b)(qs, ps) and K(b, nb)(ps, qs), column-major
        auto coupling_nb_b = [&](const Box& b, const NearF11& nr, const std::vector<int64_t>& qs, const std::vector<int64_t>& ps) {
            const int64_t m = static_cast<int64_t>(qs.size()), n = static_cast<int64_t>(ps.size());
            std::vector<DataType> out(static_cast<size_t>(m * n));
            const auto* v = nr.nb_b;
            if (v == nullptr) {
                if (m > 0 && n > 0) {
                    const auto ri = nb_ids_of(nr, qs), ci = b_ids_of(b, ps);
                    kernel->evaluate_block_by_index(ri.data(), m, ci.data(), n, out.data(), m);
                }
                return out;
            }
            for (int64_t c = 0; c < n; ++c)
                for (int64_t i = 0; i < m; ++i)
                    out[static_cast<size_t>(i + c * m)] = nr.nb_b_own ? v->a_ns(qs[static_cast<size_t>(i)], ps[static_cast<size_t>(c)])
                                                                      : v->A_SN(qs[static_cast<size_t>(i)], ps[static_cast<size_t>(c)]);
            return out;
        };
        auto coupling_b_nb = [&](const Box& b, const NearF11& nr, const std::vector<int64_t>& ps, const std::vector<int64_t>& qs) {
            const int64_t m = static_cast<int64_t>(ps.size()), n = static_cast<int64_t>(qs.size());
            std::vector<DataType> out(static_cast<size_t>(m * n));
            const auto* v = nr.b_nb;
            if (v == nullptr) {
                if (m > 0 && n > 0) {
                    const auto ri = b_ids_of(b, ps), ci = nb_ids_of(nr, qs);
                    kernel->evaluate_block_by_index(ri.data(), m, ci.data(), n, out.data(), m);
                }
                return out;
            }
            for (int64_t c = 0; c < n; ++c)
                for (int64_t i = 0; i < m; ++i)
                    out[static_cast<size_t>(i + c * m)] = nr.b_nb_own ? v->A_SN(ps[static_cast<size_t>(i)], qs[static_cast<size_t>(c)])
                                                                      : v->a_ns(ps[static_cast<size_t>(i)], qs[static_cast<size_t>(c)]);
            return out;
        };
        // ---- pass 1: the IDs and the own records
#pragma omp parallel
    {
        fmm::FactorizationThreadScratch<CoordType, DataType> scratch;
#pragma omp for schedule(dynamic, 1) reduction(+ : sparsified_boxes, sparsified_points, kept_points, t_sketch, t_target, t_id, t_records)
        for (int64_t oi = 0; oi < static_cast<int64_t>(occupied_local.size()); ++oi) {
            try {
                Box& b = level.local_boxes[static_cast<size_t>(occupied_local[static_cast<size_t>(oi)])];
                const std::vector<int64_t>& S_old = b.skeleton_before_f22;
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
                    std::vector<int64_t> pts(static_cast<size_t>(k_old));
                    for (int64_t q = 0; q < k_old; ++q)
                        pts[static_cast<size_t>(q)] = b.point_indices[static_cast<size_t>(S_old[static_cast<size_t>(q)])];
                    MS block;
                    block.allocate(k_old, k_old, MS::FULL);
                    kernel->evaluate_block_by_index(pts.data(), k_old, pts.data(), k_old, block.data.data(), k_old);
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
                std::vector<NearF11> near(b.one_hop.size());
                int64_t near_rows = 0;
                for (size_t idx = 0; idx < b.one_hop.size(); ++idx) {
                    const int64_t nm = b.one_hop[idx];
                    const Old& o = old_of.at(nm);
                    near[idx].S_old = o.S_old;
                    near[idx].ids = o.ids;
                    {
                        const int64_t n_nb = static_cast<int64_t>(o.S_old->size());
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
                            if (auto* v = a_sn_block(*nb_local, b.morton_index); v != nullptr && v->A_SN.is_allocated()) {
                                if (v->A_SN.rows != n_nb || v->A_SN.cols != k_old) {
                                    throw std::runtime_error("F22 phase: box " + std::to_string(nm) + "'s view of " +
                                                             std::to_string(b.morton_index) + " is not over the two skeletons");
                                }
                                near[idx].nb_b = v;
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
                            if (auto* v = a_ns_block(*nb_local, b.morton_index); v != nullptr && v->a_ns_is_allocated()) {
                                if (v->a_ns_rows() != k_old || v->a_ns_cols() != n_nb) {
                                    throw std::runtime_error("F22 phase: box " + std::to_string(nm) + "'s view of " +
                                                             std::to_string(b.morton_index) + " is not over the two skeletons");
                                }
                                near[idx].b_nb = v;
                                near[idx].b_nb_own = false;
                            }
                        }
                    }
                    const auto& nS = *o.S_old;
                    for (int64_t q = 0; q < static_cast<int64_t>(nS.size()); ++q) {
                        if (tree->partial_f11((*o.ids)[static_cast<size_t>(nS[static_cast<size_t>(q)])])) near[idx].q11.push_back(q);
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
                    const int64_t m = static_cast<int64_t>(near[idx].q11.size());
                    if (m == 0) continue;
                    const auto blk = coupling_nb_b(b, near[idx], near[idx].q11, pos22);
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
                    const auto blk = coupling_b_nb(b, near[idx], pos22, near[idx].q11);  // n22 x m
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
                const int64_t k22 = static_cast<int64_t>(id.skeleton_indices.size());
                const int64_t r22 = static_cast<int64_t>(id.redundant_indices.size());
                if (fmm::h2_id_trace_enabled()) {  // (BPACK_TRACE=id; the device's F22 step writes the same)
                    fmm::h2_id_trace_write("F22 L" + std::to_string(level.level) + " m=" + std::to_string(b.morton_index) +
                                           " n22=" + std::to_string(n22) + " k22=" + std::to_string(k22) +
                                           " rows=" + std::to_string(rows) + fmm::h2_id_adaptive_trace(scratch.id_adaptive));
                }
                kept_points += k22;
                // nothing to sparsify: the box stays as it is (its live set
                // unchanged, no records)
                if (r22 == 0) continue;
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
                Step& st = steps[static_cast<size_t>(oi)];
                st.posS = pos11;
                for (int64_t i = 0; i < k22; ++i) st.posS.push_back(pos22[static_cast<size_t>(id.skeleton_indices[static_cast<size_t>(ps[static_cast<size_t>(i)])])]);
                for (int64_t j = 0; j < r22; ++j) st.posR.push_back(pos22[static_cast<size_t>(id.redundant_indices[static_cast<size_t>(pr_[static_cast<size_t>(j)])])]);
                const int64_t kp = static_cast<int64_t>(st.posS.size());
                auto& f = b.f22;
                f.sparsified = true;
                for (int64_t i = 0; i < k22; ++i) f.skeleton.push_back(S_old[static_cast<size_t>(st.posS[static_cast<size_t>(k11 + i)])]);
                for (int64_t j = 0; j < r22; ++j) f.redundant.push_back(S_old[static_cast<size_t>(st.posR[static_cast<size_t>(j)])]);
                // T over S' rows (zero on S11)
                f.T.allocate(kp, r22);
                for (int64_t j = 0; j < r22; ++j)
                    for (int64_t i = 0; i < k22; ++i)
                        f.T.data[static_cast<size_t>(k11 + i + j * kp)] = id.interpolation(ps[static_cast<size_t>(i)], pr_[static_cast<size_t>(j)]);
                t0 = MPI_Wtime();
                const auto& posS = st.posS;
                const auto& posR = st.posR;
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
                f.X_SR.assign(kp, r22, A_SR);
                partial_gemm<DataType>('N', 'N', kp, r22, kp, DataType(-1), A_SS.data(), kp, T, kp, DataType(1), f.X_SR.data.data(), kp);
                f.X_RS.assign(r22, kp, A_RS);
                partial_gemm<DataType>('T', 'N', r22, kp, kp, DataType(-1), T, kp, A_SS.data(), kp, DataType(1), f.X_RS.data.data(), r22);
                f.X_RR.assign(r22, r22, A_RR);
                // X_RR = A_RR - T^T A_SR - X_RS T  (X_RS T = A_RS T - T^T A_SS T)
                partial_gemm<DataType>('T', 'N', r22, r22, kp, DataType(-1), T, kp, A_SR.data(), kp, DataType(1), f.X_RR.data.data(), r22);
                partial_gemm<DataType>('N', 'N', r22, r22, kp, DataType(-1), f.X_RS.data.data(), r22, T, kp, DataType(1), f.X_RR.data.data(), r22);
                // ---- the neighbours' raw F22 points: X_RN = A_RN - T^T A_SN
                // (X_NR, over the ones that stay live, in the second pass)
                const size_t nn = b.one_hop.size();
                f.neighbor_f22_counts.assign(nn, 0);
                f.neighbor_f22_slots.assign(nn, {});
                int64_t total = 0;
                for (size_t idx = 0; idx < nn; ++idx) {
                    for (int64_t q : near[idx].q22) f.neighbor_f22_slots[idx].push_back((*near[idx].S_old)[static_cast<size_t>(q)]);
                    f.neighbor_f22_counts[idx] = static_cast<int64_t>(near[idx].q22.size());
                    total += f.neighbor_f22_counts[idx];
                }
                f.X_RN.allocate(r22, total);
                int64_t col0 = 0;
                for (size_t idx = 0; idx < nn; ++idx) {
                    const int64_t n = f.neighbor_f22_counts[idx];
                    if (n == 0) continue;
                    // K(b, nb raw F22) over (posR, posS)
                    std::vector<DataType> A_RN = coupling_b_nb(b, near[idx], posR, near[idx].q22);  // r22 x n
                    const auto A_SN = coupling_b_nb(b, near[idx], posS, near[idx].q22);           // kp x n
                    partial_gemm<DataType>('T', 'N', r22, n, kp, DataType(-1), T, kp, A_SN.data(), kp, DataType(1), A_RN.data(), r22);
                    std::copy_n(A_RN.data(), static_cast<size_t>(r22 * n), f.X_RN.data.data() + col0 * r22);
                    col0 += n;
                }
                if (tree->partial_compress_records) {
                    for (auto* m : {&f.X_SR, &f.X_RS, &f.X_RR, &f.X_RN}) fmm::compress_record_host(*m, tolerance);
                }
                t_records += MPI_Wtime() - t0;
                st.pos11 = std::move(pos11);
                st.near = std::move(near);
                kept[static_cast<size_t>(oi)] = st.posS;
                ++sparsified_boxes;
                sparsified_points += r22;
            } catch (...) {
#pragma omp critical(partial_f22_phase_failure)
                if (!failure) failure = std::current_exception();
            }
        }
    }
        if (failure) std::rethrow_exception(failure);
        publish_kept();
        index_live();
        for (size_t oi = 0; oi < occupied_local.size(); ++oi) {
            if (!level.local_boxes[static_cast<size_t>(occupied_local[oi])].f22.sparsified) continue;
            const auto& one_hop = level.local_boxes[static_cast<size_t>(occupied_local[oi])].one_hop;
            for (size_t idx = 0; idx < one_hop.size(); ++idx) {
                auto& nr = steps[oi].near[idx];
                nr.q22_nr = live_rows(one_hop[idx], nr.q22);
            }
        }
        // ---- pass 2: X_NR = A_NR - A_NS T over the neighbours' raw F22 points
        // that stay live
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : t_records)
        for (int64_t oi = 0; oi < static_cast<int64_t>(occupied_local.size()); ++oi) {
            try {
                Box& b = level.local_boxes[static_cast<size_t>(occupied_local[static_cast<size_t>(oi)])];
                auto& f = b.f22;
                if (!f.sparsified) continue;
                const double t0 = MPI_Wtime();
                const Step& st = steps[static_cast<size_t>(oi)];
                const int64_t kp = static_cast<int64_t>(st.posS.size()), r22 = static_cast<int64_t>(st.posR.size());
                const DataType* T = f.T.data.data();
                const size_t nn = b.one_hop.size();
                f.neighbor_nr_counts.assign(nn, 0);
                f.neighbor_nr_slots.assign(nn, {});
                for (size_t idx = 0; idx < nn; ++idx) {
                    const auto& nr = st.near[idx];
                    for (int64_t q : nr.q22_nr) f.neighbor_nr_slots[idx].push_back((*nr.S_old)[static_cast<size_t>(q)]);
                    f.neighbor_nr_counts[idx] = static_cast<int64_t>(nr.q22_nr.size());
                }
                f.X_NR.assign(nn, {});
                for (size_t idx = 0; idx < nn; ++idx) {
                    const auto& nr = st.near[idx];
                    const int64_t m = static_cast<int64_t>(nr.q22_nr.size());
                    if (m == 0) continue;
                    // K(nb live F22, b) over (posR, posS)
                    std::vector<DataType> A_NR = coupling_nb_b(b, nr, nr.q22_nr, st.posR);  // m x r22
                    const auto A_NS = coupling_nb_b(b, nr, nr.q22_nr, st.posS);           // m x kp
                    partial_gemm<DataType>('N', 'N', m, r22, kp, DataType(-1), A_NS.data(), m, T, kp, DataType(1), A_NR.data(), m);
                    f.X_NR[idx].assign(m, r22, A_NR);
                    if (tree->partial_compress_records) fmm::compress_record_host(f.X_NR[idx], tolerance);
                }
                t_records += MPI_Wtime() - t0;
            } catch (...) {
#pragma omp critical(partial_f22_phase_failure)
                if (!failure) failure = std::current_exception();
            }
        }
        if (failure) std::rethrow_exception(failure);
    }
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
            std::vector<int64_t> pts(static_cast<size_t>(kp));
            for (int64_t i = 0; i < kp; ++i)
                pts[static_cast<size_t>(i)] =
                    b.point_indices[static_cast<size_t>(b.skeleton_before_f22[static_cast<size_t>(posS[static_cast<size_t>(i)])])];
            MS block;
            block.allocate(kp, kp, MS::FULL);
            kernel->evaluate_block_by_index(pts.data(), kp, pts.data(), kp, block.data.data(), kp);
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
    // ---- every box's views of its sparsified neighbours (this rank's and,
    // published, the others'): their rows and columns to the new live sets
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
    // the records' bytes as stored (on the host, kept on the device), and uncompressed
    double record_gb[3] = {0.0, 0.0, 0.0};
    for (int64_t li : occupied_local) {
        const auto& f = level.local_boxes[static_cast<size_t>(li)].f22;
        if (!f.sparsified) continue;
        auto add = [&](const fmm::PartialRecord<DataType>& m) {
            double kept = 0.0;
            if (m.device_only()) {
                if (m.compressed()) {
                    for (uint32_t st : m.strides) kept += static_cast<double>(st);
                } else {
                    kept = static_cast<double>(m.rows) * static_cast<double>(m.cols) * sizeof(DataType);
                }
            }
            record_gb[0] += (static_cast<double>(m.bytes()) + kept) / 1e9;
            record_gb[1] += static_cast<double>(m.rows) * static_cast<double>(m.cols) * sizeof(DataType) / 1e9;
            record_gb[2] += kept / 1e9;
        };
        for (const auto* m : {&f.X_SR, &f.X_RS, &f.X_RR, &f.X_RN}) add(*m);
        for (const auto& m : f.X_NR) add(m);
    }
    int64_t counts[3] = {sparsified_boxes, sparsified_points, kept_points};
    if (tree->mpi_size > 1) {
        MPI_Allreduce(MPI_IN_PLACE, counts, 3, MPI_INT64_T, MPI_SUM, level_comm);
        MPI_Allreduce(MPI_IN_PLACE, record_gb, 3, MPI_DOUBLE, MPI_SUM, level_comm);
    }
    if (print_detail && fmm::env::trace("phase")) {
        if (device) {
            std::printf("  F22 phase parts (this rank): plans %.2f s, device %.2f s, kept sets and results %.2f s, blocks "
                        "sliced %.2f s, then %.2f s serial\n", t_target, t_sketch, t_records, t_id,
                        MPI_Wtime() - t_loop0 - t_loop);
        } else {
            std::printf("  F22 phase parts (this rank): passes %.2f s, then %.2f s serial; thread-seconds: sketch %.2f, "
                        "target %.2f, ID %.2f, records %.2f\n", t_loop, MPI_Wtime() - t_loop0 - t_loop, t_sketch,
                        t_target, t_id, t_records);
        }
    }
    if (print_detail) {
        std::printf("  F22 phase: %lld boxes, %lld points sparsified, %lld kept; records %.2f GB", static_cast<long long>(counts[0]),
                    static_cast<long long>(counts[1]), static_cast<long long>(counts[2]), record_gb[0]);
        if (tree->partial_compress_records) std::printf(" (%.2f GB uncompressed)", record_gb[1]);
        if (record_gb[2] > 0.0) std::printf(", %.2f GB of them kept on the device only", record_gb[2]);
        std::printf("\n");
    }
#ifdef H2_HAVE_GPU
    // the kept records' host copies before an eviction of the operator's device data
    if (tree->partial_records_on_device) {
        fmm::gpu::partial_device_store().spill = [tree] { partial_records_to_host(*tree, true); };
    }
#endif
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

// The root's split by point id (its F11 and F22 slots), before its numbers
// (partial_root_factorize on the host, LevelEliminator's
// partial_root_on_device on the device).
template<typename CoordType, typename DataType>
void partial_root_setup(fmm::ParallelTree<CoordType, DataType>& tree, const fmm::BoxData<CoordType, DataType>& root_box) {
    auto& pr = tree.partial_root;
    const int64_t n = root_box.num_points;
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
    partial_root_setup(tree, root_box);
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
//   one hop:   vp_i^T X_RN_bi(:, bj raw) v_j + wv_i^T X_NR_bj(bi live raw, :) vp_j
// and at the root v_i^T CB_root v_j over the F22 slots.  The chains are
// zero on F11 slots (T has no S11 rows), so F11 rows never enter a block.
// Row blocks partition the requested rows, so a level's terms run in
// parallel over them.
// =====================================================================
template<typename DataType>
struct PartialChainBlock {
    int64_t box = -1;                ///< local box index at the level
    int64_t inter = 0;               ///< the block of the request its members are in
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

// The F22 points' leaf boxes and slots (local leaf boxes; -1: another
// rank's), built once.
template<typename CoordType, typename DataType>
void partial_leaf_index(fmm::ParallelTree<CoordType, DataType>& tree) {
    auto& pr = tree.partial_root;
    const int64_t pc = tree.partial_count, n22 = tree.num_points - pc;
    const int leaf = tree.num_levels - 1;
    if (!pr.leaf_box.empty() || n22 <= 0) return;
    pr.leaf_box.assign(static_cast<size_t>(n22), -1);
    pr.leaf_slot.assign(static_cast<size_t>(n22), -1);
    const auto& level = tree.levels[static_cast<size_t>(leaf)];
    if (!level.is_process_active) return;
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

// One side of a CB request: each block of the request's distinct points (the
// same on every rank), their positions in it (pos[pos_start[m] ..
// pos_start[m + 1])), and the leaf chain blocks of this rank's boxes (with
// their identity values when `values`)
template<typename DataType>
struct PartialCbSide {
    std::vector<int64_t> points, pos_start, pos;
    std::vector<PartialChainBlock<DataType>> blocks;
};
template<typename CoordType, typename DataType>
PartialCbSide<DataType> partial_cb_side(const fmm::ParallelTree<CoordType, DataType>& tree, int64_t ninter, const int64_t* idx,
                                        const std::vector<int64_t>& off, bool values) {
    const auto& pr = tree.partial_root;
    const int64_t n22 = tree.num_points - tree.partial_count;
    int nranks = 1;
    MPI_Comm_size(tree.comm, &nranks);
    // each block of the request on its own, then numbered in the request's order
    using Side = PartialCbSide<DataType>;
    std::vector<Side> parts(static_cast<size_t>(ninter));
    std::exception_ptr failure;
#pragma omp parallel for schedule(dynamic, 1)
    for (int64_t k = 0; k < ninter; ++k) {
        try {
            Side& part = parts[static_cast<size_t>(k)];
            std::unordered_map<int64_t, int64_t> of;       // F22 index -> distinct point
            std::unordered_map<int64_t, size_t> block_of;  // leaf box -> chain block
            const int64_t q0 = off[static_cast<size_t>(k)], n = off[static_cast<size_t>(k) + 1] - q0;
            std::vector<int64_t> point_of(static_cast<size_t>(n));
            for (int64_t q = 0; q < n; ++q) {
                const int64_t f = idx[q0 + q];
                if (f < 0 || f >= n22) throw std::runtime_error("partial Schur: an index out of range");
                auto it = of.find(f);
                if (it == of.end()) {
                    const int64_t m = static_cast<int64_t>(part.points.size());
                    it = of.emplace(f, m).first;
                    part.points.push_back(f);
                    const int64_t b = pr.leaf_box[static_cast<size_t>(f)];
                    if (b >= 0) {
                        auto bt = block_of.find(b);
                        if (bt == block_of.end()) {
                            bt = block_of.emplace(b, part.blocks.size()).first;
                            part.blocks.emplace_back();
                            part.blocks.back().box = b;
                            part.blocks.back().inter = k;
                        }
                        auto& blk = part.blocks[bt->second];
                        blk.members.push_back(m);
                        blk.rows.push_back(pr.leaf_slot[static_cast<size_t>(f)]);
                    } else if (nranks == 1) {
                        throw std::runtime_error("partial Schur: a point is not in a local leaf");
                    }  // (else another rank's)
                }
                point_of[static_cast<size_t>(q)] = it->second;
            }
            const size_t np = part.points.size();
            part.pos_start.assign(np + 1, 0);
            for (int64_t m : point_of) ++part.pos_start[static_cast<size_t>(m) + 1];
            for (size_t m = 0; m < np; ++m) part.pos_start[m + 1] += part.pos_start[m];
            part.pos.resize(static_cast<size_t>(n));
            std::vector<int64_t> fill(part.pos_start.begin(), part.pos_start.end() - 1);
            for (int64_t q = 0; q < n; ++q) part.pos[static_cast<size_t>(fill[static_cast<size_t>(point_of[static_cast<size_t>(q)])]++)] = q;
            if (values) {
                for (auto& blk : part.blocks) {
                    const int64_t nb = static_cast<int64_t>(blk.members.size());
                    blk.v.assign(static_cast<size_t>(nb * nb), DataType(0));
                    for (int64_t t = 0; t < nb; ++t) blk.v[static_cast<size_t>(t + t * nb)] = DataType(1);
                }
            }
        } catch (...) {
#pragma omp critical(partial_schur_blocks_failure)
            if (!failure) failure = std::current_exception();
        }
    }
    if (failure) std::rethrow_exception(failure);
    Side side;
    side.pos_start.push_back(0);
    for (auto& part : parts) {
        const int64_t m0 = static_cast<int64_t>(side.points.size()), p0 = static_cast<int64_t>(side.pos.size());
        side.points.insert(side.points.end(), part.points.begin(), part.points.end());
        for (size_t m = 1; m < part.pos_start.size(); ++m) side.pos_start.push_back(p0 + part.pos_start[m]);
        side.pos.insert(side.pos.end(), part.pos.begin(), part.pos.end());
        for (auto& blk : part.blocks) {
            for (auto& m : blk.members) m += m0;
            side.blocks.push_back(std::move(blk));
        }
    }
    return side;
}

// H2_partial_extract auto (0): the device for a call of at least this many
// entries (when the records are on it)
inline constexpr int64_t kPartialExtractDeviceEntries = int64_t{1} << 16;

#ifdef H2_HAVE_GPU
// partial_schur_blocks on the device (H2_partial_extract, h2_gpu/
// partial_device.hpp): the same levels, blocks, pairs and terms, with the
// chain blocks' values on the device for the whole call (the leaves'
// identities, each level's T22 steps, the pair terms, the merges, the root's
// pairs); the host plans the index lists, sends those of the blocks other
// ranks read (their values device to device when MPI is GPU-aware:
// partial_values_exchange) and sums the pairs' blocks into out.  The records
// and the root block are the store's (resident).
template<typename DataType>
struct PartialDeviceChain {
    int64_t box = -1, inter = 0;
    std::vector<int64_t> members, rows, ridx, wrows;
    std::vector<int32_t> rowpos, wpos;
    DataType* v = nullptr;   // |rows| x |members|
    DataType* vp = nullptr;  // |ridx| x |members|
    DataType* wv = nullptr;  // |wrows| x |members|
};

template<typename CoordType, typename DataType>
void partial_schur_blocks_device(fmm::ParallelTree<CoordType, DataType>& tree, int64_t ninter, const int64_t* rows,
                                 const int64_t* nrows, const int64_t* cols, const int64_t* ncols, DataType* out,
                                 const std::vector<int64_t>& roff, const std::vector<int64_t>& coff,
                                 const std::vector<int64_t>& ooff) {
    namespace gpu = fmm::gpu;
    using S = typename fmm::gpu::DeviceScalar<DataType>::type;
    using Chain = PartialDeviceChain<DataType>;
    using Pair = gpu::PartialPair<DataType>;
    using Term = gpu::PartialTerm<DataType>;
    const size_t D = sizeof(DataType);
    auto& pr = tree.partial_root;
    const int dimension = tree.dimension;
    const int leaf = tree.num_levels - 1;
    int me = 0, nranks = 1;
    MPI_Comm_rank(tree.comm, &me);
    MPI_Comm_size(tree.comm, &nranks);
    gpu::PartialDeviceWork work;
    gpu::PartialDeviceStore& store = work.store();
    cudaStream_t stream = work.stream();
    // BPACK_TRACE=phase: the call's time by phase on each rank (the stream
    // synchronized at each phase's end, so a phase holds its device work)
    // and the bytes the exchanges send
    const bool trace = fmm::env::trace("phase");
    enum { kPlan, kSteps, kTermsOut, kTermsMpi, kTermsIn, kPairs, kUpOut, kUpMpi, kUpIn, kRoot, kPhases };
    double tp[kPhases] = {}, sent[2] = {0.0, 0.0};
    double tmark = MPI_Wtime();
    const double tcall = tmark;
    auto mark = [&](int phase) {
        if (!trace) return;
        work.sync();
        const double t = MPI_Wtime();
        tp[phase] += t - tmark;
        tmark = t;
    };
    partial_leaf_index(tree);
    const PartialCbSide<DataType> Rh = partial_cb_side(tree, ninter, rows, roff, false);
    const PartialCbSide<DataType> Ch = partial_cb_side(tree, ninter, cols, coff, false);
    // G (|members_i| x |members_j|) into block k of out at the requested positions (as the host's)
    auto scatter = [&](int64_t k, const std::vector<int64_t>& mi, const std::vector<int64_t>& mj, const DataType* G) {
        const int64_t ni = static_cast<int64_t>(mi.size()), ld = nrows[k];
        DataType* ok = out + ooff[static_cast<size_t>(k)];
        for (size_t c = 0; c < mj.size(); ++c) {
            const size_t m = static_cast<size_t>(mj[c]);
            for (int64_t pj = Ch.pos_start[m]; pj < Ch.pos_start[m + 1]; ++pj) {
                DataType* col = ok + Ch.pos[static_cast<size_t>(pj)] * ld;
                for (int64_t r = 0; r < ni; ++r) {
                    const DataType g = G[static_cast<size_t>(r + static_cast<int64_t>(c) * ni)];
                    const size_t mr = static_cast<size_t>(mi[static_cast<size_t>(r)]);
                    for (int64_t pi = Rh.pos_start[mr]; pi < Rh.pos_start[mr + 1]; ++pi) col[Rh.pos[static_cast<size_t>(pi)]] += g;
                }
            }
        }
    };
    // ---- the leaves' chain blocks: identities on the device
    std::vector<Chain> sides[2];
    std::vector<char*> arenas;  // the current level's blocks' values
    {
        fmm::gpu::MetaBuilder& meta = work.meta();
        std::vector<fmm::gpu::IdentityItemT<S>> ids;
        int max_n = 0;
        for (int sd = 0; sd < 2; ++sd) {
            const auto& hb = (sd == 0 ? Rh : Ch).blocks;
            size_t bytes = 0;
            for (const auto& b : hb) bytes += fmm::gpu::align_up(b.members.size() * b.members.size() * D);
            char* arena = work.alloc(bytes);
            arenas.push_back(arena);
            size_t at = 0;
            for (const auto& b : hb) {
                Chain c;
                c.box = b.box;
                c.inter = b.inter;
                c.members = b.members;
                c.rows = b.rows;
                c.v = reinterpret_cast<DataType*>(arena + at);
                const int n = static_cast<int>(b.members.size());
                at += fmm::gpu::align_up(static_cast<size_t>(n) * static_cast<size_t>(n) * D);
                ids.push_back({reinterpret_cast<S*>(c.v), n, n});
                max_n = std::max(max_n, n);
                sides[sd].push_back(std::move(c));
            }
        }
        const size_t off = meta.append(ids);
        char* md = work.upload_meta();
        fmm::gpu::launch_identity(reinterpret_cast<const fmm::gpu::IdentityItemT<S>*>(md + off), static_cast<int>(ids.size()), max_n,
                                  stream);
    }
    std::vector<Chain>& R = sides[0];
    std::vector<Chain>& C = sides[1];
    // a group of pairs on the device, then their blocks into out
    struct PairInfo { int64_t inter; const std::vector<int64_t>* mi; const std::vector<int64_t>* mj; };
    auto run_pairs = [&](std::vector<Pair>& pairs, const std::vector<PairInfo>& info) {
        const size_t budget = work.budget();
        const size_t pool = store.scratch_bytes / (std::max<size_t>(1, fmm::kRecordChunkBytes / D) * D);
        size_t i0 = 0;
        while (i0 < pairs.size()) {
            size_t bytes = 0, chunks = 0, i1 = i0;
            std::unordered_map<const fmm::gpu::PartialDeviceRecord*, std::unordered_set<int64_t>> seen;
            std::vector<int64_t> ck;
            for (; i1 < pairs.size(); ++i1) {
                const size_t pb = fmm::gpu::partial_pair_bytes(pairs[i1]);
                size_t pc_new = 0;
                for (const auto& t : pairs[i1].terms) {
                    fmm::gpu::partial_term_chunks(t, ck);
                    for (int64_t k : ck) pc_new += seen.count(t.rec) && seen[t.rec].count(k) ? 0 : 1;
                }
                if (i1 > i0 && (bytes + pb > budget || chunks + pc_new > pool)) break;
                bytes += pb;
                chunks += pc_new;
                for (const auto& t : pairs[i1].terms) {
                    fmm::gpu::partial_term_chunks(t, ck);
                    for (int64_t k : ck) seen[t.rec].insert(k);
                }
            }
            std::vector<const Pair*> group;
            size_t g_elems = 0;
            for (size_t q = i0; q < i1; ++q) {
                group.push_back(&pairs[q]);
                g_elems += static_cast<size_t>(pairs[q].ni) * static_cast<size_t>(pairs[q].nj);
            }
            std::vector<DataType> g(g_elems);
            fmm::gpu::partial_device_pairs(work, group, g.data());
            size_t at = 0;
            for (size_t q = i0; q < i1; ++q) {
                scatter(info[q].inter, *info[q].mi, *info[q].mj, g.data() + at);
                at += static_cast<size_t>(pairs[q].ni) * static_cast<size_t>(pairs[q].nj);
            }
            i0 = i1;
        }
    };
    // a block of another rank's box: its members, raw rows, w part and live set, values on the device
    struct RemoteD {
        int64_t morton = -1, num_points = 0, inter = 0;
        std::vector<int64_t> members, rows, wrows, skeleton;
        std::vector<int32_t> rowpos, wpos;
        DataType* v = nullptr;
        DataType* wv = nullptr;
        size_t voff = 0;  // (its values in the received ones)
    };
    for (int l = leaf; l >= 0; --l) {
        auto& level = tree.levels[static_cast<size_t>(l)];
        const bool active = level.is_process_active;
        if (l == 0) {
            // the root: v_i^T CB_root v_j over the root's F22 slots
            if (pr.active && store.root_cb != nullptr) {
                const auto& root = level.local_boxes[0];
                std::vector<int64_t> pos22(static_cast<size_t>(root.num_points), -1);
                for (int64_t a = 0; a < pr.n22; ++a) pos22[static_cast<size_t>(pr.slots22[static_cast<size_t>(a)])] = a;
                struct Part22 { std::vector<int64_t> a, t; };
                std::vector<Part22> rp(R.size()), cp(C.size());
                auto f22_part = [&](const Chain& b, Part22& p) {
                    for (size_t q = 0; q < b.rows.size(); ++q) {
                        const int64_t x = pos22[static_cast<size_t>(b.rows[q])];
                        if (x >= 0) { p.a.push_back(x); p.t.push_back(static_cast<int64_t>(q)); }
                    }
                };
                for (size_t i = 0; i < R.size(); ++i) f22_part(R[i], rp[i]);
                for (size_t j = 0; j < C.size(); ++j) f22_part(C[j], cp[j]);
                std::unordered_map<int64_t, std::vector<size_t>> cols_of;
                for (size_t j = 0; j < C.size(); ++j)
                    if (!cp[j].a.empty()) cols_of[C[j].inter].push_back(j);
                std::vector<Pair> pairs;
                std::vector<PairInfo> info;
                for (size_t i = 0; i < R.size(); ++i) {
                    if (rp[i].a.empty()) continue;
                    const auto it = cols_of.find(R[i].inter);
                    if (it == cols_of.end()) continue;
                    for (size_t j : it->second) {
                        Pair P;
                        P.ni = static_cast<int>(R[i].members.size());
                        P.nj = static_cast<int>(C[j].members.size());
                        Term T;
                        T.kind = 5;
                        T.raw = reinterpret_cast<const DataType*>(store.root_cb);
                        T.ld = store.root_n;
                        T.rows_ext = &rp[i].a;
                        T.cols_ext = &cp[j].a;
                        T.a.dev = R[i].v;
                        T.a.ld = static_cast<int64_t>(R[i].rows.size());
                        T.a.gathered = true;
                        T.a.rows_ext = &rp[i].t;
                        T.b.dev = C[j].v;
                        T.b.ld = static_cast<int64_t>(C[j].rows.size());
                        T.b.gathered = true;
                        T.b.rows_ext = &cp[j].t;
                        T.k1 = static_cast<int>(rp[i].a.size());
                        T.k2 = static_cast<int>(cp[j].a.size());
                        P.terms.push_back(std::move(T));
                        pairs.push_back(std::move(P));
                        info.push_back({R[i].inter, &R[i].members, &C[j].members});
                    }
                }
                run_pairs(pairs, info);
            }
            mark(kRoot);
            break;
        }
        const auto& SL = store.levels[static_cast<size_t>(l)];
        auto sbox = [&](int64_t bi) -> const fmm::gpu::PartialDeviceBox* {
            if (static_cast<size_t>(bi) >= SL.box_of.size() || SL.box_of[static_cast<size_t>(bi)] < 0) return nullptr;
            return &SL.boxes[static_cast<size_t>(SL.box_of[static_cast<size_t>(bi)])];
        };
        // the slot maps of the boxes the blocks touch (as the host's)
        std::unordered_map<int64_t, std::pair<std::vector<int32_t>, std::vector<int32_t>>> maps;
        for (const auto* side : {&R, &C})
            for (const auto& b : *side) maps.emplace(b.box, std::pair<std::vector<int32_t>, std::vector<int32_t>>{});
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
        // ---- the chain steps (partial_chain_parts): plans here, values there
        struct PartPlan { std::vector<int64_t> rp, sp, dst; bool dense = false; int64_t k11 = 0, k22 = 0; };
        std::vector<PartPlan> plans[2];
        for (int sd = 0; sd < 2; ++sd) {
            auto& blocks = sides[sd];
            plans[sd].resize(blocks.size());
#pragma omp parallel for schedule(dynamic, 16)
            for (size_t i = 0; i < blocks.size(); ++i) {
                Chain& b = blocks[i];
                PartPlan& P = plans[sd][i];
                const auto& box = level.local_boxes[static_cast<size_t>(b.box)];
                const auto& m = maps.at(b.box);
                const auto& f = box.f22;
                const bool own = f.sparsified && !f.redundant.empty();
                const int64_t kp = static_cast<int64_t>(box.skeleton_indices.size());
                std::vector<int64_t> sidx;
                for (size_t t = 0; t < b.rows.size(); ++t) {
                    const int64_t sl = b.rows[t];
                    if (own && m.second[static_cast<size_t>(sl)] >= 0) {
                        P.rp.push_back(static_cast<int64_t>(t));
                        b.ridx.push_back(m.second[static_cast<size_t>(sl)]);
                    } else if (m.first[static_cast<size_t>(sl)] >= 0) {
                        P.sp.push_back(static_cast<int64_t>(t));
                        sidx.push_back(m.first[static_cast<size_t>(sl)]);
                    }
                }
                if (!P.rp.empty()) {
                    P.dense = true;
                    P.k22 = static_cast<int64_t>(f.skeleton.size());
                    P.k11 = kp - P.k22;
                    b.wrows.resize(static_cast<size_t>(P.k22));
                    for (int64_t q = 0; q < P.k22; ++q) b.wrows[static_cast<size_t>(q)] = P.k11 + q;
                    std::vector<int64_t> src;
                    for (size_t q = 0; q < P.sp.size(); ++q) {
                        const int64_t x = sidx[q] - P.k11;
                        if (x < 0) continue;  // an S11 row: zero in an F22 chain
                        P.dst.push_back(x);
                        src.push_back(P.sp[q]);
                    }
                    P.sp = std::move(src);
                } else {
                    b.wrows = sidx;
                }
                b.rowpos.assign(static_cast<size_t>(box.num_points), -1);
                for (size_t t = 0; t < b.rows.size(); ++t) b.rowpos[static_cast<size_t>(b.rows[t])] = static_cast<int32_t>(t);
                b.wpos.assign(static_cast<size_t>(kp), -1);
                for (size_t t = 0; t < b.wrows.size(); ++t) b.wpos[static_cast<size_t>(b.wrows[t])] = static_cast<int32_t>(t);
            }
        }
        mark(kPlan);
        {
            size_t bytes = 0;
            for (int sd = 0; sd < 2; ++sd) {
                for (size_t i = 0; i < sides[sd].size(); ++i) {
                    const Chain& b = sides[sd][i];
                    const size_t n = b.members.size();
                    bytes += fmm::gpu::align_up(b.ridx.size() * n * D) + fmm::gpu::align_up(b.wrows.size() * n * D);
                    if (plans[sd][i].dense) bytes += fmm::gpu::align_up(static_cast<size_t>(plans[sd][i].k22) * b.ridx.size() * D);
                }
            }
            char* arena = work.alloc(bytes);
            arenas.push_back(arena);
            work.zero(arena, bytes);
            fmm::gpu::MetaBuilder& meta = work.meta();
            std::vector<fmm::gpu::GatherItemT<S>> gathers;
            std::vector<fmm::gpu::IndexCopyItemT<S>> copies;
            fmm::gpu::VBatch<DataType> steps;
            int gm = 0, gn = 0, cm = 0, cn = 0;
            auto list = [&](const std::vector<int64_t>& v) {
                std::vector<int> ints(v.begin(), v.end());
                fmm::gpu::IndexList l;
                l.offset = static_cast<int64_t>(meta.append(ints));
                return l;
            };
            size_t at = 0;
            for (int sd = 0; sd < 2; ++sd) {
                for (size_t i = 0; i < sides[sd].size(); ++i) {
                    Chain& b = sides[sd][i];
                    const PartPlan& P = plans[sd][i];
                    const int n = static_cast<int>(b.members.size()), nr = static_cast<int>(b.rows.size());
                    const int nrp = static_cast<int>(b.ridx.size()), nw = static_cast<int>(b.wrows.size());
                    b.vp = reinterpret_cast<DataType*>(arena + at);
                    at += fmm::gpu::align_up(static_cast<size_t>(nrp) * static_cast<size_t>(n) * D);
                    b.wv = reinterpret_cast<DataType*>(arena + at);
                    at += fmm::gpu::align_up(static_cast<size_t>(nw) * static_cast<size_t>(n) * D);
                    if (n == 0) continue;
                    if (nrp > 0) {
                        fmm::gpu::GatherItemT<S> g{};
                        g.out = reinterpret_cast<S*>(b.vp);
                        g.ld = nrp;
                        g.m = nrp;
                        g.n = n;
                        g.src = reinterpret_cast<const S*>(b.v);
                        g.rs = 1;
                        g.cs = nr;
                        g.rows = list(P.rp);
                        gathers.push_back(g);
                        gm = std::max(gm, nrp);
                        gn = std::max(gn, n);
                    }
                    if (P.dense) {
                        // wv = (the S22 rows of v, zero elsewhere) + T22(S22, ridx) vp
                        if (!P.dst.empty()) {
                            fmm::gpu::IndexCopyItemT<S> c{};
                            c.out = reinterpret_cast<S*>(b.wv);
                            c.ors = 1;
                            c.ocs = nw;
                            c.out_rows = list(P.dst);
                            c.src = reinterpret_cast<const S*>(b.v);
                            c.rs = 1;
                            c.cs = nr;
                            c.src_rows = list(P.sp);
                            c.m = static_cast<int>(P.dst.size());
                            c.n = n;
                            c.add = 0;
                            copies.push_back(c);
                            cm = std::max(cm, c.m);
                            cn = std::max(cn, n);
                        }
                        const auto* B = sbox(b.box);
                        if (B == nullptr || B->t22.dev == nullptr) throw std::logic_error("GPU partial CB entries: no T22 on the device");
                        DataType* tsub = reinterpret_cast<DataType*>(arena + at);
                        at += fmm::gpu::align_up(static_cast<size_t>(P.k22) * static_cast<size_t>(nrp) * D);
                        fmm::gpu::GatherItemT<S> g{};
                        g.out = reinterpret_cast<S*>(tsub);
                        g.ld = static_cast<int>(P.k22);
                        g.m = static_cast<int>(P.k22);
                        g.n = nrp;
                        g.src = reinterpret_cast<const S*>(B->t22.dev);
                        g.rs = 1;
                        g.cs = B->t22.rows;
                        g.rows.base = static_cast<int>(P.k11);
                        g.cols = list(b.ridx);
                        gathers.push_back(g);
                        gm = std::max(gm, static_cast<int>(P.k22));
                        gn = std::max(gn, nrp);
                        steps.entries.push_back({tsub, b.vp, b.wv, static_cast<magma_int_t>(P.k22), n, nrp,
                                                 static_cast<magma_int_t>(P.k22), nrp, static_cast<magma_int_t>(P.k22)});
                    } else if (nw > 0) {
                        fmm::gpu::GatherItemT<S> g{};
                        g.out = reinterpret_cast<S*>(b.wv);
                        g.ld = nw;
                        g.m = nw;
                        g.n = n;
                        g.src = reinterpret_cast<const S*>(b.v);
                        g.rs = 1;
                        g.cs = nr;
                        g.rows = list(P.sp);
                        gathers.push_back(g);
                        gm = std::max(gm, nw);
                        gn = std::max(gn, n);
                    }
                }
            }
            const size_t goff = meta.append(gathers), coff2 = meta.append(copies);
            steps.stage(meta);
            char* md = work.upload_meta();
            fmm::gpu::launch_gather(reinterpret_cast<const fmm::gpu::GatherItemT<S>*>(md + goff), static_cast<int>(gathers.size()), gm,
                                    gn, md, stream);
            fmm::gpu::launch_index_copy(reinterpret_cast<const fmm::gpu::IndexCopyItemT<S>*>(md + coff2),
                                        static_cast<int>(copies.size()), cm, cn, md, stream);
            if (steps.count() > 0) steps.gemm(md, MagmaNoTrans, MagmaNoTrans, DataType(1), DataType(1), work.queue());
        }
        mark(kSteps);
        // ---- the level's terms
        if (l >= 2) {
            // this rank's blocks to the owners of their boxes' remote neighbours (row
            // blocks their w part and the box's live set, column blocks their raw
            // values): the index lists through the host, the values device to device
            std::vector<std::vector<char>> send(static_cast<size_t>(nranks));
            std::vector<std::vector<std::pair<const DataType*, size_t>>> send_vals(static_cast<size_t>(nranks));
            if (active) {
                for (int side_id = 0; side_id < 2; ++side_id) {
                    const bool row = side_id == 0;
                    for (const auto& b : sides[side_id]) {
                        const auto& box = level.local_boxes[static_cast<size_t>(b.box)];
                        std::set<int> dest;
                        for (int64_t nm : box.one_hop)
                            if (!level.is_box_on_process(nm)) dest.insert(level.owner_rank_of(nm));
                        if (dest.empty()) continue;
                        std::vector<char> msg;
                        const int64_t head[8] = {side_id, box.morton_index, box.num_points, static_cast<int64_t>(b.members.size()),
                                                 row ? 0 : static_cast<int64_t>(b.rows.size()), row ? static_cast<int64_t>(b.wrows.size()) : 0,
                                                 row ? static_cast<int64_t>(box.skeleton_indices.size()) : 0, b.inter};
                        partial_put(msg, head, 8);
                        partial_put(msg, b.members.data(), b.members.size());
                        if (row) {
                            partial_put(msg, b.wrows.data(), b.wrows.size());
                            partial_put(msg, box.skeleton_indices.data(), box.skeleton_indices.size());
                        } else {
                            partial_put(msg, b.rows.data(), b.rows.size());
                        }
                        const size_t e = (row ? b.wrows.size() : b.rows.size()) * b.members.size();
                        for (int r : dest) {
                            send[static_cast<size_t>(r)].insert(send[static_cast<size_t>(r)].end(), msg.begin(), msg.end());
                            send_vals[static_cast<size_t>(r)].emplace_back(row ? b.wv : b.v, e);
                            sent[0] += static_cast<double>(e * D);
                        }
                    }
                }
            }
            for (const auto& s : send) sent[0] += static_cast<double>(s.size());
            mark(kTermsOut);
            const auto recv = partial_exchange(tree.comm, send, 620);
            mark(kTermsMpi);
            std::vector<RemoteD> remote_rows, remote_cols;
            std::vector<std::pair<int, size_t>> remote_at;  // (side, index) of each received block, in order
            std::vector<int64_t> rcount(static_cast<size_t>(nranks), 0);
            size_t nvals = 0;
            for (size_t from = 0; from < recv.size(); ++from) {
                const auto& buf = recv[from];
                for (const char* p = buf.data(); p < buf.data() + buf.size();) {
                    int64_t head[8];
                    p = partial_get(p, head, 8);
                    RemoteD rb;
                    rb.morton = head[1];
                    rb.num_points = head[2];
                    rb.inter = head[7];
                    const int64_t nm = head[3], nr = head[4], nw = head[5], kp = head[6];
                    rb.members.resize(static_cast<size_t>(nm));
                    p = partial_get(p, rb.members.data(), rb.members.size());
                    int64_t e = 0;
                    if (head[0] == 0) {
                        rb.wrows.resize(static_cast<size_t>(nw));
                        p = partial_get(p, rb.wrows.data(), rb.wrows.size());
                        rb.skeleton.resize(static_cast<size_t>(kp));
                        p = partial_get(p, rb.skeleton.data(), rb.skeleton.size());
                        rb.wpos.assign(static_cast<size_t>(kp), -1);
                        for (int64_t t = 0; t < nw; ++t) rb.wpos[static_cast<size_t>(rb.wrows[static_cast<size_t>(t)])] = static_cast<int32_t>(t);
                        e = nw * nm;
                    } else {
                        rb.rows.resize(static_cast<size_t>(nr));
                        p = partial_get(p, rb.rows.data(), rb.rows.size());
                        rb.rowpos.assign(static_cast<size_t>(rb.num_points), -1);
                        for (int64_t t = 0; t < nr; ++t) rb.rowpos[static_cast<size_t>(rb.rows[static_cast<size_t>(t)])] = static_cast<int32_t>(t);
                        e = nr * nm;
                    }
                    // (the values arrive in rank order, a rank's in its lists' order)
                    rb.voff = nvals;
                    nvals += static_cast<size_t>(e);
                    rcount[from] += e;
                    if (head[0] == 0) {
                        remote_at.emplace_back(0, remote_rows.size());
                        remote_rows.push_back(std::move(rb));
                    } else {
                        remote_at.emplace_back(1, remote_cols.size());
                        remote_cols.push_back(std::move(rb));
                    }
                }
            }
            mark(kTermsIn);
            std::vector<int64_t> roff;
            if (char* vals = gpu::partial_values_exchange(work, tree.comm, send_vals, rcount, 622, roff)) {
                arenas.push_back(vals);
                for (const auto& [sd, i] : remote_at) {
                    RemoteD& rb = sd == 0 ? remote_rows[i] : remote_cols[i];
                    (sd == 0 ? rb.wv : rb.v) = reinterpret_cast<DataType*>(vals) + rb.voff;
                }
            }
            mark(kTermsMpi);

            if (active) {
                // the column blocks by box (Morton): this rank's and received (as the host's)
                std::unordered_map<int64_t, std::vector<size_t>> col_by_box, remote_col_by_box;
                for (size_t j = 0; j < C.size(); ++j)
                    col_by_box[level.local_boxes[static_cast<size_t>(C[j].box)].morton_index].push_back(j);
                for (size_t j = 0; j < remote_cols.size(); ++j) remote_col_by_box[remote_cols[j].morton].push_back(j);
                // an X_RN term: X_RN_bi(ri, bj raw slots in the column block's rows) and the block's rows of them
                auto rn_term = [&](const Chain& bi, const typename fmm::BoxData<CoordType, DataType>::F22Part& fi,
                                   const fmm::gpu::PartialDeviceBox& SB, int64_t q, const std::vector<int32_t>& rowpos,
                                   const DataType* vj, int64_t nrj, Pair& P) {
                    const int64_t ri = static_cast<int64_t>(bi.ridx.size());
                    const int64_t m = fi.neighbor_f22_counts[static_cast<size_t>(q)];
                    if (m == 0 || ri == 0) return false;
                    int64_t col0 = 0;
                    for (int64_t t = 0; t < q; ++t) col0 += fi.neighbor_f22_counts[static_cast<size_t>(t)];
                    const auto& slots = fi.neighbor_f22_slots[static_cast<size_t>(q)];
                    Term T;
                    for (int64_t c = 0; c < m; ++c) {
                        const int32_t t = rowpos[static_cast<size_t>(slots[static_cast<size_t>(c)])];
                        if (t >= 0) { T.cols_own.push_back(col0 + c); T.b.rows_own.push_back(t); }
                    }
                    if (T.cols_own.empty()) return false;
                    T.kind = 3;
                    T.rec = &SB.rn;
                    T.ld = SB.r;
                    T.rows_ext = &bi.ridx;
                    T.a.dev = bi.vp;
                    T.a.ld = ri;
                    T.b.dev = vj;
                    T.b.ld = nrj;
                    T.b.gathered = true;
                    T.k1 = static_cast<int>(ri);
                    T.k2 = static_cast<int>(T.cols_own.size());
                    P.terms.push_back(std::move(T));
                    return true;
                };
                // an X_NR term: X_NR_bj(bi raw slots in the row block's w rows, rj) and the row block's rows of wv
                auto nr_term = [&](const std::vector<int32_t>& s_of_i, const std::vector<int32_t>& wpos_i, const DataType* wv_i,
                                   int64_t nwi, int64_t morton_i, const fmm::BoxData<CoordType, DataType>& boxJ, const Chain& bj,
                                   int64_t bj_box, Pair& P) {
                    const auto& fj = boxJ.f22;
                    const auto itb = std::find(boxJ.one_hop.begin(), boxJ.one_hop.end(), morton_i);
                    if (itb == boxJ.one_hop.end()) return false;
                    const size_t idi = static_cast<size_t>(itb - boxJ.one_hop.begin());
                    const int64_t m = fj.neighbor_nr_counts[idi];
                    const auto& slots = fj.neighbor_nr_slots[idi];
                    Term T;
                    for (int64_t c = 0; c < m; ++c) {
                        const int32_t sp = s_of_i[static_cast<size_t>(slots[static_cast<size_t>(c)])];
                        if (sp < 0) continue;
                        const int32_t t = wpos_i[static_cast<size_t>(sp)];
                        if (t >= 0) { T.rows_own.push_back(c); T.a.rows_own.push_back(t); }
                    }
                    if (T.rows_own.empty()) return false;
                    const auto* SBj = sbox(bj_box);
                    if (SBj == nullptr) throw std::logic_error("GPU partial CB entries: no records of a sparsified box");
                    T.kind = 4;
                    T.rec = &SBj->nr[idi];
                    T.ld = m;
                    T.cols_ext = &bj.ridx;
                    T.a.dev = wv_i;
                    T.a.ld = nwi;
                    T.a.gathered = true;
                    T.b.dev = bj.vp;
                    T.b.ld = static_cast<int64_t>(bj.ridx.size());
                    T.k1 = static_cast<int>(T.rows_own.size());
                    T.k2 = static_cast<int>(bj.ridx.size());
                    P.terms.push_back(std::move(T));
                    return true;
                };
                // the row blocks, this rank's and received, in the host's order (by row box)
                struct Item { int64_t group; bool remote; size_t index; };
                std::vector<Item> items;
                for (size_t i = 0; i < R.size(); ++i)
                    items.push_back({level.local_boxes[static_cast<size_t>(R[i].box)].morton_index, false, i});
                for (size_t i = 0; i < remote_rows.size(); ++i)
                    if (!remote_rows[i].wrows.empty()) items.push_back({remote_rows[i].morton, true, i});
                std::stable_sort(items.begin(), items.end(), [](const Item& a, const Item& b) { return a.group < b.group; });
                std::vector<std::vector<Pair>> item_pairs(items.size());
                std::vector<std::vector<PairInfo>> item_info(items.size());
                std::exception_ptr failure;
#pragma omp parallel for schedule(dynamic, 1)
                for (size_t t = 0; t < items.size(); ++t) {
                    try {
                        auto& pairs = item_pairs[t];
                        auto& info = item_info[t];
                        if (!items[t].remote) {
                            const Chain& bi = R[items[t].index];
                            const auto& boxI = level.local_boxes[static_cast<size_t>(bi.box)];
                            const auto& fi = boxI.f22;
                            const auto* SBi = sbox(bi.box);
                            const bool own_i = fi.sparsified && !fi.redundant.empty() && !bi.ridx.empty() && SBi != nullptr;
                            const int ni = static_cast<int>(bi.members.size());
                            const int64_t ri = static_cast<int64_t>(bi.ridx.size());
                            const auto& s_of_i = maps.at(bi.box).first;
                            // the same box
                            if (const auto it = col_by_box.find(boxI.morton_index); it != col_by_box.end()) {
                                const int64_t r = static_cast<int64_t>(fi.redundant.size()), kp = static_cast<int64_t>(boxI.skeleton_indices.size());
                                for (size_t jb : it->second) {
                                    const Chain& bj = C[jb];
                                    if (bj.inter != bi.inter) continue;
                                    Pair P;
                                    P.ni = ni;
                                    P.nj = static_cast<int>(bj.members.size());
                                    const int64_t rj = static_cast<int64_t>(bj.ridx.size()), wj = static_cast<int64_t>(bj.wrows.size());
                                    bool any = false;
                                    if (own_i) {
                                        if (rj > 0) {
                                            Term T;
                                            T.kind = 0;
                                            T.rec = &SBi->rr;
                                            T.ld = r;
                                            T.rows_ext = &bi.ridx;
                                            T.cols_ext = &bj.ridx;
                                            T.a.dev = bi.vp;
                                            T.a.ld = ri;
                                            T.b.dev = bj.vp;
                                            T.b.ld = rj;
                                            T.k1 = static_cast<int>(ri);
                                            T.k2 = static_cast<int>(rj);
                                            P.terms.push_back(std::move(T));
                                        }
                                        if (wj > 0) {
                                            Term T;
                                            T.kind = 1;
                                            T.rec = &SBi->rs;
                                            T.ld = r;
                                            T.rows_ext = &bi.ridx;
                                            T.cols_ext = &bj.wrows;
                                            T.a.dev = bi.vp;
                                            T.a.ld = ri;
                                            T.b.dev = bj.wv;
                                            T.b.ld = wj;
                                            T.k1 = static_cast<int>(ri);
                                            T.k2 = static_cast<int>(wj);
                                            P.terms.push_back(std::move(T));
                                        }
                                        any = true;
                                    }
                                    const int64_t wi = static_cast<int64_t>(bi.wrows.size());
                                    if (fi.sparsified && !fi.redundant.empty() && SBi != nullptr && rj > 0 && wi > 0) {
                                        Term T;
                                        T.kind = 2;
                                        T.rec = &SBi->sr;
                                        T.ld = kp;
                                        T.rows_ext = &bi.wrows;
                                        T.cols_ext = &bj.ridx;
                                        T.a.dev = bi.wv;
                                        T.a.ld = wi;
                                        T.b.dev = bj.vp;
                                        T.b.ld = rj;
                                        T.k1 = static_cast<int>(wi);
                                        T.k2 = static_cast<int>(rj);
                                        P.terms.push_back(std::move(T));
                                        any = true;
                                    }
                                    if (any) {
                                        pairs.push_back(std::move(P));
                                        info.push_back({bi.inter, &bi.members, &bj.members});
                                    }
                                }
                            }
                            // the one-hop neighbours
                            for (size_t q = 0; q < boxI.one_hop.size(); ++q) {
                                const int64_t nm = boxI.one_hop[q];
                                if (level.is_box_on_process(nm)) {
                                    const auto it = col_by_box.find(nm);
                                    if (it == col_by_box.end()) continue;
                                    const int64_t nl = nm - level.local_morton_start;
                                    const auto& boxJ = level.local_boxes[static_cast<size_t>(nl)];
                                    const bool sp_j = boxJ.f22.sparsified && !boxJ.f22.redundant.empty();
                                    for (size_t jb : it->second) {
                                        const Chain& bj = C[jb];
                                        if (bj.inter != bi.inter) continue;
                                        Pair P;
                                        P.ni = ni;
                                        P.nj = static_cast<int>(bj.members.size());
                                        bool any = false;
                                        if (own_i)
                                            any |= rn_term(bi, fi, *SBi, static_cast<int64_t>(q), bj.rowpos, bj.v,
                                                           static_cast<int64_t>(bj.rows.size()), P);
                                        if (sp_j && !bj.ridx.empty() && !bi.wrows.empty())
                                            any |= nr_term(s_of_i, bi.wpos, bi.wv, static_cast<int64_t>(bi.wrows.size()), boxI.morton_index,
                                                           boxJ, bj, nl, P);
                                        if (any) {
                                            pairs.push_back(std::move(P));
                                            info.push_back({bi.inter, &bi.members, &bj.members});
                                        }
                                    }
                                } else if (own_i) {
                                    // a remote column box: its X_NR terms are its owner's
                                    const auto it = remote_col_by_box.find(nm);
                                    if (it == remote_col_by_box.end()) continue;
                                    for (size_t jb : it->second) {
                                        const RemoteD& bj = remote_cols[jb];
                                        if (bj.inter != bi.inter) continue;
                                        Pair P;
                                        P.ni = ni;
                                        P.nj = static_cast<int>(bj.members.size());
                                        if (rn_term(bi, fi, *SBi, static_cast<int64_t>(q), bj.rowpos, bj.v,
                                                    static_cast<int64_t>(bj.rows.size()), P)) {
                                            pairs.push_back(std::move(P));
                                            info.push_back({bi.inter, &bi.members, &bj.members});
                                        }
                                    }
                                }
                            }
                        } else {
                            // a received row block: its X_NR terms with the column blocks of this
                            // rank's boxes among its box's neighbours (they list it)
                            const RemoteD& bi = remote_rows[items[t].index];
                            std::vector<int32_t> s_of_i(static_cast<size_t>(bi.num_points), -1);
                            for (size_t q = 0; q < bi.skeleton.size(); ++q) s_of_i[static_cast<size_t>(bi.skeleton[q])] = static_cast<int32_t>(q);
                            for (const auto& kv : col_by_box) {
                                const int64_t nl = kv.first - level.local_morton_start;
                                const auto& boxJ = level.local_boxes[static_cast<size_t>(nl)];
                                if (!(boxJ.f22.sparsified && !boxJ.f22.redundant.empty())) continue;
                                if (std::find(boxJ.one_hop.begin(), boxJ.one_hop.end(), bi.morton) == boxJ.one_hop.end()) continue;
                                for (size_t jb : kv.second) {
                                    const Chain& bj = C[jb];
                                    if (bj.ridx.empty() || bj.inter != bi.inter) continue;
                                    Pair P;
                                    P.ni = static_cast<int>(bi.members.size());
                                    P.nj = static_cast<int>(bj.members.size());
                                    if (nr_term(s_of_i, bi.wpos, bi.wv, static_cast<int64_t>(bi.wrows.size()), bi.morton, boxJ, bj, nl, P)) {
                                        pairs.push_back(std::move(P));
                                        info.push_back({bi.inter, &bi.members, &bj.members});
                                    }
                                }
                            }
                        }
                    } catch (...) {
#pragma omp critical(partial_schur_blocks_failure)
                        if (!failure) failure = std::current_exception();
                    }
                }
                if (failure) std::rethrow_exception(failure);
                std::vector<Pair> pairs;
                std::vector<PairInfo> info;
                for (size_t t = 0; t < items.size(); ++t) {
                    for (auto& P : item_pairs[t]) pairs.push_back(std::move(P));
                    info.insert(info.end(), item_info[t].begin(), item_info[t].end());
                }
                run_pairs(pairs, info);
            }
            mark(kPairs);
        }
        // ---- the parent blocks (as the host's): each block's w part at its box's
        // offset, the dense parts of a box's blocks merged, blocks whose parent is
        // another rank's to it; the values copied on the device
        auto& parent_level = tree.levels[static_cast<size_t>(l) - 1];
        std::vector<std::vector<char>> send_up(static_cast<size_t>(nranks));
        std::vector<std::vector<std::pair<const DataType*, size_t>>> send_up_vals(static_cast<size_t>(nranks));
        std::vector<char*> next_arenas;
        std::vector<Chain> next[2];
        for (int side_id = 0; side_id < 2; ++side_id) {
            auto& blocks = sides[side_id];
            std::vector<Chain> up;
            std::vector<int64_t> up_parent, width;
            std::vector<int> sources;
            std::vector<std::pair<int64_t, int64_t>> place(blocks.size(), {-1, 0});
            std::unordered_map<int64_t, size_t> dense_of;
            for (size_t s2 = 0; s2 < blocks.size(); ++s2) {
                const auto& b = blocks[s2];
                const auto& box = level.local_boxes[static_cast<size_t>(b.box)];
                if (b.wrows.empty()) continue;
                const int64_t morton = box.morton_index;
                const int64_t nmem = static_cast<int64_t>(b.members.size());
                if (!b.ridx.empty()) {
                    auto it = dense_of.find(b.box * ninter + b.inter);
                    if (it != dense_of.end()) {
                        place[s2] = {static_cast<int64_t>(it->second), width[it->second]};
                        width[it->second] += nmem;
                        ++sources[it->second];
                        continue;
                    }
                    dense_of.emplace(b.box * ninter + b.inter, up.size());
                }
                const int64_t off = partial_child_offset(level, morton, dimension);
                Chain u;
                u.inter = b.inter;
                u.rows.resize(b.wrows.size());
                for (size_t t = 0; t < b.wrows.size(); ++t) u.rows[t] = off + b.wrows[t];
                place[s2] = {static_cast<int64_t>(up.size()), 0};
                up.push_back(std::move(u));
                up_parent.push_back(morton >> dimension);
                width.push_back(nmem);
                sources.push_back(1);
            }
            size_t bytes = 0;
            for (size_t k = 0; k < up.size(); ++k) bytes += fmm::gpu::align_up(up[k].rows.size() * static_cast<size_t>(width[k]) * D);
            char* arena = work.alloc(bytes);
            next_arenas.push_back(arena);
            size_t at = 0;
            for (size_t k = 0; k < up.size(); ++k) {
                up[k].v = reinterpret_cast<DataType*>(arena + at);
                up[k].members.resize(static_cast<size_t>(width[k]));
                at += fmm::gpu::align_up(up[k].rows.size() * static_cast<size_t>(width[k]) * D);
            }
            fmm::gpu::MetaBuilder& meta = work.meta();
            std::vector<fmm::gpu::CopySpan> spans;
            int64_t max_words = 0;
            for (size_t s2 = 0; s2 < blocks.size(); ++s2) {
                if (place[s2].first < 0) continue;
                const auto& b = blocks[s2];
                auto& u = up[static_cast<size_t>(place[s2].first)];
                const int64_t c0 = place[s2].second;
                std::copy(b.members.begin(), b.members.end(), u.members.begin() + c0);
                const int64_t words = static_cast<int64_t>(b.wrows.size() * b.members.size() * D / 4);
                if (words == 0) continue;
                spans.push_back({b.wv, u.v + c0 * static_cast<int64_t>(u.rows.size()), words});
                max_words = std::max(max_words, words);
            }
            const size_t soff = meta.append(spans);
            char* md = work.upload_meta();
            fmm::gpu::launch_copy_spans(reinterpret_cast<const fmm::gpu::CopySpan*>(md + soff), static_cast<int>(spans.size()), max_words,
                                        stream);
            for (size_t k = 0; k < up.size(); ++k) {
                const int64_t parent = up_parent[k];
                if (parent_level.is_process_active && parent_level.is_box_on_process(parent)) {
                    up[k].box = parent - parent_level.local_morton_start;
                    next[side_id].push_back(std::move(up[k]));
                    continue;
                }
                const auto& u = up[k];
                const int dest = level.parent_level_owner;
                if (dest < 0) throw std::runtime_error("partial Schur: no owner of a parent box");
                auto& buf = send_up[static_cast<size_t>(dest)];
                const int64_t head[5] = {side_id, parent, static_cast<int64_t>(u.members.size()), static_cast<int64_t>(u.rows.size()),
                                         u.inter};
                partial_put(buf, head, 5);
                partial_put(buf, u.members.data(), u.members.size());
                partial_put(buf, u.rows.data(), u.rows.size());
                const size_t e = u.rows.size() * u.members.size();
                send_up_vals[static_cast<size_t>(dest)].emplace_back(u.v, e);
                sent[1] += static_cast<double>(e * D);
            }
        }
        for (const auto& s : send_up) sent[1] += static_cast<double>(s.size());
        mark(kUpOut);
        const auto recv_up = partial_exchange(tree.comm, send_up, 621);
        mark(kUpMpi);
        {
            std::vector<std::pair<int, size_t>> got;  // (side, index in next) of each, values at their offsets
            std::vector<size_t> voff;
            std::vector<int64_t> rcount(static_cast<size_t>(nranks), 0);
            size_t nvals = 0;
            for (size_t from = 0; from < recv_up.size(); ++from) {
                const auto& buf = recv_up[from];
                for (const char* p = buf.data(); p < buf.data() + buf.size();) {
                    int64_t head[5];
                    p = partial_get(p, head, 5);
                    Chain u;
                    u.box = head[1] - parent_level.local_morton_start;
                    u.inter = head[4];
                    u.members.resize(static_cast<size_t>(head[2]));
                    u.rows.resize(static_cast<size_t>(head[3]));
                    p = partial_get(p, u.members.data(), u.members.size());
                    p = partial_get(p, u.rows.data(), u.rows.size());
                    voff.push_back(nvals);
                    nvals += static_cast<size_t>(head[2] * head[3]);
                    rcount[from] += head[2] * head[3];
                    got.emplace_back(static_cast<int>(head[0]), next[head[0]].size());
                    next[head[0]].push_back(std::move(u));
                }
            }
            mark(kUpIn);
            std::vector<int64_t> roff;
            if (char* vals = gpu::partial_values_exchange(work, tree.comm, send_up_vals, rcount, 623, roff)) {
                next_arenas.push_back(vals);
                for (size_t g = 0; g < got.size(); ++g)
                    next[got[g].first][got[g].second].v = reinterpret_cast<DataType*>(vals) + voff[g];
            }
        }
        work.sync();
        mark(kUpMpi);
        for (char* a : arenas) work.release(a);
        arenas = std::move(next_arenas);
        R = std::move(next[0]);
        C = std::move(next[1]);
    }
    if (trace) {
        constexpr int kRow = kPhases + 3;
        double row[kRow];
        std::copy(tp, tp + kPhases, row);
        row[kPhases] = sent[0] / 1e6;
        row[kPhases + 1] = sent[1] / 1e6;
        row[kPhases + 2] = MPI_Wtime() - tcall;
        std::vector<double> all(me == 0 ? static_cast<size_t>(nranks) * kRow : 0);
        MPI_Gather(row, kRow, MPI_DOUBLE, all.data(), kRow, MPI_DOUBLE, 0, tree.comm);
        if (me == 0) {
            for (int r = 0; r < nranks; ++r) {
                const double* x = all.data() + static_cast<size_t>(r) * kRow;
                std::printf("[partial CB entries] rank %d: %.4f s - plans %.4f, chain steps %.4f, terms' exchange (out %.4f, "
                            "MPI %.4f, in %.4f; %.1f MB sent), pairs %.4f, parents' exchange (out %.4f, MPI %.4f, in %.4f; "
                            "%.1f MB sent), root %.4f\n",
                            r, x[kPhases + 2], x[kPlan], x[kSteps], x[kTermsOut], x[kTermsMpi], x[kTermsIn], x[kPhases],
                            x[kPairs], x[kUpOut], x[kUpMpi], x[kUpIn], x[kPhases + 1], x[kRoot]);
            }
            std::fflush(stdout);
        }
    }
}
#endif  // H2_HAVE_GPU

// The blocks CB(I_n, J_n), n < ninter, of a request: I_n is nrows[n] F22
// indices (0-based) after the earlier blocks' in rows, J_n likewise in cols,
// and block n goes column major into out after the earlier blocks.
// Collective over the tree's communicator, every rank passing the same
// request; `replicate`: the result on every rank, else each rank's out holds
// its contributions, which sum to it over the ranks.
template<typename CoordType, typename DataType>
void partial_schur_blocks(fmm::ParallelTree<CoordType, DataType>& tree, int64_t ninter, const int64_t* rows,
                          const int64_t* nrows, const int64_t* cols, const int64_t* ncols, DataType* out,
                          bool replicate, int extract_mode = 0) {
    // A rank holds the chain blocks of its boxes, one a box and block of the
    // request (blocks of different ones never pair); a level's terms are
    // computed once: a box's own terms and its X_RN terms by its owner (the
    // neighbours' column blocks received), its X_NR terms by the owner of the
    // column box (the neighbours' row blocks received), the root's term by
    // the root's owner; the blocks go to the parents' owner between levels,
    // and a sum assembles the result.
    using Block = PartialChainBlock<DataType>;
    auto& pr = tree.partial_root;
    if (pr.cb_released) throw std::runtime_error("partial Schur: released (c_bpack_partial_schur_release)");
    const int64_t pc = tree.partial_count, n22 = tree.num_points - pc;
    const int dimension = tree.dimension;
    const int leaf = tree.num_levels - 1;
    int me = 0, nranks = 1;
    MPI_Comm_rank(tree.comm, &me);
    MPI_Comm_size(tree.comm, &nranks);
    std::vector<int64_t> roff(static_cast<size_t>(ninter) + 1, 0), coff(roff), ooff(roff);
    for (int64_t k = 0; k < ninter; ++k) {
        roff[static_cast<size_t>(k) + 1] = roff[static_cast<size_t>(k)] + nrows[k];
        coff[static_cast<size_t>(k) + 1] = coff[static_cast<size_t>(k)] + ncols[k];
        ooff[static_cast<size_t>(k) + 1] = ooff[static_cast<size_t>(k)] + nrows[k] * ncols[k];
    }
    std::fill(out, out + ooff.back(), DataType(0));
    if (ooff.back() == 0) return;
    auto reduce = [&] {
        if (!replicate || nranks <= 1) return;
        for (int64_t at = 0; at < ooff.back();) {  // (counts are ints)
            const int64_t n = std::min<int64_t>(ooff.back() - at, int64_t{1} << 30);
            MPI_Allreduce(MPI_IN_PLACE, out + at, static_cast<int>(n), partial_mpi_type<DataType>(), MPI_SUM, tree.comm);
            at += n;
        }
    };
#ifdef H2_HAVE_GPU
    // H2_partial_extract: 2 the host; 1 the device when the records are on it;
    // 0 (auto) also when the call is large enough (the same on every rank:
    // the request and the option are).  BPACK_CHECK=solve: the host's too.
    if constexpr (fmm::gpu::gpu_data_type<DataType>) {
        if (extract_mode != 2 &&
            (extract_mode == 1 || ooff.back() >= kPartialExtractDeviceEntries || tree.partial_records_on_device) &&
            fmm::gpu::partial_device_ready(&tree, 0)) {
            partial_schur_blocks_device(tree, ninter, rows, nrows, cols, ncols, out, roff, coff, ooff);
            if (fmm::gpu::device_solve_check()) {
                std::vector<DataType> host(static_cast<size_t>(ooff.back()));
                partial_schur_blocks(tree, ninter, rows, nrows, cols, ncols, host.data(), false, 2);
                double sums[2] = {0.0, 0.0};
                for (size_t i = 0; i < host.size(); ++i) {
                    sums[0] += std::norm(out[i] - host[i]);
                    sums[1] += std::norm(host[i]);
                }
                MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_DOUBLE, MPI_SUM, tree.comm);
                if (me == 0) {
                    std::printf("GPU partial CB entries check: |gpu - host| / |host| = %.3e (%lld entries)\n",
                                sums[1] > 0.0 ? std::sqrt(sums[0] / sums[1]) : 0.0, static_cast<long long>(ooff.back()));
                    std::fflush(stdout);
                }
            }
            reduce();
            return;
        }
    }
#endif
    partial_records_to_host(tree, false);  // (the host reads the records)
    // G (ni x nj) += A^T X B, A k1 x ni, X k1 x k2, B k2 x nj (each its rows the leading dimension), in the
    // cheaper order: a block's few members against another's many take A^T X first
    auto triple = [&](int64_t ni, int64_t nj, int64_t k1, int64_t k2, const DataType* A, const DataType* X,
                      const DataType* B, DataType* G, std::vector<DataType>& tmp) {
        if (ni == 0 || nj == 0 || k1 == 0 || k2 == 0) return;
        const double xb = static_cast<double>(k1) * static_cast<double>(nj) * static_cast<double>(k2 + ni);
        const double ax = static_cast<double>(ni) * static_cast<double>(k2) * static_cast<double>(k1 + nj);
        if (xb <= ax) {
            tmp.resize(static_cast<size_t>(k1 * nj));
            partial_gemm<DataType>('N', 'N', k1, nj, k2, DataType(1), X, k1, B, k2, DataType(0), tmp.data(), k1);
            partial_gemm<DataType>('T', 'N', ni, nj, k1, DataType(1), A, k1, tmp.data(), k1, DataType(1), G, ni);
        } else {
            tmp.resize(static_cast<size_t>(ni * k2));
            partial_gemm<DataType>('T', 'N', ni, k2, k1, DataType(1), A, k1, X, k1, DataType(0), tmp.data(), ni);
            partial_gemm<DataType>('N', 'N', ni, nj, k2, DataType(1), tmp.data(), ni, B, k2, DataType(1), G, ni);
        }
    };
    // X (nr x nc) = rec(ri, ci), ld the record's rows
    fmm::RecordChunkCache<DataType> chunks;
    using Reader = typename fmm::RecordChunkCache<DataType>::Reader;
    auto gather = [&](const fmm::PartialRecord<DataType>& rec, int64_t ld, const int64_t* ri, int64_t nr, const int64_t* ci,
                      int64_t nc, std::vector<DataType>& X, int64_t group, Reader& reader) {
        X.resize(static_cast<size_t>(nr * nc));
        if (!rec.compressed()) {
            const DataType* d = rec.data.data();
            for (int64_t c = 0; c < nc; ++c) {
                const DataType* col = d + ci[c] * ld;
                DataType* x = X.data() + c * nr;
                for (int64_t a = 0; a < nr; ++a) x[a] = col[ri[a]];
            }
        } else {
            for (int64_t c = 0; c < nc; ++c)
                for (int64_t a = 0; a < nr; ++a) X[static_cast<size_t>(a + c * nr)] = chunks.at(rec, static_cast<size_t>(ri[a] + ci[c] * ld), group, reader);
        }
    };
    partial_leaf_index(tree);
    // the two sides (partial_cb_side): the chain blocks of this rank's boxes at the current level
    using Side = PartialCbSide<DataType>;
    Side R = partial_cb_side(tree, ninter, rows, roff, true), C = partial_cb_side(tree, ninter, cols, coff, true);
    // G (|members_i| x |members_j|) into block k of out at the requested positions
    auto scatter = [&](int64_t k, const std::vector<int64_t>& mi, const std::vector<int64_t>& mj, const std::vector<DataType>& G) {
        const int64_t ni = static_cast<int64_t>(mi.size()), ld = nrows[k];
        DataType* ok = out + ooff[static_cast<size_t>(k)];
        for (size_t c = 0; c < mj.size(); ++c) {
            const size_t m = static_cast<size_t>(mj[c]);
            for (int64_t pj = C.pos_start[m]; pj < C.pos_start[m + 1]; ++pj) {
                DataType* col = ok + C.pos[static_cast<size_t>(pj)] * ld;
                for (int64_t r = 0; r < ni; ++r) {
                    const DataType g = G[static_cast<size_t>(r + static_cast<int64_t>(c) * ni)];
                    const size_t mr = static_cast<size_t>(mi[static_cast<size_t>(r)]);
                    for (int64_t pi = R.pos_start[mr]; pi < R.pos_start[mr + 1]; ++pi) col[R.pos[static_cast<size_t>(pi)]] += g;
                }
            }
        }
    };
    // a block of another rank's box: its members, raw rows and values, w part, and the box's live set
    struct Remote {
        int64_t morton = -1, num_points = 0, inter = 0;
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
                std::vector<std::vector<int64_t>> ra(R.blocks.size()), ca(C.blocks.size());
                std::vector<std::vector<DataType>> rv(R.blocks.size()), cv(C.blocks.size());
#pragma omp parallel for schedule(dynamic, 1)
                for (size_t t = 0; t < R.blocks.size() + C.blocks.size(); ++t) {
                    if (t < R.blocks.size()) {
                        f22_part(R.blocks[t], ra[t], rv[t]);
                    } else {
                        const size_t j = t - R.blocks.size();
                        f22_part(C.blocks[j], ca[j], cv[j]);
                    }
                }
                // the pairs of a row and a column block of one block of the request (their entries are disjoint)
                std::vector<std::pair<size_t, size_t>> pairs;
                {
                    std::unordered_map<int64_t, std::vector<size_t>> cols_of;
                    for (size_t j = 0; j < C.blocks.size(); ++j)
                        if (!ca[j].empty()) cols_of[C.blocks[j].inter].push_back(j);
                    for (size_t i = 0; i < R.blocks.size(); ++i) {
                        if (ra[i].empty()) continue;
                        const auto it = cols_of.find(R.blocks[i].inter);
                        if (it != cols_of.end())
                            for (size_t j : it->second) pairs.emplace_back(i, j);
                    }
                }
#pragma omp parallel for schedule(dynamic, 1)
                for (size_t p = 0; p < pairs.size(); ++p) {
                    const size_t i = pairs[p].first, j = pairs[p].second;
                    const int64_t na = static_cast<int64_t>(ra[i].size()), nb = static_cast<int64_t>(ca[j].size());
                    const int64_t ni = static_cast<int64_t>(R.blocks[i].members.size()), nj = static_cast<int64_t>(C.blocks[j].members.size());
                    std::vector<DataType> M(static_cast<size_t>(na * nb)), tmp, G(static_cast<size_t>(ni * nj), DataType(0));
                    for (int64_t c = 0; c < nb; ++c) {
                        const DataType* col = pr.cb.data.data() + ca[j][static_cast<size_t>(c)] * pr.cb.lda;
                        for (int64_t r = 0; r < na; ++r) M[static_cast<size_t>(r + c * na)] = col[ra[i][static_cast<size_t>(r)]];
                    }
                    triple(ni, nj, na, nb, rv[i].data(), M.data(), cv[j].data(), G.data(), tmp);
                    scatter(R.blocks[i].inter, R.blocks[i].members, C.blocks[j].members, G);
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
                        // a row block's X_NR terms read its w part and the box's live set, a
                        // column block's X_RN terms its raw rows and values
                        std::vector<char> msg;
                        const bool row = side_id == 0;
                        const int64_t head[8] = {side_id, box.morton_index, box.num_points, static_cast<int64_t>(b.members.size()),
                                                 row ? 0 : static_cast<int64_t>(b.rows.size()), row ? static_cast<int64_t>(b.wrows.size()) : 0,
                                                 row ? static_cast<int64_t>(box.skeleton_indices.size()) : 0, b.inter};
                        partial_put(msg, head, 8);
                        partial_put(msg, b.members.data(), b.members.size());
                        if (row) {
                            partial_put(msg, b.wrows.data(), b.wrows.size());
                            partial_put(msg, b.wv.data(), b.wv.size());
                            partial_put(msg, box.skeleton_indices.data(), box.skeleton_indices.size());
                        } else {
                            partial_put(msg, b.rows.data(), b.rows.size());
                            partial_put(msg, b.v.data(), b.v.size());
                        }
                        for (int r : dest) send[static_cast<size_t>(r)].insert(send[static_cast<size_t>(r)].end(), msg.begin(), msg.end());
                    }
                }
            }
            const auto recv = partial_exchange(tree.comm, send, 620);
            std::vector<Remote> remote_rows, remote_cols;
            for (const auto& buf : recv) {
                for (const char* p = buf.data(); p < buf.data() + buf.size();) {
                    int64_t head[8];
                    p = partial_get(p, head, 8);
                    Remote rb;
                    rb.morton = head[1];
                    rb.num_points = head[2];
                    rb.inter = head[7];
                    const int64_t nm = head[3], nr = head[4], nw = head[5], kp = head[6];
                    rb.members.resize(static_cast<size_t>(nm));
                    rb.rows.resize(static_cast<size_t>(nr));
                    rb.v.resize(static_cast<size_t>(nr * nm));
                    rb.wrows.resize(static_cast<size_t>(nw));
                    rb.wv.resize(static_cast<size_t>(nw * nm));
                    rb.skeleton.resize(static_cast<size_t>(kp));
                    p = partial_get(p, rb.members.data(), rb.members.size());
                    if (head[0] == 0) {
                        p = partial_get(p, rb.wrows.data(), rb.wrows.size());
                        p = partial_get(p, rb.wv.data(), rb.wv.size());
                        p = partial_get(p, rb.skeleton.data(), rb.skeleton.size());
                        rb.wpos.assign(static_cast<size_t>(kp), -1);
                        for (int64_t t = 0; t < nw; ++t) rb.wpos[static_cast<size_t>(rb.wrows[static_cast<size_t>(t)])] = static_cast<int32_t>(t);
                    } else {
                        p = partial_get(p, rb.rows.data(), rb.rows.size());
                        p = partial_get(p, rb.v.data(), rb.v.size());
                        rb.rowpos.assign(static_cast<size_t>(rb.num_points), -1);
                        for (int64_t t = 0; t < nr; ++t) rb.rowpos[static_cast<size_t>(rb.rows[static_cast<size_t>(t)])] = static_cast<int32_t>(t);
                    }
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
                                   const std::vector<DataType>& vj, int64_t nrj, int64_t nj, std::vector<DataType>& G,
                                   std::vector<DataType>& X, std::vector<DataType>& V, std::vector<DataType>& tmp, int64_t group,
                                   Reader& reader) {
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
                    V.resize(static_cast<size_t>(nc * nj));
                    gather(fi.X_RN, r, bi.ridx.data(), ri, cpos.data(), nc, X, group, reader);
                    for (int64_t c = 0; c < nj; ++c)
                        for (int64_t q = 0; q < nc; ++q) V[static_cast<size_t>(q + c * nc)] = vj[static_cast<size_t>(tpos[static_cast<size_t>(q)] + c * nrj)];
                    triple(ni, nj, ri, nc, bi.vp.data(), X.data(), V.data(), G.data(), tmp);
                    return true;
                };
                // X_NR_bj(bi raw slots in the row block's w rows, rj) times vp_j: G += wv_i^T (X vp_j)
                auto nr_term = [&](const std::vector<int32_t>& s_of_i, const std::vector<int32_t>& wpos_i, const std::vector<DataType>& wv_i,
                                   int64_t nwi, int64_t ni, int64_t morton_i, const fmm::BoxData<CoordType, DataType>& boxJ, const Block& bj,
                                   std::vector<DataType>& G, std::vector<DataType>& X, std::vector<DataType>& V,
                                   std::vector<DataType>& tmp, int64_t group, Reader& reader) {
                    const auto& fj = boxJ.f22;
                    const int64_t nj = static_cast<int64_t>(bj.members.size()), rj = static_cast<int64_t>(bj.ridx.size());
                    const auto itb = std::find(boxJ.one_hop.begin(), boxJ.one_hop.end(), morton_i);
                    if (itb == boxJ.one_hop.end()) return false;
                    const size_t idi = static_cast<size_t>(itb - boxJ.one_hop.begin());
                    const int64_t m = fj.neighbor_nr_counts[idi];
                    const auto& slots = fj.neighbor_nr_slots[idi];
                    std::vector<int64_t> crow, wposl;
                    for (int64_t c = 0; c < m; ++c) {
                        const int32_t sp = s_of_i[static_cast<size_t>(slots[static_cast<size_t>(c)])];
                        if (sp < 0) continue;  // (X_NR keeps bi's live rows only: R_bi is the X_RN term)
                        const int32_t t = wpos_i[static_cast<size_t>(sp)];
                        if (t >= 0) { crow.push_back(c); wposl.push_back(t); }
                    }
                    const int64_t nc = static_cast<int64_t>(crow.size());
                    if (nc == 0) return false;
                    V.resize(static_cast<size_t>(nc * ni));
                    gather(fj.X_NR[idi], m, crow.data(), nc, bj.ridx.data(), rj, X, group, reader);
                    for (int64_t c = 0; c < ni; ++c)
                        for (int64_t q = 0; q < nc; ++q) V[static_cast<size_t>(q + c * nc)] = wv_i[static_cast<size_t>(wposl[static_cast<size_t>(q)] + c * nwi)];
                    triple(ni, nj, nc, rj, V.data(), X.data(), bj.vp.data(), G.data(), tmp);
                    return true;
                };
                // the row blocks, this rank's and received, by box: their own terms, X_RN
                // terms with any column block and X_NR terms with this rank's.  A box's
                // records are read by its row blocks only: their chunks go once the
                // last of them is done.
                struct Item {
                    int64_t group = 0;  // the row box's Morton index
                    bool remote = false;
                    size_t index = 0;
                };
                std::vector<Item> items;
                for (size_t i = 0; i < R.blocks.size(); ++i)
                    items.push_back({level.local_boxes[static_cast<size_t>(R.blocks[i].box)].morton_index, false, i});
                for (size_t i = 0; i < remote_rows.size(); ++i)
                    if (!remote_rows[i].wrows.empty()) items.push_back({remote_rows[i].morton, true, i});
                std::stable_sort(items.begin(), items.end(), [](const Item& a, const Item& b) { return a.group < b.group; });
                std::vector<int64_t> group_of(items.size()), left;
                for (size_t t = 0; t < items.size(); ++t) {
                    if (t == 0 || items[t].group != items[t - 1].group) left.push_back(0);
                    group_of[t] = static_cast<int64_t>(left.size()) - 1;
                    ++left.back();
                }
#pragma omp parallel for schedule(dynamic, 1)
                for (size_t t = 0; t < items.size(); ++t) {
                    const int64_t group = items[t].group;
                    Reader reader;
                    std::vector<DataType> X, V, tmp, G;
                    if (!items[t].remote) {
                        const Block& bi = R.blocks[items[t].index];
                        const auto& boxI = level.local_boxes[static_cast<size_t>(bi.box)];
                        const auto& fi = boxI.f22;
                        const bool own_i = fi.sparsified && !fi.redundant.empty() && !bi.ridx.empty();
                        const int64_t ni = static_cast<int64_t>(bi.members.size()), ri = static_cast<int64_t>(bi.ridx.size());
                        const auto& s_of_i = maps.at(bi.box).first;
                        // the same box
                        if (const auto it = col_by_box.find(boxI.morton_index); it != col_by_box.end()) {
                            const int64_t r = static_cast<int64_t>(fi.redundant.size()), kp = static_cast<int64_t>(boxI.skeleton_indices.size());
                            for (size_t jb : it->second) {
                                const Block& bj = C.blocks[jb];
                                if (bj.inter != bi.inter) continue;
                                const int64_t nj = static_cast<int64_t>(bj.members.size()), rj = static_cast<int64_t>(bj.ridx.size());
                                G.assign(static_cast<size_t>(ni * nj), DataType(0));
                                bool any = false;
                                if (own_i) {
                                    // G += vp_i^T (X_RR(ri, rj) vp_j + X_RS(ri, wj) wv_j)
                                    if (rj > 0) {
                                        gather(fi.X_RR, r, bi.ridx.data(), ri, bj.ridx.data(), rj, X, group, reader);
                                        triple(ni, nj, ri, rj, bi.vp.data(), X.data(), bj.vp.data(), G.data(), tmp);
                                    }
                                    const int64_t wj = static_cast<int64_t>(bj.wrows.size());
                                    if (wj > 0) {
                                        gather(fi.X_RS, r, bi.ridx.data(), ri, bj.wrows.data(), wj, X, group, reader);
                                        triple(ni, nj, ri, wj, bi.vp.data(), X.data(), bj.wv.data(), G.data(), tmp);
                                    }
                                    any = true;
                                }
                                const int64_t wi = static_cast<int64_t>(bi.wrows.size());
                                if (fi.sparsified && !fi.redundant.empty() && rj > 0 && wi > 0) {
                                    // G += wv_i^T X_SR(wi, rj) vp_j
                                    gather(fi.X_SR, kp, bi.wrows.data(), wi, bj.ridx.data(), rj, X, group, reader);
                                    triple(ni, nj, wi, rj, bi.wv.data(), X.data(), bj.vp.data(), G.data(), tmp);
                                    any = true;
                                }
                                if (any) scatter(bi.inter, bi.members, bj.members, G);
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
                                    if (bj.inter != bi.inter) continue;
                                    const int64_t nj = static_cast<int64_t>(bj.members.size());
                                    G.assign(static_cast<size_t>(ni * nj), DataType(0));
                                    bool any = false;
                                    if (own_i)
                                        any |= rn_term(bi, fi, static_cast<int64_t>(q), bj.rowpos, bj.v, static_cast<int64_t>(bj.rows.size()), nj, G,
                                                       X, V, tmp, group, reader);
                                    if (sp_j && !bj.ridx.empty() && !bi.wrows.empty())
                                        any |= nr_term(s_of_i, bi.wpos, bi.wv, static_cast<int64_t>(bi.wrows.size()), ni, boxI.morton_index, boxJ, bj,
                                                       G, X, V, tmp, group, reader);
                                    if (any) scatter(bi.inter, bi.members, bj.members, G);
                                }
                            } else if (own_i) {
                                // a remote column box: its X_NR terms are its owner's
                                const auto it = remote_col_by_box.find(nm);
                                if (it == remote_col_by_box.end()) continue;
                                for (size_t jb : it->second) {
                                    const Remote& bj = remote_cols[jb];
                                    if (bj.inter != bi.inter) continue;
                                    const int64_t nj = static_cast<int64_t>(bj.members.size());
                                    G.assign(static_cast<size_t>(ni * nj), DataType(0));
                                    if (rn_term(bi, fi, static_cast<int64_t>(q), bj.rowpos, bj.v, static_cast<int64_t>(bj.rows.size()), nj, G, X, V,
                                                tmp, group, reader))
                                        scatter(bi.inter, bi.members, bj.members, G);
                                }
                            }
                        }
                    } else {
                        // a received row block: its X_NR terms with the column blocks of this
                        // rank's boxes among its box's neighbours (they list it)
                        const Remote& bi = remote_rows[items[t].index];
                        std::vector<int32_t> s_of_i(static_cast<size_t>(bi.num_points), -1);
                        for (size_t q = 0; q < bi.skeleton.size(); ++q) s_of_i[static_cast<size_t>(bi.skeleton[q])] = static_cast<int32_t>(q);
                        const int64_t ni = static_cast<int64_t>(bi.members.size());
                        for (const auto& kv : col_by_box) {
                            const auto& boxJ = level.local_boxes[static_cast<size_t>(kv.first - level.local_morton_start)];
                            if (!(boxJ.f22.sparsified && !boxJ.f22.redundant.empty())) continue;
                            if (std::find(boxJ.one_hop.begin(), boxJ.one_hop.end(), bi.morton) == boxJ.one_hop.end()) continue;
                            for (size_t jb : kv.second) {
                                const Block& bj = C.blocks[jb];
                                if (bj.ridx.empty() || bj.inter != bi.inter) continue;
                                G.assign(static_cast<size_t>(ni * static_cast<int64_t>(bj.members.size())), DataType(0));
                                if (nr_term(s_of_i, bi.wpos, bi.wv, static_cast<int64_t>(bi.wrows.size()), ni, bi.morton, boxJ, bj, G, X, V, tmp, group,
                                            reader))
                                    scatter(bi.inter, bi.members, bj.members, G);
                            }
                        }
                    }
                    int64_t rest = 0;
#pragma omp atomic capture
                    rest = --left[static_cast<size_t>(group_of[t])];
                    if (rest == 0) chunks.release(group);
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
            auto& blocks = side->blocks;
            // each block's place: its parent block and its first column there
            std::vector<Block> up;
            std::vector<int64_t> up_parent;  // the parents' Morton indices
            std::vector<int64_t> width;      // the parent blocks' members
            std::vector<int> sources;        // and how many blocks make each
            std::vector<std::pair<int64_t, int64_t>> place(blocks.size(), {-1, 0});
            std::unordered_map<int64_t, size_t> dense_of;  // (box at this level, block of the request) -> merged block at the parent
            for (size_t s = 0; s < blocks.size(); ++s) {
                const auto& b = blocks[s];
                const auto& box = level.local_boxes[static_cast<size_t>(b.box)];
                if (b.wrows.empty()) continue;  // nothing left of these chains above
                const int64_t morton = box.morton_index;
                const int64_t nmem = static_cast<int64_t>(b.members.size());
                if (!b.ridx.empty()) {
                    auto it = dense_of.find(b.box * ninter + b.inter);
                    if (it != dense_of.end()) {
                        place[s] = {static_cast<int64_t>(it->second), width[it->second]};
                        width[it->second] += nmem;
                        ++sources[it->second];
                        continue;
                    }
                    dense_of.emplace(b.box * ninter + b.inter, up.size());
                }
                const int64_t off = partial_child_offset(level, morton, dimension);
                const int64_t nw = static_cast<int64_t>(b.wrows.size());
                Block u;
                u.inter = b.inter;
                u.rows.resize(static_cast<size_t>(nw));
                for (int64_t t = 0; t < nw; ++t) u.rows[static_cast<size_t>(t)] = off + b.wrows[static_cast<size_t>(t)];
                place[s] = {static_cast<int64_t>(up.size()), 0};
                up.push_back(std::move(u));
                up_parent.push_back(morton >> dimension);
                width.push_back(nmem);
                sources.push_back(1);
            }
            // a parent block of one block takes its w part, the merged ones copies of theirs
#pragma omp parallel for schedule(dynamic, 1)
            for (size_t k = 0; k < up.size(); ++k) {
                if (sources[k] == 1) continue;
                up[k].members.resize(static_cast<size_t>(width[k]));
                up[k].v.resize(up[k].rows.size() * static_cast<size_t>(width[k]));
            }
#pragma omp parallel for schedule(dynamic, 1)
            for (size_t s = 0; s < blocks.size(); ++s) {
                if (place[s].first < 0) continue;
                auto& b = blocks[s];
                auto& u = up[static_cast<size_t>(place[s].first)];
                if (sources[static_cast<size_t>(place[s].first)] == 1) {
                    u.members = std::move(b.members);
                    u.v = std::move(b.wv);
                } else {
                    const int64_t c0 = place[s].second;
                    std::copy(b.members.begin(), b.members.end(), u.members.begin() + c0);
                    std::copy(b.wv.begin(), b.wv.end(), u.v.begin() + c0 * static_cast<int64_t>(u.rows.size()));
                }
                decltype(b.wv)().swap(b.wv);
                decltype(b.v)().swap(b.v);
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
                const int64_t head[5] = {side_id, parent, static_cast<int64_t>(up[k].members.size()), static_cast<int64_t>(up[k].rows.size()),
                                         up[k].inter};
                partial_put(buf, head, 5);
                partial_put(buf, up[k].members.data(), up[k].members.size());
                partial_put(buf, up[k].rows.data(), up[k].rows.size());
                partial_put(buf, up[k].v.data(), up[k].v.size());
            }
        }
        const auto recv_up = partial_exchange(tree.comm, send_up, 621);
        for (const auto& buf : recv_up) {
            for (const char* p = buf.data(); p < buf.data() + buf.size();) {
                int64_t head[5];
                p = partial_get(p, head, 5);
                Block u;
                u.box = head[1] - parent_level.local_morton_start;
                u.inter = head[4];
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
    reduce();
}

// One block CB(rows, cols) (c_bpack_partial_schur_extract): rows and cols
// F22 indices, 0-based, out nrows x ncols column major, replicated.
template<typename CoordType, typename DataType>
void partial_schur_entries(fmm::ParallelTree<CoordType, DataType>& tree, const int64_t* rows, int64_t nrows,
                           const int64_t* cols, int64_t ncols, DataType* out, int extract_mode = 0) {
    partial_schur_blocks(tree, 1, rows, &nrows, cols, &ncols, out, true, extract_mode);
}

// c_bpack_extractelement(_bridge) on a partial factorization (format 7): the
// Schur complement's blocks CB(I_n, J_n), n < ninter, the CB as an N22 x N22
// matrix (STRUMPACK reads a child's update block so).  allrows / allcols: the
// intersections' 1-based F22 indices (1..N22), stacked, rowidx[n] and
// colidx[n] of them, every caller rank passing all of them; block n to the
// caller's rank owner[n], its blocks in alldat_loc in intersection order,
// each column major.  tree: the factorization's, null on the caller's ranks
// outside it, whose ranks must be a subset of caller_comm's.  Collective over
// caller_comm.  The source ranks compute every block in one pass
// (partial_schur_blocks), each block's sum lands on one of them (its owner
// when the owner holds the factorization), and goes on to the owner.
template<typename CoordType, typename DataType>
void extract_partial_blocks(fmm::ParallelTree<CoordType, DataType>* tree, MPI_Comm caller_comm, int64_t ninter,
                            const int* allrows, const int* allcols, const int* rowidx, const int* colidx,
                            const std::vector<int>& owner, DataType* alldat_loc, int extract_mode = 0) {
    int crank = 0, csize = 1, srank = -1, ssize = 0;
    MPI_Comm_rank(caller_comm, &crank);
    MPI_Comm_size(caller_comm, &csize);
    if (tree != nullptr) {
        MPI_Comm_rank(tree->comm, &srank);
        MPI_Comm_size(tree->comm, &ssize);
    }
    // the source's ranks among the caller's
    std::vector<int> source_of(static_cast<size_t>(csize));  // caller rank -> source rank (-1)
    MPI_Allgather(&srank, 1, MPI_INT, source_of.data(), 1, MPI_INT, caller_comm);
    int nsource = 0;
    for (int r : source_of) nsource += r >= 0 ? 1 : 0;
    if (nsource == 0) throw std::invalid_argument("extract_partial_blocks: no rank holds the factorization");
    std::vector<int> caller_of(static_cast<size_t>(nsource), -1);  // source rank -> caller rank
    for (int c = 0; c < csize; ++c) {
        const int s = source_of[static_cast<size_t>(c)];
        if (s < 0) continue;
        if (s >= nsource || caller_of[static_cast<size_t>(s)] != -1) {
            throw std::invalid_argument("extract_partial_blocks: the factorization's ranks are not a subset of the caller's");
        }
        caller_of[static_cast<size_t>(s)] = c;
    }
    if (tree != nullptr && ssize != nsource) throw std::invalid_argument("extract_partial_blocks: a source rank outside the caller");
    // the source rank each block's sum lands on, and the blocks in that order
    std::vector<int> lands(static_cast<size_t>(ninter));
    std::vector<int64_t> r0(static_cast<size_t>(ninter) + 1, 0), c0(r0);
    for (int64_t n = 0; n < ninter; ++n) {
        if (rowidx[n] < 0 || colidx[n] < 0) throw std::invalid_argument("extract_partial_blocks: a negative block size");
        const int o = owner[static_cast<size_t>(n)];
        if (o < 0 || o >= csize) throw std::invalid_argument("extract_partial_blocks: an owner outside the caller's ranks");
        lands[static_cast<size_t>(n)] = source_of[static_cast<size_t>(o)] >= 0 ? source_of[static_cast<size_t>(o)]
                                                                              : static_cast<int>(n % nsource);
        r0[static_cast<size_t>(n) + 1] = r0[static_cast<size_t>(n)] + rowidx[n];
        c0[static_cast<size_t>(n) + 1] = c0[static_cast<size_t>(n)] + colidx[n];
    }
    std::vector<int64_t> order(static_cast<size_t>(ninter));
    std::iota(order.begin(), order.end(), int64_t{0});
    std::stable_sort(order.begin(), order.end(), [&](int64_t x, int64_t y) { return lands[x] < lands[y]; });
    // where each owned block goes in alldat_loc
    std::vector<int64_t> base(static_cast<size_t>(ninter), -1);
    int64_t mine_total = 0;
    for (int64_t n = 0; n < ninter; ++n) {
        if (owner[static_cast<size_t>(n)] != crank) continue;
        base[static_cast<size_t>(n)] = mine_total;
        mine_total += static_cast<int64_t>(rowidx[n]) * colidx[n];
    }
    std::vector<std::vector<char>> send(static_cast<size_t>(csize));
    if (tree != nullptr) {
        // ---- the blocks in landing order, 0-based F22 indices
        std::vector<int64_t> rows, cols, nr(static_cast<size_t>(ninter)), nc(static_cast<size_t>(ninter));
        std::vector<int64_t> boff(static_cast<size_t>(ninter) + 1, 0);
        for (size_t q = 0; q < order.size(); ++q) {
            const int64_t n = order[q];
            for (int64_t i = r0[static_cast<size_t>(n)]; i < r0[static_cast<size_t>(n) + 1]; ++i) rows.push_back(int64_t{allrows[i]} - 1);
            for (int64_t j = c0[static_cast<size_t>(n)]; j < c0[static_cast<size_t>(n) + 1]; ++j) cols.push_back(int64_t{allcols[j]} - 1);
            nr[q] = rowidx[n];
            nc[q] = colidx[n];
            boff[q + 1] = boff[q] + nr[q] * nc[q];
        }
        std::vector<DataType> part(static_cast<size_t>(std::max<int64_t>(boff.back(), 1)));
        partial_schur_blocks(*tree, ninter, rows.data(), nr.data(), cols.data(), nc.data(), part.data(), false, extract_mode);
        // ---- each block's sum on its landing rank
        std::vector<int> counts(static_cast<size_t>(ssize), 0);
        for (size_t q = 0; q < order.size(); ++q) {
            const int64_t c = counts[static_cast<size_t>(lands[order[q]])] + nr[q] * nc[q];
            if (c > INT_MAX) throw std::runtime_error("extract_partial_blocks: a source rank's blocks over 2^31 entries");
            counts[static_cast<size_t>(lands[order[q]])] = static_cast<int>(c);
        }
        std::vector<DataType> landed(static_cast<size_t>(std::max(counts[static_cast<size_t>(srank)], 1)));
        if (ssize > 1) {
            MPI_Reduce_scatter(part.data(), landed.data(), counts.data(), partial_mpi_type<DataType>(), MPI_SUM, tree->comm);
        } else {
            std::copy(part.begin(), part.begin() + counts[0], landed.begin());
        }
        // ---- on to the owners (their own: copied)
        int64_t at = 0;
        for (size_t q = 0; q < order.size(); ++q) {
            const int64_t n = order[q];
            if (lands[static_cast<size_t>(n)] != srank) continue;
            const int64_t len = nr[q] * nc[q];
            const int to = owner[static_cast<size_t>(n)];
            if (to == crank) {
                std::copy_n(landed.data() + at, len, alldat_loc + base[static_cast<size_t>(n)]);
            } else {
                partial_put(send[static_cast<size_t>(to)], &n, 1);
                partial_put(send[static_cast<size_t>(to)], landed.data() + at, static_cast<size_t>(len));
            }
            at += len;
        }
    }
    const auto recv = partial_exchange(caller_comm, send, 630);
    for (const auto& buf : recv) {
        for (const char* p = buf.data(); p < buf.data() + buf.size();) {
            int64_t n = -1;
            p = partial_get(p, &n, 1);
            if (n < 0 || n >= ninter || base[static_cast<size_t>(n)] < 0) {
                throw std::logic_error("extract_partial_blocks: a block for another rank");
            }
            p = partial_get(p, alldat_loc + base[static_cast<size_t>(n)], static_cast<size_t>(int64_t{rowidx[n]} * colidx[n]));
        }
    }
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
    // The record and root terms on the device (h2_gpu/partial_device.hpp:
    // the records and the root block resident there), a pass at most
    // max_columns columns; with BPACK_CHECK=solve the host product too, and
    // rank 0 prints their difference.
    bool device = false;
#ifdef H2_HAVE_GPU
    if constexpr (fmm::gpu::gpu_data_type<DataType>) {
        device = fmm::gpu::partial_device_ready(&tree, 0);
        if (device) {
            const int64_t cmax = fmm::gpu::partial_device_store().max_columns;
            if (ncol > cmax) {
                for (int64_t c0 = 0; c0 < ncol; c0 += cmax)
                    partial_schur_multiply(tree, x2 + c0 * n22, y2 + c0 * n22, std::min(cmax, ncol - c0));
                return;
            }
        }
    }
#endif
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
    if (!device) partial_records_to_host(tree, false);  // (the host reads the records)
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
        if (device) {
#ifdef H2_HAVE_GPU
            if constexpr (fmm::gpu::gpu_data_type<DataType>) fmm::gpu::partial_device_root(v.data(), ncol, g.data());
#endif
        } else {
            constexpr int64_t chunk = 256;
            const int64_t nchunks = (m + chunk - 1) / chunk;
#pragma omp parallel for schedule(static)
            for (int64_t q = 0; q < nchunks; ++q) {
                const int64_t a0 = q * chunk, rows = std::min(chunk, m - a0);
                partial_gemm<DataType>('N', 'N', rows, ncol, m, DataType(1), pr.cb.data.data() + a0, pr.cb.lda, v.data(), m,
                                       DataType(0), g.data() + a0, m);
            }
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
#ifdef H2_HAVE_GPU
        if constexpr (fmm::gpu::gpu_data_type<DataType>) {
            if (device) {
                // the sparsified boxes' p and w here, their products with the
                // records on the device: per box gR, gS and t_q = X_NR[q] p
                const auto& plan = fmm::gpu::partial_device_store().levels[static_cast<size_t>(l)];
                std::vector<DataType> outv;
                if (active) {
                    std::vector<DataType> pw(static_cast<size_t>(plan.pw_units * ncol));
#pragma omp parallel for schedule(dynamic, 1)
                    for (int64_t q = 0; q < static_cast<int64_t>(plan.boxes.size()); ++q) {
                        const auto& B = plan.boxes[static_cast<size_t>(q)];
                        const auto& box = level.local_boxes[static_cast<size_t>(B.bi)];
                        const auto& p = pv[static_cast<size_t>(l)][static_cast<size_t>(B.bi)];
                        const auto w = w_of(box, raw[static_cast<size_t>(l)][static_cast<size_t>(B.bi)].left_side, &p);
                        std::copy(p.begin(), p.end(), pw.begin() + B.pw_off * ncol);
                        std::copy(w.begin(), w.end(), pw.begin() + (B.pw_off + B.r) * ncol);
                    }
                    std::vector<const DataType*> xl(static_cast<size_t>(level.num_boxes_local));
                    for (int64_t bi = 0; bi < level.num_boxes_local; ++bi)
                        xl[static_cast<size_t>(bi)] = raw[static_cast<size_t>(l)][static_cast<size_t>(bi)].left_side.data();
                    std::unordered_map<int64_t, std::pair<int64_t, const DataType*>> xr;
                    for (const auto& kv : remote_raw) xr.emplace(kv.first, std::make_pair(kv.second.first, kv.second.second.data()));
                    outv.resize(static_cast<size_t>(plan.out_units * ncol));
                    fmm::gpu::partial_device_terms(&tree, l, ncol, xl, xr, pw.data(), outv.data());
                }
                auto t_of = [&](int q, size_t slot) {
                    const auto& B = plan.boxes[static_cast<size_t>(q)];
                    return outv.data() + (B.out_off + B.t_off[slot]) * ncol;
                };
                // (2) the products of other ranks' neighbours to their owners
                std::vector<std::vector<char>> send_nr(static_cast<size_t>(nranks));
                if (active) {
                    for (size_t q = 0; q < plan.boxes.size(); ++q) {
                        const auto& nb = level.local_boxes[static_cast<size_t>(plan.boxes[q].bi)];
                        const auto& fn = nb.f22;
                        for (size_t idx = 0; idx < nb.one_hop.size(); ++idx) {
                            const int64_t m = fn.neighbor_nr_counts[idx];
                            const int64_t target = nb.one_hop[idx];
                            if (m == 0 || level.is_box_on_process(target)) continue;
                            auto& buf = send_nr[static_cast<size_t>(level.owner_rank_of(target))];
                            const int64_t head[2] = {target, m};
                            partial_put(buf, head, 2);
                            partial_put(buf, fn.neighbor_nr_slots[idx].data(), static_cast<size_t>(m));
                            partial_put(buf, t_of(static_cast<int>(q), idx), static_cast<size_t>(m * ncol));
                        }
                    }
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
                // (3) every local box: gS with the neighbours' products added (local
                // ones in one-hop order, then the received), into the accumulators
                std::exception_ptr failure;
#pragma omp parallel for schedule(dynamic, 1)
                for (int64_t bi = 0; bi < level.num_boxes_local; ++bi) {
                    try {
                        const auto& box = level.local_boxes[static_cast<size_t>(bi)];
                        if (box.skeleton_indices.empty()) continue;
                        const int qi = plan.box_of[static_cast<size_t>(bi)];
                        const int64_t n = box.num_points, kp = static_cast<int64_t>(box.skeleton_indices.size());
                        const int64_t r = qi >= 0 ? plan.boxes[static_cast<size_t>(qi)].r : 0;
                        std::vector<DataType> gS(static_cast<size_t>(kp * ncol), DataType(0));
                        const DataType* gR = nullptr;
                        if (qi >= 0) {
                            const auto& B = plan.boxes[static_cast<size_t>(qi)];
                            gR = outv.data() + B.out_off * ncol;
                            std::copy_n(gR + r * ncol, kp * ncol, gS.begin());
                        }
                        std::vector<int32_t> s_of(static_cast<size_t>(n), -1);
                        for (int64_t i = 0; i < kp; ++i) s_of[static_cast<size_t>(box.skeleton_indices[static_cast<size_t>(i)])] = static_cast<int32_t>(i);
                        auto add_rows = [&](const std::vector<int64_t>& slots, const DataType* t, int64_t m) {
                            for (int64_t i = 0; i < m; ++i) {
                                const int32_t si = s_of[static_cast<size_t>(slots[static_cast<size_t>(i)])];
                                if (si < 0) continue;
                                for (int64_t c = 0; c < ncol; ++c) gS[static_cast<size_t>(si + c * kp)] += t[static_cast<size_t>(i + c * m)];
                            }
                        };
                        for (size_t idx = 0; idx < box.one_hop.size(); ++idx) {
                            const int64_t nm = box.one_hop[idx];
                            if (!is_local(nm)) continue;
                            const int64_t nl = nm - level.local_morton_start;
                            const auto& nb = level.local_boxes[static_cast<size_t>(nl)];
                            if (!own(nb)) continue;
                            const auto it = std::find(nb.one_hop.begin(), nb.one_hop.end(), box.morton_index);
                            if (it == nb.one_hop.end()) continue;
                            const size_t idn = static_cast<size_t>(it - nb.one_hop.begin());
                            const int64_t m = nb.f22.neighbor_nr_counts[idn];
                            if (m == 0) continue;
                            add_rows(nb.f22.neighbor_nr_slots[idn], t_of(plan.box_of[static_cast<size_t>(nl)], idn), m);
                        }
                        if (const auto it = pushed.find(box.morton_index); it != pushed.end())
                            for (const auto& [slots, tv] : it->second) add_rows(slots, tv.data(), static_cast<int64_t>(slots.size()));
                        auto& out = h[static_cast<size_t>(l)][static_cast<size_t>(bi)].left_side;
                        const auto& f = box.f22;
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
                continue;
            }
        }
#endif
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
                std::vector<DataType> scratch;  // (a compressed record's entries)
                for (size_t idx = 0; idx < nb.one_hop.size(); ++idx) {
                    const int64_t m = fn.neighbor_nr_counts[idx];
                    const int64_t target = nb.one_hop[idx];
                    if (m > 0 && !level.is_box_on_process(target)) {
                        std::vector<DataType> t(static_cast<size_t>(m * ncol));
                        partial_gemm<DataType>('N', 'N', m, ncol, rn, DataType(1), fmm::record_values(fn.X_NR[idx], scratch), m,
                                               pn.data(), rn, DataType(0), t.data(), m);
                        std::vector<char> msg;
                        const int64_t head[2] = {target, m};
                        partial_put(msg, head, 2);
                        partial_put(msg, fn.neighbor_nr_slots[idx].data(), static_cast<size_t>(m));
                        partial_put(msg, t.data(), t.size());
                        parts[static_cast<size_t>(ni)].emplace_back(level.owner_rank_of(target), std::move(msg));
                    }
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
                std::vector<DataType> scratch;  // (a compressed record's entries)
                if (mine) {
                    const auto& p = pv[static_cast<size_t>(l)][static_cast<size_t>(bi)];
                    const auto w = w_of(box, x, &p);
                    partial_gemm<DataType>('N', 'N', r, ncol, r, DataType(1), fmm::record_values(f.X_RR, scratch), r, p.data(), r,
                                           DataType(1), gR.data(), r);
                    partial_gemm<DataType>('N', 'N', r, ncol, kp, DataType(1), fmm::record_values(f.X_RS, scratch), r, w.data(), kp,
                                           DataType(1), gR.data(), r);
                    partial_gemm<DataType>('N', 'N', kp, ncol, r, DataType(1), fmm::record_values(f.X_SR, scratch), kp, p.data(), r,
                                           DataType(1), gS.data(), kp);
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
                        partial_gemm<DataType>('N', 'N', r, ncol, total, DataType(1), fmm::record_values(f.X_RN, scratch), r, rn.data(),
                                               total, DataType(1), gR.data(), r);
                    }
                }
                // the sparsified neighbours' p against their X_NR rows over this box's S' (local: here; remote: received)
                std::vector<int32_t> s_of(static_cast<size_t>(n), -1);
                for (int64_t i = 0; i < kp; ++i) s_of[static_cast<size_t>(box.skeleton_indices[static_cast<size_t>(i)])] = static_cast<int32_t>(i);
                auto add_rows = [&](const std::vector<int64_t>& slots, const DataType* t, int64_t m) {
                    for (int64_t i = 0; i < m; ++i) {
                        const int32_t si = s_of[static_cast<size_t>(slots[static_cast<size_t>(i)])];
                        if (si < 0) continue;  // (X_NR keeps this box's live rows only: its R rows are in X_RN above)
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
                    const int64_t m = fn.neighbor_nr_counts[idn], rn_ = static_cast<int64_t>(fn.redundant.size());
                    if (m == 0) continue;
                    const auto& pn = pv[static_cast<size_t>(l)][static_cast<size_t>(nl)];
                    t.resize(static_cast<size_t>(m * ncol));
                    partial_gemm<DataType>('N', 'N', m, ncol, rn_, DataType(1), fmm::record_values(fn.X_NR[idn], scratch), m, pn.data(), rn_,
                                           DataType(0), t.data(), m);
                    add_rows(fn.neighbor_nr_slots[idn], t.data(), m);
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
#ifdef H2_HAVE_GPU
    if constexpr (fmm::gpu::gpu_data_type<DataType>) {
        if (device && fmm::gpu::device_solve_check()) {
            std::vector<DataType> yh(static_cast<size_t>(n22 * ncol));
            fmm::gpu::device_solve_suspended() = true;
            partial_schur_multiply(tree, x2, yh.data(), ncol);
            fmm::gpu::device_solve_suspended() = false;
            double sums[2] = {0.0, 0.0};
            for (size_t i = 0; i < yh.size(); ++i) {
                sums[0] += std::norm(y2[i] - yh[i]);
                sums[1] += std::norm(yh[i]);
            }
            if (me == 0) {
                std::printf("GPU partial CB product check: |y_gpu - y_host| / |y_host| = %.3e\n",
                            sums[1] > 0.0 ? std::sqrt(sums[0] / sums[1]) : 0.0);
                std::fflush(stdout);
            }
        }
    }
#endif
}

// The release of the Schur complement (c_bpack_partial_schur_release): the
// root's CB block and every sparsified box's records and neighbour lists,
// which only the CB entries and products read, go.  What the partial solves
// read stays: the F11 factors, the boxes' T, the root's LU, W and A21.
// Returns the bytes this rank freed (BPACK_TRACE=memory prints them).
template<typename CoordType, typename DataType>
size_t partial_schur_release(fmm::ParallelTree<CoordType, DataType>& tree) {
    auto& pr = tree.partial_root;
    const size_t rss0 = fmm::env::trace("memory") ? h2_diag_current_rss_bytes() : 0;
    size_t records = 0, lists = 0;
    const size_t root = pr.cb.data.capacity() * sizeof(DataType) +
                        (pr.leaf_box.capacity() + pr.leaf_slot.capacity()) * sizeof(int64_t);
    std::vector<DataType>().swap(pr.cb.data);
    pr.cb.rows = pr.cb.cols = pr.cb.lda = 0;
    std::vector<int64_t>().swap(pr.leaf_box);
    std::vector<int64_t>().swap(pr.leaf_slot);
    pr.cb_released = true;
    for (auto& level : tree.levels) {
        for (auto& box : level.local_boxes) {
            auto& f = box.f22;
            if (!f.sparsified) continue;
            for (auto* m : {&f.X_RR, &f.X_RS, &f.X_SR, &f.X_RN}) {
                records += m->bytes();
                *m = fmm::PartialRecord<DataType>{};
            }
            for (const auto& m : f.X_NR) records += m.bytes();
            decltype(f.X_NR)().swap(f.X_NR);
            for (auto* v : {&f.neighbor_f22_counts, &f.neighbor_nr_counts}) {
                lists += v->capacity() * sizeof(int64_t);
                std::vector<int64_t>().swap(*v);
            }
            for (auto* vv : {&f.neighbor_f22_slots, &f.neighbor_nr_slots}) {
                for (const auto& v : *vv) lists += v.capacity() * sizeof(int64_t);
                std::vector<std::vector<int64_t>>().swap(*vv);
            }
            f.released = true;
        }
    }
#ifdef H2_HAVE_GPU
    if constexpr (fmm::gpu::gpu_data_type<DataType>) fmm::gpu::partial_device_release(&tree);  // (the kept records too)
#endif
    tree.partial_records_on_device = false;
    const size_t freed = records + root + lists;
    if (fmm::env::trace("memory")) {
        std::printf("H2_PARTIAL_RELEASE rank=%d freed=%.6f records=%.6f root=%.6f lists=%.6f rss_before=%.6f rss_after=%.6f\n",
                    tree.mpi_rank, freed / 1e9, records / 1e9, root / 1e9, lists / 1e9, rss0 / 1e9,
                    h2_diag_current_rss_bytes() / 1e9);
        std::fflush(stdout);
    }
    return freed;
}

}  // namespace color_unstructured
}  // namespace butterfly
