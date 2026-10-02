// ===========================================================================
// Multiply (A x = W D V x) on the component-owner engine — design §11.3
// ===========================================================================
//
// The verification multiply applied the factorization with the replicated
// sweeps, which fetch the FACTORS of every halo box before each colour: on
// an engine level those live only on the component owners, so at 64 ranks
// the multiply spent 54 s of its 58 s gathering factors (ledger 29/32).
// Here the boundary sweeps run once, on the owners, exactly as the solve
// does (owner_solve.hpp); only vectors cross the wire.
//
// The two sweeps have opposite data flow, so each borrows the machinery of
// the solve's OTHER sweep:
//   * W (forward, leaf -> root, blue -> purple -> green -> interior) is a
//     PULL: a box reads its neighbours' current values and writes only its
//     own entries (apply_mul_forward_W).  So, per colour, the owners step
//     their boxes, then send each stepped box's entries (REPLACE) to the
//     owners of adjacent higher-colour components — the solve's backward
//     structure, ascending.  Neighbours not yet stepped are read at their
//     level-start value, which the owner's ghost entry holds.
//   * V (backward, root -> leaf, interior -> green -> purple -> blue) is a
//     PUSH: a box adds into its neighbours' entries (apply_mul_backward_V).
//     Contributions go as ADD records to the target's holder and are applied
//     in canonical order — the solve's forward structure, descending.  The
//     canonical key is the NEGATED elimination sequence, so the ascending
//     sort of owner_solve_apply replays reverse elimination order.
// The interior is stepped on the homes, replicated, between the two owner
// sweeps of a level as in the solve; the owners write their non-local boxes'
// final values back to the homes at the end of each sweep.
#pragma once

#include "owner_solve.hpp"
#include "apply_mul.hpp"

namespace fmm {

/// V step of one owned box on its owner: own skeleton and redundant updates
/// in place, the neighbour pushes as records (same arithmetic and order as
/// apply_mul_backward_V_with_pending, whose pending map they replace).
template<typename CoordType, typename DataType>
void owner_mul_backward_step(TreeLevel<CoordType, DataType>& level,
                             std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
                             SolveDataRequest<CoordType, DataType>& e, const BoxData<CoordType, DataType>& box,
                             int32_t seq, std::vector<OwnerSolveRec<DataType>>& out) {
    using owner_solve_detail::gemm;
    const auto& S = box.skeleton_indices;
    const auto& R = box.redundant_indices;
    if (S.empty() || R.empty()) return;
    const int64_t k = static_cast<int64_t>(S.size()), r = static_cast<int64_t>(R.size());
    const int64_t nrhs = e.nrhs;
    if (nrhs <= 0 || e.num_points * nrhs != static_cast<int64_t>(e.left_side.size()))
        throw std::runtime_error("owner_mul_backward_step: invalid batched vector layout");
    std::vector<DataType> x_R(static_cast<size_t>(r * nrhs));
    for (int64_t column = 0; column < nrhs; ++column)
        for (int64_t i = 0; i < r; ++i)
            x_R[static_cast<size_t>(i + column * r)] =
                e.left_side[static_cast<size_t>(R[static_cast<size_t>(i)] +
                                                column * e.num_points)];
    // x_S -= X_SR x_R
    if (box.X_SR.is_allocated()) {
        std::vector<DataType> t(static_cast<size_t>(k * nrhs), DataType{0});
        gemm<DataType>('N', static_cast<int>(k), static_cast<int>(nrhs),
                       static_cast<int>(r), box.X_SR.data.data(),
                       static_cast<int>(box.X_SR.lda), x_R.data(), static_cast<int>(r),
                       DataType{0}, t.data(), static_cast<int>(k));
        for (int64_t column = 0; column < nrhs; ++column)
            for (int64_t i = 0; i < k; ++i)
                e.left_side[static_cast<size_t>(S[static_cast<size_t>(i)] +
                                                column * e.num_points)] -=
                    t[static_cast<size_t>(i + column * k)];
    }
    // neighbours: x_N -= X_NR x_R, one record per 1-hop slot
    if (box.X_NR.is_allocated() && !box.one_hop.empty()) {
        if (box.use_full_set.size() != box.one_hop.size())
            throw std::runtime_error("owner_mul_backward_step: use_full_set missing for box " + std::to_string(box.morton_index));
        const int64_t N = box.X_NR.rows;
        std::vector<DataType> u(static_cast<size_t>(N * nrhs), DataType{0});
        gemm<DataType>('N', static_cast<int>(N), static_cast<int>(nrhs),
                       static_cast<int>(r), box.X_NR.data.data(),
                       static_cast<int>(box.X_NR.lda), x_R.data(), static_cast<int>(r),
                       DataType{0}, u.data(), static_cast<int>(N));
        int64_t off = 0;
        for (size_t i = 0; i < box.one_hop.size(); ++i) {
            const int64_t nb = box.one_hop[i];
            const bool use_full = (box.use_full_set[i] == 1);
            const SolveDataRequest<CoordType, DataType>* ne = resolve_any_solve_data_for_morton(level, level_solve_data, nb);
            if (ne == nullptr)
                throw std::runtime_error("owner_mul_backward_step: neighbour " + std::to_string(nb) + " of box " +
                                         std::to_string(box.morton_index) + " has no solve entry here");
            if (ne->nrhs != nrhs)
                throw std::runtime_error("owner_mul_backward_step: neighbour RHS count mismatch");
            const int64_t n_i = use_full
                ? ne->num_points
                : static_cast<int64_t>(ne->skeleton_indices.size());
            if (off + n_i > N) throw std::runtime_error("owner_mul_backward_step: X_NR rows short for box " + std::to_string(box.morton_index));
            OwnerSolveRec<DataType> rec;
            rec.target = nb;
            rec.source = box.morton_index;
            rec.seq = seq;
            rec.kind = use_full ? OSR_ADD_FULL : OSR_ADD_SKEL;
            rec.v.resize(static_cast<size_t>(n_i * nrhs));
            for (int64_t column = 0; column < nrhs; ++column)
                for (int64_t j = 0; j < n_i; ++j)
                    rec.v[static_cast<size_t>(j + column * n_i)] =
                        -u[static_cast<size_t>(off + j + column * N)];
            out.push_back(std::move(rec));
            off += n_i;
        }
        if (off != N)
            throw std::runtime_error("owner_mul_backward_step: X_NR rows (" + std::to_string(N) + ") != slot sum (" +
                                     std::to_string(off) + ") for box " + std::to_string(box.morton_index));
    }
    // x_R += T^T x_S (with the updated x_S)
    if (box.interpolation_matrix.is_allocated()) {
        std::vector<DataType> x_S(static_cast<size_t>(k * nrhs));
        for (int64_t column = 0; column < nrhs; ++column)
            for (int64_t i = 0; i < k; ++i)
                x_S[static_cast<size_t>(i + column * k)] =
                    e.left_side[static_cast<size_t>(S[static_cast<size_t>(i)] +
                                                    column * e.num_points)];
        std::vector<DataType> t(static_cast<size_t>(r * nrhs), DataType{0});
        gemm<DataType>('T', static_cast<int>(r), static_cast<int>(nrhs),
                       static_cast<int>(k), box.interpolation_matrix.data.data(),
                       static_cast<int>(box.interpolation_matrix.lda), x_S.data(),
                       static_cast<int>(k), DataType{0}, t.data(), static_cast<int>(r));
        for (int64_t column = 0; column < nrhs; ++column)
            for (int64_t i = 0; i < r; ++i)
                e.left_side[static_cast<size_t>(R[static_cast<size_t>(i)] +
                                                column * e.num_points)] +=
                    t[static_cast<size_t>(i + column * r)];
    }
}

namespace owner_mul_detail {

/// REPLACE records of every owned box of colour `ci` to the owners of its
/// adjacent components of higher (`higher` = true) or lower colour.
template<typename CoordType, typename DataType>
void queue_final_values(TreeLevel<CoordType, DataType>& level,
                        std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data, const OwnerSolveRecord& rec,
                        int ci, bool higher, std::map<int, std::vector<char>>& out) {
    if (!rec.active) return;
    for (int32_t id : rec.owned) {
        if (rec.color_index(id) != ci) continue;
        for (int64_t m : rec.graph.comps[static_cast<size_t>(id)].boxes) {
            SolveDataRequest<CoordType, DataType>* e = nullptr;
            BoxData<CoordType, DataType>* box = nullptr;
            owner_solve_detail::owned_box_and_entry(level, level_solve_data, m, e, box);
            std::set<int> dests;
            for (int64_t nb : box->one_hop) {
                const int32_t cn = rec.comp_of(nb);
                if (cn < 0) continue;
                const int cc = rec.color_index(cn);
                if (higher ? cc > ci : cc < ci)
                    dests.insert(rec.rank_of_pid(static_cast<uint32_t>(rec.graph.comps[static_cast<size_t>(cn)].owner)));
            }
            for (int d : dests) {
                if (d == rec.my_rank) continue;
                OwnerSolveRec<DataType> r;
                r.target = m; r.source = m; r.seq = 0; r.kind = OSR_REPLACE; r.v = e->left_side;
                owner_solve_put(out[d], r);
            }
        }
    }
}

/// REPLACE records of every owned box that is not local, to its home.
template<typename CoordType, typename DataType>
void queue_write_back(TreeLevel<CoordType, DataType>& level,
                      std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data, const OwnerSolveRecord& rec,
                      std::map<int, std::vector<char>>& out) {
    if (!rec.active) return;
    for (int32_t id : rec.owned)
        for (int64_t m : rec.graph.comps[static_cast<size_t>(id)].boxes) {
            if (level.is_box_on_process(m)) continue;
            SolveDataRequest<CoordType, DataType>* e = nullptr;
            BoxData<CoordType, DataType>* box = nullptr;
            owner_solve_detail::owned_box_and_entry(level, level_solve_data, m, e, box);
            OwnerSolveRec<DataType> r;
            r.target = m; r.source = m; r.seq = 0; r.kind = OSR_REPLACE; r.v = e->left_side;
            owner_solve_put(out[rec.home_rank(m)], r);
        }
}

}  // namespace owner_mul_detail

/// W sweep of the boundary of one engine level (collective: every rank calls it).
/// Precondition: the level-start vectors gather has run, so every owner's
/// entries (own boxes and their neighbours) hold the level's input values.
/// Postcondition: the homes hold the stepped values of every boundary box.
template<typename CoordType, typename DataType>
void owner_mul_forward_boundary(TreeLevel<CoordType, DataType>& level,
                                std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
                                const OwnerSolveRecord& rec, OwnerSolveStats& stats,
                                MPI_Comm comm) {
    using namespace owner_solve_detail;
    using clock = std::chrono::high_resolution_clock;
    const auto groups = owned_by_color_wave(rec);
    std::vector<OwnerSolveRec<DataType>> inbox;
    P2PSends sends;
    int round = 0;
    for (int ci = 0; ci < 3; ++ci) {   // blue, purple, green
        for (int w = 0; w < rec.num_waves; ++w) {
            const auto& wave = groups[static_cast<size_t>(ci)][static_cast<size_t>(w)];
            if (wave.empty()) continue;
            const auto t0 = clock::now();
            std::exception_ptr ex;
            std::mutex ex_mutex;
            std::atomic<bool> failed{false};
            #pragma omp parallel for schedule(dynamic)
            for (int64_t i = 0; i < static_cast<int64_t>(wave.size()); ++i) {
                if (failed.load(std::memory_order_relaxed)) continue;
                try {
                    const int64_t m = wave[static_cast<size_t>(i)];
                    SolveDataRequest<CoordType, DataType>* e = nullptr;
                    BoxData<CoordType, DataType>* box = nullptr;
                    owned_box_and_entry(level, level_solve_data, m, e, box);
                    // the kernel takes the factors from the ghost BoxData when the
                    // entry carries none (vectors-only gather), as on the owner
                    apply_mul_forward_W(level, *e, level_solve_data, MatrixProperty::SYMMETRIC, !level.is_box_on_process(m));
                } catch (...) {
                    if (!failed.exchange(true)) { std::lock_guard<std::mutex> l(ex_mutex); ex = std::current_exception(); }
                }
            }
            if (ex) std::rethrow_exception(ex);
            stats.boxes += static_cast<int64_t>(wave.size());
            stats.t_step_ms += std::chrono::duration<double, std::milli>(clock::now() - t0).count();
        }
        if (ci == 2) break;   // green's readers are the interior, served by the write-back
        std::map<int, std::vector<char>> out;
        owner_mul_detail::queue_final_values(level, level_solve_data, rec, ci, /*higher=*/true, out);
        std::set<int> dests, senders;
        round_peers(rec, ci, RoundKind::HANDOFF_UP, dests, senders);
        const auto t2 = clock::now();
        p2p_send(sends, out, dests, TAG_MUL_FWD + round, stats, comm);
        p2p_receive<DataType>(senders, TAG_MUL_FWD + round, inbox, stats, comm);
        ++round;
        stats.t_round_ms += std::chrono::duration<double, std::milli>(clock::now() - t2).count();
        deliver(level, level_solve_data, static_cast<OwnerSlotInbox<DataType>*>(nullptr), inbox, stats);
    }
    {
        std::map<int, std::vector<char>> out;
        owner_mul_detail::queue_write_back(level, level_solve_data, rec, out);
        std::set<int> dests, senders;
        write_back_peers(level, rec, dests, senders);
        const auto t2 = clock::now();
        p2p_send(sends, out, dests, TAG_MUL_FWD + round, stats, comm);
        p2p_receive<DataType>(senders, TAG_MUL_FWD + round, inbox, stats, comm);
        stats.t_round_ms += std::chrono::duration<double, std::milli>(clock::now() - t2).count();
        deliver(level, level_solve_data, static_cast<OwnerSlotInbox<DataType>*>(nullptr), inbox, stats);
    }
    sends.finish();
}

/// V sweep of the boundary of one engine level (collective).  Precondition:
/// the interior has been stepped on the homes and the after-interior vectors
/// gather has run, so the owners' entries carry the interior's pushes.
/// Postcondition: the homes hold the final values of every boundary box, and
/// the interior entries carry the boundary's pushes.
template<typename CoordType, typename DataType>
void owner_mul_backward_boundary(TreeLevel<CoordType, DataType>& level,
                                 std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
                                 const OwnerSolveRecord& rec, OwnerSolveStats& stats,
                                 MPI_Comm comm) {
    using namespace owner_solve_detail;
    using clock = std::chrono::high_resolution_clock;
    const auto groups = owned_by_color_wave(rec);
    // fold order: decreasing sequence — the V sweep runs in reverse elimination order
    OwnerSlotInbox<DataType> slots = build_slot_inbox(level, level_solve_data, rec, /*descending=*/true);
    const int nthreads = std::max(1, omp_get_max_threads());
    std::vector<OwnerSolveRec<DataType>> inbox;
    P2PSends sends;
    int round = 0;
    for (int ci = 2; ci >= 0; --ci) {   // green, purple, blue
        std::map<int, std::vector<char>> out;
        for (int w = rec.num_waves - 1; w >= 0; --w) {
            const auto& wave = groups[static_cast<size_t>(ci)][static_cast<size_t>(w)];
            if (wave.empty()) continue;
            const auto t0 = clock::now();
            std::vector<std::vector<OwnerSolveRec<DataType>>> per_thread(static_cast<size_t>(nthreads));
            std::exception_ptr ex;
            std::mutex ex_mutex;
            std::atomic<bool> failed{false};
            #pragma omp parallel for schedule(dynamic)
            for (int64_t i = 0; i < static_cast<int64_t>(wave.size()); ++i) {
                if (failed.load(std::memory_order_relaxed)) continue;
                try {
                    const int64_t m = wave[static_cast<size_t>(i)];
                    SolveDataRequest<CoordType, DataType>* e = nullptr;
                    BoxData<CoordType, DataType>* box = nullptr;
                    owned_box_and_entry(level, level_solve_data, m, e, box);
                    // everything eliminated after this box (earlier in V) has contributed
                    auto* t = slots.find(m);
                    if (t == nullptr) throw std::runtime_error("owner_mul_backward_boundary: no slots for owned box " + std::to_string(m));
                    const int32_t seq = rec.seq_of(m);
                    slots.template fold<CoordType>(*t, *e, [seq](int32_t s) { return s > seq; });
                    owner_mul_backward_step(level, level_solve_data, *e, *box, seq, per_thread[static_cast<size_t>(omp_get_thread_num())]);
                } catch (...) {
                    if (!failed.exchange(true)) { std::lock_guard<std::mutex> l(ex_mutex); ex = std::current_exception(); }
                }
            }
            if (ex) std::rethrow_exception(ex);
            stats.boxes += static_cast<int64_t>(wave.size());
            for (auto& v : per_thread) route(rec, slots, v, out, stats);
            stats.t_step_ms += std::chrono::duration<double, std::milli>(clock::now() - t0).count();
        }
        std::set<int> dests, senders;
        round_peers(rec, ci, RoundKind::PUSH_ALL, dests, senders);
        const auto t2 = clock::now();
        p2p_send(sends, out, dests, TAG_MUL_BWD + round, stats, comm);
        p2p_receive<DataType>(senders, TAG_MUL_BWD + round, inbox, stats, comm);
        ++round;
        stats.t_round_ms += std::chrono::duration<double, std::milli>(clock::now() - t2).count();
        const auto t3 = clock::now();
        deliver(level, level_solve_data, &slots, inbox, stats);
        stats.t_apply_ms += std::chrono::duration<double, std::milli>(clock::now() - t3).count();
    }
    {
        const auto t3 = clock::now();
        fold_all(level, level_solve_data, slots);   // later (lower-colour) sources onto stepped boxes, and the interior's inputs
        stats.t_apply_ms += std::chrono::duration<double, std::milli>(clock::now() - t3).count();
    }
    {
        std::map<int, std::vector<char>> out;
        owner_mul_detail::queue_write_back(level, level_solve_data, rec, out);
        std::set<int> dests, senders;
        write_back_peers(level, rec, dests, senders);
        const auto t2 = clock::now();
        p2p_send(sends, out, dests, TAG_MUL_BWD + round, stats, comm);
        p2p_receive<DataType>(senders, TAG_MUL_BWD + round, inbox, stats, comm);
        stats.t_round_ms += std::chrono::duration<double, std::milli>(clock::now() - t2).count();
        deliver(level, level_solve_data, static_cast<OwnerSlotInbox<DataType>*>(nullptr), inbox, stats);
    }
    sends.finish();
}

}  // namespace fmm
