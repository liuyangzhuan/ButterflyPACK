#pragma once
// ===========================================================================
// owner_solve.hpp — the solve on the component-owner engine
// ===========================================================================
// Design: concise-algorithm-ref/component-owner-dataflow-design.md §11.1/§11.2.
//
// A level factorized by the engine (H2_CA_owner_component=3, more than one active
// rank) is swept the same way it was eliminated: every boundary box's forward
// and backward step runs ONCE, on the owner of its colour component, and only
// vectors cross ranks.  Levels factorized replicated keep the replicated sweep
// in distributed_routine_all.cpp; the choice is per level, no flag.
//
// Forward (per level >= 2), after the level-start vector gather (which hands
// every owner the current x of its owned remote boxes and the sizes of all
// their neighbours):
//   blue step on owners  -> round: contributions to their holders
//   purple step          -> round
//   green step           -> round, then WRITE-BACK of owned remote boxes to
//                           their homes; the interior then runs on the homes
//                           with the replicated code (all its neighbours are
//                           local) and gather_skeleton_to_parent is unchanged.
// A contribution x[N] += X̃_NR·x[R] from source S to neighbour G is applied by
// G's HOLDER during the boundary phase: the owner of G's component if G is in
// one, else G's home.  Contributions inside a component (earlier wave -> later
// wave) are applied on the owner at the end of the wave; every other one is
// buffered and applied after the colour's round, all in canonical (sequence,
// source, target) order — so the sums do not depend on timing, thread count
// or rank count (fast == serial bitwise; the replicated sweep sums in thread
// order and agrees only up to rounding).
//
// Backward (per level >= 2), after the parent scatter, the diagonal solve on
// the homes, the interior step on the homes (it reads only skeleton entries,
// final since the scatter) and the after-interior vector gather (which hands
// every owner the final x of the interior boxes and the scattered/diagonal
// x of its owned boxes and their neighbours):
//   green step on owners  -> round: final x of green boxes to the owners of
//                            adjacent lower-colour components
//   purple step           -> round: to blue owners
//   blue step             -> WRITE-BACK of every owned remote box to its home
// The backward step REPLACES x[R] and subtracts from x[S], so it runs exactly
// once and its result is shipped; a box's step needs the final x of every
// neighbour eliminated after it (higher colour, later wave, interior) — the
// rounds above deliver exactly those.
//
// Ownership information outlives the factorization through OwnerSolveRecord
// (the OwnerScheduleState dies with its level).
// ===========================================================================

#include "owner_schedule.hpp"
#include "solver.hpp"
#include <map>

namespace fmm {

// ---------------------------------------------------------------------------
// Per-level record of the factorization schedule
// ---------------------------------------------------------------------------
struct OwnerSolveRecord {
    bool engine_level = false;      ///< the level was factorized by the engine (same on every rank)
    bool active = false;            ///< this rank took part (has the graph)
    int dimension = 3;
    int num_waves = 8;
    int num_boundary_colors = 3;
    dataflow::ProcessGrid grid;
    dataflow::ComponentGraph graph;
    uint32_t my_pid = 0;
    int my_rank = -1;
    std::unordered_map<int, int> pid_to_rank;
    std::vector<int32_t> comp_seq_color;
    std::unordered_map<int64_t, int32_t> box_wave;   ///< component-aware wave of every boundary box (as factorized)
    std::vector<int32_t> owned;     ///< components owned by this rank

    int32_t comp_of(int64_t m) const { return graph.comp_id_of_box(m); }
    int rank_of_pid(uint32_t pid) const {
        auto it = pid_to_rank.find(static_cast<int>(pid));
        if (it == pid_to_rank.end()) throw std::runtime_error("owner_solve: no rank for pid " + std::to_string(pid));
        return it->second;
    }
    int home_rank(int64_t m) const { return rank_of_pid(grid.proc_of_box(m)); }
    /// Rank holding the authoritative x of box m during the boundary phase.
    int holder_rank(int64_t m) const {
        const int32_t c = comp_of(m);
        return c < 0 ? home_rank(m) : rank_of_pid(static_cast<uint32_t>(graph.comps[static_cast<size_t>(c)].owner));
    }
    /// Elimination sequence of a boundary box: (colour index, wave), as in the factorization.
    int32_t seq_of(int64_t m) const {
        const int32_t c = comp_of(m);
        const int32_t ci = c >= 0 ? comp_seq_color[static_cast<size_t>(c)] : num_boundary_colors;
        return ci * num_waves + wave_of(m);
    }
    /// Wave of a boundary box, exactly as the factorization eliminated it.
    int32_t wave_of(int64_t m) const {
        auto it = box_wave.find(m);
        if (it == box_wave.end()) throw std::runtime_error("owner_solve: box " + std::to_string(m) + " has no wave");
        return it->second;
    }
    int color_index(int32_t comp) const { return static_cast<int>(graph.comps[static_cast<size_t>(comp)].color); }
};

inline std::map<int32_t, OwnerSolveRecord>& owner_solve_records() {
    static std::map<int32_t, OwnerSolveRecord> records;
    return records;
}

/// Wave grouping for the replicated sweeps (the multiply) on a level: the
/// factorization's own order must be replayed, so on an engine level a
/// boundary box takes its recorded component-aware wave, while interior
/// boxes (and every box of a non-engine level) keep the Morton parity
/// classes.  `parity_waves` is 4 in 2D, 8 in 3D.
inline int owner_level_num_waves(int32_t level_index, int parity_waves) {
    const auto& recs = owner_solve_records();
    auto it = recs.find(level_index);
    if (it == recs.end() || !it->second.engine_level || !it->second.active) return parity_waves;
    return std::max(parity_waves, it->second.num_waves);
}
inline int owner_level_wave_of(int32_t level_index, int64_t m, int parity_waves) {
    const auto& recs = owner_solve_records();
    auto it = recs.find(level_index);
    if (it != recs.end() && it->second.engine_level && it->second.active) {
        auto w = it->second.box_wave.find(m);
        if (w != it->second.box_wave.end()) return w->second;
    }
    return static_cast<int>(m % parity_waves);
}

/// Record the schedule of a level for the solve (called by the factorization
/// driver on every rank at the end of the level, active or not).
template<typename CoordType, typename DataType>
void owner_solve_record_level(const OwnerScheduleState<CoordType, DataType>& st, int32_t level_index, bool engine_level) {
    OwnerSolveRecord& rec = owner_solve_records()[level_index];
    rec = OwnerSolveRecord{};
    rec.engine_level = engine_level;
    rec.active = st.active;
    if (!st.active) return;
    rec.dimension = st.dimension;
    rec.num_waves = st.num_waves;
    rec.box_wave = st.box_wave;
    rec.num_boundary_colors = st.num_boundary_colors;
    rec.grid = st.grid;
    rec.graph = st.graph;
    rec.my_pid = st.my_pid;
    rec.my_rank = st.my_rank;
    rec.pid_to_rank = st.pid_to_rank;
    rec.comp_seq_color = st.comp_seq_color;
    rec.owned = st.rt.owned();
}

inline bool owner_solve_engine_level(int32_t level_index) {
    auto it = owner_solve_records().find(level_index);
    return it != owner_solve_records().end() && it->second.engine_level;
}

// ---------------------------------------------------------------------------
// Records exchanged between ranks
// ---------------------------------------------------------------------------
/// ADD_FULL: x_target[j] += v[j] (point order); ADD_SKEL: x_target[skel[j]] += v[j];
/// REPLACE: x_target = v (the box's authoritative vector moving between ranks).
enum : int8_t { OSR_ADD_FULL = 0, OSR_ADD_SKEL = 1, OSR_REPLACE = 2 };

template<typename DataType>
struct OwnerSolveRec {
    int64_t target = 0;
    int64_t source = 0;   ///< contributing box (== target for REPLACE)
    int32_t seq = 0;      ///< elimination sequence of the source (canonical order key)
    int8_t kind = OSR_ADD_FULL;
    std::vector<DataType> v;
};

template<typename DataType>
inline void owner_solve_put(std::vector<char>& buf, const OwnerSolveRec<DataType>& r) {
    const int64_t n = static_cast<int64_t>(r.v.size());
    const size_t off = buf.size();
    buf.resize(off + 3 * sizeof(int64_t) + sizeof(int32_t) + sizeof(int8_t) + static_cast<size_t>(n) * sizeof(DataType));
    char* p = buf.data() + off;
    std::memcpy(p, &r.target, sizeof(int64_t)); p += sizeof(int64_t);
    std::memcpy(p, &r.source, sizeof(int64_t)); p += sizeof(int64_t);
    std::memcpy(p, &r.seq, sizeof(int32_t)); p += sizeof(int32_t);
    std::memcpy(p, &r.kind, sizeof(int8_t)); p += sizeof(int8_t);
    std::memcpy(p, &n, sizeof(int64_t)); p += sizeof(int64_t);
    if (n) std::memcpy(p, r.v.data(), static_cast<size_t>(n) * sizeof(DataType));
}

template<typename DataType>
inline void owner_solve_parse(const char* p, const char* end, std::vector<OwnerSolveRec<DataType>>& out) {
    while (p < end) {
        OwnerSolveRec<DataType> r;
        int64_t n = 0;
        std::memcpy(&r.target, p, sizeof(int64_t)); p += sizeof(int64_t);
        std::memcpy(&r.source, p, sizeof(int64_t)); p += sizeof(int64_t);
        std::memcpy(&r.seq, p, sizeof(int32_t)); p += sizeof(int32_t);
        std::memcpy(&r.kind, p, sizeof(int8_t)); p += sizeof(int8_t);
        std::memcpy(&n, p, sizeof(int64_t)); p += sizeof(int64_t);
        r.v.resize(static_cast<size_t>(n));
        if (n) { std::memcpy(r.v.data(), p, static_cast<size_t>(n) * sizeof(DataType)); p += static_cast<size_t>(n) * sizeof(DataType); }
        out.push_back(std::move(r));
    }
    if (p != end) throw std::runtime_error("owner_solve_parse: trailing bytes");
}

struct OwnerSolveStats {
    double t_step_ms = 0, t_round_ms = 0, t_apply_ms = 0;
    size_t bytes_sent = 0, bytes_recv = 0;
    int rounds = 0;
    int64_t boxes = 0, records = 0;
};

/// One round: every rank sends its per-destination buffers and receives the
/// records addressed to it (collective over the active level communicator).
/// Records for this rank are parsed
/// directly.  The order of arrival is irrelevant: the caller sorts.
template<typename DataType>
void owner_solve_exchange(std::map<int, std::vector<char>>& out, std::vector<OwnerSolveRec<DataType>>& inbox,
                          OwnerSolveStats& stats, MPI_Comm comm) {
    static constexpr int TAG = 3401;
    int size = 1, rank = 0;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    std::vector<int64_t> send_sizes(static_cast<size_t>(size), 0), recv_sizes(static_cast<size_t>(size), 0);
    for (const auto& kv : out)
        if (kv.first != rank) send_sizes[static_cast<size_t>(kv.first)] = static_cast<int64_t>(kv.second.size());
    MPI_Alltoall(send_sizes.data(), 1, MPI_INT64_T, recv_sizes.data(), 1, MPI_INT64_T, comm);
    std::vector<std::vector<char>> rbuf(static_cast<size_t>(size));
    std::vector<MPI_Request> reqs;
    for (int s = 0; s < size; ++s) {
        const int64_t n = recv_sizes[static_cast<size_t>(s)];
        if (n <= 0) continue;
        if (n > std::numeric_limits<int>::max()) throw std::runtime_error("owner_solve_exchange: message too large");
        rbuf[static_cast<size_t>(s)] = OwnerBufferPool::instance().acquire(static_cast<size_t>(n));
        MPI_Request r;
        owner_mpi_check(MPI_Irecv(rbuf[static_cast<size_t>(s)].data(), static_cast<int>(n), MPI_CHAR, s, TAG, comm, &r),
                        "MPI_Irecv(owner_solve)");
        reqs.push_back(r);
        stats.bytes_recv += static_cast<size_t>(n);
    }
    for (auto& kv : out) {
        if (kv.first == rank || kv.second.empty()) continue;
        if (kv.second.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
            throw std::runtime_error("owner_solve_exchange: message too large");
        MPI_Request r;
        owner_mpi_check(MPI_Isend(kv.second.data(), static_cast<int>(kv.second.size()), MPI_CHAR, kv.first, TAG,
                                  comm, &r), "MPI_Isend(owner_solve)");
        reqs.push_back(r);
        stats.bytes_sent += kv.second.size();
    }
    if (!reqs.empty()) owner_mpi_check(MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE), "MPI_Waitall(owner_solve)");
    auto self = out.find(rank);
    if (self != out.end()) owner_solve_parse<DataType>(self->second.data(), self->second.data() + self->second.size(), inbox);
    for (int s = 0; s < size; ++s) {
        auto& b = rbuf[static_cast<size_t>(s)];
        if (b.empty()) continue;
        owner_solve_parse<DataType>(b.data(), b.data() + b.size(), inbox);
        OwnerBufferPool::instance().release(std::move(b));
    }
    out.clear();
    ++stats.rounds;
}

/// Apply records in canonical order (sequence, source, target, kind).
template<typename CoordType, typename DataType>
void owner_solve_apply(TreeLevel<CoordType, DataType>& level,
                       std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
                       std::vector<OwnerSolveRec<DataType>>& recs, OwnerSolveStats& stats) {
    std::sort(recs.begin(), recs.end(), [](const OwnerSolveRec<DataType>& a, const OwnerSolveRec<DataType>& b) {
        if (a.seq != b.seq) return a.seq < b.seq;
        if (a.source != b.source) return a.source < b.source;
        if (a.target != b.target) return a.target < b.target;
        return a.kind < b.kind;
    });
    for (const auto& r : recs) {
        SolveDataRequest<CoordType, DataType>* e = resolve_any_solve_data_for_morton(level, level_solve_data, r.target);
        if (e == nullptr)
            throw std::runtime_error("owner_solve_apply: no solve entry for box " + std::to_string(r.target) +
                                     " (record from " + std::to_string(r.source) + ")");
        const size_t n = r.v.size();
        if (r.kind == OSR_ADD_FULL || r.kind == OSR_REPLACE) {
            if (e->left_side.size() != n)
                throw std::runtime_error("owner_solve_apply: size mismatch for box " + std::to_string(r.target) + ": entry " +
                                         std::to_string(e->left_side.size()) + ", record " + std::to_string(n));
            if (r.kind == OSR_REPLACE) std::copy(r.v.begin(), r.v.end(), e->left_side.begin());
            else for (size_t j = 0; j < n; ++j) e->left_side[j] += r.v[j];
        } else {
            const size_t skeleton_values =
                e->skeleton_indices.size() * static_cast<size_t>(e->nrhs);
            if (skeleton_values != n)
                throw std::runtime_error("owner_solve_apply: skeleton size mismatch for box " + std::to_string(r.target) +
                                         ": entry " + std::to_string(skeleton_values) + ", record " + std::to_string(n));
            const int64_t k = static_cast<int64_t>(e->skeleton_indices.size());
            for (int64_t column = 0; column < e->nrhs; ++column)
                for (int64_t j = 0; j < k; ++j)
                    e->left_side[static_cast<size_t>(
                        e->skeleton_indices[static_cast<size_t>(j)] +
                        column * e->num_points)] +=
                        r.v[static_cast<size_t>(j + column * k)];
        }
    }
    stats.records += static_cast<int64_t>(recs.size());
    recs.clear();
}

// ---------------------------------------------------------------------------
// Owner-side steps (the replicated kernels' arithmetic, operand for operand)
// ---------------------------------------------------------------------------
namespace owner_solve_detail {

template<typename DataType>
inline void gemv(char trans, int m, int n, const DataType* A, int lda, const DataType* x, DataType beta, DataType* y) {
    int incx = 1, incy = 1;
    DataType alpha = 1.0;
    if constexpr (std::is_same_v<DataType, float>) {
        sgemv_(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
    } else if constexpr (std::is_same_v<DataType, double>) {
        dgemv_(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
    } else if constexpr (std::is_same_v<DataType, std::complex<float>>) {
        cgemv_(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
    } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
        zgemv_(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
    } else {
        static_assert(sizeof(DataType) == 0, "owner_solve: unsupported DataType");
    }
}

template<typename DataType>
inline void gemm(char trans_a, int m, int n, int k, const DataType* A,
                 int lda, const DataType* B, int ldb, DataType beta,
                 DataType* C, int ldc) {
    char trans_b = 'N';
    DataType alpha = DataType{1};
    if constexpr (std::is_same_v<DataType, float>) {
        sgemm_(&trans_a, &trans_b, &m, &n, &k, &alpha, A, &lda,
               B, &ldb, &beta, C, &ldc);
    } else if constexpr (std::is_same_v<DataType, double>) {
        dgemm_(&trans_a, &trans_b, &m, &n, &k, &alpha, A, &lda,
               B, &ldb, &beta, C, &ldc);
    } else if constexpr (std::is_same_v<DataType, std::complex<float>>) {
        cgemm_(&trans_a, &trans_b, &m, &n, &k, &alpha, A, &lda,
               B, &ldb, &beta, C, &ldc);
    } else if constexpr (std::is_same_v<DataType, std::complex<double>>) {
        zgemm_(&trans_a, &trans_b, &m, &n, &k, &alpha, A, &lda,
               B, &ldb, &beta, C, &ldc);
    } else {
        static_assert(sizeof(DataType) == 0,
                      "owner_solve: unsupported DataType");
    }
}

}  // namespace owner_solve_detail

/// Forward step of box `box` on its entry `e`; the contributions to its
/// neighbours come back as records (one per 1-hop slot, in X̃_NR's row order:
/// full point order or the neighbour's skeleton order per use_full_set).
template<typename CoordType, typename DataType>
void owner_solve_forward_step(TreeLevel<CoordType, DataType>& level,
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
        throw std::runtime_error("owner_solve_forward_step: invalid batched vector layout");
    std::vector<DataType> x_S(static_cast<size_t>(k * nrhs));
    std::vector<DataType> x_R(static_cast<size_t>(r * nrhs));
    for (int64_t column = 0; column < nrhs; ++column)
        for (int64_t i = 0; i < k; ++i)
            x_S[static_cast<size_t>(i + column * k)] =
                e.left_side[static_cast<size_t>(S[static_cast<size_t>(i)] +
                                                column * e.num_points)];
    // x[R] -= T^T x[S]
    if (box.interpolation_matrix.is_allocated()) {
        std::vector<DataType> t(static_cast<size_t>(r * nrhs), DataType{0});
        gemm<DataType>('T', static_cast<int>(r), static_cast<int>(nrhs),
                       static_cast<int>(k), box.interpolation_matrix.data.data(),
                       static_cast<int>(box.interpolation_matrix.lda), x_S.data(),
                       static_cast<int>(k), DataType{0}, t.data(), static_cast<int>(r));
        for (int64_t column = 0; column < nrhs; ++column)
            for (int64_t i = 0; i < r; ++i)
                e.left_side[static_cast<size_t>(R[static_cast<size_t>(i)] +
                                                column * e.num_points)] -=
                    t[static_cast<size_t>(i + column * r)];
    }
    for (int64_t column = 0; column < nrhs; ++column)
        for (int64_t i = 0; i < r; ++i)
            x_R[static_cast<size_t>(i + column * r)] =
                e.left_side[static_cast<size_t>(R[static_cast<size_t>(i)] +
                                                column * e.num_points)];
    // x[S] += X̃_SR x[R]   (stored X_SR = -X_SR X_RR^{-1})
    if (box.X_SR.is_allocated()) {
        std::vector<DataType> t(static_cast<size_t>(k * nrhs), DataType{0});
        gemm<DataType>('N', static_cast<int>(k), static_cast<int>(nrhs),
                       static_cast<int>(r), box.X_SR.data.data(),
                       static_cast<int>(box.X_SR.lda), x_R.data(),
                       static_cast<int>(r), DataType{0}, t.data(), static_cast<int>(k));
        for (int64_t column = 0; column < nrhs; ++column)
            for (int64_t i = 0; i < k; ++i)
                e.left_side[static_cast<size_t>(S[static_cast<size_t>(i)] +
                                                column * e.num_points)] +=
                    t[static_cast<size_t>(i + column * k)];
    }
    if (!box.X_NR.is_allocated() || box.one_hop.empty()) return;
    if (box.use_full_set.size() != box.one_hop.size())
        throw std::runtime_error("owner_solve_forward_step: use_full_set missing for box " + std::to_string(box.morton_index));
    // u = X̃_NR x[R]   (stored X_NR = -X_NR X_RR^{-1}), cut per neighbour slot
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
            throw std::runtime_error("owner_solve_forward_step: neighbour " + std::to_string(nb) + " of box " +
                                     std::to_string(box.morton_index) + " has no solve entry here");
        if (ne->nrhs != nrhs)
            throw std::runtime_error("owner_solve_forward_step: neighbour RHS count mismatch");
        const int64_t n_i = use_full
            ? ne->num_points
            : static_cast<int64_t>(ne->skeleton_indices.size());
        if (off + n_i > N)
            throw std::runtime_error("owner_solve_forward_step: X_NR rows short for box " + std::to_string(box.morton_index));
        OwnerSolveRec<DataType> rec;
        rec.target = nb;
        rec.source = box.morton_index;
        rec.seq = seq;
        rec.kind = use_full ? OSR_ADD_FULL : OSR_ADD_SKEL;
        rec.v.resize(static_cast<size_t>(n_i * nrhs));
        for (int64_t column = 0; column < nrhs; ++column)
            std::copy_n(u.begin() + off + column * N, n_i,
                        rec.v.begin() + column * n_i);
        out.push_back(std::move(rec));
        off += n_i;
    }
    if (off != N)
        throw std::runtime_error("owner_solve_forward_step: X_NR rows (" + std::to_string(N) + ") != slot sum (" +
                                 std::to_string(off) + ") for box " + std::to_string(box.morton_index));
}

/// Backward step of box `box` on its entry `e` (symmetric branch): needs the
/// final x of every neighbour eliminated after it (use_full_set == 1) and the
/// final skeleton entries of the others, read from this rank's entries.
template<typename CoordType, typename DataType>
void owner_solve_backward_step(TreeLevel<CoordType, DataType>& level,
                               std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
                               SolveDataRequest<CoordType, DataType>& e, const BoxData<CoordType, DataType>& box) {
    using owner_solve_detail::gemm;
    const auto& S = box.skeleton_indices;
    const auto& R = box.redundant_indices;
    if (S.empty() || R.empty()) return;
    const int64_t k = static_cast<int64_t>(S.size()), r = static_cast<int64_t>(R.size());
    const int64_t nrhs = e.nrhs;
    if (nrhs <= 0 || e.num_points * nrhs != static_cast<int64_t>(e.left_side.size()))
        throw std::runtime_error("owner_solve_backward_step: invalid batched vector layout");
    std::vector<DataType> x_S(static_cast<size_t>(k * nrhs));
    std::vector<DataType> x_R(static_cast<size_t>(r * nrhs));
    for (int64_t column = 0; column < nrhs; ++column) {
        for (int64_t i = 0; i < k; ++i)
            x_S[static_cast<size_t>(i + column * k)] =
                e.left_side[static_cast<size_t>(S[static_cast<size_t>(i)] +
                                                column * e.num_points)];
        for (int64_t i = 0; i < r; ++i)
            x_R[static_cast<size_t>(i + column * r)] =
                e.left_side[static_cast<size_t>(R[static_cast<size_t>(i)] +
                                                column * e.num_points)];
    }
    if (!box.X_SR.is_allocated())
        throw std::runtime_error("owner_solve_backward_step: X_SR missing for box " + std::to_string(box.morton_index));
    // x_R += X̃_SR^T x_S
    gemm<DataType>('T', static_cast<int>(r), static_cast<int>(nrhs),
                   static_cast<int>(k), box.X_SR.data.data(),
                   static_cast<int>(box.X_SR.lda), x_S.data(), static_cast<int>(k),
                   DataType{1}, x_R.data(), static_cast<int>(r));
    if (box.X_NR.is_allocated() && !box.one_hop.empty()) {
        if (box.use_full_set.size() != box.one_hop.size())
            throw std::runtime_error("owner_solve_backward_step: use_full_set missing for box " + std::to_string(box.morton_index));
        const int64_t N = box.X_NR.rows;
        std::vector<DataType> v(static_cast<size_t>(N * nrhs), DataType{0});
        int64_t off = 0;
        for (size_t i = 0; i < box.one_hop.size(); ++i) {
            const int64_t nb = box.one_hop[i];
            const SolveDataRequest<CoordType, DataType>* ne = resolve_any_solve_data_for_morton(level, level_solve_data, nb);
            if (ne == nullptr)
                throw std::runtime_error("owner_solve_backward_step: neighbour " + std::to_string(nb) + " of box " +
                                         std::to_string(box.morton_index) + " has no solve entry here");
            if (ne->nrhs != nrhs)
                throw std::runtime_error("owner_solve_backward_step: neighbour RHS count mismatch");
            if (box.use_full_set[i] == 1) {
                const int64_t n_i = ne->num_points;
                if (off + n_i > N) throw std::runtime_error("owner_solve_backward_step: X_NR rows short (full slot)");
                for (int64_t column = 0; column < nrhs; ++column)
                    std::copy_n(ne->left_side.begin() + column * ne->num_points,
                                n_i, v.begin() + off + column * N);
                off += n_i;
            } else {
                const int64_t n_i = static_cast<int64_t>(ne->skeleton_indices.size());
                if (off + n_i > N) throw std::runtime_error("owner_solve_backward_step: X_NR rows short (skeleton slot)");
                for (int64_t column = 0; column < nrhs; ++column)
                    for (int64_t j = 0; j < n_i; ++j)
                        v[static_cast<size_t>(off + j + column * N)] =
                            ne->left_side[static_cast<size_t>(
                                ne->skeleton_indices[static_cast<size_t>(j)] +
                                column * ne->num_points)];
                off += n_i;
            }
        }
        if (off != N)
            throw std::runtime_error("owner_solve_backward_step: X_NR rows (" + std::to_string(N) + ") != slot sum (" +
                                     std::to_string(off) + ") for box " + std::to_string(box.morton_index));
        // x_R += X̃_NR^T v
        gemm<DataType>('T', static_cast<int>(r), static_cast<int>(nrhs),
                       static_cast<int>(N), box.X_NR.data.data(),
                       static_cast<int>(box.X_NR.lda), v.data(), static_cast<int>(N),
                       DataType{1}, x_R.data(), static_cast<int>(r));
    }
    // x[R] = x_R (replace); x[S] -= T x_R
    for (int64_t column = 0; column < nrhs; ++column)
        for (int64_t i = 0; i < r; ++i)
            e.left_side[static_cast<size_t>(R[static_cast<size_t>(i)] +
                                            column * e.num_points)] =
                x_R[static_cast<size_t>(i + column * r)];
    if (!box.interpolation_matrix.is_allocated())
        throw std::runtime_error("owner_solve_backward_step: T missing for box " + std::to_string(box.morton_index));
    std::vector<DataType> t(static_cast<size_t>(k * nrhs), DataType{0});
    gemm<DataType>('N', static_cast<int>(k), static_cast<int>(nrhs),
                   static_cast<int>(r), box.interpolation_matrix.data.data(),
                   static_cast<int>(box.interpolation_matrix.lda), x_R.data(),
                   static_cast<int>(r), DataType{0}, t.data(), static_cast<int>(k));
    for (int64_t column = 0; column < nrhs; ++column)
        for (int64_t i = 0; i < k; ++i)
            e.left_side[static_cast<size_t>(S[static_cast<size_t>(i)] +
                                            column * e.num_points)] -=
                t[static_cast<size_t>(i + column * k)];
}

// ---------------------------------------------------------------------------
// The boundary sweeps
// ---------------------------------------------------------------------------
namespace owner_solve_detail {

/// Owned boxes grouped by colour index (dataflow::Color value) and wave.
inline std::vector<std::vector<std::vector<int64_t>>> owned_by_color_wave(const OwnerSolveRecord& rec) {
    std::vector<std::vector<std::vector<int64_t>>> g(3, std::vector<std::vector<int64_t>>(static_cast<size_t>(rec.num_waves)));
    if (!rec.active) return g;
    for (int32_t id : rec.owned) {
        const dataflow::Component& c = rec.graph.comps[static_cast<size_t>(id)];
        const int ci = static_cast<int>(c.color);
        if (ci < 0 || ci > 2) throw std::runtime_error("owner_solve: owned component " + std::to_string(id) + " has no boundary colour");
        for (int64_t m : c.boxes) g[static_cast<size_t>(ci)][static_cast<size_t>(rec.wave_of(m))].push_back(m);
    }
    for (auto& per_color : g)
        for (auto& wave : per_color) std::sort(wave.begin(), wave.end());
    return g;
}

/// The solve entry and the BoxData of an owned box on this rank (local box or
/// the owner's ghost copy — the copy it eliminated, holding T, X_SR, X̃_NR).
template<typename CoordType, typename DataType>
inline void owned_box_and_entry(TreeLevel<CoordType, DataType>& level,
                                std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data, int64_t m,
                                SolveDataRequest<CoordType, DataType>*& e, BoxData<CoordType, DataType>*& box) {
    e = resolve_local_or_ghost_solve_data_for_morton(level, level_solve_data, m);
    box = level.find_local_box(m);
    if (box == nullptr) box = level.find_ghost_box(m);
    if (e == nullptr || box == nullptr) {
        const auto solve_it =
            level.ghost_and_assisting_box_points_for_solve_map.find(m);
        const bool in_solve_map =
            solve_it != level.ghost_and_assisting_box_points_for_solve_map.end();
        const bool classified_as_ghost =
            in_solve_map && solve_it->second >= 0 &&
            static_cast<size_t>(solve_it->second) < level.is_ghost_solve.size() &&
            level.is_ghost_solve[static_cast<size_t>(solve_it->second)] != 0;
        throw std::runtime_error(
            "owner_solve: owned box " + std::to_string(m) + " has no " +
            (e == nullptr ? "solve entry" : "BoxData") + " on its owner" +
            " (solve_map=" + (in_solve_map ? "yes" : "no") +
            ", solve_ghost=" + (classified_as_ghost ? "yes" : "no") +
            ", ghost_box=" +
            (level.find_ghost_box(m) != nullptr ? "yes" : "no") + ")");
    }
}

}  // namespace owner_solve_detail

/// Forward boundary sweep of one engine level (collective: every rank calls it).
// ---------------------------------------------------------------------------
// Slot inbox: contributions land in fixed per-target slots (ledger 42)
// ---------------------------------------------------------------------------
// A push sweep (solve forward, multiply V) sends every source's contribution
// to each neighbour.  Records that were buffered, sorted for the canonical
// order and applied one by one are replaced by a slot per (target,
// neighbour), laid out once per level: a source writes its contribution
// into its own slot of the target (no race: slots are disjoint, and a
// source's neighbours are never stepped in the same wave), and the target
// folds its slots into its vector in a FIXED order — the neighbours sorted
// by elimination sequence — which is exactly the order the sorted records
// produced, so the numbers are unchanged.  A target folds twice: the
// contributions of sources eliminated before it just before its own step,
// the rest at the end of the sweep (nothing reads a target's entries in
// between: the push sweeps read only the stepping box's own entries).
// Slot length: a target eliminated after the source receives its full
// entries, otherwise its skeleton entries (the source's use_full_set).
template<typename DataType>
struct OwnerSlotInbox {
    struct Target {
        int64_t morton = 0;
        int32_t seq = 0;                          ///< elimination sequence of the target (INT32_MAX: interior)
        std::vector<int64_t> nbr;                 ///< neighbours (one_hop order)
        std::unordered_map<int64_t, int> slot_of; ///< neighbour morton -> slot
        std::vector<int64_t> off, len;            ///< slot offsets and lengths in `data`
        std::vector<int32_t> nseq;                ///< neighbour sequences
        std::vector<char> full;                   ///< slot carries full entries (else skeleton entries)
        std::vector<char> filled;
        std::vector<int> order;                   ///< slots in fold order
        std::vector<DataType> data;
    };
    std::vector<Target> targets;
    std::unordered_map<int64_t, size_t> index;
    bool descending = false;

    Target* find(int64_t m) {
        auto it = index.find(m);
        return it == index.end() ? nullptr : &targets[it->second];
    }
    /// A contribution of `n` values from `source` for `target`.
    void put(int64_t target, int64_t source, const DataType* v, size_t n) {
        Target* t = find(target);
        if (t == nullptr) throw std::runtime_error("OwnerSlotInbox::put: no slots for box " + std::to_string(target));
        auto it = t->slot_of.find(source);
        if (it == t->slot_of.end())
            throw std::runtime_error("OwnerSlotInbox::put: box " + std::to_string(source) + " is not a neighbour of " + std::to_string(target));
        const int i = it->second;
        if (static_cast<size_t>(t->len[static_cast<size_t>(i)]) != n)
            throw std::runtime_error("OwnerSlotInbox::put: length " + std::to_string(n) + " for slot of " + std::to_string(source) +
                                     " in " + std::to_string(target) + ", expected " + std::to_string(t->len[static_cast<size_t>(i)]));
        std::memcpy(t->data.data() + t->off[static_cast<size_t>(i)], v, n * sizeof(DataType));
        t->filled[static_cast<size_t>(i)] = 1;
    }
    /// Fold the filled slots whose source satisfies `take(seq)` into `e`, in fold order.
    template<typename CoordType, typename Pred>
    int64_t fold(Target& t, SolveDataRequest<CoordType, DataType>& e, Pred take) {
        int64_t n_folded = 0;
        for (int i : t.order) {
            const size_t si = static_cast<size_t>(i);
            if (!t.filled[si] || !take(t.nseq[si])) continue;
            const DataType* v = t.data.data() + t.off[si];
            const size_t n = static_cast<size_t>(t.len[si]);
            if (t.full[si]) {
                if (e.left_side.size() != n) throw std::runtime_error("OwnerSlotInbox::fold: full slot size mismatch for box " + std::to_string(t.morton));
                for (size_t j = 0; j < n; ++j) e.left_side[j] += v[j];
            } else {
                const int64_t k = static_cast<int64_t>(e.skeleton_indices.size());
                if (static_cast<size_t>(k * e.nrhs) != n)
                    throw std::runtime_error("OwnerSlotInbox::fold: skeleton slot size mismatch for box " + std::to_string(t.morton));
                for (int64_t column = 0; column < e.nrhs; ++column)
                    for (int64_t j = 0; j < k; ++j)
                        e.left_side[static_cast<size_t>(
                            e.skeleton_indices[static_cast<size_t>(j)] +
                            column * e.num_points)] +=
                            v[static_cast<size_t>(j + column * k)];
            }
            t.filled[si] = 0;
            ++n_folded;
        }
        return n_folded;
    }
};

namespace owner_solve_detail {

/// Slots for every box this rank holds as a target of a push sweep: its owned
/// boundary boxes and its local interior boxes.  `descending`: fold order by
/// decreasing sequence (multiply V); else increasing (solve forward).
template<typename CoordType, typename DataType>
OwnerSlotInbox<DataType> build_slot_inbox(TreeLevel<CoordType, DataType>& level,
                                          std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
                                          const OwnerSolveRecord& rec, bool descending) {
    OwnerSlotInbox<DataType> ib;
    ib.descending = descending;
    if (!rec.active) return ib;
    auto seq_or_interior = [&](int64_t m) -> int32_t {
        return rec.comp_of(m) < 0 ? std::numeric_limits<int32_t>::max() : rec.seq_of(m);
    };
    std::vector<int64_t> boxes;
    for (int32_t id : rec.owned)
        for (int64_t m : rec.graph.comps[static_cast<size_t>(id)].boxes) boxes.push_back(m);
    for (const auto& b : level.local_boxes)
        if (rec.comp_of(b.morton_index) < 0) boxes.push_back(b.morton_index);
    ib.targets.reserve(boxes.size());
    for (int64_t m : boxes) {
        SolveDataRequest<CoordType, DataType>* e = nullptr;
        BoxData<CoordType, DataType>* box = nullptr;
        owned_box_and_entry(level, level_solve_data, m, e, box);
        typename OwnerSlotInbox<DataType>::Target t;
        t.morton = m;
        t.seq = seq_or_interior(m);
        // an interior target receives only from boundary sources (interior
        // boxes are stepped by the replicated loop, which keeps its own
        // accumulation), so it gets slots for its boundary neighbours only
        const bool interior_target = rec.comp_of(m) < 0;
        for (int64_t s : box->one_hop)
            if (!interior_target || rec.comp_of(s) >= 0) t.nbr.push_back(s);
        const size_t nn = t.nbr.size();
        t.off.resize(nn); t.len.resize(nn); t.nseq.resize(nn); t.full.resize(nn); t.filled.assign(nn, 0);
        int64_t off = 0;
        for (size_t i = 0; i < nn; ++i) {
            const int64_t s = t.nbr[i];
            const int32_t ss = seq_or_interior(s);
            const bool full = ss < t.seq;   // target still live when the source was eliminated
            t.nseq[i] = ss;
            t.full[i] = full ? 1 : 0;
            t.len[i] = full
                ? static_cast<int64_t>(e->left_side.size())
                : static_cast<int64_t>(e->skeleton_indices.size()) * e->nrhs;
            t.off[i] = off;
            off += t.len[i];
            t.slot_of.emplace(s, static_cast<int>(i));
        }
        t.data.assign(static_cast<size_t>(off), DataType{0});
        t.order.resize(nn);
        for (size_t i = 0; i < nn; ++i) t.order[i] = static_cast<int>(i);
        std::sort(t.order.begin(), t.order.end(), [&](int a, int b) {
            const int32_t sa = t.nseq[static_cast<size_t>(a)], sb = t.nseq[static_cast<size_t>(b)];
            if (sa != sb) return descending ? sa > sb : sa < sb;
            return t.nbr[static_cast<size_t>(a)] < t.nbr[static_cast<size_t>(b)];
        });
        ib.index.emplace(m, ib.targets.size());
        ib.targets.push_back(std::move(t));
    }
    return ib;
}

// ---------------------------------------------------------------------------
// Point-to-point rounds (ledger 42)
// ---------------------------------------------------------------------------
// A round used to be an all-to-all of sizes followed by every rank waiting
// for all of its sends and receives, so one late rank held all 64.  Now a
// rank sends one message (possibly empty) to each rank that can need
// something from it in this round, and waits for exactly the messages the
// graph says it can receive.  Both sets are computed from the shared
// component graph, so sender and receiver agree without announcing sizes.
// Each round of each sweep has its own tag, so an early neighbour's next
// round cannot be mistaken for this one.

/// Which ranks a component's boxes can address in a round.
enum class RoundKind {
    PUSH_ALL,      ///< contributions to every neighbour: owners of preds and succs, homes (interior)
    HANDOFF_DOWN,  ///< final values to the owners of adjacent lower-colour components
    HANDOFF_UP     ///< final values to the owners of adjacent higher-colour components
};

inline void component_dests(const OwnerSolveRecord& rec, const dataflow::Component& k, RoundKind kind, std::set<int>& out) {
    auto owner_rank = [&](int32_t c) { return rec.rank_of_pid(static_cast<uint32_t>(rec.graph.comps[static_cast<size_t>(c)].owner)); };
    if (kind == RoundKind::PUSH_ALL || kind == RoundKind::HANDOFF_DOWN)
        for (int32_t c : k.preds) out.insert(owner_rank(c));
    if (kind == RoundKind::PUSH_ALL || kind == RoundKind::HANDOFF_UP)
        for (int32_t c : k.succs) out.insert(owner_rank(c));
    if (kind == RoundKind::PUSH_ALL)
        for (uint32_t pid : k.sharers) out.insert(rec.rank_of_pid(pid));   // interior neighbours live at home
}

/// Destinations of this rank's colour-`ci` components in a round, and the
/// ranks whose colour-`ci` components can address this rank (self excluded).
inline void round_peers(const OwnerSolveRecord& rec, int ci, RoundKind kind, std::set<int>& dests, std::set<int>& senders) {
    dests.clear();
    senders.clear();
    if (!rec.active) return;
    for (int32_t id : rec.owned)
        if (rec.color_index(id) == ci) component_dests(rec, rec.graph.comps[static_cast<size_t>(id)], kind, dests);
    dests.erase(rec.my_rank);
    for (const auto& k : rec.graph.comps) {
        if (rec.color_index(k.id) != ci || k.owner < 0) continue;
        const int r = rec.rank_of_pid(static_cast<uint32_t>(k.owner));
        if (r == rec.my_rank) continue;
        std::set<int> d;
        component_dests(rec, k, kind, d);
        if (d.count(rec.my_rank)) senders.insert(r);
    }
}

/// Write-back round: owners send their non-local boxes to the homes.
template<typename CoordType, typename DataType>
void write_back_peers(TreeLevel<CoordType, DataType>& level, const OwnerSolveRecord& rec, std::set<int>& dests, std::set<int>& senders) {
    dests.clear();
    senders.clear();
    if (!rec.active) return;
    for (int32_t id : rec.owned)
        for (int64_t m : rec.graph.comps[static_cast<size_t>(id)].boxes)
            if (!level.is_box_on_process(m)) dests.insert(rec.home_rank(m));
    for (const auto& b : level.local_boxes) {
        const int32_t c = rec.comp_of(b.morton_index);
        if (c < 0) continue;
        const int r = rec.rank_of_pid(static_cast<uint32_t>(rec.graph.comps[static_cast<size_t>(c)].owner));
        if (r != rec.my_rank) senders.insert(r);
    }
    dests.erase(rec.my_rank);
}

/// Sends of one sweep, completed at its end.
struct P2PSends {
    std::list<std::vector<char>> bufs;
    std::vector<MPI_Request> reqs;
    void finish() {
        if (!reqs.empty()) owner_mpi_check(MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE), "MPI_Waitall(p2p)");
        reqs.clear();
        bufs.clear();
    }
};

/// Post one message per destination (empty when nothing was queued for it).
inline void p2p_send(P2PSends& sends, std::map<int, std::vector<char>>& out,
                     const std::set<int>& dests, int tag,
                     OwnerSolveStats& stats, MPI_Comm comm) {
    for (int d : dests) {
        std::vector<char> buf;
        auto it = out.find(d);
        if (it != out.end()) buf = std::move(it->second);
        if (buf.size() > static_cast<size_t>(std::numeric_limits<int>::max())) throw std::runtime_error("p2p_send: message too large");
        sends.bufs.push_back(std::move(buf));
        MPI_Request r;
        owner_mpi_check(MPI_Isend(sends.bufs.back().data(), static_cast<int>(sends.bufs.back().size()), MPI_CHAR, d, tag,
                                  comm, &r), "MPI_Isend(p2p)");
        sends.reqs.push_back(r);
        stats.bytes_sent += sends.bufs.back().size();
    }
    for (auto& kv : out)
        if (kv.first != -1 && !dests.count(kv.first) && !kv.second.empty())
            throw std::runtime_error("p2p_send: contributions queued for rank " + std::to_string(kv.first) + " which is not a destination of this round");
    out.clear();
    ++stats.rounds;
}

/// Receive exactly one message from each expected sender, parsed into `inbox`.
template<typename DataType>
void p2p_receive(const std::set<int>& senders, int tag,
                 std::vector<OwnerSolveRec<DataType>>& inbox,
                 OwnerSolveStats& stats, MPI_Comm comm) {
    for (size_t i = 0; i < senders.size(); ++i) {
        MPI_Status st;
        owner_mpi_check(MPI_Probe(MPI_ANY_SOURCE, tag, comm, &st), "MPI_Probe(p2p)");
        if (!senders.count(st.MPI_SOURCE))
            throw std::runtime_error("p2p_receive: message from rank " + std::to_string(st.MPI_SOURCE) + " which is not an expected sender");
        int n = 0;
        owner_mpi_check(MPI_Get_count(&st, MPI_CHAR, &n), "MPI_Get_count(p2p)");
        std::vector<char> buf(static_cast<size_t>(n));
        owner_mpi_check(MPI_Recv(buf.data(), n, MPI_CHAR, st.MPI_SOURCE, tag, comm, MPI_STATUS_IGNORE), "MPI_Recv(p2p)");
        stats.bytes_recv += static_cast<size_t>(n);
        if (n) owner_solve_parse<DataType>(buf.data(), buf.data() + n, inbox);
    }
}

/// Deliver received records: ADDs into slots, REPLACEs straight into the entry.
template<typename CoordType, typename DataType>
void deliver(TreeLevel<CoordType, DataType>& level, std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
             OwnerSlotInbox<DataType>* slots, std::vector<OwnerSolveRec<DataType>>& recs, OwnerSolveStats& stats) {
    for (const auto& r : recs) {
        if (r.kind == OSR_REPLACE) {
            SolveDataRequest<CoordType, DataType>* e = resolve_any_solve_data_for_morton(level, level_solve_data, r.target);
            if (e == nullptr) throw std::runtime_error("owner_solve: no solve entry for box " + std::to_string(r.target));
            if (e->left_side.size() != r.v.size()) throw std::runtime_error("owner_solve: REPLACE size mismatch for box " + std::to_string(r.target));
            std::copy(r.v.begin(), r.v.end(), e->left_side.begin());
        } else {
            if (slots == nullptr) throw std::runtime_error("owner_solve: ADD record in a sweep without slots");
            slots->put(r.target, r.source, r.v.data(), r.v.size());
        }
    }
    stats.records += static_cast<int64_t>(recs.size());
    recs.clear();
}

/// Route a step's records: targets held here go to the slots now, the rest
/// are queued for their holder.
template<typename DataType>
void route(const OwnerSolveRecord& rec, OwnerSlotInbox<DataType>& slots, std::vector<OwnerSolveRec<DataType>>& recs,
           std::map<int, std::vector<char>>& out, OwnerSolveStats& stats) {
    for (auto& r : recs) {
        const int h = rec.holder_rank(r.target);
        if (h == rec.my_rank) slots.put(r.target, r.source, r.v.data(), r.v.size());
        else owner_solve_put(out[h], r);
    }
    stats.records += static_cast<int64_t>(recs.size());
    recs.clear();
}

/// Fold every remaining slot of every target (end of a push sweep).
template<typename CoordType, typename DataType>
void fold_all(TreeLevel<CoordType, DataType>& level, std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
              OwnerSlotInbox<DataType>& slots) {
    #pragma omp parallel for schedule(dynamic)
    for (int64_t i = 0; i < static_cast<int64_t>(slots.targets.size()); ++i) {
        auto& t = slots.targets[static_cast<size_t>(i)];
        SolveDataRequest<CoordType, DataType>* e = resolve_any_solve_data_for_morton(level, level_solve_data, t.morton);
        if (e == nullptr) throw std::runtime_error("owner_solve: no solve entry for slot target " + std::to_string(t.morton));
        slots.template fold<CoordType>(t, *e, [](int32_t) { return true; });
    }
}

constexpr int TAG_SOLVE_FWD = 3410, TAG_SOLVE_BWD = 3420, TAG_MUL_FWD = 3430, TAG_MUL_BWD = 3440;

}  // namespace owner_solve_detail

/// Forward boundary sweep of one engine level (collective: every rank calls it).
template<typename CoordType, typename DataType>
void owner_solve_forward_boundary(TreeLevel<CoordType, DataType>& level,
                                  std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
                                  const OwnerSolveRecord& rec, OwnerSolveStats& stats,
                                  MPI_Comm comm) {
    using namespace owner_solve_detail;
    using clock = std::chrono::high_resolution_clock;
    const auto groups = owned_by_color_wave(rec);
    OwnerSlotInbox<DataType> slots = build_slot_inbox(level, level_solve_data, rec, /*descending=*/false);
    const int nthreads = std::max(1, omp_get_max_threads());
    std::vector<OwnerSolveRec<DataType>> inbox;
    P2PSends sends;
    int round = 0;
    for (int ci = 0; ci < 3; ++ci) {   // blue, purple, green
        std::map<int, std::vector<char>> out;
        for (int w = 0; w < rec.num_waves; ++w) {
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
                    // everything eliminated before this box has contributed: fold it in
                    auto* t = slots.find(m);
                    if (t == nullptr) throw std::runtime_error("owner_solve_forward_boundary: no slots for owned box " + std::to_string(m));
                    const int32_t seq = rec.seq_of(m);
                    slots.template fold<CoordType>(*t, *e, [seq](int32_t s) { return s < seq; });
                    owner_solve_forward_step(level, level_solve_data, *e, *box, seq, per_thread[static_cast<size_t>(omp_get_thread_num())]);
                } catch (...) {
                    if (!failed.exchange(true)) { std::lock_guard<std::mutex> l(ex_mutex); ex = std::current_exception(); }
                }
            }
            if (ex) std::rethrow_exception(ex);
            stats.boxes += static_cast<int64_t>(wave.size());
            for (auto& v : per_thread) route(rec, slots, v, out, stats);
            stats.t_step_ms += std::chrono::duration<double, std::milli>(clock::now() - t0).count();
        }
        // hand this colour's contributions to their holders; take what this
        // colour's components elsewhere produced for boxes held here
        std::set<int> dests, senders;
        round_peers(rec, ci, RoundKind::PUSH_ALL, dests, senders);
        const auto t2 = clock::now();
        p2p_send(sends, out, dests, TAG_SOLVE_FWD + round, stats, comm);
        p2p_receive<DataType>(senders, TAG_SOLVE_FWD + round, inbox, stats, comm);
        ++round;
        stats.t_round_ms += std::chrono::duration<double, std::milli>(clock::now() - t2).count();
        const auto t3 = clock::now();
        deliver(level, level_solve_data, &slots, inbox, stats);
        stats.t_apply_ms += std::chrono::duration<double, std::milli>(clock::now() - t3).count();
    }
    // later sources' skeleton contributions and the interior's inputs
    {
        const auto t3 = clock::now();
        fold_all(level, level_solve_data, slots);
        stats.t_apply_ms += std::chrono::duration<double, std::milli>(clock::now() - t3).count();
    }
    // write-back of non-local owned boxes to their homes
    {
        std::map<int, std::vector<char>> out;
        if (rec.active)
            for (int32_t id : rec.owned)
                for (int64_t m : rec.graph.comps[static_cast<size_t>(id)].boxes) {
                    if (level.is_box_on_process(m)) continue;
                    SolveDataRequest<CoordType, DataType>* e = nullptr;
                    BoxData<CoordType, DataType>* box = nullptr;
                    owned_box_and_entry(level, level_solve_data, m, e, box);
                    OwnerSolveRec<DataType> r;
                    r.target = m; r.source = m; r.seq = 0; r.kind = OSR_REPLACE; r.v = e->left_side;
                    owner_solve_put(out[rec.home_rank(m)], r);
                }
        std::set<int> dests, senders;
        write_back_peers(level, rec, dests, senders);
        const auto t2 = clock::now();
        p2p_send(sends, out, dests, TAG_SOLVE_FWD + round, stats, comm);
        p2p_receive<DataType>(senders, TAG_SOLVE_FWD + round, inbox, stats, comm);
        stats.t_round_ms += std::chrono::duration<double, std::milli>(clock::now() - t2).count();
        deliver(level, level_solve_data, static_cast<OwnerSlotInbox<DataType>*>(nullptr), inbox, stats);
    }
    sends.finish();
}

/// Backward boundary sweep of one engine level (collective: every rank calls it).
template<typename CoordType, typename DataType>
void owner_solve_backward_boundary(TreeLevel<CoordType, DataType>& level,
                                   std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data,
                                   const OwnerSolveRecord& rec, OwnerSolveStats& stats,
                                   MPI_Comm comm) {
    using namespace owner_solve_detail;
    using clock = std::chrono::high_resolution_clock;
    const auto groups = owned_by_color_wave(rec);
    std::vector<OwnerSolveRec<DataType>> inbox;
    P2PSends sends;
    int round = 0;
    for (int ci = 2; ci >= 0; --ci) {   // green, purple, blue
        for (int w = rec.num_waves - 1; w >= 0; --w) {
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
                    owner_solve_backward_step(level, level_solve_data, *e, *box);
                } catch (...) {
                    if (!failed.exchange(true)) { std::lock_guard<std::mutex> l(ex_mutex); ex = std::current_exception(); }
                }
            }
            if (ex) std::rethrow_exception(ex);
            stats.boxes += static_cast<int64_t>(wave.size());
            stats.t_step_ms += std::chrono::duration<double, std::milli>(clock::now() - t0).count();
        }
        if (ci == 0) break;   // nothing is eliminated before blue
        // final values of this colour's boxes to the owners of adjacent lower-colour components
        std::map<int, std::vector<char>> out;
        if (rec.active)
            for (int32_t id : rec.owned) {
                if (rec.color_index(id) != ci) continue;
                for (int64_t m : rec.graph.comps[static_cast<size_t>(id)].boxes) {
                    SolveDataRequest<CoordType, DataType>* e = nullptr;
                    BoxData<CoordType, DataType>* box = nullptr;
                    owned_box_and_entry(level, level_solve_data, m, e, box);
                    std::set<int> dd;
                    for (int64_t nb : box->one_hop) {
                        const int32_t cn = rec.comp_of(nb);
                        if (cn >= 0 && rec.color_index(cn) < ci)
                            dd.insert(rec.rank_of_pid(static_cast<uint32_t>(rec.graph.comps[static_cast<size_t>(cn)].owner)));
                    }
                    for (int d : dd) {
                        if (d == rec.my_rank) continue;
                        OwnerSolveRec<DataType> r;
                        r.target = m; r.source = m; r.seq = 0; r.kind = OSR_REPLACE; r.v = e->left_side;
                        owner_solve_put(out[d], r);
                    }
                }
            }
        std::set<int> dests, senders;
        round_peers(rec, ci, RoundKind::HANDOFF_DOWN, dests, senders);
        const auto t2 = clock::now();
        p2p_send(sends, out, dests, TAG_SOLVE_BWD + round, stats, comm);
        p2p_receive<DataType>(senders, TAG_SOLVE_BWD + round, inbox, stats, comm);
        ++round;
        stats.t_round_ms += std::chrono::duration<double, std::milli>(clock::now() - t2).count();
        deliver(level, level_solve_data, static_cast<OwnerSlotInbox<DataType>*>(nullptr), inbox, stats);
    }
    {
        std::map<int, std::vector<char>> out;
        if (rec.active)
            for (int32_t id : rec.owned)
                for (int64_t m : rec.graph.comps[static_cast<size_t>(id)].boxes) {
                    if (level.is_box_on_process(m)) continue;
                    SolveDataRequest<CoordType, DataType>* e = nullptr;
                    BoxData<CoordType, DataType>* box = nullptr;
                    owned_box_and_entry(level, level_solve_data, m, e, box);
                    OwnerSolveRec<DataType> r;
                    r.target = m; r.source = m; r.seq = 0; r.kind = OSR_REPLACE; r.v = e->left_side;
                    owner_solve_put(out[rec.home_rank(m)], r);
                }
        std::set<int> dests, senders;
        write_back_peers(level, rec, dests, senders);
        const auto t2 = clock::now();
        p2p_send(sends, out, dests, TAG_SOLVE_BWD + round, stats, comm);
        p2p_receive<DataType>(senders, TAG_SOLVE_BWD + round, inbox, stats, comm);
        stats.t_round_ms += std::chrono::duration<double, std::milli>(clock::now() - t2).count();
        deliver(level, level_solve_data, static_cast<OwnerSlotInbox<DataType>*>(nullptr), inbox, stats);
    }
    sends.finish();
}

inline void owner_solve_print_stats(const OwnerSolveStats& s, int level_index, const char* sweep) {
    std::printf("  [owner-solve] level %d %s: %lld owned boxes, %d rounds, %lld records, sent %.2f MB, recv %.2f MB | ms: steps %.1f, rounds %.1f, apply %.1f\n",
                level_index, sweep, static_cast<long long>(s.boxes), s.rounds, static_cast<long long>(s.records),
                s.bytes_sent / 1048576.0, s.bytes_recv / 1048576.0, s.t_step_ms, s.t_round_ms, s.t_apply_ms);
    std::fflush(stdout);
}

/// Checksum of the local solution vectors of a level (Σ|x| over local boxes,
/// folded in box order; ranks summed in rank order on the print rank), under
/// FMM_BOX_CHECKSUM — the bitwise oracle for the solve.  Collective over comm.
template<typename CoordType, typename DataType>
void print_solve_checksum(const std::vector<SolveDataRequest<CoordType, DataType>>& level_solve_data, int level_index,
                          const char* sweep, MPI_Comm comm, int print_rank, int rank, bool active) {
    if (!box_checksum_enabled()) return;
    double v = 0.0;
    if (active)
        for (const auto& e : level_solve_data) {
            double s = 0.0;
            for (const auto& x : e.left_side) s += static_cast<double>(std::abs(x));
            v += s;
        }
    int size = 1;
    MPI_Comm_size(comm, &size);
    std::vector<double> all(static_cast<size_t>(size), 0.0);
    MPI_Gather(&v, 1, MPI_DOUBLE, all.data(), 1, MPI_DOUBLE, print_rank, comm);
    if (rank != print_rank) return;
    double tot = 0.0;
    for (int r = 0; r < size; ++r) tot += all[static_cast<size_t>(r)];
    std::printf("  [solve-checksum] %s level %d: %.17g\n", sweep, level_index, tot);
    std::fflush(stdout);
}

}  // namespace fmm
