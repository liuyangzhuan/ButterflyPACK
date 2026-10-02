#pragma once
// Adaptive ID training rows (H2_ID_proxy 2) of a batch of boxes on the device
// (GPU_PROXY2_PLAN.md, A3).
//
// Mirrors append_adaptive_id_training_rows (color_CA/factorization.hpp) with
// the box's sketch as the ID target: per box, the far field is sampled node by
// node of the source tree, breadth first; a node's sampled rows are tested
// against the current ID (their residual), the independent residual rows join
// the target, and a node that has not converged is subdivided.  The numbers
// are computed on the device for all boxes together; the decisions, the index
// of unused far points and the node queues stay on the host, which steers by
// one small download per round.
//
//   selector.start(targets);
//   while (selector.active()) {
//       selector.plan(meta);            // launch items of the round, into the image
//       md = meta.upload(...);
//       selector.launch(md, stream);    // (the first round: after the sketches)
//       copy selector.status() to the host; synchronize;
//       selector.finish(host copy);     // decisions: the next round of each box
//   }
//
// A round of a box is one of
//   sample   the next node: its sampled and hold-out rows, their residuals and
//            norms, the row ID of the sampled rows' residual, and those rows
//            appended to the target; before it, where due, the hold-out rows
//            of the previous node and a new ID of the target;
//   recheck  a node that converged under an ID with residual IDs: a new ID of
//            the target, then the residuals again;
//   extra    the ID of a residual that did not converge (a residual ID);
//   final    the ID of the finished target, unless it is current.
// The common case is one sample round (the root node converges and its
// hold-out rows confirm it) and the final round.
//
// The box's ID is left as launch_qrcp leaves it: rank, pivots, flag in the
// target's outputs, T in id_factor() with leading dimension id_ld().

#ifdef H2_HAVE_GPU

#include "device_adaptive.hpp"
#include "device_heap.hpp"
#include "device_kernels.hpp"
#include "gpu_runtime.hpp"

#include <omp.h>

#include <algorithm>
#include <cstdint>
#include <exception>
#include <functional>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <utility>
#include <vector>

namespace fmm {
namespace gpu {

template<typename CoordType, typename DataType>
class AdaptiveRowSelector {
public:
    using S = typename DeviceScalar<DataType>::type;
    using Tree = ParallelTree<CoordType, DataType>;
    using Box = BoxData<CoordType, DataType>;

    // One box: its points (the columns), its ID target w (cap x n: the
    // sketch in rows 0..n-1), a block f of the same shape for the in-place
    // ID, and the outputs of its ID.
    struct Target {
        const Box* box = nullptr;
        int n = 0;
        IndexList cols;
        S* w = nullptr;
        S* f = nullptr;
        int cap = 0;
        int* jpvt = nullptr;
        int* rank = nullptr;
        double* norm = nullptr;  // norm of the sketch (the first ID's input)
        int* flag = nullptr;
    };

    struct Hooks {
        // Point slots (consecutive, from the returned one) of the points with
        // these global indices, until the next call.
        std::function<int(const std::vector<int64_t>&, cudaStream_t)> place_points;
        // Kernel entries of the items, which are in `meta` at `offset`
        // (device image md): out(i, j) = K(point rows[i], point cols[j]).
        std::function<void(const std::vector<EvalItemT<S>>&, MetaBuilder&, const char*, size_t, int, int, cudaStream_t)>
            eval;
    };

    AdaptiveRowSelector(const Tree* tree, double tolerance, DeviceHeap& heap, Hooks hooks)
        : tree_(tree), tolerance_(tolerance), heap_(heap), hooks_(std::move(hooks)) {}
    AdaptiveRowSelector(const AdaptiveRowSelector&) = delete;
    AdaptiveRowSelector& operator=(const AdaptiveRowSelector&) = delete;
    // (after the stream has run the last round)
    ~AdaptiveRowSelector() {
        for (char* p : transient_) heap_.free(p);
        for (FrameBlock& b : frame_blocks_) heap_.free(b.ptr);
        for (char* p : factor_blocks_) heap_.free(p);
        heap_.free(d_status_);
    }

    static int sample_batch(const Tree* tree) { return static_cast<int>(std::max<int64_t>(1, tree->id_adaptive_batch)); }
    static int holdout_count(const Tree* tree) {
        return static_cast<int>(std::min<int64_t>(4, std::max<int64_t>(1, sample_batch(tree) / 2)));
    }
    // Rows of a box's target for its sketch and the rows its first node can add.
    static int first_capacity(const Tree* tree, int n) { return n + sample_batch(tree) + holdout_count(tree); }

    void start(std::vector<Target> targets) {
        validate_id_source_index(tree_);
        states_.clear();
        states_.resize(targets.size());
        std::exception_ptr failure;
        std::mutex failure_mutex;
        #pragma omp parallel for schedule(dynamic, 8)
        for (int64_t i = 0; i < static_cast<int64_t>(targets.size()); ++i) {
            try {
                State& st = states_[static_cast<size_t>(i)];
                st.t = targets[static_cast<size_t>(i)];
                st.rows = st.t.n;
                st.stats.active = true;
                const std::vector<IDAdaptiveFrame> far = make_initial_id_far_frames(tree_, st.t.box);
                if (!far.empty()) {
                    std::vector<std::pair<int64_t, int64_t>> ranges;
                    ranges.reserve(far.size());
                    for (const IDAdaptiveFrame& frame : far) ranges.emplace_back(frame.lo, frame.hi);
                    st.unused = std::make_unique<IDUnusedIndexTree>(tree_->num_points, std::move(ranges));
                    const auto root = id_source_point_range(tree_, 0, 0);
                    st.queue.push_back({0, 0, root.first, root.second - 1});
                }
                next_node(st);
            } catch (...) {
                std::lock_guard<std::mutex> lock(failure_mutex);
                if (!failure) failure = std::current_exception();
            }
        }
        if (failure) std::rethrow_exception(failure);
        // status of every box: doubles (residual norm, target norm, hold-out
        // residual norm, two unused QR norms), then ints (rows selected,
        // hold-out rows selected, residual ID rank, two unused QR flags)
        const size_t nb = states_.size();
        status_ints_ = align_up(kStatusDoubles * nb * sizeof(double), 16);
        status_bytes_ = align_up(status_ints_ + kStatusInts * nb * sizeof(int), 16);
        d_status_ = heap_.alloc(std::max<size_t>(status_bytes_, 16));
        rounds_ = 0;
    }

    bool active() const {
        for (const State& st : states_) {
            if (st.step != Step::Done) return true;
        }
        return false;
    }
    int rounds() const { return rounds_; }
    char* status() const { return d_status_; }
    size_t status_bytes() const { return status_bytes_; }

    const S* id_factor(size_t i) const { return states_[i].id_f; }
    int id_ld(size_t i) const { return states_[i].id_ld; }
    const IDAdaptiveStats& stats(size_t i) const { return states_[i].stats; }
    // Blocks that hold a final ID's factor outside the targets' own f: the
    // caller frees them with those (from the heap).
    std::vector<char*> take_factor_blocks() { return std::move(factor_blocks_); }

    // ---- one round: the items of every box's next step, into the image
    void plan(MetaBuilder& meta) {
        ++rounds_;
        grow_.clear(); copies_.clear(); late_copies_.clear();
        appends_.clear(); sample_appends_.clear();
        ids_.clear(); selections_.clear(); extra_ids_.clear();
        evals_.clear(); residuals_.clear(); transposes_.clear();
        max_id_n_ = max_sel_n_ = max_extra_n_ = 0;
        max_copy_m_ = max_copy_n_ = max_eval_m_ = max_eval_n_ = max_tr_m_ = max_tr_n_ = 0;
        const size_t nb = states_.size();
        const size_t D = sizeof(S), I = sizeof(int);

        // the nodes sampled in this round: their rows' global indices and the
        // round's buffers
        std::vector<int64_t> point_ids;
        std::vector<int> point_at(nb, 0);
        size_t frame_bytes = 0;
        std::vector<size_t> frame_at(nb, 0);
        const int batch = sample_batch(tree_), holdouts = holdout_count(tree_);
        for (size_t i = 0; i < nb; ++i) {
            State& st = states_[i];
            st.launched = st.step;
            if (st.step != Step::Sample) continue;
            const std::vector<int64_t> samples =
                pick_evenly_spaced_unused(*st.unused, st.frame.lo, st.frame.hi, batch, true);
            st.holdout_positions =
                pick_evenly_spaced_unused(*st.unused, st.frame.lo, st.frame.hi, holdouts, false);
            st.s = static_cast<int>(samples.size());
            st.h = static_cast<int>(st.holdout_positions.size());
            ++st.stats.frames;
            st.stats.sampled += st.s;
            point_at[i] = static_cast<int>(point_ids.size());
            for (int64_t p : samples) point_ids.push_back(tree_->id_source_point_order[static_cast<size_t>(p)]);
            for (int64_t p : st.holdout_positions) {
                point_ids.push_back(tree_->id_source_point_order[static_cast<size_t>(p)]);
            }
            const size_t sh = static_cast<size_t>(st.s + st.h), n = static_cast<size_t>(st.t.n);
            frame_at[i] = frame_bytes;
            frame_bytes += 2 * align_up(sh * n * D) + align_up(n * static_cast<size_t>(st.s) * D) +
                           align_up(n * static_cast<size_t>(st.h) * D) + align_up(static_cast<size_t>(st.s) * I) +
                           align_up(static_cast<size_t>(st.h) * I);
        }
        int frame_block = -1;
        char* frame = nullptr;
        if (frame_bytes > 0) {
            frame = heap_.alloc(frame_bytes);
            frame_block = static_cast<int>(frame_blocks_.size());
            frame_blocks_.push_back(FrameBlock{frame, 0});
        }
        point_slot0_ = point_ids.empty() ? 0 : hooks_.place_points(point_ids, Context::instance().stream());

        for (size_t i = 0; i < nb; ++i) {
            State& st = states_[i];
            const int n = st.t.n;
            double* dstat = reinterpret_cast<double*>(d_status_) + kStatusDoubles * i;
            int* istat = reinterpret_cast<int*>(d_status_ + status_ints_) + kStatusInts * i;
            // the IDs a residual is taken under: the target's, then the residual IDs
            auto id_list = [&]() {
                std::vector<IdRefT<S>> ids;
                ids.push_back(IdRefT<S>{st.id_f, st.id_ld, st.t.jpvt, st.t.rank, 0});
                for (const Extra& e : st.extras) ids.push_back(IdRefT<S>{e.f, e.ld, e.jpvt, nullptr, e.rank});
                return static_cast<int64_t>(meta.append(ids));
            };
            auto copy = [&](std::vector<GatherItemT<S>>& list, S* out, int ld, int m, const S* src, int ld_src) {
                list.push_back(GatherItemT<S>{out, ld, m, n, src, 1, ld_src, IndexList{}, IndexList{}});
                max_copy_m_ = std::max(max_copy_m_, m);
                max_copy_n_ = std::max(max_copy_n_, n);
            };
            // hold-out rows of the previous node, selected on the device
            auto append_holdout = [&]() {
                if (st.append_holdout == 0) return;
                appends_.push_back(AppendRowsItemT<S>{st.t.w, st.t.cap, st.rows, st.app_src, st.app_lds, st.app_h, n,
                                                      istat + kHoldoutSelected, st.app_sel});
                st.rows += st.append_holdout;
                st.append_holdout = 0;
                st.released_block = st.app_block;  // (after this round)
                st.app_block = -1;
            };
            // a new ID of the target (in place in f)
            auto new_id = [&]() {
                copy(copies_, st.t.f, st.t.cap, st.rows, st.t.w, st.t.cap);
                ids_.push_back(QrcpItemT<S>{st.t.f, st.rows, n, st.t.cap, st.t.jpvt, st.t.rank,
                                            st.id_rows < 0 ? st.t.norm : dstat + kNormA, st.t.flag});
                max_id_n_ = std::max(max_id_n_, n);
                st.id_f = st.t.f;
                st.id_ld = st.t.cap;
                st.id_rows = st.rows;
                st.extras.clear();
                st.extra_rank = 0;
                st.need_id = false;
            };
            // residuals of the node's rows (sampled, and hold-out when the
            // ID they are judged under is this one), with their norms
            auto residuals = [&](bool with_holdout) {
                const int sh = st.s + st.h;
                copy(late_copies_, st.res, sh, sh, st.sh, sh);
                const int64_t ids = id_list();
                const int nids = 1 + static_cast<int>(st.extras.size());
                residuals_.push_back(IdResidualItemT<S>{st.res, sh, st.s, n, ids, nids, dstat + kResidualNorm});
                if (with_holdout) {
                    residuals_.push_back(
                        IdResidualItemT<S>{st.res + st.s, sh, st.h, n, ids, nids, dstat + kHoldoutNorm});
                    transposes_.push_back(TransposeItemT<S>{st.res + st.s, sh, st.hrest, n, st.h, n});
                    selections_.push_back(QrcpItemT<S>{st.hrest, n, st.h, n, st.hsel, istat + kHoldoutSelected,
                                                       dstat + kNormB, istat + kFlagB});
                    max_sel_n_ = std::max(max_sel_n_, st.h);
                    max_tr_m_ = std::max(max_tr_m_, st.h);
                    max_tr_n_ = std::max(max_tr_n_, n);
                }
                // the target's norm (before this node's rows join it)
                residuals_.push_back(IdResidualItemT<S>{st.t.w, st.t.cap, st.rows, n, 0, 0, dstat + kTargetNorm});
                st.holdout_valid = with_holdout;
            };

            switch (st.step) {
                case Step::Sample: {
                    // room for the node's rows
                    const int need = st.rows + st.append_holdout + st.s + st.h;
                    if (need > st.t.cap) {
                        const int cap = need + 2 * (batch + holdouts);
                        const size_t bytes = align_up(static_cast<size_t>(cap) * n * D);
                        S* w = heap_.template alloc<S>(bytes);
                        S* f = heap_.template alloc<S>(bytes);
                        transient_.push_back(reinterpret_cast<char*>(w));
                        factor_blocks_.push_back(reinterpret_cast<char*>(f));
                        copy(grow_, w, cap, st.rows, st.t.w, st.t.cap);
                        st.t.w = w;
                        st.t.f = f;
                        st.t.cap = cap;
                    }
                    append_holdout();
                    if (st.need_id) new_id();
                    char* at = frame + frame_at[i];
                    const size_t sh = static_cast<size_t>(st.s + st.h);
                    st.sh = reinterpret_cast<S*>(at);    at += align_up(sh * n * D);
                    st.res = reinterpret_cast<S*>(at);   at += align_up(sh * n * D);
                    st.rest = reinterpret_cast<S*>(at);  at += align_up(static_cast<size_t>(n) * st.s * D);
                    st.hrest = reinterpret_cast<S*>(at); at += align_up(static_cast<size_t>(n) * st.h * D);
                    st.sel = reinterpret_cast<int*>(at); at += align_up(static_cast<size_t>(st.s) * I);
                    st.hsel = reinterpret_cast<int*>(at);
                    st.frame_block = frame_block;
                    ++frame_blocks_[static_cast<size_t>(frame_block)].live;
                    IndexList rows;
                    rows.base = point_slot0_ + point_at[i];
                    evals_.push_back(EvalItemT<S>{st.sh, st.s + st.h, st.s + st.h, n, rows, st.t.cols});
                    max_eval_m_ = std::max(max_eval_m_, st.s + st.h);
                    max_eval_n_ = std::max(max_eval_n_, n);
                    residuals(st.extras.empty() && st.h > 0);
                    transposes_.push_back(TransposeItemT<S>{st.res, st.s + st.h, st.rest, n, st.s, n});
                    max_tr_m_ = std::max(max_tr_m_, st.s);
                    max_tr_n_ = std::max(max_tr_n_, n);
                    selections_.push_back(QrcpItemT<S>{st.rest, n, st.s, n, st.sel, istat + kSelected, dstat + kNormA,
                                                       istat + kFlagA});
                    max_sel_n_ = std::max(max_sel_n_, st.s);
                    sample_appends_.push_back(AppendRowsItemT<S>{st.t.w, st.t.cap, st.rows, st.sh, st.s + st.h, st.s, n,
                                                                 istat + kSelected, st.sel});
                    break;
                }
                case Step::Recheck:
                    new_id();
                    residuals(st.h > 0);
                    break;
                case Step::Extra: {
                    // ID of the residual, to the tolerance that keeps the absolute threshold
                    const size_t bytes = align_up(static_cast<size_t>(st.s) * n * D);
                    char* block = heap_.alloc(bytes + align_up(static_cast<size_t>(n) * I));
                    transient_.push_back(block);
                    st.extra_f = reinterpret_cast<S*>(block);
                    st.extra_jpvt = reinterpret_cast<int*>(block + bytes);
                    copy(late_copies_, st.extra_f, st.s, st.s, st.res, st.s + st.h);
                    QrcpItemT<S> item{st.extra_f, st.s, n, st.s, st.extra_jpvt, istat + kExtraRank, dstat + kNormA,
                                      istat + kFlagA};
                    item.tol = std::min(1.0, st.ref / st.rn * tolerance_);
                    extra_ids_.push_back(item);
                    max_extra_n_ = std::max(max_extra_n_, n);
                    break;
                }
                case Step::Final:
                    append_holdout();
                    if (st.need_id || st.id_rows != st.rows) new_id();
                    break;
                case Step::Done:
                    break;
            }
        }
        off_grow_ = meta.append(grow_);
        off_appends_ = meta.append(appends_);
        off_copies_ = meta.append(copies_);
        off_ids_ = meta.append(ids_);
        off_evals_ = meta.append(evals_);
        off_late_copies_ = meta.append(late_copies_);
        off_residuals_ = meta.append(residuals_);
        off_transposes_ = meta.append(transposes_);
        off_selections_ = meta.append(selections_);
        off_sample_appends_ = meta.append(sample_appends_);
        off_extra_ids_ = meta.append(extra_ids_);
        meta_ = &meta;
    }

    // ---- ... and their launches, in the order the steps need
    void launch(const char* md, cudaStream_t stream) {
        auto count = [](const auto& v) { return static_cast<int>(v.size()); };
        auto qrcp = [&](const std::vector<QrcpItemT<S>>& items, size_t offset, int max_n) {
            if (items.empty()) return;
            char* work = nullptr;
            if (const size_t w = qrcp_work_bytes<S>(count(items), max_n)) work = heap_.alloc(w);
            launch_qrcp(reinterpret_cast<const QrcpItemT<S>*>(md + offset), count(items), max_n, tolerance_, work,
                        stream);
            heap_.free(work);  // later launches are ordered after it
        };
        launch_gather(reinterpret_cast<const GatherItemT<S>*>(md + off_grow_), count(grow_), max_copy_m_, max_copy_n_,
                      md, stream);
        launch_append_rows(reinterpret_cast<const AppendRowsItemT<S>*>(md + off_appends_), count(appends_), stream);
        launch_gather(reinterpret_cast<const GatherItemT<S>*>(md + off_copies_), count(copies_), max_copy_m_,
                      max_copy_n_, md, stream);
        qrcp(ids_, off_ids_, max_id_n_);
        if (!evals_.empty()) hooks_.eval(evals_, *meta_, md, off_evals_, max_eval_m_, max_eval_n_, stream);
        launch_gather(reinterpret_cast<const GatherItemT<S>*>(md + off_late_copies_), count(late_copies_),
                      max_copy_m_, max_copy_n_, md, stream);
        launch_id_residual(reinterpret_cast<const IdResidualItemT<S>*>(md + off_residuals_), count(residuals_), md,
                           stream);
        launch_transpose(reinterpret_cast<const TransposeItemT<S>*>(md + off_transposes_), count(transposes_),
                         max_tr_m_, max_tr_n_, stream);
        qrcp(selections_, off_selections_, max_sel_n_);
        launch_append_rows(reinterpret_cast<const AppendRowsItemT<S>*>(md + off_sample_appends_),
                           count(sample_appends_), stream);
        qrcp(extra_ids_, off_extra_ids_, max_extra_n_);
    }

    // ---- the round's results (a host copy of status()): each box's next step
    void finish(const char* status) {
        for (size_t i = 0; i < states_.size(); ++i) {
            State& st = states_[i];
            const double* dstat = reinterpret_cast<const double*>(status) + kStatusDoubles * i;
            const int* istat = reinterpret_cast<const int*>(status + status_ints_) + kStatusInts * i;
            if (st.released_block >= 0) {
                release_frame(st.released_block);
                st.released_block = -1;
            }
            switch (st.launched) {
                case Step::Sample: {
                    const int selected = istat[kSelected];
                    st.rows += selected;
                    st.stats.appended += selected;
                    st.rn = dstat[kResidualNorm];
                    st.ref = dstat[kTargetNorm];
                    const bool converged = st.rn < st.ref * tolerance_;
                    if (converged && !st.extras.empty()) {
                        ++st.stats.recomputes;
                        st.step = Step::Recheck;
                    } else {
                        judge(st, converged, dstat[kHoldoutNorm], istat[kHoldoutSelected]);
                    }
                    break;
                }
                case Step::Recheck:
                    st.rn = dstat[kResidualNorm];
                    st.ref = dstat[kTargetNorm];
                    judge(st, st.rn < st.ref * tolerance_, dstat[kHoldoutNorm], istat[kHoldoutSelected]);
                    break;
                case Step::Extra: {
                    const int rank = istat[kExtraRank];
                    if (rank > 0 && st.extra_rank + rank <= st.t.n) {
                        st.extras.push_back(Extra{st.extra_f, st.s, st.extra_jpvt, rank});
                        st.extra_rank += rank;
                        ++st.stats.extra_ids;
                    } else {
                        st.need_id = true;
                        ++st.stats.recomputes;
                    }
                    append_id_source_children(tree_, st.frame, st.queue);
                    node_done(st);
                    break;
                }
                case Step::Final:
                    st.step = Step::Done;
                    break;
                case Step::Done:
                    break;
            }
        }
    }

private:
    enum class Step { Sample, Recheck, Extra, Final, Done };
    struct Extra {
        const S* f;
        int ld;
        const int* jpvt;
        int rank;
    };
    struct FrameBlock {
        char* ptr;
        int live;  // boxes whose node still uses it
    };
    struct State {
        Target t;
        std::unique_ptr<IDUnusedIndexTree> unused;
        std::vector<IDAdaptiveFrame> queue;
        size_t head = 0;
        IDAdaptiveFrame frame{};
        Step step = Step::Final, launched = Step::Done;
        int rows = 0;      // rows of the target
        int id_rows = -1;  // ... when its ID was last computed (-1: never)
        const S* id_f = nullptr;
        int id_ld = 0;
        std::vector<Extra> extras;  // residual IDs since the last ID of the target
        int extra_rank = 0;
        bool need_id = true;  // the next step starts with a new ID of the target
        // the node being sampled: s sampled and h hold-out rows (sh, their
        // residuals res, the transposed residuals, the row IDs' pivots)
        int s = 0, h = 0;
        std::vector<int64_t> holdout_positions;
        S* sh = nullptr;
        S* res = nullptr;
        S* rest = nullptr;
        S* hrest = nullptr;
        int* sel = nullptr;
        int* hsel = nullptr;
        int frame_block = -1;
        bool holdout_valid = false;
        double rn = 0.0, ref = 0.0;
        S* extra_f = nullptr;
        int* extra_jpvt = nullptr;
        // hold-out rows that failed, to append at the box's next step
        int append_holdout = 0;
        const S* app_src = nullptr;
        int app_lds = 0, app_h = 0;
        const int* app_sel = nullptr;
        int app_block = -1, released_block = -1;
        IDAdaptiveStats stats;
    };

    enum { kResidualNorm = 0, kTargetNorm, kHoldoutNorm, kNormA, kNormB, kStatusDoubles };
    enum { kSelected = 0, kHoldoutSelected, kExtraRank, kFlagA, kFlagB, kStatusInts };

    void release_frame(int block) {
        FrameBlock& b = frame_blocks_[static_cast<size_t>(block)];
        if (--b.live == 0) {
            heap_.free(b.ptr);  // (its launches have run: finish follows a synchronization)
            b.ptr = nullptr;
        }
    }
    // The next node with unused far points, or the final ID.
    void next_node(State& st) {
        while (st.head < st.queue.size()) {
            const IDAdaptiveFrame frame = st.queue[st.head++];
            if (st.unused->count(frame.lo, frame.hi) == 0) continue;
            st.frame = frame;
            st.step = Step::Sample;
            return;
        }
        st.step = (st.need_id || st.append_holdout > 0 || st.id_rows != st.rows) ? Step::Final : Step::Done;
    }
    void node_done(State& st) {
        if (st.frame_block >= 0) {
            if (st.app_block == st.frame_block) {
                // (the pending hold-out rows live in it: released after their append)
            } else {
                release_frame(st.frame_block);
            }
            st.frame_block = -1;
        }
        next_node(st);
    }
    // After the residual of a node's sampled rows (norm rn against the
    // target's norm ref): the hold-out test of a converged node, then whether
    // the node is subdivided, as append_adaptive_id_training_rows decides.
    void judge(State& st, bool converged, double holdout_norm, int holdout_selected) {
        const double threshold = st.ref * tolerance_;
        bool holdout_failed = false;
        if (converged && st.h > 0) {
            if (!st.holdout_valid) throw std::runtime_error("AdaptiveRowSelector: hold-out residual missing");
            st.stats.holdout += st.h;
            if (holdout_norm >= threshold) {
                holdout_failed = true;
                st.append_holdout = holdout_selected;
                st.stats.appended += holdout_selected;
                st.app_src = st.sh + st.s;
                st.app_lds = st.s + st.h;
                st.app_h = st.h;
                st.app_sel = st.hsel;
                st.app_block = st.frame_block;
                for (int64_t position : st.holdout_positions) st.unused->mark_used(position);
                st.need_id = true;
                ++st.stats.recomputes;
            }
        }
        if (converged && !holdout_failed) {
            node_done(st);
            return;
        }
        if (st.unused->count(st.frame.lo, st.frame.hi) == 0) {
            node_done(st);
            return;
        }
        if (!holdout_failed && st.rn > 0.0) {
            st.step = Step::Extra;
            return;
        }
        append_id_source_children(tree_, st.frame, st.queue);
        node_done(st);
    }

    const Tree* tree_;
    double tolerance_;
    DeviceHeap& heap_;
    Hooks hooks_;
    std::vector<State> states_;
    char* d_status_ = nullptr;
    size_t status_ints_ = 0, status_bytes_ = 0;
    int rounds_ = 0;
    int point_slot0_ = 0;
    std::vector<FrameBlock> frame_blocks_;
    std::vector<char*> transient_;       // grown targets, residual IDs: until the end
    std::vector<char*> factor_blocks_;   // grown f blocks (may hold a final ID)
    // the round's items
    MetaBuilder* meta_ = nullptr;
    std::vector<GatherItemT<S>> grow_, copies_, late_copies_;
    std::vector<AppendRowsItemT<S>> appends_, sample_appends_;
    std::vector<QrcpItemT<S>> ids_, selections_, extra_ids_;
    std::vector<EvalItemT<S>> evals_;
    std::vector<IdResidualItemT<S>> residuals_;
    std::vector<TransposeItemT<S>> transposes_;
    size_t off_grow_ = 0, off_appends_ = 0, off_copies_ = 0, off_ids_ = 0, off_evals_ = 0, off_late_copies_ = 0,
           off_residuals_ = 0, off_transposes_ = 0, off_selections_ = 0, off_sample_appends_ = 0, off_extra_ids_ = 0;
    int max_id_n_ = 0, max_sel_n_ = 0, max_extra_n_ = 0;
    int max_copy_m_ = 0, max_copy_n_ = 0, max_eval_m_ = 0, max_eval_n_ = 0, max_tr_m_ = 0, max_tr_n_ = 0;
};

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
