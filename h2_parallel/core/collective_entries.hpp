#pragma once
// Entries through a collective block callback (elem_extract=1 with a block
// callback; GPU_BACKEND/h2_gpu/GPU_CA_PLAN.md, "Entries from a collective
// callback"): the callback is called by every rank together (as
// c_bpack_extractelement_bridge, which extracts them from another
// matrix), so the H2 construction cannot ask for a block when it needs it.
// Instead a phase runs its items (boxes, candidates, parents, sources) with
// the entries it has; a request for a block it does not have records the
// block and abandons the item (H2EntryMiss), which must not have changed
// anything outside itself yet; then one round fetches every rank's recorded
// blocks in one callback call, and the abandoned items run again from the
// start (their computations are deterministic: the sketches are seeded per
// box).  A phase ends when a round finds no rank with a request.
//
// The rounds are collective over the H2's communicator, so every rank takes
// part in each: the ranks running a section of the factorization (a level)
// with flag kRoundActive, the others serving (serve) until the section's
// end (end_section).  With entries through a local callback (elem_extract 0
// or 2) none of this exists: H2Kernel::collective is null.
//
// With a sketch callback too (c_bpack_set_sketch_callback) a round also
// fetches sketched blocks, Omega A(rows, columns) for a sparse Omega the
// asker gives, and with a transposed one (c_bpack_set_sketch_transpose_callback)
// Omega A(columns, rows)^T: the ID's sketch of its kernel rows asks for those instead of
// the rows' entries (a d x n block for an m x n one, m >> d).
//
// A round is the same calls wherever a rank makes it, so the ranks need not
// be at the same point of their work: a phase ends at the first round in
// which no rank asks for a block, whichever rounds the others are in.  The
// device path (H2_use_gpu; gpu/collective_evaluator.hpp) fetches the blocks
// a batch lacks in one round on the spot (fetch_now, flag kRoundBusy), and
// before any other communication between the level's ranks (a transport,
// the level's end, the transition's decisions) waits in rounds until no rank
// is busy (quiesce): a rank blocked elsewhere would leave the others' rounds
// without it.

#include "extraction_profile.hpp"
#include "tree.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <exception>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace fmm {

// A block not fetched yet: the item that asked for it is abandoned.
struct H2EntryMiss : std::exception {
    const char* what() const noexcept override {
        return "H2 entries: a block not fetched yet (collective callback)";
    }
};

// The item a thread runs in a phase (entry_item_begin .. entry_item_end).
// A request for a block not fetched yet inside an item records the block and
// goes on with zeros (discovery: the item's other requests are recorded in
// the same pass), the item doomed; at its next checkpoint (before it changes
// anything outside itself, and at its end) a doomed item is abandoned.
// Outside an item such a request throws at once, and so does one after the
// item's checkpoint (its changes have begun, so it could not run again): a
// request the item should have made before it.
struct EntryItemState {
    bool active = false;
    bool doomed = false;
    bool committed = false;  // (past its checkpoint)
};
inline EntryItemState& entry_item() {
    thread_local EntryItemState state;
    return state;
}
inline void entry_item_begin() { entry_item() = EntryItemState{true, false, false}; }
inline void entry_checkpoint() {
    if (entry_item().doomed) {
        entry_item() = EntryItemState{};
        throw H2EntryMiss();
    }
    entry_item().committed = entry_item().active;
}
// A request for a block not fetched yet (recorded by the caller): the item
// is doomed and goes on with zeros; see above.
inline void entry_item_miss() {
    EntryItemState& state = entry_item();
    if (!state.active) throw H2EntryMiss();
    if (state.committed) {
        throw std::logic_error("H2 entries: a block not fetched yet asked for after its item's checkpoint");
    }
    state.doomed = true;
}
inline void entry_item_end() {
    const bool doomed = entry_item().doomed;
    entry_item() = EntryItemState{};
    if (doomed) throw H2EntryMiss();
}
// In an item's catch-all: whether the failure is that of a doomed item
// (computing on zeros), which is abandoned, not reported.
inline bool entry_item_abandoned() {
    const bool doomed = entry_item().doomed;
    entry_item() = EntryItemState{};
    return doomed;
}

// What asked for a block (the profile of the rounds, report): the host
// phases, or the device path's ID sketch, box region, the owner pass's new
// fill targets, remote generators and parent transition, set around its
// evaluations by the thread that runs them.
enum class EntryOrigin : int { host = 0, sketch, box, fill, remote, transition, count };
inline EntryOrigin& entry_origin() {
    thread_local EntryOrigin origin = EntryOrigin::host;
    return origin;
}
struct EntryOriginScope {
    explicit EntryOriginScope(EntryOrigin origin) : saved_(entry_origin()) { entry_origin() = origin; }
    ~EntryOriginScope() { entry_origin() = saved_; }
    EntryOriginScope(const EntryOriginScope&) = delete;
    EntryOriginScope& operator=(const EntryOriginScope&) = delete;

private:
    EntryOrigin saved_;
};

template<typename DataType, typename CallbackValue>
class CollectiveEntries {
public:
    using Callback = void (*)(int*, int*, int*, int64_t*, int*, int*, CallbackValue*, int*, int*, int*, int*, int*,
                              void*);
    // (c_bpack_func_zmn_sketch, BPACK_wrapper.h)
    using SketchCallback = void (*)(int*, int*, int*, int64_t*, int64_t*, int*, int*, int*, int*, int*, int*, int*,
                                    int*, double*, CallbackValue*, int*, int*, int*, void*);
    // (c_bpack_func_zmn_sketch_device, BPACK_wrapper.h)
    using SketchDeviceCallback = void (*)(int*, int*, int*, int64_t*, int*, int*, int*, int*, int*, int*, int*, int*,
                                          double*, void**, int64_t*, int*, void*, int*, void*);
    // (c_bpack_func_zmn_block_device, BPACK_wrapper.h)
    using BlockDeviceCallback = void (*)(int*, int*, int*, int64_t*, int*, int*, int*, int*, void**, int64_t*, void*,
                                         int*, void*);
    static constexpr int kRoundRequest = 1;  // (a rank has recorded blocks)
    static constexpr int kRoundActive = 2;   // (a section is running)
    static constexpr int kRoundBusy = 4;     // (a rank may still ask, in the device path's sense)

    // symmetric: the matrix is transpose-symmetric (option%sym=1), so a block
    // also serves as its transpose's; a nonsymmetric matrix asks for both
    // orientations of a pair of boxes, and each is fetched apart.
    CollectiveEntries(MPI_Comm comm, Callback callback, void* quant, int callback_pid, bool symmetric = true)
        : callback_(callback), quant_(quant), symmetric_(symmetric) {
        MPI_Comm_dup(comm, &comm_);
        MPI_Comm_size(comm_, &size_);
        pids_.resize(static_cast<size_t>(size_));
        MPI_Allgather(&callback_pid, 1, MPI_INT, pids_.data(), 1, MPI_INT, comm_);
    }
    ~CollectiveEntries() {
        int finalized = 0;
        MPI_Finalized(&finalized);
        if (!finalized && comm_ != MPI_COMM_NULL) MPI_Comm_free(&comm_);
    }
    CollectiveEntries(const CollectiveEntries&) = delete;
    CollectiveEntries& operator=(const CollectiveEntries&) = delete;

    // A[i + j lda] = K(x_i, y_j) if this block (or its transpose, for a
    // symmetric matrix) was fetched; else false.  Thread safe.
    bool lookup(const int64_t* x, int64_t nx, const int64_t* y, int64_t ny, DataType* A, int64_t lda) const {
        if (const Block* b = find(x, nx, y, ny)) {
            for (int64_t j = 0; j < ny; ++j) {
                for (int64_t i = 0; i < nx; ++i) A[i + j * lda] = b->values[static_cast<size_t>(i + j * nx)];
            }
            return true;
        }
        if (!symmetric_) return false;
        if (const Block* b = find(y, ny, x, nx)) {
            for (int64_t j = 0; j < ny; ++j) {
                for (int64_t i = 0; i < nx; ++i) A[i + j * lda] = b->values[static_cast<size_t>(j + i * ny)];
            }
            return true;
        }
        return false;
    }

    // Drops the block (or its transpose, as lookup) once read: the device
    // path reads each block it fetched once.  Not concurrent with lookups.
    void release(const int64_t* x, int64_t nx, const int64_t* y, int64_t ny) {
        for (int t = 0; t < (symmetric_ ? 2 : 1); ++t) {
            const int64_t* a = t == 0 ? x : y;
            const int64_t* b = t == 0 ? y : x;
            const int64_t na = t == 0 ? nx : ny, nb = t == 0 ? ny : nx;
            auto range = index_.equal_range(hash(a, na, b, nb));
            for (auto it = range.first; it != range.second; ++it) {
                Block& blk = blocks_[it->second];
                if (!same(blk, a, na, b, nb)) continue;
                blk = Block{};
                index_.erase(it);
                return;
            }
        }
    }

    // A sketch callback: sketched blocks are fetched through it; the
    // transposed sketches (Omega A(columns, rows)^T) through the transposed
    // one.  A nonsymmetric matrix's ID needs both (its K(F, B) rows sketch
    // A(F, B) from the left, its K(B, F)^T rows A(B, F) from the right), and
    // fetches its kernel rows' blocks without them.
    void set_sketch(SketchCallback callback) { sketch_callback_ = callback; }
    void set_sketch_transposed(SketchCallback callback) { sketch_t_callback_ = callback; }
    // A sketch callback that adds the sketches into device memory instead
    // (one rank: fetch_sketches_device); either orientation.
    void set_sketch_device(SketchDeviceCallback callback) { sketch_device_callback_ = callback; }
    bool has_sketch_device() const { return sketch_device_callback_ != nullptr; }
    // A block callback that adds the blocks into device memory instead (one
    // rank: fetch_blocks_device).
    void set_block_device(BlockDeviceCallback callback) { block_device_callback_ = callback; }
    bool has_block_device() const { return block_device_callback_ != nullptr; }
    bool symmetric() const { return symmetric_; }
    bool has_sketch() const {
        return sketch_callback_ != nullptr && (symmetric_ || sketch_t_callback_ != nullptr);
    }

    // A sketch request: Y = Omega K(x, y), Omega (d x nx) given by its nnz
    // entries (orow: the position in x, odest: the row of Y, oval); with
    // transposed, Y = Omega K(y, x)^T.
    struct SketchKey {
        const int64_t* x;
        int64_t nx;
        const int64_t* y;
        int64_t ny;
        int d;
        int64_t nnz;
        const int* orow;
        const int* odest;
        const double* oval;
        bool transposed = false;
    };
    // Y (d x ny, leading dimension ldy) if this request was fetched; else
    // false.  Thread safe.
    bool lookup_sketch(const SketchKey& k, DataType* Y, int64_t ldy) const {
        const Sketch* s = find_sketch(k, sketch_index_, sketches_);
        if (s == nullptr) return false;
        for (int64_t j = 0; j < k.ny; ++j) {
            for (int64_t i = 0; i < k.d; ++i) Y[i + j * ldy] = s->values[static_cast<size_t>(i + j * k.d)];
        }
        return true;
    }
    // Records the request for the next round.  Thread safe.
    void record_sketch(const SketchKey& k) { record_sketch_to(k, nullptr); }
    // As record_sketch, the round writing Y (d x ny, column major, leading
    // dimension d) into dst instead of keeping it, so it is never looked up;
    // dst must stay valid until the round.  Thread safe.
    void record_sketch_into(const SketchKey& k, DataType* dst) { record_sketch_to(k, dst); }
    void record_sketch_to(const SketchKey& k, DataType* dst) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (Sketch* p = const_cast<Sketch*>(find_sketch(k, pending_sketch_index_, pending_sketches_))) {
            add_destination(*p, dst);
            return;
        }
        Sketch s;
        s.dst = dst;
        s.keep = dst == nullptr;
        s.rows.assign(k.x, k.x + k.nx);
        s.cols.assign(k.y, k.y + k.ny);
        s.d = k.d;
        s.orow.assign(k.orow, k.orow + k.nnz);
        s.odest.assign(k.odest, k.odest + k.nnz);
        s.oval.assign(k.oval, k.oval + k.nnz);
        s.transposed = k.transposed;
        const size_t o = static_cast<size_t>(entry_origin());
        origin_blocks_[o] += 1.0;
        origin_entries_[o] += static_cast<double>(k.d) * static_cast<double>(k.ny);
        pending_sketch_index_.emplace(sketch_hash(k), pending_sketches_.size());
        pending_sketches_.push_back(std::move(s));
    }
    // Drops a fetched request once read (the device path).  Not concurrent
    // with lookups.
    void release_sketch(const SketchKey& k) {
        auto range = sketch_index_.equal_range(sketch_hash(k));
        for (auto it = range.first; it != range.second; ++it) {
            if (!same_sketch(sketches_[it->second], k)) continue;
            sketches_[it->second] = Sketch{};
            sketch_index_.erase(it);
            return;
        }
    }

    // One rank (a sequential front's H2, comm MPI_COMM_SELF): the device
    // path's blocks and sketches come straight from the callbacks into its
    // buffer, packed in its order, without being recorded, gathered, stored
    // and looked up again.  Those three copies of every value and the
    // recording's hashing were a quarter of a sequential front's host time
    // (Poisson 400^3, 64 ranks, GPU); the device path reads each block once,
    // and with one rank nobody else asks.
    bool alone() const { return size_ == 1; }
    struct BlockRef {
        const int64_t* x;
        int64_t nx;
        const int64_t* y;
        int64_t ny;
    };
    // K(x_n, y_n), column major, packed in order into out (zeroed, as the
    // callback expects): one callback call.
    void fetch_blocks_direct(const std::vector<BlockRef>& blocks, DataType* out) {
        if (!alone()) throw std::logic_error("CollectiveEntries: a direct fetch over several ranks");
        if (blocks.empty()) return;
        const double t_gather = MPI_Wtime();
        DirectBlocks a = pack_blocks(blocks);
        count_blocks(blocks);
        std::vector<int> pgidx(blocks.size(), 0);
        int npmap = 1;
        int pmaps[3] = {1, 1, pids_[0]};
        const double t_callback = MPI_Wtime();
        gather_seconds_ += t_callback - t_gather;
        callback_(&a.ninter, &a.nallrows, &a.nallcols, &a.nloc, a.allrows.data(), a.allcols.data(),
                  reinterpret_cast<CallbackValue*>(out), a.rowidx.data(), a.colidx.data(), pgidx.data(), &npmap, pmaps,
                  quant_);
        const double t_end = MPI_Wtime();
        callback_seconds_ += t_end - t_callback;
        fetched_entries_ += static_cast<double>(a.nloc);
        fetched_blocks_ += a.ninter;
        ++fetch_rounds_;
        ++rounds_;
        round_seconds_ += t_end - t_gather;
    }
    // As fetch_blocks_direct, through the device block callback: block n
    // added into dout[n] (device memory, leading dimension ldo[n], zeroed by
    // the caller; distinct blocks' outputs do not overlap).  False, with
    // nothing written, when the callback declined; on true its writes are
    // complete.  stream: the caller's, for ordering.
    bool fetch_blocks_device(const std::vector<BlockRef>& blocks, void* const* dout, const int64_t* ldo, void* stream) {
        if (!alone()) throw std::logic_error("CollectiveEntries: a direct fetch over several ranks");
        if (blocks.empty()) return true;
        if (block_device_callback_ == nullptr) return false;
        const double t_gather = MPI_Wtime();
        DirectBlocks a = pack_blocks(blocks);
        std::vector<void*> outs(dout, dout + blocks.size());
        std::vector<int64_t> lds(ldo, ldo + blocks.size());
        int done = 0;
        const double t_callback = MPI_Wtime();
        gather_seconds_ += t_callback - t_gather;
        block_device_callback_(&a.ninter, &a.nallrows, &a.nallcols, &a.nloc, a.allrows.data(), a.allcols.data(),
                               a.rowidx.data(), a.colidx.data(), outs.data(), lds.data(), stream, &done, quant_);
        const double t_end = MPI_Wtime();
        callback_seconds_ += t_end - t_callback;
        if (done == 0) return false;
        count_blocks(blocks);
        fetched_entries_ += static_cast<double>(a.nloc);
        fetched_blocks_ += a.ninter;
        ++fetch_rounds_;
        ++rounds_;
        round_seconds_ += t_end - t_gather;
        return true;
    }
    // Every key's Omega K(x, y) (with `transposed`, all of them: Omega K(y,
    // x)^T), d x ny column major, packed in order into out (zeroed): one call
    // of that orientation's callback.
    void fetch_sketches_direct(const std::vector<SketchKey>& keys, bool transposed, DataType* out) {
        if (!alone()) throw std::logic_error("CollectiveEntries: a direct fetch over several ranks");
        if (keys.empty()) return;
        SketchCallback callback = transposed ? sketch_t_callback_ : sketch_callback_;
        if (callback == nullptr) throw std::logic_error("CollectiveEntries: a sketch request without its callback");
        const double t_gather = MPI_Wtime();
        DirectSketches a = pack_direct(keys, transposed);
        count_direct(keys);
        int npmap = 1;
        int pmaps[3] = {1, 1, pids_[0]};
        std::vector<int> pgidx(keys.size(), 0);
        const double t_callback = MPI_Wtime();
        gather_seconds_ += t_callback - t_gather;
        callback(&a.nsketch, &a.nallrows, &a.nallcols, &a.nallomega, &a.nloc, a.allrows.data(), a.allcols.data(),
                 a.rowidx.data(), a.colidx.data(), a.dims.data(), a.omegaidx.data(), a.omega_row.data(),
                 a.omega_dest.data(), a.omega_val.data(), reinterpret_cast<CallbackValue*>(out), pgidx.data(), &npmap,
                 pmaps, quant_);
        const double t_end = MPI_Wtime();
        callback_seconds_ += t_end - t_callback;
        fetched_entries_ += static_cast<double>(a.nloc);
        ++fetch_rounds_;
        ++rounds_;
        round_seconds_ += t_end - t_gather;
    }
    // As fetch_sketches_direct, through the device sketch callback: key n's
    // sketch added into dout[n] (device memory, d x ny column major, leading
    // dimension ldo[n], zeroed by the caller; distinct keys' outputs do not
    // overlap).  False, with nothing written, when the callback declined;
    // on true its writes are complete.  stream: the caller's, for ordering.
    bool fetch_sketches_device(const std::vector<SketchKey>& keys, bool transposed, void* const* dout,
                               const int64_t* ldo, void* stream) {
        if (!alone()) throw std::logic_error("CollectiveEntries: a direct fetch over several ranks");
        if (keys.empty()) return true;
        if (sketch_device_callback_ == nullptr) return false;
        const double t_gather = MPI_Wtime();
        DirectSketches a = pack_direct(keys, transposed);
        std::vector<void*> outs(dout, dout + keys.size());
        std::vector<int64_t> lds(ldo, ldo + keys.size());
        int t = transposed ? 1 : 0, done = 0;
        const double t_callback = MPI_Wtime();
        gather_seconds_ += t_callback - t_gather;
        sketch_device_callback_(&a.nsketch, &a.nallrows, &a.nallcols, &a.nallomega, a.allrows.data(), a.allcols.data(),
                                a.rowidx.data(), a.colidx.data(), a.dims.data(), a.omegaidx.data(), a.omega_row.data(),
                                a.omega_dest.data(), a.omega_val.data(), outs.data(), lds.data(), &t, stream, &done,
                                quant_);
        const double t_end = MPI_Wtime();
        callback_seconds_ += t_end - t_callback;
        if (done == 0) return false;
        count_direct(keys);
        fetched_entries_ += static_cast<double>(a.nloc);
        ++fetch_rounds_;
        ++rounds_;
        round_seconds_ += t_end - t_gather;
        return true;
    }

    // Records the block for the next round.  Thread safe.
    void record(const int64_t* x, int64_t nx, const int64_t* y, int64_t ny) { record_to(x, nx, y, ny, nullptr); }
    // As record, the round writing the block (nx x ny, column major, leading
    // dimension nx) into dst instead of keeping it, so it is never looked up;
    // dst must stay valid until the round.  The device path's batches read
    // their blocks so: the fetched values went from a fresh callback buffer
    // into the store and from there into the batch's image (Poisson 400^3,
    // 64 ranks, STRUMPACK on the GPU: 53.3 -> 52.4 s).  Thread safe.
    void record_into(const int64_t* x, int64_t nx, const int64_t* y, int64_t ny, DataType* dst) {
        record_to(x, nx, y, ny, dst);
    }
    void record_to(const int64_t* x, int64_t nx, const int64_t* y, int64_t ny, DataType* dst) {
        std::lock_guard<std::mutex> lock(mutex_);
        const uint64_t h = hash(x, nx, y, ny);
        auto range = pending_index_.equal_range(h);
        for (auto it = range.first; it != range.second; ++it) {
            Block& b = pending_[it->second];
            if (same(b, x, nx, y, ny)) {
                const size_t o = static_cast<size_t>(entry_origin());
                dup_blocks_[o] += 1.0;
                dup_entries_[o] += static_cast<double>(nx) * static_cast<double>(ny);
                add_destination(b, dst);
                return;
            }
        }
        Block b;
        b.dst = dst;
        b.keep = dst == nullptr;
        b.rows.assign(x, x + nx);
        b.cols.assign(y, y + ny);
        const size_t o = static_cast<size_t>(entry_origin());
        origin_blocks_[o] += 1.0;
        origin_entries_[o] += static_cast<double>(nx) * static_cast<double>(ny);
        pending_index_.emplace(h, pending_.size());
        pending_.push_back(std::move(b));
    }

    // Collective: one round.  flags | kRoundRequest where this rank has
    // recorded blocks, all ranks' maximum returned; with kRoundRequest set
    // every rank's blocks are fetched (one callback call over all ranks).
    int round(int flags) {
        const double t0 = MPI_Wtime();
        int local = flags | (pending_.empty() && pending_sketches_.empty() ? 0 : kRoundRequest), global = 0;
        MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, comm_);
        wait_seconds_ += MPI_Wtime() - t0;
        if ((global & kRoundRequest) != 0) {
            fetch();
            ++fetch_rounds_;
        }
        ++rounds_;
        round_seconds_ += MPI_Wtime() - t0;
        return global;
    }
    // The blocks recorded so far, fetched now (one round, flags kRoundActive
    // | kRoundBusy): the device path, from one thread.
    void fetch_now() { round(kRoundActive | kRoundBusy); }
    // Rounds until one in which no rank is busy: before communication between
    // the section's ranks that is not a round (see the notes at the top).
    // Every rank of the section calls it at the same points.
    void quiesce() {
        while ((round(kRoundActive) & kRoundBusy) != 0) {}
    }
    // A rank outside a running section: rounds until the section's end.
    void serve() {
        while ((round(0) & kRoundActive) != 0) {}
    }
    // A rank of a section, at its end: serves the others' rounds until every
    // rank of the section has ended it.
    void end_section() {
        serve();
        clear();
    }
    // The fetched blocks are no longer needed (a phase's end).
    void clear() {
        blocks_.clear();
        index_.clear();
        sketches_.clear();
        sketch_index_.clear();
    }
    int64_t rounds() const { return rounds_; }
    int64_t fetch_rounds() const { return fetch_rounds_; }
    double round_seconds() const { return round_seconds_; }
    // The device path's own time with the fetched blocks (the copies into
    // the device's blocks, gpu/collective_evaluator.hpp), for the report.
    void count_delivery(double seconds) { delivery_seconds_ += seconds; }
    // On rank 0 of comm: the rounds, the blocks and entries fetched (all
    // ranks), the slowest rank's time in rounds; then the rounds' time by
    // part, the blocks by what asked for them, and when the callback was a
    // batched extraction from a compressed H2 in this process, its time by
    // step (extraction_profile.hpp).  Collective.
    void report(const char* what) {
        double v[2] = {fetched_entries_, static_cast<double>(fetched_blocks_)}, t = round_seconds_;
        MPI_Allreduce(MPI_IN_PLACE, v, 2, MPI_DOUBLE, MPI_SUM, comm_);
        MPI_Allreduce(MPI_IN_PLACE, &t, 1, MPI_DOUBLE, MPI_MAX, comm_);
        const ExtractionProfile& xp = extraction_profile();
        // rounds' parts, the extraction's steps and its all-to-alls, its work
        std::vector<double> mine = {wait_seconds_, gather_seconds_, callback_seconds_, store_seconds_,
                                    delivery_seconds_};
        mine.insert(mine.end(), xp.seconds, xp.seconds + ExtractionProfile::kSteps);
        mine.insert(mine.end(), {xp.route_mpi, xp.entries_computed, xp.bytes_routed});
        const int n = static_cast<int>(mine.size());
        std::vector<double> sum(mine.size()), min(mine.size()), max(mine.size());
        MPI_Allreduce(mine.data(), sum.data(), n, MPI_DOUBLE, MPI_SUM, comm_);
        MPI_Allreduce(mine.data(), min.data(), n, MPI_DOUBLE, MPI_MIN, comm_);
        MPI_Allreduce(mine.data(), max.data(), n, MPI_DOUBLE, MPI_MAX, comm_);
        constexpr int kOrigins = static_cast<int>(EntryOrigin::count);
        double origins[4 * kOrigins];
        std::copy(origin_blocks_, origin_blocks_ + kOrigins, origins);
        std::copy(origin_entries_, origin_entries_ + kOrigins, origins + kOrigins);
        std::copy(dup_blocks_, dup_blocks_ + kOrigins, origins + 2 * kOrigins);
        std::copy(dup_entries_, dup_entries_ + kOrigins, origins + 3 * kOrigins);
        MPI_Allreduce(MPI_IN_PLACE, origins, 4 * kOrigins, MPI_DOUBLE, MPI_SUM, comm_);
        double sketch_counts[3] = {sketch_requests_, sketch_values_, sketch_covered_};
        MPI_Allreduce(MPI_IN_PLACE, sketch_counts, 3, MPI_DOUBLE, MPI_SUM, comm_);
        long long calls = static_cast<long long>(xp.calls);
        MPI_Allreduce(MPI_IN_PLACE, &calls, 1, MPI_LONG_LONG, MPI_MAX, comm_);
        int rank = 0;
        MPI_Comm_rank(comm_, &rank);
        if (rank == 0) {
            std::printf("H2 entries through the collective callback (%s): %lld rounds, %lld fetching, %.0f blocks, "
                        "%.3g entries, %.2f s in rounds (slowest rank)\n",
                        what, static_cast<long long>(rounds_), static_cast<long long>(fetch_rounds_), v[1], v[0], t);
            auto mean_max = [&](int i) {
                char b[64];
                std::snprintf(b, sizeof(b), "%.2f / %.2f / %.2f", sum[static_cast<size_t>(i)] / size_,
                              min[static_cast<size_t>(i)], max[static_cast<size_t>(i)]);
                return std::string(b);
            };
            std::printf("  rounds by part (s, mean / min / max over ranks): waiting for the other ranks %s, request lists "
                        "%s, callback %s, storing %s; the device's copies of the blocks %s\n",
                        mean_max(0).c_str(), mean_max(1).c_str(), mean_max(2).c_str(), mean_max(3).c_str(),
                        mean_max(4).c_str());
            if (sketch_counts[0] > 0) {
                std::printf("  sketched blocks: %.0f, %.3g values for %.3g entries of their rows\n", sketch_counts[0],
                            sketch_counts[1], sketch_counts[2]);
            }
            static const char* const names[kOrigins] = {"host phases", "ID sketch",         "box region",
                                                        "fill targets", "remote generators", "transition"};
            std::printf("  fetched by origin (blocks, entries):");
            for (int o = 0; o < kOrigins; ++o) {
                std::printf(" %s %.0f, %.3g%s", names[o], origins[o], origins[kOrigins + o], o + 1 < kOrigins ? ";" : "\n");
            }
            // (requests of a block already pending in the round: fetched once)
            std::printf("  requested again while pending, not fetched (blocks, entries):");
            for (int o = 0; o < kOrigins; ++o) {
                std::printf(" %s %.0f, %.3g%s", names[o], origins[2 * kOrigins + o], origins[3 * kOrigins + o],
                            o + 1 < kOrigins ? ";" : "\n");
            }
            if (calls > 0) {
                constexpr int x0 = 5, xs = ExtractionProfile::kSteps;
                std::printf("  the callback's extraction by step (s, mean / min / max over ranks, %lld calls):", calls);
                for (int s = 0; s < xs; ++s) {
                    if (s == ExtractionProfile::kRoute) {
                        std::printf(" %s %s (all-to-alls %s);", ExtractionProfile::kNames[s], mean_max(x0 + s).c_str(),
                                    mean_max(x0 + xs).c_str());
                    } else {
                        std::printf(" %s %s%s", ExtractionProfile::kNames[s], mean_max(x0 + s).c_str(),
                                    s + 1 < xs ? ";" : "");
                    }
                }
                std::printf("\n  the extraction's work: %.3g entries computed, at most %.3g on a rank; %.2f GB of "
                            "values routed, at most %.2f GB from a rank\n",
                            sum[static_cast<size_t>(x0 + xs + 1)], max[static_cast<size_t>(x0 + xs + 1)],
                            sum[static_cast<size_t>(x0 + xs + 2)] / 1e9, max[static_cast<size_t>(x0 + xs + 2)] / 1e9);
            }
            std::fflush(stdout);
        }
    }
    void reset_counts() {
        rounds_ = fetch_rounds_ = fetched_blocks_ = 0;
        fetched_entries_ = round_seconds_ = 0.0;
        wait_seconds_ = gather_seconds_ = callback_seconds_ = store_seconds_ = delivery_seconds_ = 0.0;
        sketch_requests_ = sketch_values_ = sketch_covered_ = 0.0;
        std::fill(origin_blocks_, origin_blocks_ + static_cast<int>(EntryOrigin::count), 0.0);
        std::fill(origin_entries_, origin_entries_ + static_cast<int>(EntryOrigin::count), 0.0);
        std::fill(dup_blocks_, dup_blocks_ + static_cast<int>(EntryOrigin::count), 0.0);
        std::fill(dup_entries_, dup_entries_ + static_cast<int>(EntryOrigin::count), 0.0);
        reset_extraction_profile();
    }
    int64_t fetched_blocks() const { return fetched_blocks_; }
    double fetched_entries() const { return fetched_entries_; }

private:
    // A pending request's destinations (record_into): dst and, for the same
    // request recorded again, extra; keep: also recorded without one, so its
    // values are stored for lookups too.
    struct Block {
        std::vector<int64_t> rows, cols;
        std::vector<DataType> values;
        DataType* dst = nullptr;
        std::vector<DataType*> extra;
        bool keep = true;
    };
    struct Sketch {
        std::vector<int64_t> rows, cols;
        int d = 0;
        std::vector<int> orow, odest;
        std::vector<double> oval;
        bool transposed = false;
        std::vector<DataType> values;
        DataType* dst = nullptr;
        std::vector<DataType*> extra;
        bool keep = true;
    };
    template<typename Request>
    static void add_destination(Request& r, DataType* dst) {
        if (dst == nullptr) {
            r.keep = true;
        } else if (r.dst == nullptr) {
            r.dst = dst;
        } else {
            r.extra.push_back(dst);
        }
    }
    // The callback's buffer for a round's values (kept over the rounds: a
    // fresh one a round faulted its pages in again), or, when every request
    // of this rank's has a destination and these lie back to back in its
    // order, the first of them, so that the callback writes them in place.
    // Zeroed either way, as the callback expects.
    template<typename Request, typename Tag>
    struct RoundValues {
        RoundValues(const std::vector<Request*>& mine, int64_t nloc, int64_t (*len)(const Request&))
            : in_place(direct(mine, nloc, len)),
              buffer(in_place ? 0 : static_cast<size_t>(std::max<int64_t>(nloc, 1))),
              data(in_place ? reinterpret_cast<CallbackValue*>(mine.front()->dst) : buffer.data()) {
            parallel_zero(data, static_cast<size_t>(in_place ? nloc : std::max<int64_t>(nloc, 1)));
        }
        static bool direct(const std::vector<Request*>& mine, int64_t nloc, int64_t (*len)(const Request&)) {
            if (nloc <= 0) return false;
            const DataType* next = nullptr;
            for (const Request* r : mine) {
                if (r->dst == nullptr || r->keep || (next != nullptr && r->dst != next)) return false;
                next = r->dst + len(*r);
            }
            return true;
        }
        // the requests' values to their destinations (in place: to their
        // duplicates' only), a request a thread
        void deliver(const std::vector<Request*>& mine, int64_t (*len)(const Request&)) const {
            std::vector<int64_t> at(mine.size() + 1, 0);
            for (size_t k = 0; k < mine.size(); ++k) at[k + 1] = at[k] + len(*mine[k]);
#pragma omp parallel for schedule(dynamic, 1) if (mine.size() > 1 && !omp_in_parallel())
            for (size_t k = 0; k < mine.size(); ++k) {
                const Request& r = *mine[k];
                const DataType* v = reinterpret_cast<const DataType*>(data) + at[k];
                const int64_t n = len(r);
                if (r.dst != nullptr && !in_place) std::copy(v, v + n, r.dst);
                for (DataType* e : r.extra) std::copy(v, v + n, e);
            }
        }
        const bool in_place;
        KeptBuffer<CallbackValue, Tag> buffer;
        CallbackValue* const data;
    };
    static uint64_t sketch_hash(const SketchKey& k) {
        uint64_t h = hash(k.x, k.nx, k.y, k.ny) ^ (static_cast<uint64_t>(k.d) * 0x9e3779b97f4a7c15ULL) ^
                     (k.transposed ? 0x5bd1e995ULL : 0ULL);
        for (int64_t e = 0; e < k.nnz; ++e) {
            uint64_t bits = 0;
            std::memcpy(&bits, &k.oval[e], sizeof(bits));
            h ^= (static_cast<uint64_t>(k.orow[e]) << 20 ^ static_cast<uint64_t>(k.odest[e]) ^ bits) +
                 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        }
        return h;
    }
    static bool same_sketch(const Sketch& s, const SketchKey& k) {
        if (s.d != k.d || s.transposed != k.transposed || static_cast<int64_t>(s.rows.size()) != k.nx ||
            static_cast<int64_t>(s.cols.size()) != k.ny ||
            static_cast<int64_t>(s.orow.size()) != k.nnz) {
            return false;
        }
        return std::equal(s.rows.begin(), s.rows.end(), k.x) && std::equal(s.cols.begin(), s.cols.end(), k.y) &&
               std::equal(s.orow.begin(), s.orow.end(), k.orow) && std::equal(s.odest.begin(), s.odest.end(), k.odest) &&
               std::equal(s.oval.begin(), s.oval.end(), k.oval);
    }
    static const Sketch* find_sketch(const SketchKey& k, const std::unordered_multimap<uint64_t, size_t>& index,
                                     const std::vector<Sketch>& list) {
        auto range = index.equal_range(sketch_hash(k));
        for (auto it = range.first; it != range.second; ++it) {
            if (same_sketch(list[it->second], k)) return &list[it->second];
        }
        return nullptr;
    }
    // Every rank's list, stacked by rank (the same on every rank), each
    // rank's length known everywhere (lengths: from the lists gathered
    // before it).
    template<typename T>
    std::vector<T> gatherv(const std::vector<T>& local, MPI_Datatype type, const std::vector<int64_t>& lengths) {
        std::vector<int> ns(static_cast<size_t>(size_)), displs(static_cast<size_t>(size_));
        int64_t total = 0;
        for (int p = 0; p < size_; ++p) {
            ns[static_cast<size_t>(p)] = checked(lengths[static_cast<size_t>(p)], "a request list");
            displs[static_cast<size_t>(p)] = checked(total, "the request lists");
            total += lengths[static_cast<size_t>(p)];
        }
        std::vector<T> all(static_cast<size_t>(total));
        MPI_Allgatherv(local.data(), checked(static_cast<int64_t>(local.size()), "a request list"), type, all.data(),
                       ns.data(), displs.data(), type, comm_);
        return all;
    }

    // A sketch request's Omega on the wire.  Every rank receives every
    // rank's requests, and Omega's entries as (row, destination, value) were
    // 89% of the bytes the rounds gathered (Poisson 400^3, 64 ranks, the
    // distributed fronts), although both paths draw +-one scale a request.
    // An entry is two varints: its row's step from the previous entry's
    // (zigzag), and 2 x its destination (zigzag) + its value's sign; values
    // that are not +-one scale follow as doubles.  About 3 bytes an entry
    // instead of 16, read back to the same lists in the same order.
    static void put_varint(std::vector<unsigned char>& out, uint64_t v) {
        for (; v >= 0x80; v >>= 7) out.push_back(static_cast<unsigned char>(v | 0x80));
        out.push_back(static_cast<unsigned char>(v));
    }
    static uint64_t get_varint(const unsigned char*& p) {
        uint64_t v = 0;
        for (int shift = 0;; shift += 7) {
            const unsigned char c = *p++;
            v |= static_cast<uint64_t>(c & 0x7f) << shift;
            if (c < 0x80) return v;
        }
    }
    static uint64_t zigzag(int64_t v) { return (static_cast<uint64_t>(v) << 1) ^ static_cast<uint64_t>(v >> 63); }
    static int64_t unzigzag(uint64_t v) { return static_cast<int64_t>(v >> 1) ^ -static_cast<int64_t>(v & 1); }
    // s's Omega appended to out; whether its values are +-scale (scale:
    // the first one's magnitude)
    static bool put_omega(const Sketch& s, std::vector<unsigned char>& out, double& scale) {
        scale = s.oval.empty() ? 0.0 : std::fabs(s.oval.front());
        bool uniform = true;
        for (double v : s.oval) uniform = uniform && std::fabs(v) == scale;
        int64_t prev = 0;
        for (size_t e = 0; e < s.orow.size(); ++e) {
            put_varint(out, zigzag(static_cast<int64_t>(s.orow[e]) - prev));
            prev = s.orow[e];
            put_varint(out, zigzag(s.odest[e]) << 1 | (std::signbit(s.oval[e]) ? 1 : 0));
        }
        if (!uniform) {
            const unsigned char* v = reinterpret_cast<const unsigned char*>(s.oval.data());
            out.insert(out.end(), v, v + s.oval.size() * sizeof(double));
        }
        return uniform;
    }
    // a request's nnz entries of Omega from p, 1-based as the callbacks take
    // them
    static void get_omega(const unsigned char* p, int64_t nnz, bool uniform, double scale, int* row, int* dest,
                          double* val) {
        int64_t prev = 0;
        for (int64_t e = 0; e < nnz; ++e) {
            prev += unzigzag(get_varint(p));
            const uint64_t d = get_varint(p);
            row[e] = static_cast<int>(prev + 1);
            dest[e] = static_cast<int>(unzigzag(d >> 1) + 1);
            val[e] = (d & 1) != 0 ? -scale : scale;
        }
        if (!uniform) std::memcpy(val, p, static_cast<size_t>(nnz) * sizeof(double));
    }

    static uint64_t hash(const int64_t* x, int64_t nx, const int64_t* y, int64_t ny) {
        uint64_t h = 1469598103934665603ULL ^ static_cast<uint64_t>(nx * 1000003 + ny);
        auto mix = [&](int64_t v) {
            h ^= static_cast<uint64_t>(v) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        };
        for (int64_t i = 0; i < nx; ++i) mix(x[i]);
        mix(-1);
        for (int64_t j = 0; j < ny; ++j) mix(y[j]);
        return h;
    }
    static bool same(const Block& b, const int64_t* x, int64_t nx, const int64_t* y, int64_t ny) {
        if (static_cast<int64_t>(b.rows.size()) != nx || static_cast<int64_t>(b.cols.size()) != ny) return false;
        for (int64_t i = 0; i < nx; ++i) {
            if (b.rows[static_cast<size_t>(i)] != x[i]) return false;
        }
        for (int64_t j = 0; j < ny; ++j) {
            if (b.cols[static_cast<size_t>(j)] != y[j]) return false;
        }
        return true;
    }
    const Block* find(const int64_t* x, int64_t nx, const int64_t* y, int64_t ny) const {
        auto range = index_.equal_range(hash(x, nx, y, ny));
        for (auto it = range.first; it != range.second; ++it) {
            const Block& b = blocks_[it->second];
            if (same(b, x, nx, y, ny)) return &b;
        }
        return nullptr;
    }
    // A direct block fetch's lists, 1-based as the callbacks take them.
    struct DirectBlocks {
        int ninter = 0, nallrows = 0, nallcols = 0;
        int64_t nloc = 0;
        std::vector<int> rowidx, colidx, allrows, allcols;
    };
    DirectBlocks pack_blocks(const std::vector<BlockRef>& blocks) {
        DirectBlocks a;
        a.ninter = checked(static_cast<int64_t>(blocks.size()), "blocks");
        a.rowidx.resize(blocks.size());
        a.colidx.resize(blocks.size());
        int64_t nr = 0, nc = 0;
        for (const BlockRef& b : blocks) {
            nr += b.nx;
            nc += b.ny;
        }
        a.allrows.reserve(static_cast<size_t>(nr));
        a.allcols.reserve(static_cast<size_t>(nc));
        for (size_t n = 0; n < blocks.size(); ++n) {
            const BlockRef& b = blocks[n];
            a.rowidx[n] = checked(b.nx, "a block's rows");
            a.colidx[n] = checked(b.ny, "a block's columns");
            for (int64_t i = 0; i < b.nx; ++i) a.allrows.push_back(checked(b.x[i] + 1, "an index"));
            for (int64_t j = 0; j < b.ny; ++j) a.allcols.push_back(checked(b.y[j] + 1, "an index"));
            a.nloc += b.nx * b.ny;
        }
        a.nallrows = checked(nr, "the rows");
        a.nallcols = checked(nc, "the columns");
        return a;
    }
    void count_blocks(const std::vector<BlockRef>& blocks) {
        const size_t o = static_cast<size_t>(entry_origin());
        for (const BlockRef& b : blocks) {
            origin_blocks_[o] += 1.0;
            origin_entries_[o] += static_cast<double>(b.nx) * static_cast<double>(b.ny);
        }
    }
    // A direct fetch's request lists, 1-based as the callbacks take them (and
    // the keys counted in the report).
    struct DirectSketches {
        int nsketch = 0, nallrows = 0, nallcols = 0;
        int64_t nallomega = 0, nloc = 0;
        std::vector<int> rowidx, colidx, dims, omegaidx, allrows, allcols, omega_row, omega_dest;
        std::vector<double> omega_val;
    };
    DirectSketches pack_direct(const std::vector<SketchKey>& keys, bool transposed) {
        DirectSketches a;
        a.nsketch = checked(static_cast<int64_t>(keys.size()), "sketches");
        const size_t ns = keys.size();
        a.rowidx.resize(ns);
        a.colidx.resize(ns);
        a.dims.resize(ns);
        a.omegaidx.resize(ns);
        int64_t nr = 0, nc = 0, ne = 0;
        for (const SketchKey& k : keys) {
            nr += k.nx;
            nc += k.ny;
            ne += k.nnz;
        }
        a.allrows.reserve(static_cast<size_t>(nr));
        a.allcols.reserve(static_cast<size_t>(nc));
        a.omega_row.reserve(static_cast<size_t>(ne));
        a.omega_dest.reserve(static_cast<size_t>(ne));
        a.omega_val.reserve(static_cast<size_t>(ne));
        for (size_t q = 0; q < ns; ++q) {
            const SketchKey& k = keys[q];
            if (k.transposed != transposed) throw std::logic_error("CollectiveEntries: a direct sketch of the other orientation");
            a.rowidx[q] = checked(k.nx, "a sketch's rows");
            a.colidx[q] = checked(k.ny, "a sketch's columns");
            a.dims[q] = k.d;
            a.omegaidx[q] = checked(k.nnz, "a sketch's entries");
            for (int64_t i = 0; i < k.nx; ++i) a.allrows.push_back(checked(k.x[i] + 1, "an index"));
            for (int64_t j = 0; j < k.ny; ++j) a.allcols.push_back(checked(k.y[j] + 1, "an index"));
            for (int64_t e = 0; e < k.nnz; ++e) {
                a.omega_row.push_back(k.orow[e] + 1);
                a.omega_dest.push_back(k.odest[e] + 1);
                a.omega_val.push_back(k.oval[e]);
            }
            a.nloc += static_cast<int64_t>(k.d) * k.ny;
        }
        a.nallrows = checked(nr, "the rows");
        a.nallcols = checked(nc, "the columns");
        a.nallomega = ne;
        return a;
    }
    void count_direct(const std::vector<SketchKey>& keys) {
        const size_t o = static_cast<size_t>(entry_origin());
        for (const SketchKey& k : keys) {
            origin_blocks_[o] += 1.0;
            origin_entries_[o] += static_cast<double>(k.d) * static_cast<double>(k.ny);
            sketch_requests_ += 1.0;
            sketch_values_ += static_cast<double>(k.d) * static_cast<double>(k.ny);
            sketch_covered_ += static_cast<double>(k.nx) * static_cast<double>(k.ny);
        }
    }
    static int checked(int64_t v, const char* what) {
        if (v < 0 || v > std::numeric_limits<int>::max()) {
            throw std::overflow_error(std::string("CollectiveEntries: ") + what + " exceeds INT_MAX");
        }
        return static_cast<int>(v);
    }

    // Every rank's recorded blocks, in one callback call: the list stacked by
    // rank, each block owned (a 1 x 1 process group) by the rank that asked.
    void fetch() {
        const double t_gather = MPI_Wtime();
        int64_t transposed = 0;
        for (const Sketch& s : pending_sketches_) transposed += s.transposed ? 1 : 0;
        int mine[3] = {checked(static_cast<int64_t>(pending_.size()), "blocks"),
                       checked(static_cast<int64_t>(pending_sketches_.size()) - transposed, "sketches"),
                       checked(transposed, "transposed sketches")};
        std::vector<int> all(static_cast<size_t>(3 * size_)), blocks(static_cast<size_t>(size_)),
            sketches(static_cast<size_t>(size_)), sketches_t(static_cast<size_t>(size_));
        MPI_Allgather(mine, 3, MPI_INT, all.data(), 3, MPI_INT, comm_);
        int64_t nblocks = 0, nsketches = 0, nsketches_t = 0;
        for (int p = 0; p < size_; ++p) {
            blocks[static_cast<size_t>(p)] = all[static_cast<size_t>(3 * p)];
            sketches[static_cast<size_t>(p)] = all[static_cast<size_t>(3 * p + 1)];
            sketches_t[static_cast<size_t>(p)] = all[static_cast<size_t>(3 * p + 2)];
            nblocks += blocks[static_cast<size_t>(p)];
            nsketches += sketches[static_cast<size_t>(p)];
            nsketches_t += sketches_t[static_cast<size_t>(p)];
        }
        gather_seconds_ += MPI_Wtime() - t_gather;
        // (every rank skips the same calls: the counts are everyone's)
        if (nblocks > 0) fetch_blocks(blocks);
        if (nsketches > 0) fetch_sketches(sketches, false);
        if (nsketches_t > 0) fetch_sketches(sketches_t, true);
        pending_sketches_.clear();
        pending_sketch_index_.clear();
    }
    void fetch_blocks(const std::vector<int>& counts) {
        const double t_gather = MPI_Wtime();
        // rows, cols per block, then the 1-based rows and the columns, each
        // gathered straight into the callback's list (a rank's lengths from
        // its sizes)
        std::vector<int> my_sizes, my_rows, my_cols;
        for (const Block& b : pending_) {
            my_sizes.push_back(checked(static_cast<int64_t>(b.rows.size()), "a block's rows"));
            my_sizes.push_back(checked(static_cast<int64_t>(b.cols.size()), "a block's columns"));
            for (int64_t r : b.rows) my_rows.push_back(checked(r + 1, "an index"));
            for (int64_t c : b.cols) my_cols.push_back(checked(c + 1, "an index"));
        }
        std::vector<int64_t> nsizes(static_cast<size_t>(size_)), nrows(static_cast<size_t>(size_)),
            ncols(static_cast<size_t>(size_));
        for (int p = 0; p < size_; ++p) nsizes[static_cast<size_t>(p)] = 2 * int64_t{counts[static_cast<size_t>(p)]};
        const std::vector<int> sizes = gatherv(my_sizes, MPI_INT, nsizes);
        int ninter = 0;
        for (int c : counts) ninter += c;
        std::vector<int> rowidx(static_cast<size_t>(ninter)), colidx(static_cast<size_t>(ninter)),
            pgidx(static_cast<size_t>(ninter));
        for (int p = 0, n = 0; p < size_; ++p) {
            for (int b = 0; b < counts[static_cast<size_t>(p)]; ++b, ++n) {
                rowidx[static_cast<size_t>(n)] = sizes[static_cast<size_t>(2 * n)];
                colidx[static_cast<size_t>(n)] = sizes[static_cast<size_t>(2 * n + 1)];
                pgidx[static_cast<size_t>(n)] = p;
                nrows[static_cast<size_t>(p)] += rowidx[static_cast<size_t>(n)];
                ncols[static_cast<size_t>(p)] += colidx[static_cast<size_t>(n)];
            }
        }
        std::vector<int> allrows = gatherv(my_rows, MPI_INT, nrows);
        std::vector<int> allcols = gatherv(my_cols, MPI_INT, ncols);
        std::vector<int> pmaps(static_cast<size_t>(3 * size_));
        for (int p = 0; p < size_; ++p) {
            pmaps[static_cast<size_t>(p)] = 1;
            pmaps[static_cast<size_t>(size_ + p)] = 1;
            pmaps[static_cast<size_t>(2 * size_ + p)] = pids_[static_cast<size_t>(p)];
        }
        int64_t nloc = 0;
        for (const Block& b : pending_) nloc += static_cast<int64_t>(b.rows.size() * b.cols.size());
        std::vector<Block*> mine;
        mine.reserve(pending_.size());
        for (Block& b : pending_) mine.push_back(&b);
        auto len = [](const Block& b) { return static_cast<int64_t>(b.rows.size() * b.cols.size()); };
        struct BlockValuesTag;
        const RoundValues<Block, BlockValuesTag> values(mine, nloc, len);
        int nallrows = checked(static_cast<int64_t>(allrows.size()), "the rows"),
            nallcols = checked(static_cast<int64_t>(allcols.size()), "the columns"), npmap = size_;
        const double t_callback = MPI_Wtime();
        gather_seconds_ += t_callback - t_gather;
        callback_(&ninter, &nallrows, &nallcols, &nloc, allrows.data(), allcols.data(), values.data, rowidx.data(),
                  colidx.data(), pgidx.data(), &npmap, pmaps.data(), quant_);
        const double t_store = MPI_Wtime();
        callback_seconds_ += t_store - t_callback;
        // (the blocks with destinations there, the others stored)
        values.deliver(mine, len);
        int64_t at = 0;
        for (Block& b : pending_) {
            const int64_t n = len(b);
            const DataType* v = reinterpret_cast<const DataType*>(values.data) + at;
            at += n;
            fetched_entries_ += static_cast<double>(n);
            ++fetched_blocks_;
            if (!b.keep) continue;
            b.values.assign(v, v + n);
            b.dst = nullptr;
            b.extra.clear();
            index_.emplace(hash(b.rows.data(), static_cast<int64_t>(b.rows.size()), b.cols.data(),
                                static_cast<int64_t>(b.cols.size())),
                           blocks_.size());
            blocks_.push_back(std::move(b));
        }
        pending_.clear();
        pending_index_.clear();
        store_seconds_ += MPI_Wtime() - t_store;
    }
    // Every rank's recorded sketch requests of one kind (transposed or not),
    // in one call of that kind's callback (as fetch_blocks; Omega's entries
    // 1-based, as the indices).
    void fetch_sketches(const std::vector<int>& counts, bool transposed) {
        const double t_gather = MPI_Wtime();
        SketchCallback callback = transposed ? sketch_t_callback_ : sketch_callback_;
        if (callback == nullptr) throw std::logic_error("CollectiveEntries: a sketch request without its callback");
        std::vector<Sketch*> mine;
        for (Sketch& s : pending_sketches_) {
            if (s.transposed == transposed) mine.push_back(&s);
        }
        // per request rows, cols, d, nnz, the bytes of its Omega, whether its
        // values are +-one scale and that scale's bits; then the 1-based rows,
        // the columns and the Omegas (put_omega, padded to whole words a
        // rank), each rank's lengths from its heads
        constexpr size_t kHead = 7;
        std::vector<int64_t> my_heads;
        std::vector<int> my_rows, my_cols;
        std::vector<unsigned char> my_omega;
        for (const Sketch* sp : mine) {
            const Sketch& s = *sp;
            const size_t start = my_omega.size();
            double scale = 0.0;
            const bool uniform = put_omega(s, my_omega, scale);
            int64_t scale_bits = 0;
            std::memcpy(&scale_bits, &scale, sizeof(scale_bits));
            my_heads.insert(my_heads.end(), {checked(static_cast<int64_t>(s.rows.size()), "a sketch's rows"),
                                             checked(static_cast<int64_t>(s.cols.size()), "a sketch's columns"), s.d,
                                             checked(static_cast<int64_t>(s.orow.size()), "a sketch's entries"),
                                             static_cast<int64_t>(my_omega.size() - start), uniform ? 1 : 0,
                                             scale_bits});
            for (int64_t r : s.rows) my_rows.push_back(checked(r + 1, "an index"));
            for (int64_t c : s.cols) my_cols.push_back(checked(c + 1, "an index"));
        }
        std::vector<uint64_t> my_words((my_omega.size() + sizeof(uint64_t) - 1) / sizeof(uint64_t), 0);
        if (!my_omega.empty()) std::memcpy(my_words.data(), my_omega.data(), my_omega.size());
        std::vector<int64_t> nheads(static_cast<size_t>(size_)), nrows(static_cast<size_t>(size_)),
            ncols(static_cast<size_t>(size_)), nwords(static_cast<size_t>(size_));
        for (int p = 0; p < size_; ++p) nheads[static_cast<size_t>(p)] = kHead * counts[static_cast<size_t>(p)];
        const std::vector<int64_t> heads = gatherv(my_heads, MPI_INT64_T, nheads);
        int nsketch = 0;
        for (int c : counts) nsketch += c;
        std::vector<int> rowidx(static_cast<size_t>(nsketch)), colidx(static_cast<size_t>(nsketch)),
            dims(static_cast<size_t>(nsketch)), omegaidx(static_cast<size_t>(nsketch)), pgidx(static_cast<size_t>(nsketch));
        std::vector<int64_t> entry_at(static_cast<size_t>(nsketch) + 1, 0), byte_at(static_cast<size_t>(nsketch));
        int64_t word_at = 0;  // (a rank's first word)
        for (int p = 0, n = 0; p < size_; ++p) {
            int64_t bytes = 0;
            for (int b = 0; b < counts[static_cast<size_t>(p)]; ++b, ++n) {
                const size_t q = static_cast<size_t>(n);
                const int64_t* h = heads.data() + kHead * q;
                rowidx[q] = static_cast<int>(h[0]);
                colidx[q] = static_cast<int>(h[1]);
                dims[q] = static_cast<int>(h[2]);
                omegaidx[q] = static_cast<int>(h[3]);
                pgidx[q] = p;
                nrows[static_cast<size_t>(p)] += h[0];
                ncols[static_cast<size_t>(p)] += h[1];
                entry_at[q + 1] = entry_at[q] + h[3];
                byte_at[q] = word_at * static_cast<int64_t>(sizeof(uint64_t)) + bytes;
                bytes += h[4];
            }
            nwords[static_cast<size_t>(p)] = (bytes + static_cast<int64_t>(sizeof(uint64_t)) - 1) /
                                             static_cast<int64_t>(sizeof(uint64_t));
            word_at += nwords[static_cast<size_t>(p)];
        }
        std::vector<int> allrows = gatherv(my_rows, MPI_INT, nrows);
        std::vector<int> allcols = gatherv(my_cols, MPI_INT, ncols);
        const std::vector<uint64_t> words = gatherv(my_words, MPI_UINT64_T, nwords);
        const size_t nomega = static_cast<size_t>(entry_at.back());
        std::vector<int> omega_row(nomega), omega_dest(nomega);
        std::vector<double> omega_val(nomega);
        {
            const unsigned char* bytes = reinterpret_cast<const unsigned char*>(words.data());
#pragma omp parallel for schedule(dynamic, 1) if (nsketch > 1 && !omp_in_parallel())
            for (int n = 0; n < nsketch; ++n) {
                const size_t q = static_cast<size_t>(n);
                const int64_t* h = heads.data() + kHead * q;
                double scale = 0.0;
                std::memcpy(&scale, &h[6], sizeof(scale));
                const size_t e = static_cast<size_t>(entry_at[q]);
                get_omega(bytes + byte_at[q], h[3], h[5] != 0, scale, omega_row.data() + e, omega_dest.data() + e,
                          omega_val.data() + e);
            }
        }
        std::vector<int> pmaps(static_cast<size_t>(3 * size_));
        for (int p = 0; p < size_; ++p) {
            pmaps[static_cast<size_t>(p)] = 1;
            pmaps[static_cast<size_t>(size_ + p)] = 1;
            pmaps[static_cast<size_t>(2 * size_ + p)] = pids_[static_cast<size_t>(p)];
        }
        int64_t nloc = 0;
        for (const Sketch* s : mine) nloc += static_cast<int64_t>(s->d) * static_cast<int64_t>(s->cols.size());
        auto len = [](const Sketch& s) { return static_cast<int64_t>(s.d) * static_cast<int64_t>(s.cols.size()); };
        struct SketchValuesTag;
        const RoundValues<Sketch, SketchValuesTag> values(mine, nloc, len);
        int nallrows = checked(static_cast<int64_t>(allrows.size()), "the rows"),
            nallcols = checked(static_cast<int64_t>(allcols.size()), "the columns"), npmap = size_;
        int64_t nallomega = static_cast<int64_t>(nomega);
        const double t_callback = MPI_Wtime();
        gather_seconds_ += t_callback - t_gather;
        callback(&nsketch, &nallrows, &nallcols, &nallomega, &nloc, allrows.data(), allcols.data(),
                         rowidx.data(), colidx.data(), dims.data(), omegaidx.data(), omega_row.data(), omega_dest.data(),
                         omega_val.data(), values.data, pgidx.data(), &npmap, pmaps.data(), quant_);
        const double t_store = MPI_Wtime();
        callback_seconds_ += t_store - t_callback;
        // (the requests with destinations there, the others stored)
        values.deliver(mine, len);
        int64_t at = 0;
        for (Sketch* sp : mine) {
            Sketch& s = *sp;
            const int64_t n = len(s);
            const DataType* v = reinterpret_cast<const DataType*>(values.data) + at;
            at += n;
            fetched_entries_ += static_cast<double>(n);
            sketch_requests_ += 1.0;
            sketch_values_ += static_cast<double>(n);
            sketch_covered_ += static_cast<double>(s.rows.size()) * static_cast<double>(s.cols.size());
            if (!s.keep) continue;
            s.values.assign(v, v + n);
            s.dst = nullptr;
            s.extra.clear();
            const SketchKey k{s.rows.data(), static_cast<int64_t>(s.rows.size()), s.cols.data(),
                              static_cast<int64_t>(s.cols.size()), s.d, static_cast<int64_t>(s.orow.size()),
                              s.orow.data(), s.odest.data(), s.oval.data(), s.transposed};
            sketch_index_.emplace(sketch_hash(k), sketches_.size());
            sketches_.push_back(std::move(s));
        }
        store_seconds_ += MPI_Wtime() - t_store;
    }

    MPI_Comm comm_ = MPI_COMM_NULL;
    int size_ = 1;
    std::vector<int> pids_;
    Callback callback_;
    void* quant_;
    bool symmetric_ = true;
    std::mutex mutex_;
    std::vector<Block> pending_, blocks_;
    SketchCallback sketch_callback_ = nullptr, sketch_t_callback_ = nullptr;
    SketchDeviceCallback sketch_device_callback_ = nullptr;
    BlockDeviceCallback block_device_callback_ = nullptr;
    std::vector<Sketch> pending_sketches_, sketches_;
    std::unordered_multimap<uint64_t, size_t> pending_sketch_index_, sketch_index_;
    double sketch_requests_ = 0.0, sketch_values_ = 0.0, sketch_covered_ = 0.0;
    std::unordered_multimap<uint64_t, size_t> pending_index_, index_;
    int64_t rounds_ = 0, fetch_rounds_ = 0, fetched_blocks_ = 0;
    double fetched_entries_ = 0.0, round_seconds_ = 0.0;
    // (report) the rounds' parts, the device's copies, the blocks by origin
    double wait_seconds_ = 0.0, gather_seconds_ = 0.0, callback_seconds_ = 0.0, store_seconds_ = 0.0,
           delivery_seconds_ = 0.0;
    double origin_blocks_[static_cast<int>(EntryOrigin::count)] = {};
    double origin_entries_[static_cast<int>(EntryOrigin::count)] = {};
    double dup_blocks_[static_cast<int>(EntryOrigin::count)] = {};
    double dup_entries_[static_cast<int>(EntryOrigin::count)] = {};
};

// A phase's items until none is abandoned: run(items) runs the given items
// and returns those abandoned (H2EntryMiss), each of which changed nothing
// outside itself; a round then fetches the blocks they asked for, and they
// run again.  Every rank of the section calls this for the same phases (with
// no items where it has none), since the rounds are collective.  The blocks
// fetched stay until the next phase starts (a fetch-only phase's are what the
// code after it takes).  Without collective entries (entries null), one pass.
template<typename Entries, typename Run>
void run_entry_phase(Entries* entries, std::vector<int64_t> items, Run&& run) {
    if (entries != nullptr) entries->clear();
    for (;;) {
        std::vector<int64_t> abandoned = run(items);
        if (entries == nullptr) {
            if (!abandoned.empty()) throw std::logic_error("run_entry_phase: an item abandoned without collective entries");
            return;
        }
        const int global = entries->round(Entries::kRoundActive | Entries::kRoundBusy);
        if ((global & Entries::kRoundRequest) == 0) {
            if (!abandoned.empty()) throw std::logic_error("run_entry_phase: an item abandoned with nothing to fetch");
            return;
        }
        items = std::move(abandoned);
    }
}

template<typename Run>
void run_entry_phase(std::nullptr_t, std::vector<int64_t> items, Run&& run) {
    if (!run(items).empty()) throw std::logic_error("run_entry_phase: an item abandoned without collective entries");
}

// A kernel's collective entries, or null (also for kernels without them).
template<typename Kernel>
auto kernel_entries(Kernel* kernel, int) -> decltype(kernel->collective.get()) {
    return kernel != nullptr ? kernel->collective.get() : nullptr;
}
template<typename Kernel>
std::nullptr_t kernel_entries(Kernel*, long) {
    return nullptr;
}

// Whether a kernel takes its entries through a collective callback.
template<typename Kernel>
auto kernel_entries_collective(const Kernel* kernel, int) -> decltype(kernel->collective != nullptr) {
    return kernel != nullptr && kernel->collective != nullptr;
}
template<typename Kernel>
bool kernel_entries_collective(const Kernel*, long) {
    return false;
}

// Whether a kernel's entries are transpose-symmetric (H2Kernel::symmetric;
// true for kernels without the member).
template<typename Kernel>
auto kernel_symmetric(const Kernel* kernel, int) -> decltype(static_cast<bool>(kernel->symmetric)) {
    return kernel == nullptr || kernel->symmetric;
}
template<typename Kernel>
bool kernel_symmetric(const Kernel*, long) {
    return true;
}

// Whether a kernel fetches the ID's sketches of kernel rows sketched (entries
// through a collective callback with a sketch callback), and that request
// (H2Kernel::sketch_block_by_index; false for kernels without one).
template<typename Kernel>
auto kernel_sketches(const Kernel* kernel, int) -> decltype(kernel->sketches_kernel_rows()) {
    return kernel != nullptr && kernel->sketches_kernel_rows();
}
template<typename Kernel>
bool kernel_sketches(const Kernel*, long) {
    return false;
}
template<typename Kernel, typename... Args>
auto kernel_sketch_block(const Kernel* kernel, int, Args... args) -> decltype(kernel->sketch_block_by_index(args...)) {
    return kernel->sketch_block_by_index(args...);
}
template<typename Kernel, typename... Args>
bool kernel_sketch_block(const Kernel*, long, Args...) {
    return false;
}

}  // namespace fmm
