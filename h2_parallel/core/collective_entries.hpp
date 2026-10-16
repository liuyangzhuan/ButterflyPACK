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
// asker gives: the ID's sketch of its kernel rows asks for those instead of
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

#include <mpi.h>

#include <algorithm>
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
// Outside an item such a request throws at once.
struct EntryItemState {
    bool active = false;
    bool doomed = false;
};
inline EntryItemState& entry_item() {
    thread_local EntryItemState state;
    return state;
}
inline void entry_item_begin() { entry_item() = EntryItemState{true, false}; }
inline void entry_checkpoint() {
    if (entry_item().doomed) {
        entry_item() = EntryItemState{};
        throw H2EntryMiss();
    }
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
    static constexpr int kRoundRequest = 1;  // (a rank has recorded blocks)
    static constexpr int kRoundActive = 2;   // (a section is running)
    static constexpr int kRoundBusy = 4;     // (a rank may still ask, in the device path's sense)

    CollectiveEntries(MPI_Comm comm, Callback callback, void* quant, int callback_pid)
        : callback_(callback), quant_(quant) {
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

    // A[i + j lda] = K(x_i, y_j) if this block (or its transpose: the matrix
    // is symmetric) was fetched; else false.  Thread safe.
    bool lookup(const int64_t* x, int64_t nx, const int64_t* y, int64_t ny, DataType* A, int64_t lda) const {
        if (const Block* b = find(x, nx, y, ny)) {
            for (int64_t j = 0; j < ny; ++j) {
                for (int64_t i = 0; i < nx; ++i) A[i + j * lda] = b->values[static_cast<size_t>(i + j * nx)];
            }
            return true;
        }
        if (const Block* b = find(y, ny, x, nx)) {
            for (int64_t j = 0; j < ny; ++j) {
                for (int64_t i = 0; i < nx; ++i) A[i + j * lda] = b->values[static_cast<size_t>(j + i * ny)];
            }
            return true;
        }
        return false;
    }

    // Drops the block (or its transpose) once read: the device path reads
    // each block it fetched once.  Not concurrent with lookups.
    void release(const int64_t* x, int64_t nx, const int64_t* y, int64_t ny) {
        for (int t = 0; t < 2; ++t) {
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

    // A sketch callback: sketched blocks are fetched through it.
    void set_sketch(SketchCallback callback) { sketch_callback_ = callback; }
    bool has_sketch() const { return sketch_callback_ != nullptr; }

    // A sketch request: Y = Omega K(x, y), Omega (d x nx) given by its nnz
    // entries (orow: the position in x, odest: the row of Y, oval).
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
    void record_sketch(const SketchKey& k) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (find_sketch(k, pending_sketch_index_, pending_sketches_) != nullptr) return;
        Sketch s;
        s.rows.assign(k.x, k.x + k.nx);
        s.cols.assign(k.y, k.y + k.ny);
        s.d = k.d;
        s.orow.assign(k.orow, k.orow + k.nnz);
        s.odest.assign(k.odest, k.odest + k.nnz);
        s.oval.assign(k.oval, k.oval + k.nnz);
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

    // Records the block for the next round.  Thread safe.
    void record(const int64_t* x, int64_t nx, const int64_t* y, int64_t ny) {
        std::lock_guard<std::mutex> lock(mutex_);
        const uint64_t h = hash(x, nx, y, ny);
        auto range = pending_index_.equal_range(h);
        for (auto it = range.first; it != range.second; ++it) {
            const Block& b = pending_[it->second];
            if (same(b, x, nx, y, ny)) return;
        }
        Block b;
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
        double origins[2 * kOrigins];
        std::copy(origin_blocks_, origin_blocks_ + kOrigins, origins);
        std::copy(origin_entries_, origin_entries_ + kOrigins, origins + kOrigins);
        MPI_Allreduce(MPI_IN_PLACE, origins, 2 * kOrigins, MPI_DOUBLE, MPI_SUM, comm_);
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
        extraction_profile() = ExtractionProfile{};
    }
    int64_t fetched_blocks() const { return fetched_blocks_; }
    double fetched_entries() const { return fetched_entries_; }

private:
    struct Block {
        std::vector<int64_t> rows, cols;
        std::vector<DataType> values;
    };
    struct Sketch {
        std::vector<int64_t> rows, cols;
        int d = 0;
        std::vector<int> orow, odest;
        std::vector<double> oval;
        std::vector<DataType> values;
    };
    static uint64_t sketch_hash(const SketchKey& k) {
        uint64_t h = hash(k.x, k.nx, k.y, k.ny) ^ (static_cast<uint64_t>(k.d) * 0x9e3779b97f4a7c15ULL);
        for (int64_t e = 0; e < k.nnz; ++e) {
            uint64_t bits = 0;
            std::memcpy(&bits, &k.oval[e], sizeof(bits));
            h ^= (static_cast<uint64_t>(k.orow[e]) << 20 ^ static_cast<uint64_t>(k.odest[e]) ^ bits) +
                 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        }
        return h;
    }
    static bool same_sketch(const Sketch& s, const SketchKey& k) {
        if (s.d != k.d || static_cast<int64_t>(s.rows.size()) != k.nx || static_cast<int64_t>(s.cols.size()) != k.ny ||
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
    // Every rank's list, stacked by rank (the same on every rank).
    template<typename T>
    std::vector<T> gatherv(const std::vector<T>& local, MPI_Datatype type) {
        int n = checked(static_cast<int64_t>(local.size()), "a request list");
        std::vector<int> ns(static_cast<size_t>(size_)), displs(static_cast<size_t>(size_));
        MPI_Allgather(&n, 1, MPI_INT, ns.data(), 1, MPI_INT, comm_);
        int64_t total = 0;
        for (int p = 0; p < size_; ++p) {
            displs[static_cast<size_t>(p)] = checked(total, "the request lists");
            total += ns[static_cast<size_t>(p)];
        }
        std::vector<T> all(static_cast<size_t>(total));
        MPI_Allgatherv(local.data(), n, type, all.data(), ns.data(), displs.data(), type, comm_);
        return all;
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
        int mine[2] = {checked(static_cast<int64_t>(pending_.size()), "blocks"),
                       checked(static_cast<int64_t>(pending_sketches_.size()), "sketches")};
        std::vector<int> both(static_cast<size_t>(2 * size_)), blocks(static_cast<size_t>(size_)),
            sketches(static_cast<size_t>(size_));
        MPI_Allgather(mine, 2, MPI_INT, both.data(), 2, MPI_INT, comm_);
        int64_t nblocks = 0, nsketches = 0;
        for (int p = 0; p < size_; ++p) {
            blocks[static_cast<size_t>(p)] = both[static_cast<size_t>(2 * p)];
            sketches[static_cast<size_t>(p)] = both[static_cast<size_t>(2 * p + 1)];
            nblocks += blocks[static_cast<size_t>(p)];
            nsketches += sketches[static_cast<size_t>(p)];
        }
        gather_seconds_ += MPI_Wtime() - t_gather;
        // (every rank skips the same calls: the counts are everyone's)
        if (nblocks > 0) fetch_blocks(blocks);
        if (nsketches > 0) fetch_sketches(sketches);
    }
    void fetch_blocks(const std::vector<int>& counts) {
        const double t_gather = MPI_Wtime();
        std::vector<int> my_sizes;  // rows, cols per block
        std::vector<int> my_ids;    // 1-based rows then columns per block
        for (const Block& b : pending_) {
            my_sizes.push_back(checked(static_cast<int64_t>(b.rows.size()), "a block's rows"));
            my_sizes.push_back(checked(static_cast<int64_t>(b.cols.size()), "a block's columns"));
            for (int64_t r : b.rows) my_ids.push_back(checked(r + 1, "an index"));
            for (int64_t c : b.cols) my_ids.push_back(checked(c + 1, "an index"));
        }
        const std::vector<int> sizes = gatherv(my_sizes, MPI_INT);
        const std::vector<int> ids = gatherv(my_ids, MPI_INT);
        int ninter = 0;
        for (int c : counts) ninter += c;
        std::vector<int> rowidx(static_cast<size_t>(ninter)), colidx(static_cast<size_t>(ninter)),
            pgidx(static_cast<size_t>(ninter)), allrows, allcols;
        {
            int n = 0;
            size_t at = 0;
            for (int p = 0; p < size_; ++p) {
                for (int b = 0; b < counts[static_cast<size_t>(p)]; ++b, ++n) {
                    rowidx[static_cast<size_t>(n)] = sizes[static_cast<size_t>(2 * n)];
                    colidx[static_cast<size_t>(n)] = sizes[static_cast<size_t>(2 * n + 1)];
                    pgidx[static_cast<size_t>(n)] = p;
                    allrows.insert(allrows.end(), ids.begin() + static_cast<std::ptrdiff_t>(at),
                                   ids.begin() + static_cast<std::ptrdiff_t>(at + rowidx[static_cast<size_t>(n)]));
                    at += static_cast<size_t>(rowidx[static_cast<size_t>(n)]);
                    allcols.insert(allcols.end(), ids.begin() + static_cast<std::ptrdiff_t>(at),
                                   ids.begin() + static_cast<std::ptrdiff_t>(at + colidx[static_cast<size_t>(n)]));
                    at += static_cast<size_t>(colidx[static_cast<size_t>(n)]);
                }
            }
        }
        std::vector<int> pmaps(static_cast<size_t>(3 * size_));
        for (int p = 0; p < size_; ++p) {
            pmaps[static_cast<size_t>(p)] = 1;
            pmaps[static_cast<size_t>(size_ + p)] = 1;
            pmaps[static_cast<size_t>(2 * size_ + p)] = pids_[static_cast<size_t>(p)];
        }
        int64_t nloc = 0;
        for (const Block& b : pending_) nloc += static_cast<int64_t>(b.rows.size() * b.cols.size());
        std::vector<CallbackValue> values(static_cast<size_t>(std::max<int64_t>(nloc, 1)));
        int nallrows = checked(static_cast<int64_t>(allrows.size()), "the rows"),
            nallcols = checked(static_cast<int64_t>(allcols.size()), "the columns"), npmap = size_;
        const double t_callback = MPI_Wtime();
        gather_seconds_ += t_callback - t_gather;
        callback_(&ninter, &nallrows, &nallcols, &nloc, allrows.data(), allcols.data(), values.data(), rowidx.data(),
                  colidx.data(), pgidx.data(), &npmap, pmaps.data(), quant_);
        const double t_store = MPI_Wtime();
        callback_seconds_ += t_store - t_callback;
        int64_t at = 0;
        for (Block& b : pending_) {
            const int64_t n = static_cast<int64_t>(b.rows.size() * b.cols.size());
            const DataType* v = reinterpret_cast<const DataType*>(values.data()) + at;
            b.values.assign(v, v + n);
            at += n;
            fetched_entries_ += static_cast<double>(n);
            ++fetched_blocks_;
            index_.emplace(hash(b.rows.data(), static_cast<int64_t>(b.rows.size()), b.cols.data(),
                                static_cast<int64_t>(b.cols.size())),
                           blocks_.size());
            blocks_.push_back(std::move(b));
        }
        pending_.clear();
        pending_index_.clear();
        store_seconds_ += MPI_Wtime() - t_store;
    }
    // Every rank's recorded sketch requests, in one sketch callback call (as
    // fetch_blocks; Omega's entries 1-based, as the indices).
    void fetch_sketches(const std::vector<int>& counts) {
        const double t_gather = MPI_Wtime();
        std::vector<int> my_sizes;  // rows, cols, d, nnz per request
        std::vector<int> my_ids;    // 1-based rows then columns per request
        std::vector<int> my_omega;  // 1-based (row, destination) per entry
        std::vector<double> my_values;
        for (const Sketch& s : pending_sketches_) {
            my_sizes.insert(my_sizes.end(), {checked(static_cast<int64_t>(s.rows.size()), "a sketch's rows"),
                                             checked(static_cast<int64_t>(s.cols.size()), "a sketch's columns"), s.d,
                                             checked(static_cast<int64_t>(s.orow.size()), "a sketch's entries")});
            for (int64_t r : s.rows) my_ids.push_back(checked(r + 1, "an index"));
            for (int64_t c : s.cols) my_ids.push_back(checked(c + 1, "an index"));
            for (size_t e = 0; e < s.orow.size(); ++e) my_omega.insert(my_omega.end(), {s.orow[e] + 1, s.odest[e] + 1});
            my_values.insert(my_values.end(), s.oval.begin(), s.oval.end());
        }
        const std::vector<int> sizes = gatherv(my_sizes, MPI_INT);
        const std::vector<int> ids = gatherv(my_ids, MPI_INT);
        const std::vector<int> omega = gatherv(my_omega, MPI_INT);
        const std::vector<double> omega_val = gatherv(my_values, MPI_DOUBLE);
        int nsketch = 0;
        for (int c : counts) nsketch += c;
        std::vector<int> rowidx(static_cast<size_t>(nsketch)), colidx(static_cast<size_t>(nsketch)),
            dims(static_cast<size_t>(nsketch)), omegaidx(static_cast<size_t>(nsketch)), pgidx(static_cast<size_t>(nsketch)),
            allrows, allcols, omega_row(omega_val.size()), omega_dest(omega_val.size());
        {
            int n = 0;
            size_t at = 0;
            for (int p = 0; p < size_; ++p) {
                for (int b = 0; b < counts[static_cast<size_t>(p)]; ++b, ++n) {
                    const size_t q = static_cast<size_t>(n);
                    rowidx[q] = sizes[4 * q];
                    colidx[q] = sizes[4 * q + 1];
                    dims[q] = sizes[4 * q + 2];
                    omegaidx[q] = sizes[4 * q + 3];
                    pgidx[q] = p;
                    allrows.insert(allrows.end(), ids.begin() + static_cast<std::ptrdiff_t>(at),
                                   ids.begin() + static_cast<std::ptrdiff_t>(at + rowidx[q]));
                    at += static_cast<size_t>(rowidx[q]);
                    allcols.insert(allcols.end(), ids.begin() + static_cast<std::ptrdiff_t>(at),
                                   ids.begin() + static_cast<std::ptrdiff_t>(at + colidx[q]));
                    at += static_cast<size_t>(colidx[q]);
                }
            }
            for (size_t e = 0; e < omega_val.size(); ++e) {
                omega_row[e] = omega[2 * e];
                omega_dest[e] = omega[2 * e + 1];
            }
        }
        std::vector<int> pmaps(static_cast<size_t>(3 * size_));
        for (int p = 0; p < size_; ++p) {
            pmaps[static_cast<size_t>(p)] = 1;
            pmaps[static_cast<size_t>(size_ + p)] = 1;
            pmaps[static_cast<size_t>(2 * size_ + p)] = pids_[static_cast<size_t>(p)];
        }
        int64_t nloc = 0;
        for (const Sketch& s : pending_sketches_) nloc += static_cast<int64_t>(s.d) * static_cast<int64_t>(s.cols.size());
        std::vector<CallbackValue> values(static_cast<size_t>(std::max<int64_t>(nloc, 1)));
        int nallrows = checked(static_cast<int64_t>(allrows.size()), "the rows"),
            nallcols = checked(static_cast<int64_t>(allcols.size()), "the columns"), npmap = size_;
        int64_t nallomega = static_cast<int64_t>(omega_val.size());
        const double t_callback = MPI_Wtime();
        gather_seconds_ += t_callback - t_gather;
        sketch_callback_(&nsketch, &nallrows, &nallcols, &nallomega, &nloc, allrows.data(), allcols.data(),
                         rowidx.data(), colidx.data(), dims.data(), omegaidx.data(), omega_row.data(), omega_dest.data(),
                         const_cast<double*>(omega_val.data()), values.data(), pgidx.data(), &npmap, pmaps.data(),
                         quant_);
        const double t_store = MPI_Wtime();
        callback_seconds_ += t_store - t_callback;
        int64_t at = 0;
        for (Sketch& s : pending_sketches_) {
            const int64_t n = static_cast<int64_t>(s.d) * static_cast<int64_t>(s.cols.size());
            const DataType* v = reinterpret_cast<const DataType*>(values.data()) + at;
            s.values.assign(v, v + n);
            at += n;
            fetched_entries_ += static_cast<double>(n);
            sketch_requests_ += 1.0;
            sketch_values_ += static_cast<double>(n);
            sketch_covered_ += static_cast<double>(s.rows.size()) * static_cast<double>(s.cols.size());
            const SketchKey k{s.rows.data(), static_cast<int64_t>(s.rows.size()), s.cols.data(),
                              static_cast<int64_t>(s.cols.size()), s.d, static_cast<int64_t>(s.orow.size()),
                              s.orow.data(), s.odest.data(), s.oval.data()};
            sketch_index_.emplace(sketch_hash(k), sketches_.size());
            sketches_.push_back(std::move(s));
        }
        pending_sketches_.clear();
        pending_sketch_index_.clear();
        store_seconds_ += MPI_Wtime() - t_store;
    }

    MPI_Comm comm_ = MPI_COMM_NULL;
    int size_ = 1;
    std::vector<int> pids_;
    Callback callback_;
    void* quant_;
    std::mutex mutex_;
    std::vector<Block> pending_, blocks_;
    SketchCallback sketch_callback_ = nullptr;
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
