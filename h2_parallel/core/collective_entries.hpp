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

#include <mpi.h>

#include <cstdint>
#include <cstdio>
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

template<typename DataType, typename CallbackValue>
class CollectiveEntries {
public:
    using Callback = void (*)(int*, int*, int*, int64_t*, int*, int*, CallbackValue*, int*, int*, int*, int*, int*,
                              void*);
    static constexpr int kRoundRequest = 1;  // (a rank has recorded blocks)
    static constexpr int kRoundActive = 2;   // (a section is running)

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
        pending_index_.emplace(h, pending_.size());
        pending_.push_back(std::move(b));
    }

    // Collective: one round.  flags | kRoundRequest where this rank has
    // recorded blocks, all ranks' maximum returned; with kRoundRequest set
    // every rank's blocks are fetched (one callback call over all ranks).
    int round(int flags) {
        const double t0 = MPI_Wtime();
        int local = flags | (pending_.empty() ? 0 : kRoundRequest), global = 0;
        MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, comm_);
        if ((global & kRoundRequest) != 0) {
            fetch();
            ++fetch_rounds_;
        }
        ++rounds_;
        round_seconds_ += MPI_Wtime() - t0;
        return global;
    }
    // A rank outside a running section: rounds until the section's end.
    void serve() {
        while ((round(0) & kRoundActive) != 0) {}
    }
    // The ranks of a section, at its end (all of them, together).
    void end_section() {
        if (round(0) != 0) throw std::logic_error("CollectiveEntries::end_section: a section still running");
        clear();
    }
    // The fetched blocks are no longer needed (a phase's end).
    void clear() {
        blocks_.clear();
        index_.clear();
    }
    int64_t rounds() const { return rounds_; }
    int64_t fetch_rounds() const { return fetch_rounds_; }
    double round_seconds() const { return round_seconds_; }
    // One line on rank 0 of comm: the rounds, the blocks and entries fetched
    // (all ranks), the slowest rank's time in rounds.
    void report(const char* what) {
        double v[2] = {fetched_entries_, static_cast<double>(fetched_blocks_)}, t = round_seconds_;
        MPI_Allreduce(MPI_IN_PLACE, v, 2, MPI_DOUBLE, MPI_SUM, comm_);
        MPI_Allreduce(MPI_IN_PLACE, &t, 1, MPI_DOUBLE, MPI_MAX, comm_);
        int rank = 0;
        MPI_Comm_rank(comm_, &rank);
        if (rank == 0) {
            std::printf("H2 entries through the collective callback (%s): %lld rounds, %lld fetching, %.0f blocks, "
                        "%.3g entries, %.2f s in rounds (slowest rank)\n",
                        what, static_cast<long long>(rounds_), static_cast<long long>(fetch_rounds_), v[1], v[0], t);
            std::fflush(stdout);
        }
    }
    void reset_counts() {
        rounds_ = fetch_rounds_ = fetched_blocks_ = 0;
        fetched_entries_ = round_seconds_ = 0.0;
    }
    int64_t fetched_blocks() const { return fetched_blocks_; }
    double fetched_entries() const { return fetched_entries_; }

private:
    struct Block {
        std::vector<int64_t> rows, cols;
        std::vector<DataType> values;
    };

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
        std::vector<int> mine{checked(static_cast<int64_t>(pending_.size()), "blocks")}, counts(static_cast<size_t>(size_));
        MPI_Allgather(mine.data(), 1, MPI_INT, counts.data(), 1, MPI_INT, comm_);
        std::vector<int> my_sizes;  // rows, cols per block
        std::vector<int> my_ids;    // 1-based rows then columns per block
        for (const Block& b : pending_) {
            my_sizes.push_back(checked(static_cast<int64_t>(b.rows.size()), "a block's rows"));
            my_sizes.push_back(checked(static_cast<int64_t>(b.cols.size()), "a block's columns"));
            for (int64_t r : b.rows) my_ids.push_back(checked(r + 1, "an index"));
            for (int64_t c : b.cols) my_ids.push_back(checked(c + 1, "an index"));
        }
        auto gatherv = [&](const std::vector<int>& local) {
            int n = checked(static_cast<int64_t>(local.size()), "a request list");
            std::vector<int> ns(static_cast<size_t>(size_)), displs(static_cast<size_t>(size_));
            MPI_Allgather(&n, 1, MPI_INT, ns.data(), 1, MPI_INT, comm_);
            int64_t total = 0;
            for (int p = 0; p < size_; ++p) {
                displs[static_cast<size_t>(p)] = checked(total, "the request lists");
                total += ns[static_cast<size_t>(p)];
            }
            std::vector<int> all(static_cast<size_t>(total));
            MPI_Allgatherv(local.data(), n, MPI_INT, all.data(), ns.data(), displs.data(), MPI_INT, comm_);
            return all;
        };
        const std::vector<int> sizes = gatherv(my_sizes);
        const std::vector<int> ids = gatherv(my_ids);
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
        callback_(&ninter, &nallrows, &nallcols, &nloc, allrows.data(), allcols.data(), values.data(), rowidx.data(),
                  colidx.data(), pgidx.data(), &npmap, pmaps.data(), quant_);
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
    }

    MPI_Comm comm_ = MPI_COMM_NULL;
    int size_ = 1;
    std::vector<int> pids_;
    Callback callback_;
    void* quant_;
    std::mutex mutex_;
    std::vector<Block> pending_, blocks_;
    std::unordered_multimap<uint64_t, size_t> pending_index_, index_;
    int64_t rounds_ = 0, fetch_rounds_ = 0, fetched_blocks_ = 0;
    double fetched_entries_ = 0.0, round_seconds_ = 0.0;
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
        const int global = entries->round(Entries::kRoundActive);
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

}  // namespace fmm
