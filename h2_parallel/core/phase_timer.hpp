#pragma once
// ---------------------------------------------------------------------------
// Elimination phase timers (diagnostic; removable).
//
// Build with -DFMM_PHASE_TIMERS=0 to compile every probe to nothing. To remove
// entirely: delete this header, its #include in factorization.hpp, and every
// FMM_PHASE_* line (grep FMM_PHASE_).
//
// Two kinds of buckets, all accumulated per thread and summed at report time:
//  - WALL_*  : recorded by the master thread around each per-wave parallel
//              region (so they are rank wall-clock phases).
//  - BOX_*   : inclusive time inside one box's compress+factor, laps in order:
//              sketch -> id -> factor -> near.
//  - SK_*    : sub-buckets of BOX_SKETCH (kernel evaluation, sketch row
//              consumption, lazy fill: temp2 gathers, sparse axpys, GEMMs).
// WALL_BOX_MAXSUM is the per-wave slowest-box time summed over waves: the
// critical-path floor of the per-box phase under the wave schedule.
// Box time by neighbourhood (FMM_PHASE_BOX_SCOPE_NB): each box's time, binned
// by the boxes within two of it (one_hop + two_hop: 124 inside a 3D domain,
// fewer on and near its surface), summed over the level's ranks.
// ---------------------------------------------------------------------------
#ifndef FMM_PHASE_TIMERS
#define FMM_PHASE_TIMERS 1
#endif

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <mutex>
#include <string>
#include <vector>
#include <mpi.h>
#include <sys/resource.h>

namespace fmm {
namespace prof {

enum Phase : int {
    WALL_BOX = 0,
    WALL_BOX_MAXSUM,
    WALL_CANDIDATES,
    WALL_OWNER,
    WALL_MIRROR,
    WALL_FINALIZE,
    WALL_STAGED,
    BOX_SKETCH,
    BOX_ID,
    BOX_FACTOR,
    BOX_NEAR,
    SK_KERNEL,
    SK_CONSUME,
    SK_FILL_GATHER,
    SK_FILL_AXPY,
    SK_FILL_GEMM,
    NEAR_KERNEL,
    NEAR_ASSEMBLE,
    S2_TEMP,
    S2_ANS,
    S2_REST,
    OWN_MATERIALIZE,
    OWN_ALLOC,
    OWN_GEMM,
    OWN_ACCUM,
    OWN_FLUSH,
    NUM_PHASES
};

/// Bins of the box time by neighbourhood: boxes within two of a box (3D at
/// most 124, 2D 24).
constexpr int kNbBins = 128;

struct ThreadRecord {
    double acc[NUM_PHASES];
    double nb_time[kNbBins];     ///< box time by neighbourhood size (s)
    int64_t nb_count[kNbBins];   ///< ... and the boxes
    double box_max;
    int64_t waves, boxes_min, boxes_max;
    int64_t boxes, minflt;   ///< boxes processed; minor page faults inside them
    ThreadRecord();
    // Nested-team threads (split parallelism) are destroyed with their
    // region; their record must leave the registry and keep its totals.
    ~ThreadRecord();
    void clear() {
        std::fill(acc, acc + NUM_PHASES, 0.0);
        std::fill(nb_time, nb_time + kNbBins, 0.0);
        std::fill(nb_count, nb_count + kNbBins, int64_t{0});
        box_max = 0.0;
        waves = 0;
        boxes_min = INT64_MAX;
        boxes_max = 0;
        boxes = 0;
        minflt = 0;
    }
};

inline int64_t thread_minflt() {
    struct rusage ru;
    getrusage(RUSAGE_THREAD, &ru);
    return static_cast<int64_t>(ru.ru_minflt);
}

inline std::mutex& registry_mutex() { static std::mutex m; return m; }
inline std::vector<ThreadRecord*>& registry() { static std::vector<ThreadRecord*> v; return v; }

/// Totals folded in from records of threads that have exited (nested teams).
struct RetiredTotals {
    double acc[NUM_PHASES] = {};
    double nb_time[kNbBins] = {};
    int64_t nb_count[kNbBins] = {};
    int64_t boxes = 0;
    int64_t minflt = 0;
};
inline RetiredTotals& retired_totals() { static RetiredTotals t; return t; }

inline ThreadRecord::ThreadRecord() {
    clear();
    std::lock_guard<std::mutex> lock(registry_mutex());
    registry().push_back(this);
}

inline ThreadRecord::~ThreadRecord() {
    std::lock_guard<std::mutex> lock(registry_mutex());
    auto& v = registry();
    v.erase(std::remove(v.begin(), v.end(), this), v.end());
    auto& R = retired_totals();
    for (int p = 0; p < NUM_PHASES; ++p) R.acc[p] += acc[p];
    for (int b = 0; b < kNbBins; ++b) { R.nb_time[b] += nb_time[b]; R.nb_count[b] += nb_count[b]; }
    R.boxes += boxes;
    R.minflt += minflt;
}

inline ThreadRecord& rec() {
    static thread_local ThreadRecord r;
    return r;
}

inline double now() {
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
}

/// RAII: adds the scope's elapsed time to one bucket of the calling thread.
struct Scope {
    int phase; double t0;
    explicit Scope(int p) : phase(p), t0(now()) {}
    ~Scope() { rec().acc[phase] += now() - t0; }
};

/// RAII: tracks the calling thread's slowest box in the current wave, the
/// box count, and the minor page faults taken inside boxes.
/// With `nb` >= 0 (the boxes within two of this one) it also bins the box's
/// time by that count; -1 leaves the box out of the bins (a box split over
/// several threads, whose time a taskloop shares with other boxes).
struct BoxScope {
    double t0; int64_t f0; int nb;
    explicit BoxScope(int neighbourhood = -1) : t0(now()), f0(thread_minflt()), nb(neighbourhood) {}
    ~BoxScope() {
        auto& r = rec();
        const double dt = now() - t0;
        r.box_max = std::max(r.box_max, dt);
        r.boxes += 1;
        r.minflt += thread_minflt() - f0;
        if (nb >= 0 && nb < kNbBins) { r.nb_time[nb] += dt; r.nb_count[nb] += 1; }
    }
};

/// Sequential laps: each lap(p) charges the time since the previous lap to p.
struct Lap {
    double t;
    Lap() : t(now()) {}
    void operator()(int p) { const double n = now(); rec().acc[p] += n - t; t = n; }
};

inline void note_wave(int64_t boxes) {
    auto& r = rec();
    r.waves += 1;
    r.boxes_min = std::min(r.boxes_min, boxes);
    r.boxes_max = std::max(r.boxes_max, boxes);
}

/// Master thread, after a wave's per-box region: fold the slowest box across
/// threads into WALL_BOX_MAXSUM and reset the per-thread maxima.
inline void wave_box_done() {
    double m = 0.0;
    {
        std::lock_guard<std::mutex> lock(registry_mutex());
        for (ThreadRecord* r : registry()) { m = std::max(m, r->box_max); r->box_max = 0.0; }
    }
    rec().acc[WALL_BOX_MAXSUM] += m;
}

inline void reset() {
    std::lock_guard<std::mutex> lock(registry_mutex());
    for (ThreadRecord* r : registry()) r->clear();
    retired_totals() = RetiredTotals{};
}

/// Collective over comm (every rank calls it; inactive ranks contribute 0).
/// Prints two lines on root: root's own values and the max over active ranks.
inline void report(MPI_Comm comm, int root, bool active, int rank, int level,
                   int threads, bool verbose) {
    constexpr int NV = NUM_PHASES + 5;
    double local[NV] = {0.0};
    if (active) {
        std::lock_guard<std::mutex> lock(registry_mutex());
        int64_t waves = 0, bmin = INT64_MAX, bmax = 0, boxes = 0, minflt = 0;
        for (ThreadRecord* r : registry()) {
            for (int p = 0; p < NUM_PHASES; ++p) local[p] += r->acc[p];
            waves += r->waves;
            bmin = std::min(bmin, r->boxes_min);
            bmax = std::max(bmax, r->boxes_max);
            boxes += r->boxes;
            minflt += r->minflt;
        }
        const auto& R = retired_totals();
        for (int p = 0; p < NUM_PHASES; ++p) local[p] += R.acc[p];
        boxes += R.boxes;
        minflt += R.minflt;
        local[NUM_PHASES] = static_cast<double>(waves);
        local[NUM_PHASES + 1] = (bmin == INT64_MAX) ? 0.0 : -static_cast<double>(bmin);  // MAX(-min)
        local[NUM_PHASES + 2] = static_cast<double>(bmax);
        local[NUM_PHASES + 3] = static_cast<double>(boxes);
        local[NUM_PHASES + 4] = static_cast<double>(minflt);
    }
    double gmax[NV] = {0.0};
    MPI_Reduce(local, gmax, NV, MPI_DOUBLE, MPI_MAX, root, comm);
    // the box time by neighbourhood, summed over the ranks
    double nb_local[2 * kNbBins] = {0.0}, nb_sum[2 * kNbBins] = {0.0};
    if (active) {
        std::lock_guard<std::mutex> lock(registry_mutex());
        for (ThreadRecord* r : registry())
            for (int b = 0; b < kNbBins; ++b) { nb_local[b] += r->nb_time[b]; nb_local[kNbBins + b] += static_cast<double>(r->nb_count[b]); }
        const auto& R = retired_totals();
        for (int b = 0; b < kNbBins; ++b) { nb_local[b] += R.nb_time[b]; nb_local[kNbBins + b] += static_cast<double>(R.nb_count[b]); }
    }
    MPI_Reduce(nb_local, nb_sum, 2 * kNbBins, MPI_DOUBLE, MPI_SUM, root, comm);
    if (!(verbose && rank == root)) return;

    auto line = [&](const char* tag, const double* v) {
        const double th = threads > 0 ? static_cast<double>(threads) : 1.0;
        const double boxes = v[NUM_PHASES + 3] > 0 ? v[NUM_PHASES + 3] : 1.0;
        std::printf(
            "  [phase] level %d %s: waves %.0f (%.0f..%.0f boxes), threads %d, boxes %.0f, "
            "minflt %.3gM (%.0f/box) | wall(s): "
            "box %.1f (slowest-box sum %.1f) cand %.1f owner %.1f mirror %.1f fin %.1f staged %.1f | "
            "per-thread avg(s): sketch %.1f [kernel %.1f consume %.1f gather %.1f axpy %.1f gemm %.1f] "
            "id %.1f factor %.1f near %.1f [kernel %.1f assemble %.1f s2-temp %.1f s2-ans %.1f s2-rest %.1f] "
            "| owner per-thread avg(s): mat %.1f alloc %.1f gemm %.1f accum %.1f flush %.1f\n",
            level, tag, v[NUM_PHASES], -v[NUM_PHASES + 1], v[NUM_PHASES + 2], threads,
            v[NUM_PHASES + 3], v[NUM_PHASES + 4] / 1e6, v[NUM_PHASES + 4] / boxes,
            v[WALL_BOX], v[WALL_BOX_MAXSUM], v[WALL_CANDIDATES], v[WALL_OWNER],
            v[WALL_MIRROR], v[WALL_FINALIZE], v[WALL_STAGED],
            v[BOX_SKETCH] / th, v[SK_KERNEL] / th, v[SK_CONSUME] / th,
            v[SK_FILL_GATHER] / th, v[SK_FILL_AXPY] / th, v[SK_FILL_GEMM] / th,
            v[BOX_ID] / th, v[BOX_FACTOR] / th, v[BOX_NEAR] / th,
            v[NEAR_KERNEL] / th, v[NEAR_ASSEMBLE] / th,
            v[S2_TEMP] / th, v[S2_ANS] / th, v[S2_REST] / th,
            v[OWN_MATERIALIZE] / th, v[OWN_ALLOC] / th, v[OWN_GEMM] / th,
            v[OWN_ACCUM] / th, v[OWN_FLUSH] / th);
    };
    line("rank-root", local);
    line("max-ranks", gmax);
    {
        // boxes within two of a box: count of boxes, ms each (all ranks, boxes
        // eliminated one per thread)
        std::string bins;
        for (int b = kNbBins - 1; b >= 0; --b) {
            if (nb_sum[kNbBins + b] <= 0) continue;
            char buf[64];
            std::snprintf(buf, sizeof(buf), " %d: %.0f, %.2f ms |", b, nb_sum[kNbBins + b], 1e3 * nb_sum[b] / nb_sum[kNbBins + b]);
            bins += buf;
        }
        if (!bins.empty())
            std::printf("  [phase] level %d box time by neighbourhood (boxes within two: boxes, time each; all ranks):%s\n",
                        level, bins.c_str());
    }
    std::fflush(stdout);
}

}  // namespace prof
}  // namespace fmm

#define FMM_PHASE_CAT2(a, b) a##b
#define FMM_PHASE_CAT(a, b) FMM_PHASE_CAT2(a, b)

#if FMM_PHASE_TIMERS
#define FMM_PHASE_SCOPE(P) ::fmm::prof::Scope FMM_PHASE_CAT(fmm_phase_scope_, __LINE__)(::fmm::prof::P)
#define FMM_PHASE_BOX_SCOPE() ::fmm::prof::BoxScope FMM_PHASE_CAT(fmm_phase_box_, __LINE__)
#define FMM_PHASE_BOX_SCOPE_NB(n) ::fmm::prof::BoxScope FMM_PHASE_CAT(fmm_phase_box_, __LINE__)(n)
#define FMM_PHASE_LAP_BEGIN(name) ::fmm::prof::Lap name
#define FMM_PHASE_LAP(name, P) name(::fmm::prof::P)
#define FMM_PHASE_NOTE_WAVE(n) ::fmm::prof::note_wave(static_cast<int64_t>(n))
#define FMM_PHASE_WAVE_BOX_DONE() ::fmm::prof::wave_box_done()
#define FMM_PHASE_RESET() ::fmm::prof::reset()
#define FMM_PHASE_REPORT(comm, root, active, rank, level, threads, verbose) \
    ::fmm::prof::report(comm, root, active, rank, level, threads, verbose)
#else
#define FMM_PHASE_SCOPE(P) ((void)0)
#define FMM_PHASE_BOX_SCOPE() ((void)0)
#define FMM_PHASE_BOX_SCOPE_NB(n) ((void)0)
#define FMM_PHASE_LAP_BEGIN(name) ((void)0)
#define FMM_PHASE_LAP(name, P) ((void)0)
#define FMM_PHASE_NOTE_WAVE(n) ((void)0)
#define FMM_PHASE_WAVE_BOX_DONE() ((void)0)
#define FMM_PHASE_RESET() ((void)0)
#define FMM_PHASE_REPORT(comm, root, active, rank, level, threads, verbose) ((void)0)
#endif
