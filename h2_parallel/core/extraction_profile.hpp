#pragma once
// The batched extraction from a compressed H2 (unstructured/extraction.hpp)
// by step: this process's seconds and work, summed over its calls since the
// last reset.  An H2 that takes its entries through a collective callback
// prints them with its rounds (CollectiveEntries::report, verbosity >= 1)
// when its callback ran such an extraction (the bridge of the example).

#include <cstdint>

namespace fmm {

struct ExtractionProfile {
    enum Step { kResolve, kTasks, kUp, kColumns, kCompute, kPack, kRoute, kUnpack, kSteps };
    static constexpr const char* kNames[kSteps] = {"resolve", "tasks",   "vectors up", "column side",
                                                   "compute", "pack",    "route",      "unpack"};
    double seconds[kSteps] = {};
    double route_mpi = 0.0;         // (of route: the all-to-alls themselves)
    double entries_computed = 0.0;  // the tasks' entries this rank computed
    double bytes_routed = 0.0;      // the values this rank sent to the blocks' owners
    int64_t calls = 0;
};

inline ExtractionProfile& extraction_profile() {
    static ExtractionProfile profile;
    return profile;
}

// What the resets (CollectiveEntries::reset_counts, at every factorization)
// took out of extraction_profile(): with it, the process's totals since the
// start (c_bpack_get_extraction_profile).
inline ExtractionProfile& extraction_profile_reset() {
    static ExtractionProfile profile;
    return profile;
}

// extraction_profile() zeroed, its values kept in extraction_profile_reset().
inline void reset_extraction_profile() {
    ExtractionProfile& p = extraction_profile();
    ExtractionProfile& r = extraction_profile_reset();
    for (int s = 0; s < ExtractionProfile::kSteps; ++s) r.seconds[s] += p.seconds[s];
    r.route_mpi += p.route_mpi;
    r.entries_computed += p.entries_computed;
    r.bytes_routed += p.bytes_routed;
    r.calls += p.calls;
    p = ExtractionProfile{};
}

}  // namespace fmm
