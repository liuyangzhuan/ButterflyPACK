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

}  // namespace fmm
