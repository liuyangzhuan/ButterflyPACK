#pragma once
// ---------------------------------------------------------------------------
// OpenMP parallel region under a condition, without a one-thread region when
// the condition is false.
//
// `#pragma omp parallel if (false)` still opens a region, an inactive one, and
// omp_in_parallel() is false inside it.  An OpenMP-threaded BLAS called there
// (OpenBLAS's OpenMP build, MKL's OpenMP layer) then opens its own team, which
// is nested, and libgomp starts fresh threads for every nested team: ~0.2 s for
// one 287x3624 matrix-vector product, so the H2 matvec and solve of a 64K-point
// block ran 25x slower (the single-box colour waves of its top level).  Without
// the region, that BLAS team and any parallel loop inside use the thread pool.
// ---------------------------------------------------------------------------

#include <omp.h>

namespace fmm {

// `#pragma omp parallel default(shared) if (parallel)` around region(), except
// that with parallel false outside any parallel region region() runs on the
// calling thread with no construct at all.  region() may hold orphaned
// worksharing constructs (omp for, single, barrier): they bind to the region
// when there is one and run on the calling thread alone otherwise, as on a
// one-thread team.  Inside an enclosing active region the construct is kept,
// so they still bind to it and not to the enclosing team.
template<typename Region>
inline void omp_parallel_if(bool parallel, Region&& region) {
    if (parallel || omp_in_parallel()) {
        #pragma omp parallel default(shared) if (parallel)
        region();
    } else {
        region();
    }
}

} // namespace fmm
