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

#include <algorithm>
#include <cstdint>
#include <cstring>

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

// Column pieces of H2_GEMM_PIECE columns as OpenMP tasks, for the threads a
// small wave leaves idle: for a function that is one gemm call on its
// columns, the same entries bit for bit as the call cut into them (OpenBLAS's
// zgemm packs op(B) in blocks of 6 columns and computes a call's last block
// differently, so a call cut at multiples of 6 from its first column gives the
// same bits, dgemm at multiples of 4: h2_c1_work/topcpu/bench
// zgemm_probe*.cpp, 2026-11-01).
constexpr int64_t H2_GEMM_PIECE = 48;
// one gemm call's columns [0, columns) in pieces of `width` columns (a
// multiple of 6)
template <typename Function>
inline void gemm_column_pieces(int64_t columns, const Function& function, int64_t width = H2_GEMM_PIECE) {
    const int64_t count = (columns + width - 1) / width;
    #pragma omp taskloop grainsize(1) default(shared)
    for (int64_t i = 0; i < count; ++i) {
        function(i * width, std::min(columns, (i + 1) * width));
    }
}

// memcpy of a large block over the threads: 8 MB pieces as tasks inside an
// active region (threads waiting at its barriers take them), else in a region
// of its own; a block under 32 MB, or one thread, one memcpy.  (The top
// levels' generators are 0.1-1 GB: one thread copied and first-touched them
// on the path from one wave's box to the next one's.)
inline void copy_bytes_parallel(void* dst, const void* src, size_t bytes) {
    constexpr size_t piece = size_t{8} << 20;
    if (bytes < 4 * piece || omp_get_max_threads() == 1) {
        std::memcpy(dst, src, bytes);
        return;
    }
    char* d = static_cast<char*>(dst);
    const char* s = static_cast<const char*>(src);
    const int64_t count = static_cast<int64_t>((bytes + piece - 1) / piece);
    auto copy_piece = [&](int64_t i) {
        const size_t offset = static_cast<size_t>(i) * piece;
        std::memcpy(d + offset, s + offset, std::min(piece, bytes - offset));
    };
    if (omp_in_parallel()) {
        #pragma omp taskloop grainsize(1) default(shared)
        for (int64_t i = 0; i < count; ++i) copy_piece(i);
    } else {
        #pragma omp parallel for schedule(dynamic) default(shared)
        for (int64_t i = 0; i < count; ++i) copy_piece(i);
    }
}

} // namespace fmm
