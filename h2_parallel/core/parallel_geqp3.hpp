#pragma once
// ---------------------------------------------------------------------------
// QR with column pivoting (LAPACK's geqp3) with the large BLAS calls of its
// blocked steps split over the threads a small wave leaves idle, bitwise
// identical to the reference.  A box's ID factors its whole sketch (about
// n + 4 rows by n columns, n ~1700 at the PPW-10 Helmholtz top levels) with
// geqp3 twice under H2_ID_proxy 2, on the box's thread alone: ~1 s each, 97%
// of it the per-column gemv over the trailing matrix and the panels' trailing
// gemm, while a one-box wave leaves the process's other threads idle.
//
// The port follows ?geqp3/?laqps of the bundled reference LAPACK (OpenBLAS
// 0.3.30) call for call, with the same BLAS routines on the same operands; a
// call is split only where its pieces give the same entries bit for bit
// (measured on this OpenBLAS, h2_c1_work/topcpu/bench/split_test.cpp,
// zgemm_probe.cpp): zgemv 'C' and the one-row zgemm by columns anywhere;
// gemv 'N' by rows, dgemv 'T' and dgemm by columns at multiples of 4 from
// the call's first entry; zgemm by columns at multiples of 6.  The column
// norms are recomputed by column (independent ?nrm2 calls), and the final
// unblocked columns go to the reference ?laqp2.  Everything else runs on the
// calling thread in the reference's order.
// ---------------------------------------------------------------------------
#include <algorithm>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
#include <omp.h>
#include "blas_declare.hpp"

extern "C" {
void dlarfg_(const int* n, double* alpha, double* x, const int* incx, double* tau);
void zlarfg_(const int* n, std::complex<double>* alpha, std::complex<double>* x, const int* incx,
             std::complex<double>* tau);
int idamax_(const int* n, const double* x, const int* incx);
void dswap_(const int* n, double* x, const int* incx, double* y, const int* incy);
void zswap_(const int* n, std::complex<double>* x, const int* incx, std::complex<double>* y,
            const int* incy);
void dlaqp2_(const int* m, const int* n, const int* offset, double* a, const int* lda, int* jpvt,
             double* tau, double* vn1, double* vn2, double* work);
void zlaqp2_(const int* m, const int* n, const int* offset, std::complex<double>* a, const int* lda,
             int* jpvt, std::complex<double>* tau, double* vn1, double* vn2,
             std::complex<double>* work);
double dlamch_(const char* cmach, std::size_t cmach_len);
int ilaenv_(const int* ispec, const char* name, const char* opts, const int* n1, const int* n2,
            const int* n3, const int* n4, std::size_t name_len, std::size_t opts_len);
}

namespace fmm {

// Box work under way on this thread: such a thread reaches task scheduling
// points (waiting in its own task loops) without being idle, and a QRCP
// helper task taken there would hold it for that whole factorization while
// its own box waited, so helpers return at once on it.
inline thread_local int box_work_depth = 0;
struct BoxWorkScope {
    BoxWorkScope() { ++box_work_depth; }
    ~BoxWorkScope() { --box_work_depth; }
    BoxWorkScope(const BoxWorkScope&) = delete;
    BoxWorkScope& operator=(const BoxWorkScope&) = delete;
};

namespace geqp3_detail {

// a step over fewer trailing columns runs on the owner alone
constexpr int min_chunk = 32;
// pieces of the per-column products start at multiples of row_step from
// their call's first entry (gemv 'N' and real gemv 'T' split bitwise at
// multiples of 4); the trailing gemm's at multiples of gemm_step (zgemm only
// at multiples of 6: its driver packs B in blocks of 3 x UNROLL_N = 6
// columns and computes a call's last ones differently)
constexpr int row_step = 4;
constexpr int gemm_step = 6;
// below this many columns the reference runs (the split's task overhead wins)
constexpr int min_columns = 256;
// A team only with this many tasks and while no other team runs in the
// process: a team's workers make many small BLAS calls a step, and this
// sequential OpenBLAS (USE_LOCKING) serializes concurrent small calls, so
// several teams at once (a wave of several boxes) ran 2-7x slower than the
// reference (7 boxes of 850 columns, 16 tasks: 375 vs 137 ms a QRCP; 32 of
// 420, 3 tasks: 415 vs 60 ms), while one team of 16 runs 3.5x faster.
constexpr int min_tasks = 8;
inline std::atomic<int> active_teams{0};

template <typename DataType>
constexpr bool is_complex = std::is_same_v<DataType, std::complex<double>>;

template <typename DataType>
void reference_geqp3(int m, int n, DataType* a, int lda, int* jpvt, DataType* tau) {
    int info = 0, lwork = -1;
    DataType query{};
    if constexpr (is_complex<DataType>) {
        std::vector<double> rwork(static_cast<size_t>(2 * std::max(n, 1)));
        zgeqp3_(&m, &n, a, &lda, jpvt, tau, &query, &lwork, rwork.data(), &info);
        lwork = static_cast<int>(query.real());
        std::vector<DataType> work(static_cast<size_t>(std::max(lwork, 1)));
        zgeqp3_(&m, &n, a, &lda, jpvt, tau, work.data(), &lwork, rwork.data(), &info);
    } else {
        dgeqp3_(&m, &n, a, &lda, jpvt, tau, &query, &lwork, &info);
        lwork = static_cast<int>(query);
        std::vector<DataType> work(static_cast<size_t>(std::max(lwork, 1)));
        dgeqp3_(&m, &n, a, &lda, jpvt, tau, work.data(), &lwork, &info);
    }
    if (info != 0) {
        throw std::runtime_error("RRQR factorization failed with INFO = " + std::to_string(info));
    }
}

template <typename DataType>
void gemv(const char* trans, int m, int n, DataType alpha, const DataType* a, int lda,
          const DataType* x, int incx, DataType beta, DataType* y, int incy) {
    if constexpr (is_complex<DataType>) {
        zgemv_(trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
    } else {
        dgemv_(trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
    }
}

template <typename DataType>
void swap(int n, DataType* x, int incx, DataType* y, int incy) {
    if constexpr (is_complex<DataType>) {
        zswap_(&n, x, &incx, y, &incy);
    } else {
        dswap_(&n, x, &incx, y, &incy);
    }
}

template <typename DataType>
void larfg(int n, DataType* alpha, DataType* x, int incx, DataType* tau) {
    if constexpr (is_complex<DataType>) {
        zlarfg_(&n, alpha, x, &incx, tau);
    } else {
        dlarfg_(&n, alpha, x, &incx, tau);
    }
}

template <typename DataType>
double nrm2(int n, const DataType* x) {
    const int one = 1;
    if constexpr (is_complex<DataType>) {
        return dznrm2_(&n, x, &one);
    } else {
        return dnrm2_(&n, x, &one);
    }
}

template <typename DataType>
double magnitude(const DataType& value) {
    // gfortran's ABS of a COMPLEX*16 is cabs, which glibc computes as hypot
    if constexpr (is_complex<DataType>) {
        return std::hypot(value.real(), value.imag());
    } else {
        return std::fabs(value);
    }
}

template <typename DataType>
DataType conj_if(const DataType& value) {
    if constexpr (is_complex<DataType>) {
        return std::conj(value);
    } else {
        return value;
    }
}

inline void cpu_relax() {
#if defined(__x86_64__) || defined(__i386__)
    asm volatile("pause" ::: "memory");
#elif defined(__aarch64__)
    asm volatile("yield" ::: "memory");
#endif
}

// The threads a small wave leaves idle, held for one factorization: helper
// tasks register as workers (slots 1..), then spin on a step counter; every
// step runs a function on each worker, which takes the pieces it owns.  A
// step costs microseconds where a task loop cost hundreds (a 128-thread
// team, ~2000 steps a factorization: a 1708 x 1704 complex QRCP takes 0.92 s
// on one thread, 1.3 s over task loops, 0.25 s on 16 workers here; CPU
// node, 2026-11-01), and a worker keeps its columns for the whole
// factorization (blocks dealt round-robin, see `owned`), so its share of the
// trailing matrix stays in its cache: dealing the pieces of each step to
// whichever thread came first moved the matrix between caches every step
// (slower than one thread past 8 workers, and on low-rank matrices, whose
// panels end after a column or two).  Helpers that have not started within
// `wait_us` (no idle thread) stay out; one the runtime runs on the owner (an
// undeferred task) or on a thread doing box work (BoxWorkScope) returns at
// once.
class SpinTeam {
public:
    // body(team) with up to `helpers` helper tasks; on return the helpers
    // that joined have left, and a helper task not yet started holds its own
    // reference to the team and leaves when it runs (waiting for it, as a
    // taskgroup would, made the owner run other boxes' queued tasks first:
    // 27-37 boxes a wave, the boxes' ID time grew 6x)
    template <typename Body>
    static void with(int helpers, const Body& body) {
        std::shared_ptr<SpinTeam> team(new SpinTeam());
        if (!omp_in_parallel() || omp_get_num_threads() <= 1) helpers = 0;
        for (int h = 0; h < helpers; ++h) {
            std::shared_ptr<SpinTeam> reference = team;
            #pragma omp task default(none) firstprivate(reference)
            reference->help();
        }
        team->close_registration(helpers);
        body(*team);
        team->stop_.store(true, std::memory_order_release);
        while (team->left_.load(std::memory_order_acquire) < team->workers_ - 1) cpu_relax();
    }

    int workers() const { return workers_; }

    // function(slot) on every worker (the owner: slot 0), all done on return
    template <typename Function>
    void run(const Function& function) {
        if (workers_ == 1) {
            function(0);
            return;
        }
        context_ = &function;
        invoke_ = [](const void* context, int slot) {
            (*static_cast<const Function*>(context))(slot);
        };
        done_.store(0, std::memory_order_relaxed);
        step_.fetch_add(1, std::memory_order_release);
        function(0);
        while (done_.load(std::memory_order_acquire) < workers_ - 1) cpu_relax();
    }

private:
    static constexpr uint32_t closed = 0x80000000u;
    static constexpr double wait_us = 200.0;

    SpinTeam() : owner_(omp_get_thread_num()) {}

    void close_registration(int helpers) {
        const double start = omp_get_wtime();
        while (static_cast<int>(registration_.load(std::memory_order_acquire) & ~closed) < helpers &&
               (omp_get_wtime() - start) * 1e6 < wait_us) {
            cpu_relax();
        }
        workers_ = static_cast<int>(registration_.fetch_or(closed, std::memory_order_acq_rel) & ~closed) + 1;
    }

    void help() {
        if (omp_get_thread_num() == owner_ || box_work_depth > 0) return;
        uint32_t word = registration_.load(std::memory_order_acquire);
        do {
            if (word & closed) return;
        } while (!registration_.compare_exchange_weak(word, word + 1, std::memory_order_acq_rel,
                                                      std::memory_order_acquire));
        const int slot = static_cast<int>(word) + 1;
        uint64_t seen = 0;
        while (true) {
            const uint64_t step = step_.load(std::memory_order_acquire);
            if (step != seen) {
                seen = step;
                invoke_(context_, slot);
                done_.fetch_add(1, std::memory_order_release);
            } else if (stop_.load(std::memory_order_acquire)) {
                left_.fetch_add(1, std::memory_order_release);
                return;
            } else {
                cpu_relax();
            }
        }
    }

    const int owner_;
    int workers_ = 1;
    const void* context_ = nullptr;
    void (*invoke_)(const void*, int) = nullptr;
    alignas(64) std::atomic<uint32_t> registration_{0};
    alignas(64) std::atomic<uint64_t> step_{0};
    alignas(64) std::atomic<int> done_{0};
    alignas(64) std::atomic<bool> stop_{false};
    std::atomic<int> left_{0};
};

// The pieces of local columns [lo, hi) that worker `slot` owns: the matrix's
// columns in blocks of `block`, dealt round-robin over the workers by their
// index in the whole matrix (`offset`: the index of local column 1, 0-based),
// each piece cut where a block ends, rounded up to the grid base + q step the
// call allows (lo on the grid).
constexpr int block = 32;
template <typename Body>
void owned(int lo, int hi, int base, int step, int offset, int slot, int workers, const Body& body) {
    int c = lo;
    while (c < hi) {
        const int global = offset + c - 1;
        const int owner = (global / block) % workers;
        int end = (global / block + 1) * block - offset + 1;
        end = base + ((end - base + step - 1) / step) * step;
        end = std::min(end, hi);
        if (owner == slot) body(c, end);
        c = end;
    }
}

// ?laqps: one panel of at most nb columns, 1-based as the reference.  Each
// column's products over the trailing columns (F(k+1:n, k), its correction by
// the panel's earlier columns, the row update of A(rk, k+1:n) and the partial
// column norms) run over the team, each worker on the pieces it owns:
//   complex: one step over row pieces of F starting at multiples of
//     `row_step` from row 1 (gemv 'C' and the one-row gemm split anywhere,
//     gemv 'N' at multiples of 4 from its first row);
//   real: the gemv 'T' in pieces from column k+1, the correction whole, then
//     the row update (a gemv 'N' of the rows k+1:n) and the norms in pieces
//     from k+1 (every split at multiples of 4 from its call's first entry).
template <typename DataType>
void laqps(int m, int n, int offset, int nb, int& kb, DataType* a, int lda, int* jpvt,
           DataType* tau, double* vn1, double* vn2, DataType* auxv, DataType* f, int ldf,
           SpinTeam& team, double tol3z) {
    auto A = [&](int i, int j) { return a + (i - 1) + static_cast<size_t>(j - 1) * lda; };
    auto F = [&](int i, int j) { return f + (i - 1) + static_cast<size_t>(j - 1) * ldf; };
    const DataType one(1.0), zero(0.0), minus_one(-1.0);
    const char* conj_trans = is_complex<DataType> ? "C" : "T";
    const int lastrk = std::min(m, n + offset);
    std::vector<int> flagged;  // the reference's LSTICC list (order irrelevant)
    const int workers = team.workers();
    const int offset_col = offset;  // (the panel's first column, 0-based: owned()'s offset)
    std::vector<std::vector<int>> worker_flagged(static_cast<size_t>(workers));
    int k = 0;
    while (k < nb && flagged.empty()) {
        ++k;
        const int rk = offset + k;
        {
            const int len = n - k + 1, inc = 1;
            const int pvt = (k - 1) + idamax_(&len, vn1 + (k - 1), &inc);
            if (pvt != k) {
                swap(m, A(1, pvt), 1, A(1, k), 1);
                swap(k - 1, F(pvt, 1), ldf, F(k, 1), ldf);
                std::swap(jpvt[pvt - 1], jpvt[k - 1]);
                vn1[pvt - 1] = vn1[k - 1];
                vn2[pvt - 1] = vn2[k - 1];
            }
        }
        if (k > 1) {
            for (int j = 1; j < k; ++j) *F(k, j) = conj_if(*F(k, j));
            gemv("N", m - rk + 1, k - 1, minus_one, A(rk, 1), lda, F(k, 1), ldf, one, A(rk, k), 1);
            for (int j = 1; j < k; ++j) *F(k, j) = conj_if(*F(k, j));
        }
        if (rk < m) {
            larfg(m - rk + 1, A(rk, k), A(rk + 1, k), 1, tau + (k - 1));
        } else {
            larfg(1, A(rk, k), A(rk, k), 1, tau + (k - 1));
        }
        const DataType akk = *A(rk, k);
        *A(rk, k) = one;
        const DataType tk = tau[k - 1];
        const int rows = m - rk + 1;
        if (k > 1) {
            gemv(conj_trans, rows, k - 1, -tk, A(rk, 1), lda, A(rk, k), 1, zero, auxv, 1);
        }
        const bool update_norms = rk < lastrk;

        // F(c0:c1-1, k) = tau A(rk:m, c0:c1-1)^H v
        auto products = [&](int c0, int c1) {
            if (c1 > c0) {
                gemv(conj_trans, rows, c1 - c0, tk, A(rk, c0), lda, A(rk, k), 1, zero, F(c0, k), 1);
            }
        };
        // F(r0:r1-1, k) += F(r0:r1-1, 1:k-1) auxv (rows <= k: zero first)
        auto correction = [&](int r0, int r1) {
            for (int r = r0; r < std::min(r1, k + 1); ++r) *F(r, k) = zero;
            if (k > 1 && r1 > r0) {
                gemv("N", r1 - r0, k - 1, one, F(r0, 1), ldf, auxv, 1, one, F(r0, k), 1);
            }
        };
        // A(rk, c0:c1-1) -= A(rk, 1:k) F(c0:c1-1, 1:k)^H, then the columns' norms
        auto row_and_norms = [&](int c0, int c1, std::vector<int>& marks) {
            if (c1 <= c0) return;
            const int cols = c1 - c0;
            if constexpr (is_complex<DataType>) {
                const int one_row = 1;
                zgemm_("N", "C", &one_row, &cols, &k, &minus_one, A(rk, 1), &lda, F(c0, 1), &ldf,
                       &one, A(rk, c0), &lda);
            } else {
                gemv("N", cols, k, minus_one, F(c0, 1), ldf, A(rk, 1), lda, one, A(rk, c0), lda);
            }
            if (!update_norms) return;
            for (int j = c0; j < c1; ++j) {
                double& v1 = vn1[j - 1];
                if (v1 != 0.0) {
                    double temp = magnitude(*A(rk, j)) / v1;
                    temp = std::max(0.0, (1.0 + temp) * (1.0 - temp));
                    const double ratio = v1 / vn2[j - 1];
                    const double temp2 = temp * (ratio * ratio);
                    if (temp2 <= tol3z) {
                        marks.push_back(j);
                    } else {
                        v1 = v1 * std::sqrt(temp);
                    }
                }
            }
        };

        const int trailing = n - k;
        if (workers == 1 || trailing < min_chunk) {
            products(k + 1, n + 1);
            correction(1, n + 1);
            if (trailing > 0) row_and_norms(k + 1, n + 1, flagged);
        } else if constexpr (is_complex<DataType>) {
            // pieces of the rows of F on the grid 1 + 4q (gemv 'N'), the same
            // columns' products, row update and norms
            team.run([&](int slot) {
                auto& marks = worker_flagged[static_cast<size_t>(slot)];
                marks.clear();
                owned(1, n + 1, 1, row_step, offset_col, slot, workers, [&](int r0, int r1) {
                    const int c0 = std::max(r0, k + 1);
                    products(c0, r1);
                    correction(r0, r1);
                    row_and_norms(c0, r1, marks);
                });
            });
        } else {
            // dgemv 'T' and the row update's dgemv 'N' cut on the grid k+1 + 4q
            team.run([&](int slot) {
                owned(k + 1, n + 1, k + 1, row_step, offset_col, slot, workers,
                      [&](int c0, int c1) { products(c0, c1); });
            });
            correction(1, n + 1);
            team.run([&](int slot) {
                auto& marks = worker_flagged[static_cast<size_t>(slot)];
                marks.clear();
                owned(k + 1, n + 1, k + 1, row_step, offset_col, slot, workers,
                      [&](int c0, int c1) { row_and_norms(c0, c1, marks); });
            });
        }
        if (workers > 1 && trailing >= min_chunk) {
            for (const auto& marks : worker_flagged) {
                flagged.insert(flagged.end(), marks.begin(), marks.end());
            }
        }
        *A(rk, k) = akk;
    }
    kb = k;
    const int rk = offset + kb;

    // the trailing block A(rk+1:m, kb+1:n) -= A(rk+1:m, 1:kb) F(kb+1:n, 1:kb)^H
    if (kb < std::min(n, m - offset)) {
        const int rows = m - rk, cols = n - kb;
        const char* tb = is_complex<DataType> ? "C" : "T";
        auto update = [&](int c0, int c1) {
            const int w = c1 - c0;
            gemm_("N", tb, &rows, &w, &kb, &minus_one, A(rk + 1, 1), &lda, F(c0, 1), &ldf,
                  &one, A(rk + 1, c0), &lda);
        };
        if (workers == 1 || cols < min_chunk) {
            update(kb + 1, n + 1);
        } else {
            team.run([&](int slot) {
                owned(kb + 1, n + 1, kb + 1, gemm_step, offset_col, slot, workers, update);
            });
        }
    }
    // the flagged columns' norms recomputed
    const int rows = m - rk;
    auto recompute = [&](int j) {
        vn1[j - 1] = nrm2(rows, A(rk + 1, j));
        vn2[j - 1] = vn1[j - 1];
    };
    if (workers == 1 || flagged.size() < 8) {
        for (const int j : flagged) recompute(j);
    } else {
        team.run([&](int slot) {
            for (const int j : flagged) {
                if (((offset_col + j - 1) / block) % workers == slot) recompute(j);
            }
        });
    }
}

}  // namespace geqp3_detail

// geqp3 of A (m x n, column-major): A, jpvt (1-based) and tau as the
// reference's, bit for bit.  `tasks`: the products of each step over the
// calling thread and up to tasks - 1 helper tasks (below min_tasks, or with
// another team running: the reference itself).
// jpvt is input as the reference's (nonzero entries: columns moved to the
// front); only an all-zero jpvt takes the split path.
template <typename DataType>
void geqp3_tasks(int m, int n, DataType* a, int lda, int* jpvt, DataType* tau, int tasks) {
    using namespace geqp3_detail;
    static_assert(is_complex<DataType> || std::is_same_v<DataType, double>,
                  "geqp3_tasks: double or complex<double>");
    const int minmn = std::min(m, n);
    bool free_columns = true;
    for (int j = 0; j < n && free_columns; ++j) free_columns = jpvt[j] == 0;
    const bool eligible = tasks >= min_tasks && minmn >= min_columns && free_columns &&
                          lda >= std::max(1, m);
    if (!eligible || active_teams.fetch_add(1, std::memory_order_acq_rel) != 0) {
        if (eligible) active_teams.fetch_sub(1, std::memory_order_acq_rel);
        reference_geqp3(m, n, a, lda, jpvt, tau);
        return;
    }
    struct TeamSlot {
        ~TeamSlot() { active_teams.fetch_sub(1, std::memory_order_acq_rel); }
    } slot;
    // (the reference with jpvt all zero: no fixed columns, NFXD = 0)
    for (int j = 0; j < n; ++j) jpvt[j] = j + 1;
    const char* geqrf = is_complex<DataType> ? "ZGEQRF" : "DGEQRF";
    const int ispec_nb = 1, ispec_nx = 3, unused = -1;
    int nb = ilaenv_(&ispec_nb, geqrf, " ", &m, &n, &unused, &unused, 6, 1);
    int nx = 0;
    if (nb > 1 && nb < minmn) {
        nx = std::max(0, ilaenv_(&ispec_nx, geqrf, " ", &m, &n, &unused, &unused, 6, 1));
    }
    const int nbmin = 2;
    std::vector<double> vn1(static_cast<size_t>(n)), vn2(static_cast<size_t>(n));
    int j = 1;
    SpinTeam::with(tasks - 1, [&](SpinTeam& team) {
        team.run([&](int slot) {
            owned(1, n + 1, 1, 1, 0, slot, team.workers(), [&](int c0, int c1) {
                for (int c = c0 - 1; c < c1 - 1; ++c) {
                    vn1[static_cast<size_t>(c)] = nrm2(m, a + static_cast<size_t>(c) * lda);
                    vn2[static_cast<size_t>(c)] = vn1[static_cast<size_t>(c)];
                }
            });
        });
        if (nb >= nbmin && nb < minmn && nx < minmn) {
            const double tol3z = std::sqrt(dlamch_("Epsilon", 7));
            std::vector<DataType> auxv(static_cast<size_t>(nb));
            std::vector<DataType> f(static_cast<size_t>(n) * nb);
            const int topbmn = minmn - nx;
            while (j <= topbmn) {
                const int jb = std::min(nb, topbmn - j + 1);
                int fjb = 0;
                laqps(m, n - j + 1, j - 1, jb, fjb, a + static_cast<size_t>(j - 1) * lda, lda,
                      jpvt + (j - 1), tau + (j - 1), vn1.data() + (j - 1), vn2.data() + (j - 1),
                      auxv.data(), f.data(), n - j + 1, team, tol3z);
                j += fjb;
            }
        }
    });
    if (j <= minmn) {
        const int cols = n - j + 1, offset = j - 1;
        std::vector<DataType> work(static_cast<size_t>(n));
        if constexpr (is_complex<DataType>) {
            zlaqp2_(&m, &cols, &offset, a + static_cast<size_t>(j - 1) * lda, &lda, jpvt + (j - 1),
                    tau + (j - 1), vn1.data() + (j - 1), vn2.data() + (j - 1), work.data());
        } else {
            dlaqp2_(&m, &cols, &offset, a + static_cast<size_t>(j - 1) * lda, &lda, jpvt + (j - 1),
                    tau + (j - 1), vn1.data() + (j - 1), vn2.data() + (j - 1), work.data());
        }
    }
}

}  // namespace fmm
