#pragma once

#ifdef H2_HAVE_GPU

// A butterfly from a black-box matvec, on the device (GPU_BACKEND/bf_gpu):
// BF_randomized with less_adapt 1 (Bplus_randomized.f90), for a butterfly on
// one process.  The algorithm is the host's, level by level:
//
//   the outer levels, V (level 0) and U (level L+1), with a rank guess that
//   grows (rank0 * rankrate^(tt-1) + 1) until neither reaches it;
//   then the left half's levels 1 .. level_half and the right half's
//   L .. level_half+1, their sample counts from the levels already built.
//
// A left level l takes 2^l bundles of random vectors, each nonzero on the
// rows of one block row of that level; applies the black box transposed;
// sweeps the half-built butterfly over the result up to level l-1 (its V and
// Kerl(1 .. l-1), the 'N' direction); and finds the range of each block's
// panel - the two vector blocks its pair of kernel blocks maps, over the
// bundle's columns - by a column-pivoted QR (launch_qrcp, tol_Rdetect), whose
// Q^T are the new kernel blocks.  The right levels mirror it over block
// columns with the black box itself and the 'T' sweep from U; at the level
// where the halves meet (level_half+1) the kernel blocks instead solve a
// least-squares problem against the left half's sweep of the same random
// vectors.  The factors stay on the device as they are built, which is what
// the sweeps of the later levels need.
//
// Two departures from the host, neither in what is computed: the solve is a
// pivoted QR's basic solution where the host's LinearSolve truncates an SVD
// (the same answer unless the left half's sketch is rank deficient to
// tol_LS), and the black box's norms behind the outer levels' negligible-block
// test are taken once rather than at every try of the rank guess.
//
// The host side is kept to one arena allocation and one metadata upload per
// batch of panels, and one upload per batch of GEMMs (their nine arrays
// packed): measured at 262144^2 (L = 12), a first version that allocated
// every block and uploaded every array on its own spent 1.5 s of its 2.3 s of
// range finding on the host.
//
// On several ranks (the butterfly's process group) the levels are the host's
// too, over the blocks this rank owns, and so are the sweeps: the
// distributed matvec's plans stopped at a level (bf_sweep_n.hpp,
// BF_block_MVP_partial), with their exchanges and, where the halves meet, the
// midpoint redistribution.  The rank guess, the sample counts and the norms
// are the maxima and sums over the group.  The random vectors are hashed by
// the global row, so a run on several ranks draws the numbers one process
// draws.  Not the masked sweeps of one process: the host computes every
// block of every level there, and so does this.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <mpi.h>

#include "bf_device.hpp"
#include "bf_sweep.hpp"
#include "bf_sweep_n.hpp"
#include "../h2_gpu/device_kernels.hpp"
#include "../h2_gpu/gpu_runtime.hpp"

namespace bf {
namespace gpu {

// A panel of the range finder or of the solve: m0 rows of s0 then m1 of s1,
// columns col0 .. col0 + ncols, into dst (ld ldd), or its transpose.
template<typename T>
struct PanelItem {
    using DataType = T;
    const DataType* s0 = nullptr;
    int ld0 = 0, m0 = 0;
    const DataType* s1 = nullptr;
    int ld1 = 0, m1 = 0;
    int col0 = 0, ncols = 0;
    DataType* dst = nullptr;
    int ldd = 0;
    int transposed = 0;
};

// Q of a factorized panel (m x n, ld m), its first `rank` columns, split:
// rows 0 .. m0 into out0 and m0 .. m into out1, as they are (rows x rank) or
// transposed (rank x rows).  zero: Q is zero (a panel the host's range finder
// calls empty or negligible).  q: scratch, m x rank.
template<typename T>
struct RangeItem {
    using DataType = T;
    DataType* a = nullptr;
    const DataType* tau = nullptr;
    DataType* q = nullptr;
    int m = 0, rank = 0, m0 = 0;
    DataType* out0 = nullptr;
    DataType* out1 = nullptr;
    int transposed = 0;
    int zero = 0;
};

// X = A^+ B for A (mm x r, ld mm) factorized by launch_qrcp (rank *k, pivots
// jpvt, scalars tau, normalized by *scale) and B (mm x n, ld mm, destroyed);
// X(:, 0 .. m0)^T into out0 (m0 x r) and X(:, m0 .. n)^T into out1, both
// zeroed beforehand.
template<typename T>
struct SolveItem {
    using DataType = T;
    DataType* a = nullptr;
    const DataType* tau = nullptr;
    const int* jpvt = nullptr;
    const int* k = nullptr;
    const double* scale = nullptr;
    int mm = 0, r = 0;
    DataType* b = nullptr;
    int n = 0, m0 = 0;
    DataType* out0 = nullptr;
    DataType* out1 = nullptr;
};

// n elements from src to dst
template<typename T>
struct CopyItem {
    using DataType = T;
    const DataType* src = nullptr;
    DataType* dst = nullptr;
    int64_t n = 0;
};

// (each defined for double and fmm::gpu::dcomplex in bf_randomized.cu)
template<typename T>
void launch_fill_uniform(T* x, int64_t ld, int64_t row0, int64_t nrows, int64_t col0, int64_t ncols, uint64_t seed,
                         int64_t hash_row0, cudaStream_t stream);
template<typename T>
void launch_fill_value(T* x, int64_t n, T v, cudaStream_t stream);
template<typename T>
void launch_panels(const PanelItem<T>* items, int count, double* amax, cudaStream_t stream);
template<typename T>
void launch_range(const RangeItem<T>* items, int count, cudaStream_t stream);
template<typename T>
void launch_solve(const SolveItem<T>* items, int count, cudaStream_t stream);
template<typename T>
void launch_copies(const CopyItem<T>* items, int count, cudaStream_t stream);

// The pinned host buffers of the constructions, kept for the process (as
// fmm::gpu::pinned_pool: pinning costs tens of milliseconds, and freeing from
// a static destructor would race the CUDA runtime's shutdown).  Two metadata
// images, because a MetaBuilder waits for its own last upload before it is
// refilled: the host fills one while the device still reads the other.  Two
// download chunks, copied into while the last one is unpacked.
struct RandomizedHost {
    static constexpr int kImages = 2;
    fmm::gpu::MetaBuilder meta[kImages];
    fmm::gpu::DeviceBuffer device[kImages];
    int next = 0;
    static constexpr size_t kChunk = size_t{32} << 20;
    fmm::gpu::PinnedBuffer chunk[2];
    cudaEvent_t copied[2] = {nullptr, nullptr};
};
inline RandomizedHost& randomized_host() {
    static RandomizedHost* host = new RandomizedHost;
    return *host;
}

// One level of the butterfly being built: an nr x nc grid of column-major
// blocks (ld = rows), row-major in b.  On several ranks the grid is this
// rank's blocks, the global (one-based) index of local block (i, j) being
// idx_r + i inc_r, idx_c + j inc_c, as the Fortran's descriptors give it.
template<typename T>
struct RBlock {
    using DataType = T;
    DataType* ptr = nullptr;
    int rows = 0, cols = 0;
};
template<typename T>
struct RLevel {
    int nr = 0, nc = 0;
    int idx_r = 1, inc_r = 1, idx_c = 1, inc_c = 1;
    std::vector<RBlock<T>> b;
    std::vector<char*> owned;
    RBlock<T>& at(int i, int j) { return b[static_cast<size_t>(i) * nc + j]; }
    const RBlock<T>& at(int i, int j) const { return b[static_cast<size_t>(i) * nc + j]; }
};

// The black box: y = op(A) x, all on the device; op 'N' takes x N x nvec to
// y M x nvec, 'T' x M x nvec to y N x nvec (ld = rows).  It writes all of y,
// and is complete on return.
template<typename T>
using BlackBox = std::function<void(char op, int nvec, const T* x, T* y)>;

template<typename T>
struct RandomizedConfig {
    int L = 0, level_half = 0;
    int64_t M = 0, N = 0;
    // the leaves' ends, from 0: of every leaf (2^L + 1 each) on one process,
    // of this rank's U and V blocks in its rows and columns on several
    std::vector<int64_t> row_off, col_off;
    // several ranks: the block's communicator and size, this rank's rows and
    // columns and the global index of the first of each, the butterfly's
    // local grids and its sweeps' ranges (layout: levels and ranges, no
    // blocks), and the Fortran's rules for its exchanges (bf_ranges.hpp)
    int nproc = 1;
    MPI_Comm comm = MPI_COMM_NULL;
    int64_t M_loc = 0, N_loc = 0, row0 = 0, col0 = 0;
    BfDevice<T> layout;
    BfBridge bridge;
    int rank0 = 32;
    double rankrate = 1.2;
    int itermax = 10, nbundle = 1;
    double tol_rdetect = 1e-5, tol_ls = 1e-12, abs_tol = 0.0;
    uint64_t seed = 20261025;
};

// Seconds: the black box, the sweeps of the half-built butterfly, the range
// finding and solves, the closing error test, all of the construction; the
// black box's calls and columns, the tries of the rank guess, and the test's
// relative error (BF_Test_Reconstruction_Error's).
struct RandomizedProfile {
    double matvec = 0, sweep = 0, range = 0, test = 0, total = 0;
    int64_t calls = 0, columns = 0;
    int tries = 0;
    double error = 0;
};

template<typename T>
class Randomized {
public:
    // the scalar's types, by the names the members use
    using DataType = T;
    using PanelItem = bf::gpu::PanelItem<T>;
    using RangeItem = bf::gpu::RangeItem<T>;
    using SolveItem = bf::gpu::SolveItem<T>;
    using CopyItem = bf::gpu::CopyItem<T>;
    using RBlock = bf::gpu::RBlock<T>;
    using RLevel = bf::gpu::RLevel<T>;
    using BlackBox = bf::gpu::BlackBox<T>;
    using RandomizedConfig = bf::gpu::RandomizedConfig<T>;
    using BfDevice = bf::gpu::BfDevice<T>;
    using SweepPlan = bf::gpu::SweepPlan<T>;
    using VecLevel = bf::gpu::VecLevel<T>;
    using GemmItem = bf::gpu::GemmItem<T>;
    using SweepScratch = bf::gpu::SweepScratch<T>;

    Randomized(const RandomizedConfig& config, BlackBox blackbox) : c_(config), bb_(std::move(blackbox)) {}
    ~Randomized() { release(); }
    Randomized(const Randomized&) = delete;
    Randomized& operator=(const Randomized&) = delete;

    // false: no room on the device, or the outer levels' rank guess never
    // settled (the host's construction stops there too); nothing kept
    bool run() {
        const auto t0 = clock::now();
        try {
            const int L = c_.L;
            lv_.assign(static_cast<size_t>(L) + 2, RLevel());
            for (int l = 0; l <= L + 1; ++l) {
                RLevel& v = lv_[l];
                if (dist()) {
                    const LevelDesc& d = c_.layout.level(l);
                    v.nr = d.nr;
                    v.nc = d.nc;
                    v.idx_r = d.idx_r;
                    v.inc_r = d.inc_r;
                    v.idx_c = d.idx_c;
                    v.inc_c = d.inc_c;
                } else {
                    v.nr = l == 0 ? 1 : (l == L + 1 ? (1 << L) : (1 << l));
                    v.nc = l == 0 ? (1 << L) : (l == L + 1 ? 1 : (1 << (L - l + 1)));
                }
                v.b.assign(static_cast<size_t>(v.nr) * v.nc, RBlock());
            }
            bool converged = false;
            for (int tt = 1; tt <= c_.itermax && !converged; ++tt) {
                const int cap = static_cast<int>(std::ceil(c_.rank0 * std::pow(c_.rankrate, tt - 1))) + 1;
                drop_level(0);
                drop_level(L + 1);
                left_level(0, cap);
                right_level(L + 1, cap);
                int64_t r0 = 0, r1 = 0;
                for (const RBlock& b : lv_[0].b) r0 = std::max<int64_t>(r0, b.cols);
                for (const RBlock& b : lv_[L + 1].b) r1 = std::max<int64_t>(r1, b.cols);
                converged = all_max(r0) < cap && all_max(r1) < cap;
                prof.tries = tt;
            }
            if (!converged) {
                release();
                return false;
            }
            for (int l = 1; l <= c_.level_half; ++l) left_level(l, -1);
            for (int l = L; l >= c_.level_half + 1; --l) right_level(l, -1);
            test_error();
        } catch (const NoRoom&) {
            release();
            return false;
        }
        prof.total = seconds_since(t0);
        return true;
    }

    // Every block's shape, V's, Kerl(1 .. L)'s row by row, then U's.
    void shapes(std::vector<int>& rows, std::vector<int>& cols) const {
        rows.clear();
        cols.clear();
        for (const RLevel& v : lv_) {
            for (const RBlock& b : v.b) {
                rows.push_back(b.rows);
                cols.push_back(b.cols);
            }
        }
    }

    // The blocks, in that order, into dst (column major; null or empty:
    // skipped): one gather on the device, then chunks through pinned memory,
    // each unpacked while the next one is copied.
    void download(const std::vector<DataType*>& dst) {
        std::vector<const RBlock*> blocks;
        for (const RLevel& v : lv_)
            for (const RBlock& b : v.b) blocks.push_back(&b);
        if (dst.size() != blocks.size()) throw std::runtime_error("bf::gpu::Randomized::download: wrong block count");
        std::vector<int64_t> at(blocks.size() + 1, 0);
        for (size_t k = 0; k < blocks.size(); ++k)
            at[k + 1] = at[k] + static_cast<int64_t>(blocks[k]->rows) * blocks[k]->cols;
        const int64_t total = at.back();
        if (total == 0) return;
        std::vector<char*> own;
        DataType* packed = take(total, own);
        std::vector<CopyItem> items;
        for (size_t k = 0; k < blocks.size(); ++k) {
            if (at[k + 1] > at[k]) items.push_back(CopyItem{blocks[k]->ptr, packed + at[k], at[k + 1] - at[k]});
        }
        Image img = begin_image();
        const size_t oi = img.meta->append(items);
        char* d = img.meta->upload(*img.device, stream());
        launch_copies(reinterpret_cast<const CopyItem*>(d + oi), static_cast<int>(items.size()), stream());

        RandomizedHost& h = randomized_host();
        const int64_t step = static_cast<int64_t>(RandomizedHost::kChunk / sizeof(DataType));
        DataType* buf[2];
        for (int s = 0; s < 2; ++s) {
            buf[s] = static_cast<DataType*>(h.chunk[s].reserve(RandomizedHost::kChunk));
            if (h.copied[s] == nullptr) {
                fmm::gpu::check_cuda(cudaEventCreateWithFlags(&h.copied[s], cudaEventDisableTiming), "bf randomized download");
            }
        }
        auto copy = [&](int64_t first, int s) {
            const int64_t n = std::min(step, total - first);
            fmm::gpu::check_cuda(cudaMemcpyAsync(buf[s], packed + first, static_cast<size_t>(n) * sizeof(DataType),
                                                 cudaMemcpyDeviceToHost, stream()),
                                 "bf randomized download");
            fmm::gpu::check_cuda(cudaEventRecord(h.copied[s], stream()), "bf randomized download");
        };
        copy(0, 0);
        size_t k0 = 0;  // the first block the chunk may hold
        for (int64_t first = 0, c = 0; first < total; first += step, ++c) {
            const int s = static_cast<int>(c & 1);
            if (first + step < total) copy(first + step, 1 - s);
            fmm::gpu::check_cuda(cudaEventSynchronize(h.copied[s]), "bf randomized download");
            const int64_t last = std::min(first + step, total);
            while (k0 < blocks.size() && at[k0 + 1] <= first) ++k0;
            size_t k1 = k0;
            while (k1 < blocks.size() && at[k1] < last) ++k1;
            const DataType* src = buf[s];
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
            for (int64_t k = static_cast<int64_t>(k0); k < static_cast<int64_t>(k1); ++k) {
                const int64_t lo = std::max(at[k], first), hi = std::min(at[k + 1], last);
                if (dst[k] != nullptr && hi > lo) {
                    std::memcpy(dst[k] + (lo - at[k]), src + (lo - first), static_cast<size_t>(hi - lo) * sizeof(DataType));
                }
            }
        }
        give_back(own);
    }

    void release() {
        drop_plans();
        for (size_t l = 0; l < lv_.size(); ++l) drop_level(static_cast<int>(l));
        lv_.clear();
    }

    RandomizedProfile prof;

private:
    using clock = std::chrono::steady_clock;
    struct NoRoom {};

    static double seconds_since(clock::time_point t0) {
        return std::chrono::duration<double>(clock::now() - t0).count();
    }

    // ---- several ranks
    bool dist() const { return c_.nproc > 1; }
    int64_t all_max(int64_t v) const {
        if (dist()) MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_INT64_T, MPI_MAX, c_.comm);
        return v;
    }
    double all_sum(double v) const {
        if (dist()) MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_SUM, c_.comm);
        return v;
    }
    // The bundle's vectors, in and out, and whether every rank has room for
    // them and twice as much again for the sweeps and the range finding: declining
    // must happen on all ranks at once, since the host's construction, which
    // takes over, is collective too.  Past this point a rank out of room
    // stops the run instead (take_bytes).
    void bundle_room(int64_t nin, int64_t nout, DataType*& in, DataType*& out, std::vector<char*>& own) {
        auto& heap = fmm::gpu::DeviceHeap::instance();
        heap.ensure_initialized();
        const size_t bin = static_cast<size_t>(std::max<int64_t>(nin, 1)) * sizeof(DataType);
        const size_t bout = static_cast<size_t>(std::max<int64_t>(nout, 1)) * sizeof(DataType);
        char* a = heap.try_alloc(bin, fmm::gpu::ResidentKind::Matvec);
        char* b = a != nullptr ? heap.try_alloc(bout, fmm::gpu::ResidentKind::Matvec) : nullptr;
        char* spare = b != nullptr ? heap.try_alloc(2 * (bin + bout), fmm::gpu::ResidentKind::Matvec) : nullptr;
        int room = spare != nullptr ? 1 : 0;
        if (spare != nullptr) heap.free(spare);
        if (dist()) MPI_Allreduce(MPI_IN_PLACE, &room, 1, MPI_INT, MPI_MIN, c_.comm);
        if (!room) {
            if (b != nullptr) heap.free(b);
            if (a != nullptr) heap.free(a);
            throw NoRoom();
        }
        own.push_back(a);
        own.push_back(b);
        in = reinterpret_cast<DataType*>(a);
        out = reinterpret_cast<DataType*>(b);
    }

    // ---- device memory, from the shared arena
    char* take_bytes(size_t bytes, std::vector<char*>& owner) {
        auto& heap = fmm::gpu::DeviceHeap::instance();
        heap.ensure_initialized();
        char* p = heap.try_alloc(std::max<size_t>(bytes, 1), fmm::gpu::ResidentKind::Matvec);
        if (p == nullptr) {
            // one rank alone cannot decline: the others would wait for it in
            // the next collective
            if (dist()) throw std::runtime_error("bf::gpu::Randomized: no room on the device");
            throw NoRoom();
        }
        owner.push_back(p);
        return p;
    }
    DataType* take(int64_t elements, std::vector<char*>& owner) {
        return reinterpret_cast<DataType*>(take_bytes(static_cast<size_t>(std::max<int64_t>(elements, 1)) * sizeof(DataType), owner));
    }
    static void give_back(std::vector<char*>& owner) {
        auto& heap = fmm::gpu::DeviceHeap::instance();
        for (char* p : owner) heap.free(p);
        owner.clear();
    }
    void drop_level(int l) {
        if (l < 0 || static_cast<size_t>(l) >= lv_.size()) return;
        give_back(lv_[l].owned);
        for (RBlock& b : lv_[l].b) b = RBlock();
    }

    // ---- several ranks: the half-built butterfly as the sweep plans read
    // it, its blocks where they are, and the plans of its partial sweeps,
    // each built once (collectively) when its levels are final
    void refresh_layout() {
        dev_ = c_.layout;
        dev_.blocks.clear();
        dev_.values = nullptr;
        dev_.elements = 0;
        for (size_t l = 0; l < lv_.size(); ++l) {
            dev_.levels[l].first = static_cast<int64_t>(dev_.blocks.size());
            for (const RBlock& b : lv_[l].b) {
                BlockDesc d;
                d.rows = b.rows;
                d.cols = b.cols;
                d.lda = std::max(b.rows, 1);
                d.ptr = b.ptr;
                dev_.blocks.push_back(d);
            }
        }
    }
    // the 'N' sweep to level stop (the midpoint redistribution when it is
    // level_half), or the 'T' sweep from U down to level stop
    SweepPlan& partial_plan(char op, int stop) {
        std::map<int, SweepPlan>& plans = op == 'N' ? plans_n_ : plans_t_;
        auto it = plans.find(stop);
        if (it != plans.end()) return it->second;
        refresh_layout();
        SweepPlan& plan = plans[stop];
        build_sweep_plan(c_.bridge, dev_, op, 1, plan, stop);
        return plan;
    }
    void drop_plans() {
        for (auto& e : plans_n_) e.second.release();
        for (auto& e : plans_t_) e.second.release();
        plans_n_.clear();
        plans_t_.clear();
    }
    // The zero-based local index of a one-based global index on a grid, or -1
    static int local_of(int index, int idx, int inc, int n) {
        const int d = index - idx;
        if (inc <= 0 || d < 0 || d % inc != 0 || d / inc >= n) return -1;
        return d / inc;
    }
    // byte offsets of several arrays in one allocation
    struct Layout {
        size_t size = 0;
        template<typename E>
        size_t add(int64_t count) {
            const size_t at = size;
            size = fmm::gpu::align_up(at + static_cast<size_t>(std::max<int64_t>(count, 1)) * sizeof(E));
            return at;
        }
    };

    // ---- the metadata of the launches: an image of the ring, filled, sent
    struct Image {
        fmm::gpu::MetaBuilder* meta;
        fmm::gpu::DeviceBuffer* device;
    };
    static Image begin_image() {
        RandomizedHost& h = randomized_host();
        const int s = h.next;
        h.next = (s + 1) % RandomizedHost::kImages;
        h.meta[s].clear();
        return Image{&h.meta[s], &h.device[s]};
    }
    static cudaStream_t stream() { return fmm::gpu::Context::instance().stream(); }

    // ---- a vector level of a sweep: an nr x nc grid of blocks h x nvec
    struct VLev {
        int nr = 0, nc = 0, nvec = 0;
        std::vector<int> h;
        std::vector<int64_t> off;
        DataType* data = nullptr;
        int64_t elements = 0;
        std::vector<char*> owned;
        DataType* at(int i, int j) const { return data + off[static_cast<size_t>(i) * nc + j]; }
        int rows(int i, int j) const { return h[static_cast<size_t>(i) * nc + j]; }
    };
    void shape(VLev& v, int nr, int nc, int nvec) {
        v.nr = nr;
        v.nc = nc;
        v.nvec = nvec;
        v.h.assign(static_cast<size_t>(nr) * nc, 0);
    }
    void place(VLev& v) {
        v.off.assign(v.h.size(), 0);
        int64_t at = 0;
        for (size_t b = 0; b < v.h.size(); ++b) {
            v.off[b] = at;
            at += static_cast<int64_t>(v.h[b]) * v.nvec;
        }
        v.elements = at;
        v.data = take(at, v.owned);
    }
    // a level some of whose blocks no product writes with beta 0
    void zero(VLev& v) {
        fmm::gpu::check_cuda(cudaMemsetAsync(v.data, 0, static_cast<size_t>(std::max<int64_t>(v.elements, 1)) * sizeof(DataType),
                                             stream()),
                             "bf randomized vectors");
    }

    // A batch of GEMMs, C = op(A) B + beta C, its nine arrays in one upload,
    // in chunks MAGMA takes
    void gemms(std::vector<GemmItem>& items, DataType beta, bool trans_a) {
        constexpr size_t kChunk = 32768;
        for (size_t first = 0; first < items.size(); first += kChunk) {
            const size_t count = std::min(items.size() - first, kChunk);
            Image img = begin_image();
            const size_t oa = img.meta->alloc(count * sizeof(DataType*));
            const size_t ob = img.meta->alloc(count * sizeof(DataType*));
            const size_t oc = img.meta->alloc(count * sizeof(DataType*));
            const size_t oi = img.meta->alloc(6 * (count + 1) * sizeof(magma_int_t));  // MAGMA reads batch + 1
            const DataType** ha = reinterpret_cast<const DataType**>(img.meta->host(oa));
            const DataType** hb = reinterpret_cast<const DataType**>(img.meta->host(ob));
            DataType** hc = reinterpret_cast<DataType**>(img.meta->host(oc));
            magma_int_t* hi = reinterpret_cast<magma_int_t*>(img.meta->host(oi));
            magma_int_t *hm = hi, *hn = hi + (count + 1), *hk = hi + 2 * (count + 1);
            magma_int_t *hlda = hi + 3 * (count + 1), *hldb = hi + 4 * (count + 1), *hldc = hi + 5 * (count + 1);
            magma_int_t max_m = 0, max_n = 0, max_k = 0;
            for (size_t t = 0; t < count; ++t) {
                const GemmItem& it = items[first + t];
                ha[t] = it.a;
                hb[t] = it.b;
                hc[t] = it.c;
                hm[t] = it.m;
                hn[t] = it.n;
                hk[t] = it.k;
                hlda[t] = it.lda;
                hldb[t] = it.ldb;
                hldc[t] = it.ldc;
                max_m = std::max<magma_int_t>(max_m, it.m);
                max_n = std::max<magma_int_t>(max_n, it.n);
                max_k = std::max<magma_int_t>(max_k, it.k);
            }
            hm[count] = hn[count] = hk[count] = hlda[count] = hldb[count] = hldc[count] = 0;
            char* d = img.meta->upload(*img.device, stream());
            magma_int_t* di = reinterpret_cast<magma_int_t*>(d + oi);
            fmm::gpu::gemm_vbatched<DataType>(
                trans_a ? MagmaTrans : MagmaNoTrans, MagmaNoTrans, di, di + (count + 1), di + 2 * (count + 1),
                DataType(1), reinterpret_cast<const DataType* const*>(d + oa), di + 3 * (count + 1),
                reinterpret_cast<const DataType* const*>(d + ob), di + 4 * (count + 1), beta,
                reinterpret_cast<DataType**>(d + oc), di + 5 * (count + 1), static_cast<magma_int_t>(count),
                fmm::gpu::Context::instance().queue(), max_m, max_n, max_k);
        }
        items.clear();
    }

    // ---- the sweeps of the half-built butterfly.  A bundle of random vectors
    // reads one block row of the left half's last vector level, or one block
    // column of the right half's, or (where the halves meet) one column of the
    // left half's; a block depends on one row of each level below it (left),
    // or one column (right), or a cone of columns.  So a sweep computes the
    // blocks its caller reads and what they depend on, from `want` - a mask
    // over the last level's grid, row-major, empty for all of it - where the
    // host's partial products compute every block of every level.
    using Mask = std::vector<char>;

    // vec(e + 1) of the 'N' sweep of x (N x nvec): V^T, then Kerl(1 .. e).
    // vec(k) is a 2^(k-1) x 2^(L-k+1) grid; block (i, j) of vec(k + 1) is
    // Kerl(k)(i, 2j) vec(k)(i/2, 2j) + Kerl(k)(i, 2j+1) vec(k)(i/2, 2j+1).
    VLev sweep_n(const DataType* x, int nvec, int e, const Mask& want) {
        const auto t0 = clock::now();
        const int L = c_.L;
        std::vector<Mask> need(static_cast<size_t>(e) + 2);
        need[e + 1] = want.empty() ? Mask(static_cast<size_t>(1) << L, 1) : want;
        for (int k = e; k >= 1; --k) {
            const int nc_next = 1 << (L - k), nc = 1 << (L - k + 1);
            need[k].assign(static_cast<size_t>(1) << L, 0);
            for (size_t b = 0; b < need[k + 1].size(); ++b) {
                if (!need[k + 1][b]) continue;
                const int i = static_cast<int>(b / nc_next), j = static_cast<int>(b % nc_next);
                need[k][static_cast<size_t>(i >> 1) * nc + 2 * j] = 1;
                need[k][static_cast<size_t>(i >> 1) * nc + 2 * j + 1] = 1;
            }
        }
        std::vector<GemmItem> items, items1;
        VLev cur;
        shape(cur, 1, 1 << L, nvec);
        for (int b = 0; b < (1 << L); ++b) cur.h[b] = lv_[0].at(0, b).cols;
        place(cur);
        bool partial = false;
        for (int b = 0; b < (1 << L); ++b) {
            if (!need[1][b]) continue;
            const RBlock& v = lv_[0].at(0, b);
            if (v.rows == 0 || v.cols == 0) {
                partial = true;
                continue;
            }
            items.push_back(GemmItem{v.ptr, x + c_.col_off[b], cur.at(0, b), v.cols, nvec, v.rows, v.rows,
                                     static_cast<int>(c_.N), v.cols});
        }
        if (partial) zero(cur);
        gemms(items, DataType(0), true);
        for (int l = 1; l <= e; ++l) {
            VLev next;
            shape(next, 1 << l, 1 << (L - l), nvec);
            for (int i = 0; i < next.nr; ++i)
                for (int j = 0; j < next.nc; ++j) next.h[static_cast<size_t>(i) * next.nc + j] = lv_[l].at(i, 2 * j).rows;
            place(next);
            partial = false;
            for (int w = 0; w < 2; ++w) {
                std::vector<GemmItem>& list = w == 0 ? items : items1;
                for (int i = 0; i < next.nr; ++i) {
                    for (int j = 0; j < next.nc; ++j) {
                        if (!need[l + 1][static_cast<size_t>(i) * next.nc + j]) continue;
                        const RBlock& k = lv_[l].at(i, 2 * j + w);
                        const int hin = cur.rows(i >> 1, 2 * j + w);
                        if (k.rows == 0 || k.cols == 0 || hin == 0) {
                            if (w == 0) partial = true;
                            continue;
                        }
                        list.push_back(GemmItem{k.ptr, cur.at(i >> 1, 2 * j + w), next.at(i, j), k.rows, nvec,
                                                k.cols, k.rows, hin, k.rows});
                    }
                }
            }
            if (partial) zero(next);
            gemms(items, DataType(0), false);
            gemms(items1, DataType(1), false);
            give_back(cur.owned);
            cur = std::move(next);
        }
        prof.sweep += seconds_since(t0);
        return cur;
    }

    // vec(t) of the 'T' sweep of y (M x nvec): U^T, then Kerl(L .. t)^T.
    // After Kerl(k)^T the grid is 2^(k-1) x 2^(L-k+1) (after U^T, 2^L x 1);
    // its block (I, J) is the sum over w of Kerl(k)(2I+w, J)^T times block
    // (2I+w, J/2) of the level before.
    VLev sweep_t(const DataType* y, int nvec, int t, const Mask& want) {
        const auto t0 = clock::now();
        const int L = c_.L;
        // need[k]: the grid after Kerl(k)^T, k = t .. L; need[L + 1]: after U^T
        std::vector<Mask> need(static_cast<size_t>(L) + 2);
        need[t] = want.empty() ? Mask(static_cast<size_t>(1) << L, 1) : want;
        for (int k = t; k <= L; ++k) {
            const int nc = 1 << (L - k + 1), nc_prev = k == L ? 1 : (1 << (L - k));
            need[k + 1].assign(static_cast<size_t>(1) << L, 0);
            for (size_t b = 0; b < need[k].size(); ++b) {
                if (!need[k][b]) continue;
                const int I = static_cast<int>(b / nc), J = static_cast<int>(b % nc);
                need[k + 1][static_cast<size_t>(2 * I) * nc_prev + (J >> 1)] = 1;
                need[k + 1][static_cast<size_t>(2 * I + 1) * nc_prev + (J >> 1)] = 1;
            }
        }
        std::vector<GemmItem> items, items1;
        VLev cur;
        shape(cur, 1 << L, 1, nvec);
        for (int a = 0; a < (1 << L); ++a) cur.h[a] = lv_[L + 1].at(a, 0).cols;
        place(cur);
        bool partial = false;
        for (int a = 0; a < (1 << L); ++a) {
            if (!need[L + 1][a]) continue;
            const RBlock& u = lv_[L + 1].at(a, 0);
            if (u.rows == 0 || u.cols == 0) {
                partial = true;
                continue;
            }
            items.push_back(GemmItem{u.ptr, y + c_.row_off[a], cur.at(a, 0), u.cols, nvec, u.rows, u.rows,
                                     static_cast<int>(c_.M), u.cols});
        }
        if (partial) zero(cur);
        gemms(items, DataType(0), true);
        for (int l = L; l >= t; --l) {
            VLev next;
            shape(next, 1 << (l - 1), 1 << (L - l + 1), nvec);
            for (int i = 0; i < next.nr; ++i)
                for (int j = 0; j < next.nc; ++j) next.h[static_cast<size_t>(i) * next.nc + j] = lv_[l].at(2 * i, j).cols;
            place(next);
            partial = false;
            for (int w = 0; w < 2; ++w) {
                std::vector<GemmItem>& list = w == 0 ? items : items1;
                for (int i = 0; i < next.nr; ++i) {
                    for (int j = 0; j < next.nc; ++j) {
                        if (!need[l][static_cast<size_t>(i) * next.nc + j]) continue;
                        const RBlock& k = lv_[l].at(2 * i + w, j);
                        const int hin = cur.rows(2 * i + w, j >> 1);
                        if (k.rows == 0 || k.cols == 0 || hin == 0) {
                            if (w == 0) partial = true;
                            continue;
                        }
                        list.push_back(GemmItem{k.ptr, cur.at(2 * i + w, j >> 1), next.at(i, j), k.cols, nvec,
                                                k.rows, k.rows, hin, k.cols});
                    }
                }
            }
            if (partial) zero(next);
            gemms(items, DataType(0), true);
            gemms(items1, DataType(1), true);
            give_back(cur.owned);
            cur = std::move(next);
        }
        prof.sweep += seconds_since(t0);
        return cur;
    }

    // ---- the black box, timed
    void apply(char op, int nvec, const DataType* x, DataType* y) {
        const auto t0 = clock::now();
        bb_(op, nvec, x, y);
        prof.matvec += seconds_since(t0);
        ++prof.calls;
        prof.columns += nvec;
    }
    // ||op(A) 1/sqrt(n)||: the scale of the host's norm_tol (on several
    // ranks, over this rank's rows, then summed)
    double norm_of(char op) {
        double& cached = op == 'T' ? norm_t_ : norm_n_;
        if (cached >= 0.0) return cached;
        std::vector<char*> own;
        const int64_t m = dist() ? c_.M_loc : c_.M, n = dist() ? c_.N_loc : c_.N;
        const int64_t nin = op == 'T' ? m : n, nout = op == 'T' ? n : m;
        DataType* x = take(nin, own);
        DataType* y = take(nout, own);
        launch_fill_value(x, nin, DataType(1.0 / std::sqrt(static_cast<double>(op == 'T' ? c_.M : c_.N))), stream());
        apply(op, 1, x, y);
        std::vector<DataType> h(static_cast<size_t>(nout));
        fmm::gpu::check_cuda(cudaMemcpyAsync(h.data(), y, h.size() * sizeof(DataType), cudaMemcpyDeviceToHost,
                                             stream()), "bf randomized norm");
        fmm::gpu::check_cuda(cudaStreamSynchronize(stream()), "bf randomized norm");
        give_back(own);
        double s = 0.0;
        for (DataType v : h) s += v * v;
        cached = std::sqrt(all_sum(s));
        return cached;
    }

    // BF_Test_Reconstruction_Error: the black box and the new butterfly on a
    // vector of ones, ||BF x - A x|| / ||A x||
    void test_error() {
        if (dist()) {
            test_error_dist();
            return;
        }
        const auto t0 = clock::now();
        const int L = c_.L;
        std::vector<char*> own;
        DataType* x = take(c_.N, own);
        DataType* yref = take(c_.M, own);
        DataType* y = take(c_.M, own);
        launch_fill_value(x, c_.N, DataType(1), stream());
        fmm::gpu::check_cuda(cudaMemsetAsync(y, 0, static_cast<size_t>(c_.M) * sizeof(DataType), stream()),
                             "bf randomized test");
        apply('N', 1, x, yref);
        VLev v = sweep_n(x, 1, L, Mask());
        std::vector<GemmItem> items;
        for (int a = 0; a < (1 << L); ++a) {
            const RBlock& u = lv_[L + 1].at(a, 0);
            if (u.rows == 0 || u.cols == 0 || v.rows(a, 0) == 0) continue;
            items.push_back(GemmItem{u.ptr, v.at(a, 0), y + c_.row_off[a], u.rows, 1, u.cols, u.rows, v.rows(a, 0),
                                     static_cast<int>(c_.M)});
        }
        gemms(items, DataType(0), false);
        std::vector<DataType> h(static_cast<size_t>(c_.M)), href(h.size());
        fmm::gpu::check_cuda(cudaMemcpyAsync(h.data(), y, h.size() * sizeof(DataType), cudaMemcpyDeviceToHost,
                                             stream()), "bf randomized test");
        fmm::gpu::check_cuda(cudaMemcpyAsync(href.data(), yref, href.size() * sizeof(DataType),
                                             cudaMemcpyDeviceToHost, stream()), "bf randomized test");
        fmm::gpu::check_cuda(cudaStreamSynchronize(stream()), "bf randomized test");
        give_back(v.owned);
        give_back(own);
        double d = 0.0, r = 0.0;
        for (size_t i = 0; i < h.size(); ++i) {
            d += (h[i] - href[i]) * (h[i] - href[i]);
            r += href[i] * href[i];
        }
        constexpr double kSafeUnderflow = 2.2250738585072014e-308;
        prof.error = (d < kSafeUnderflow && r < kSafeUnderflow) ? 0.0 : std::sqrt(d) / std::sqrt(r);
        prof.test += seconds_since(t0);
    }

    // ---- where a range finder's Q goes
    struct Target {
        int level = 0;
        int i0 = 0, j0 = 0, i1 = -1, j1 = -1;  // the second block (i1 < 0: none)
        bool transposed = false;              // the blocks are Q^T (rank x rows)
        int cap = -1;
    };

    // The range of each panel (tol_Rdetect), into its target's blocks.  All
    // targets lie in one level, whose blocks share one allocation.
    void range_batch(std::vector<PanelItem>& panels, const std::vector<Target>& targets, double norm_tol) {
        const auto t0 = clock::now();
        const int count = static_cast<int>(panels.size());
        if (count == 0) return;
        int64_t elems = 0, piv = 0, taus = 0;
        int max_n = 0;
        std::vector<int64_t> at(static_cast<size_t>(count)), pat(at.size()), tat(at.size());
        for (int p = 0; p < count; ++p) {
            const int m = panels[p].m0 + panels[p].m1, n = panels[p].ncols;
            at[p] = elems;
            pat[p] = piv;
            tat[p] = taus;
            elems += static_cast<int64_t>(m) * n;
            piv += n;
            taus += std::max(1, std::min(m, n));
            max_n = std::max(max_n, n);
        }
        // one block of scratch: the norms and ranks first (one copy back)
        const size_t wbytes = fmm::gpu::qrcp_work_bytes<DataType>(count, max_n);
        Layout lay;
        const size_t o_norm = lay.template add<double>(count), o_rank = lay.template add<int>(count), o_flag = lay.template add<int>(count);
        const size_t o_mats = lay.template add<DataType>(elems), o_jpvt = lay.template add<int>(piv), o_tau = lay.template add<DataType>(taus);
        const size_t o_work = lay.template add<char>(static_cast<int64_t>(wbytes));
        std::vector<char*> own;
        char* base = take_bytes(lay.size, own);
        double* norm = reinterpret_cast<double*>(base + o_norm);
        int* rank = reinterpret_cast<int*>(base + o_rank);
        int* flag = reinterpret_cast<int*>(base + o_flag);
        DataType* mats = reinterpret_cast<DataType*>(base + o_mats);
        int* jpvt = reinterpret_cast<int*>(base + o_jpvt);
        DataType* tau = reinterpret_cast<DataType*>(base + o_tau);
        void* work = wbytes > 0 ? static_cast<void*>(base + o_work) : nullptr;
        for (int p = 0; p < count; ++p) {
            panels[p].dst = mats + at[p];
            panels[p].ldd = std::max(panels[p].m0 + panels[p].m1, 1);
            panels[p].transposed = 0;
        }
        std::vector<fmm::gpu::QrcpItemT<DataType>> q(static_cast<size_t>(count));
        for (int p = 0; p < count; ++p) {
            const int m = panels[p].m0 + panels[p].m1;
            q[p] = fmm::gpu::QrcpItemT<DataType>{mats + at[p], m, panels[p].ncols, std::max(m, 1), jpvt + pat[p],
                                                 rank + p, norm + p, flag + p};
            q[p].tau = tau + tat[p];
        }
        Image img = begin_image();
        const size_t op = img.meta->append(panels), oq = img.meta->append(q);
        char* d = img.meta->upload(*img.device, stream());
        launch_panels(reinterpret_cast<const PanelItem*>(d + op), count, nullptr, stream());
        fmm::gpu::launch_qrcp(reinterpret_cast<const fmm::gpu::QrcpItemT<DataType>*>(d + oq), count, max_n,
                              c_.tol_rdetect, c_.abs_tol, work, stream());
        std::vector<char> back(o_rank + static_cast<size_t>(count) * sizeof(int));
        fmm::gpu::check_cuda(cudaMemcpyAsync(back.data(), base, back.size(), cudaMemcpyDeviceToHost, stream()),
                             "bf randomized ranks");
        fmm::gpu::check_cuda(cudaStreamSynchronize(stream()), "bf randomized ranks");
        const double* hnorm = reinterpret_cast<const double*>(back.data() + o_norm);
        const int* hrank = reinterpret_cast<const int*>(back.data() + o_rank);
        // the host's ComputeRange: a negligible panel (norm_tol) or rank 0
        // keeps one zero column, and the rank guess caps the outer levels
        std::vector<RangeItem> items(static_cast<size_t>(count));
        int64_t qsize = 0, osize = 0;
        for (int p = 0; p < count; ++p) {
            const int m = panels[p].m0 + panels[p].m1, n = panels[p].ncols;
            int r = hrank[p];
            bool zero = false;
            if (norm_tol > 0.0 && hnorm[p] / std::sqrt(static_cast<double>(n)) < norm_tol) {
                zero = true;
                r = 1;
            }
            if (r == 0) {
                zero = true;
                r = 1;
            }
            if (targets[p].cap > 0 && r > targets[p].cap) r = targets[p].cap;
            RangeItem& it = items[p];
            it.a = mats + at[p];
            it.tau = tau + tat[p];
            it.m = m;
            it.rank = r;
            it.m0 = panels[p].m0;
            it.transposed = targets[p].transposed ? 1 : 0;
            it.zero = zero ? 1 : 0;
            qsize += static_cast<int64_t>(m) * r;
            osize += static_cast<int64_t>(m) * r;
        }
        DataType* qs = take(qsize, own);
        RLevel& v = lv_[targets[0].level];
        DataType* out = take(osize, v.owned);
        int64_t qat = 0, oat = 0;
        for (int p = 0; p < count; ++p) {
            RangeItem& it = items[p];
            const Target& t = targets[p];
            it.q = qs + qat;
            qat += static_cast<int64_t>(it.m) * it.rank;
            const int m1 = it.m - it.m0;
            RBlock b0, b1;
            b0.rows = t.transposed ? it.rank : it.m0;
            b0.cols = t.transposed ? it.m0 : it.rank;
            b0.ptr = out + oat;
            oat += static_cast<int64_t>(it.m0) * it.rank;
            it.out0 = b0.ptr;
            v.at(t.i0, t.j0) = b0;
            if (t.i1 >= 0) {
                b1.rows = t.transposed ? it.rank : m1;
                b1.cols = t.transposed ? m1 : it.rank;
                b1.ptr = out + oat;
                it.out1 = b1.ptr;
                v.at(t.i1, t.j1) = b1;
            }
            oat += static_cast<int64_t>(m1) * it.rank;
        }
        Image img2 = begin_image();
        const size_t oi = img2.meta->append(items);
        char* d2 = img2.meta->upload(*img2.device, stream());
        launch_range(reinterpret_cast<const RangeItem*>(d2 + oi), count, stream());
        // the scratch goes back with the range kernel still queued: whatever
        // reuses it runs after, on the same stream
        give_back(own);
        prof.range += seconds_since(t0);
    }

    // The least-squares blocks of the level where the halves meet:
    // a_panels give A^T's windows (r x nvs, transposed into nvs x r),
    // b_panels B's (the 'T' sweep's pair, transposed into nvs x (m0 + m1)).
    void solve_batch(std::vector<PanelItem>& a_panels, std::vector<PanelItem>& b_panels,
                     const std::vector<Target>& targets) {
        const auto t0 = clock::now();
        const int count = static_cast<int>(a_panels.size());
        if (count == 0) return;
        int64_t ae = 0, be = 0, piv = 0, taus = 0, osize = 0;
        int max_n = 0;
        std::vector<int64_t> aat(static_cast<size_t>(count)), bat(aat.size()), pat(aat.size()), tat(aat.size());
        for (int p = 0; p < count; ++p) {
            const int mm = a_panels[p].ncols, r = a_panels[p].m0, n = b_panels[p].m0 + b_panels[p].m1;
            aat[p] = ae;
            bat[p] = be;
            pat[p] = piv;
            tat[p] = taus;
            ae += static_cast<int64_t>(mm) * r;
            be += static_cast<int64_t>(mm) * n;
            piv += r;
            taus += std::max(1, std::min(mm, r));
            osize += static_cast<int64_t>(n) * r;
            max_n = std::max(max_n, r);
        }
        const size_t wbytes = fmm::gpu::qrcp_work_bytes<DataType>(count, max_n);
        Layout lay;
        const size_t o_a = lay.template add<DataType>(ae), o_b = lay.template add<DataType>(be), o_amax = lay.template add<double>(count);
        const size_t o_jpvt = lay.template add<int>(piv), o_rank = lay.template add<int>(count), o_flag = lay.template add<int>(count);
        const size_t o_norm = lay.template add<double>(count), o_tau = lay.template add<DataType>(taus);
        const size_t o_work = lay.template add<char>(static_cast<int64_t>(wbytes));
        std::vector<char*> own;
        char* base = take_bytes(lay.size, own);
        DataType* A = reinterpret_cast<DataType*>(base + o_a);
        DataType* B = reinterpret_cast<DataType*>(base + o_b);
        double* amax = reinterpret_cast<double*>(base + o_amax);
        int* jpvt = reinterpret_cast<int*>(base + o_jpvt);
        int* rank = reinterpret_cast<int*>(base + o_rank);
        int* flag = reinterpret_cast<int*>(base + o_flag);
        double* norm = reinterpret_cast<double*>(base + o_norm);
        DataType* tau = reinterpret_cast<DataType*>(base + o_tau);
        void* work = wbytes > 0 ? static_cast<void*>(base + o_work) : nullptr;
        for (int p = 0; p < count; ++p) {
            a_panels[p].dst = A + aat[p];
            a_panels[p].ldd = std::max(a_panels[p].ncols, 1);
            a_panels[p].transposed = 1;
            b_panels[p].dst = B + bat[p];
            b_panels[p].ldd = std::max(b_panels[p].ncols, 1);
            b_panels[p].transposed = 1;
        }
        std::vector<fmm::gpu::QrcpItemT<DataType>> q(static_cast<size_t>(count));
        for (int p = 0; p < count; ++p) {
            q[p] = fmm::gpu::QrcpItemT<DataType>{A + aat[p], a_panels[p].ncols, a_panels[p].m0,
                                                 std::max(a_panels[p].ncols, 1), jpvt + pat[p], rank + p, norm + p,
                                                 flag + p};
            q[p].tau = tau + tat[p];
        }
        // the blocks, zeroed: the solve writes only the rows of the columns it keeps
        RLevel& v = lv_[targets[0].level];
        DataType* out = take(osize, v.owned);
        fmm::gpu::check_cuda(cudaMemsetAsync(out, 0, static_cast<size_t>(std::max<int64_t>(osize, 1)) * sizeof(DataType), stream()),
                             "bf randomized solve");
        std::vector<SolveItem> items(static_cast<size_t>(count));
        int64_t oat = 0;
        for (int p = 0; p < count; ++p) {
            const Target& t = targets[p];
            const int r = a_panels[p].m0, m0 = b_panels[p].m0, m1 = b_panels[p].m1;
            RBlock b0{out + oat, m0, r};
            oat += static_cast<int64_t>(m0) * r;
            RBlock b1{out + oat, m1, r};
            oat += static_cast<int64_t>(m1) * r;
            v.at(t.i0, t.j0) = b0;
            v.at(t.i1, t.j1) = b1;
            SolveItem& it = items[p];
            it.a = A + aat[p];
            it.tau = tau + tat[p];
            it.jpvt = jpvt + pat[p];
            it.k = rank + p;
            it.scale = amax + p;
            it.mm = a_panels[p].ncols;
            it.r = r;
            it.b = B + bat[p];
            it.n = m0 + m1;
            it.m0 = m0;
            it.out0 = b0.ptr;
            it.out1 = b1.ptr;
        }
        Image img = begin_image();
        const size_t oa = img.meta->append(a_panels), ob = img.meta->append(b_panels);
        const size_t oq = img.meta->append(q), os = img.meta->append(items);
        char* d = img.meta->upload(*img.device, stream());
        launch_panels(reinterpret_cast<const PanelItem*>(d + oa), count, amax, stream());
        launch_panels(reinterpret_cast<const PanelItem*>(d + ob), count, nullptr, stream());
        fmm::gpu::launch_qrcp(reinterpret_cast<const fmm::gpu::QrcpItemT<DataType>*>(d + oq), count, max_n,
                              c_.tol_ls, 0.0, work, stream());
        launch_solve(reinterpret_cast<const SolveItem*>(d + os), count, stream());
        give_back(own);
        prof.range += seconds_since(t0);
    }

    // ---- the levels.  nth counts from 1, as the host's bundles do.
    static int nint(double x) { return static_cast<int>(std::lround(x)); }

    void left_level(int l, int cap) {
        if (dist()) {
            left_level_dist(l, cap);
            return;
        }
        const int L = c_.L, lh = c_.level_half;
        const int nsub = nint(std::pow(2.0, std::ceil((L - 1) / 2.0)) / std::pow(2.0, lh - l));
        const int ng = (1 << L) / nsub;
        const int nbind = std::min(c_.nbundle, nsub);
        const int num_col = 1 << (L - l + 1);
        double norm_tol = 0.0;
        if (l == 0) {
            const double norm = norm_of('T');
            norm_tol = std::sqrt(norm * norm / static_cast<double>(1 << L)) * c_.tol_rdetect / 1e2;
        }
        for (int s = 0; s < nsub / nbind; ++s) {
            const int nth_s = s * nbind + 1, nth_e = (s + 1) * nbind;
            auto block_row = [&](int nth) {
                return l == 0 ? 0 : static_cast<int>(std::ceil(nth * ng * 2.0 / num_col)) - 1;
            };
            // samples per nth: the rank guess at V, else the blocks' panel heights
            int nvs = 0;
            if (l == 0) {
                nvs = cap + 5;
            } else {
                for (int nth = nth_s; nth <= nth_e; ++nth) {
                    const int i = block_row(nth);
                    for (int j = 0; j < (1 << (L - l)); ++j) {
                        nvs = std::max(nvs, vec_height_n(l, i >> 1, 2 * j) + vec_height_n(l, i >> 1, 2 * j + 1));
                    }
                }
            }
            const int nvec = nvs * (nth_e - nth_s + 1);
            std::vector<char*> own;
            DataType* in = take(c_.M * nvec, own);
            DataType* out = take(c_.N * nvec, own);
            fmm::gpu::check_cuda(cudaMemsetAsync(in, 0, static_cast<size_t>(c_.M) * nvec * sizeof(DataType), stream()),
                                 "bf randomized vectors");
            for (int nth = nth_s; nth <= nth_e; ++nth) {
                const int64_t r0 = c_.row_off[static_cast<size_t>((nth - 1) * ng)];
                const int64_t r1 = c_.row_off[static_cast<size_t>(nth * ng)];
                launch_fill_uniform(in, c_.M, r0, r1 - r0, static_cast<int64_t>(nth - nth_s) * nvs, nvs,
                                    next_seed(), 0, stream());
            }
            apply('T', nvec, in, out);
            std::vector<PanelItem> panels;
            std::vector<Target> targets;
            if (l == 0) {
                for (int b = 0; b < (1 << L); ++b) {
                    PanelItem p;
                    p.s0 = out + c_.col_off[b];
                    p.ld0 = static_cast<int>(c_.N);
                    p.m0 = static_cast<int>(c_.col_off[b + 1] - c_.col_off[b]);
                    p.col0 = 0;
                    p.ncols = nvs;
                    panels.push_back(p);
                    targets.push_back(Target{0, 0, b, -1, -1, false, cap});
                }
                range_batch(panels, targets, norm_tol);
            } else {
                // the bundle reads the block rows i/2 of vec(l), 2^(l-1) x 2^(L-l+1)
                const int nc = 1 << (L - l + 1);
                Mask want(static_cast<size_t>(1) << L, 0);
                for (int nth = nth_s; nth <= nth_e; ++nth) {
                    const int i = block_row(nth);
                    std::fill(want.begin() + static_cast<int64_t>(i >> 1) * nc, want.begin() + static_cast<int64_t>((i >> 1) + 1) * nc, 1);
                }
                VLev v = sweep_n(out, nvec, l - 1, want);
                for (int nth = nth_s; nth <= nth_e; ++nth) {
                    const int i = block_row(nth);
                    for (int j = 0; j < (1 << (L - l)); ++j) {
                        PanelItem p;
                        p.s0 = v.at(i >> 1, 2 * j);
                        p.ld0 = v.rows(i >> 1, 2 * j);
                        p.m0 = p.ld0;
                        p.s1 = v.at(i >> 1, 2 * j + 1);
                        p.ld1 = v.rows(i >> 1, 2 * j + 1);
                        p.m1 = p.ld1;
                        p.col0 = (nth - nth_s) * nvs;
                        p.ncols = nvs;
                        panels.push_back(p);
                        targets.push_back(Target{l, i, 2 * j, i, 2 * j + 1, true, -1});
                    }
                }
                range_batch(panels, targets, 0.0);
                give_back(v.owned);
            }
            give_back(own);
        }
    }

    void right_level(int l, int cap) {
        if (dist()) {
            right_level_dist(l, cap);
            return;
        }
        const int L = c_.L, lh = c_.level_half, lls = lh + 1;
        const int nsub = nint(std::pow(2.0, std::ceil(L / 2.0)) / std::pow(2.0, l - lls));
        const int ng = (1 << L) / nsub;
        const int nbind = std::min(c_.nbundle, nsub);
        double norm_tol = 0.0;
        if (l == L + 1) {
            const double norm = norm_of('N');
            norm_tol = std::sqrt(norm * norm / static_cast<double>(1 << L)) * c_.tol_rdetect / 1e2;
        }
        for (int s = 0; s < nsub / nbind; ++s) {
            const int nth_s = s * nbind + 1, nth_e = (s + 1) * nbind;
            auto block_col = [&](int nth) {
                return l == L + 1 ? 0 : static_cast<int>(std::ceil(nth * ng * 2.0 / (1 << l))) - 1;
            };
            // samples per nth: the rank guess at U, else the blocks' panel
            // heights, and where the halves meet the left half's ranks too
            int nvs = 0;
            if (l == L + 1) {
                nvs = cap + 5;
            } else {
                for (int nth = nth_s; nth <= nth_e; ++nth) {
                    const int J = block_col(nth);
                    for (int I = 0; I < (1 << (l - 1)); ++I) {
                        nvs = std::max(nvs, vec_height_t(l + 1, 2 * I, J >> 1) + vec_height_t(l + 1, 2 * I + 1, J >> 1));
                        if (l == lls) nvs = std::max(nvs, vec_height_n(l, I, J));
                    }
                }
            }
            const int nvec = nvs * (nth_e - nth_s + 1);
            std::vector<char*> own;
            DataType* in = take(c_.N * nvec, own);
            DataType* out = take(c_.M * nvec, own);
            fmm::gpu::check_cuda(cudaMemsetAsync(in, 0, static_cast<size_t>(c_.N) * nvec * sizeof(DataType), stream()),
                                 "bf randomized vectors");
            for (int nth = nth_s; nth <= nth_e; ++nth) {
                const int64_t c0 = c_.col_off[static_cast<size_t>((nth - 1) * ng)];
                const int64_t c1 = c_.col_off[static_cast<size_t>(nth * ng)];
                launch_fill_uniform(in, c_.N, c0, c1 - c0, static_cast<int64_t>(nth - nth_s) * nvs, nvs, next_seed(),
                                    0, stream());
            }
            apply('N', nvec, in, out);
            std::vector<PanelItem> panels;
            std::vector<Target> targets;
            if (l == L + 1) {
                for (int a = 0; a < (1 << L); ++a) {
                    PanelItem p;
                    p.s0 = out + c_.row_off[a];
                    p.ld0 = static_cast<int>(c_.M);
                    p.m0 = static_cast<int>(c_.row_off[a + 1] - c_.row_off[a]);
                    p.col0 = 0;
                    p.ncols = nvs;
                    panels.push_back(p);
                    targets.push_back(Target{L + 1, a, 0, -1, -1, false, cap});
                }
                range_batch(panels, targets, norm_tol);
                give_back(own);
                continue;
            }
            // the bundle reads the block columns J/2 of the grid after
            // Kerl(l+1)^T, 2^l x 2^(L-l) (after U^T, 2^L x 1)
            Mask want_t(static_cast<size_t>(1) << L, 0);
            {
                const int nc = l == L ? 1 : (1 << (L - l));
                for (int nth = nth_s; nth <= nth_e; ++nth) {
                    const int J = block_col(nth);
                    for (int I = 0; I < (1 << l); ++I) want_t[static_cast<size_t>(I) * nc + (J >> 1)] = 1;
                }
            }
            VLev vt = sweep_t(out, nvec, l + 1, want_t);
            if (l == lls) {
                // the halves meet: the left half's sweep of the same random
                // vectors gives A - its block columns J of vec(l), 2^(l-1) x
                // 2^(L-l+1), a cone of the leaves the bundle's vectors live on
                // - and the right half's gives B
                Mask want_n(static_cast<size_t>(1) << L, 0);
                const int nc = 1 << (L - l + 1);
                for (int nth = nth_s; nth <= nth_e; ++nth) {
                    const int J = block_col(nth);
                    for (int I = 0; I < (1 << (l - 1)); ++I) want_n[static_cast<size_t>(I) * nc + J] = 1;
                }
                VLev vn = sweep_n(in, nvec, l - 1, want_n);
                std::vector<PanelItem> a_panels;
                for (int nth = nth_s; nth <= nth_e; ++nth) {
                    const int J = block_col(nth);
                    for (int I = 0; I < (1 << (l - 1)); ++I) {
                        PanelItem a;
                        a.s0 = vn.at(I, J);
                        a.ld0 = vn.rows(I, J);
                        a.m0 = a.ld0;
                        a.col0 = (nth - nth_s) * nvs;
                        a.ncols = nvs;
                        a_panels.push_back(a);
                        PanelItem b;
                        b.s0 = vt.at(2 * I, J >> 1);
                        b.ld0 = vt.rows(2 * I, J >> 1);
                        b.m0 = b.ld0;
                        b.s1 = vt.at(2 * I + 1, J >> 1);
                        b.ld1 = vt.rows(2 * I + 1, J >> 1);
                        b.m1 = b.ld1;
                        b.col0 = a.col0;
                        b.ncols = nvs;
                        panels.push_back(b);
                        targets.push_back(Target{l, 2 * I, J, 2 * I + 1, J, false, -1});
                    }
                }
                solve_batch(a_panels, panels, targets);
                give_back(vn.owned);
            } else {
                for (int nth = nth_s; nth <= nth_e; ++nth) {
                    const int J = block_col(nth);
                    for (int I = 0; I < (1 << (l - 1)); ++I) {
                        PanelItem p;
                        p.s0 = vt.at(2 * I, J >> 1);
                        p.ld0 = vt.rows(2 * I, J >> 1);
                        p.m0 = p.ld0;
                        p.s1 = vt.at(2 * I + 1, J >> 1);
                        p.ld1 = vt.rows(2 * I + 1, J >> 1);
                        p.m1 = p.ld1;
                        p.col0 = (nth - nth_s) * nvs;
                        p.ncols = nvs;
                        panels.push_back(p);
                        targets.push_back(Target{l, 2 * I, J, 2 * I + 1, J, false, -1});
                    }
                }
                range_batch(panels, targets, 0.0);
            }
            give_back(vt.owned);
            give_back(own);
        }
    }

    // ---- the levels on several ranks: BF_Reconstruction_LL / RR and the
    // host's resolving routines (BF_Resolving_Butterfly_LL_dat,
    // BF_OneBlock_LL and their RR mirrors) over this rank's blocks.  Global
    // block indices are one based here, as there.

    // The random vectors of bundles nth_s .. nth_e, nvs columns each, into x
    // (rows x nvec): every bundle's rows of this rank's leaves among the
    // bundle's ng.  off: the local leaves' ends; first: the global index,
    // zero based, of the first local leaf; hash0: the global index of x's
    // first row.  A seed per bundle on every rank, as on one process.
    void fill_bundles(DataType* x, int64_t rows, int nvec, const std::vector<int64_t>& off, int first,
                      int64_t hash0, int nth_s, int nth_e, int ng, int nvs) {
        fmm::gpu::check_cuda(cudaMemsetAsync(x, 0, static_cast<size_t>(std::max<int64_t>(rows * nvec, 1)) * sizeof(DataType),
                                             stream()),
                             "bf randomized vectors");
        const int nleaf = static_cast<int>(off.size()) - 1;
        for (int nth = nth_s; nth <= nth_e; ++nth) {
            const uint64_t seed = next_seed();
            const int a0 = std::max((nth - 1) * ng - first, 0), a1 = std::min(nth * ng - first, nleaf);
            if (a1 <= a0) continue;
            launch_fill_uniform(x, std::max<int64_t>(rows, 1), off[static_cast<size_t>(a0)],
                                off[static_cast<size_t>(a1)] - off[static_cast<size_t>(a0)],
                                static_cast<int64_t>(nth - nth_s) * nvs, nvs, seed, hash0, stream());
        }
    }

    // a partial sweep, timed: the level it ends in, its data the caller's to free
    VecLevel& sweep_part(SweepPlan& plan, int nvec, const DataType* x, int64_t ldx) {
        const auto t0 = clock::now();
        VecLevel& v = sweep_partial(dev_, plan, nvec, x, static_cast<int>(std::max<int64_t>(ldx, 1)), nullptr, stream());
        prof.sweep += seconds_since(t0);
        return v;
    }
    static void free_level(VecLevel& v) {
        if (v.data != nullptr) fmm::gpu::DeviceHeap::instance().free(v.data);
        v.data = nullptr;
    }

    void left_level_dist(int l, int cap) {
        const int L = c_.L, lh = c_.level_half;
        const int nsub = nint(std::pow(2.0, std::ceil((L - 1) / 2.0)) / std::pow(2.0, lh - l));
        const int ng = (1 << L) / nsub;
        const int nbind = std::min(c_.nbundle, nsub);
        const int num_col = 1 << (L - l + 1);
        double norm_tol = 0.0;
        if (l == 0) {
            const double norm = norm_of('T');
            norm_tol = std::sqrt(norm * norm / static_cast<double>(1 << L)) * c_.tol_rdetect / 1e2;
        }
        // vec(l): the 'N' sweep to level l-1, with its pair-wise broadcasts
        SweepPlan* plan = l >= 1 ? &partial_plan('N', l - 1) : nullptr;
        const SweepRanges& g = c_.layout.kernel_r;  // the host's loop grid ('R' ordering)
        const size_t lu = static_cast<size_t>(l);
        const RLevel& K = lv_[static_cast<size_t>(l)];
        // block row index_i of Kerl(l), if this rank resolves it
        auto owns_row = [&](int index_i) { return g.nc[lu] > 0 && local_of(index_i, g.idx_r[lu], g.inc_r[lu], g.nr[lu]) >= 0; };
        // the pair of vec(l)'s blocks the loop grid's column jl reads for row
        // index_i: (index_i+1)/2, 2 index_j - 1, and the next local block
        auto pair_of = [&](const VecLevel& v, int index_i, int jl) {
            const int index_j = g.idx_c[lu] + jl * g.inc_c[lu];
            const size_t b = static_cast<size_t>(v.row_local((index_i + 1) / 2)) * v.nc + v.col_local(2 * index_j - 1);
            return std::make_pair(b, b + 1);
        };
        for (int s = 0; s < nsub / nbind; ++s) {
            const int nth_s = s * nbind + 1, nth_e = (s + 1) * nbind;
            auto block_row = [&](int nth) { return l == 0 ? 1 : static_cast<int>(std::ceil(nth * ng * 2.0 / num_col)); };
            int64_t nvs = 0;
            if (l == 0) {
                nvs = cap + 5;
            } else {
                const VecLevel& v = plan->vec[static_cast<size_t>(plan->result)];
                for (int nth = nth_s; nth <= nth_e; ++nth) {
                    const int i = block_row(nth);
                    if (!owns_row(i)) continue;
                    for (int jl = 0; jl < g.nc[lu]; ++jl) {
                        const auto b = pair_of(v, i, jl);
                        nvs = std::max<int64_t>(nvs, v.rows[b.first] + v.rows[b.second]);
                    }
                }
                nvs = all_max(nvs);
            }
            const int nvec = static_cast<int>(nvs) * (nth_e - nth_s + 1);
            std::vector<char*> own;
            DataType *in = nullptr, *out = nullptr;
            bundle_room(c_.M_loc * nvec, c_.N_loc * nvec, in, out, own);
            fill_bundles(in, c_.M_loc, nvec, c_.row_off, lv_[static_cast<size_t>(L) + 1].idx_r - 1, c_.row0, nth_s,
                         nth_e, ng, static_cast<int>(nvs));
            apply('T', nvec, in, out);
            std::vector<PanelItem> panels;
            std::vector<Target> targets;
            if (l == 0) {
                for (int b = 0; b < K.nc; ++b) {
                    PanelItem p;
                    p.s0 = out + c_.col_off[static_cast<size_t>(b)];
                    p.ld0 = static_cast<int>(std::max<int64_t>(c_.N_loc, 1));
                    p.m0 = static_cast<int>(c_.col_off[static_cast<size_t>(b) + 1] - c_.col_off[static_cast<size_t>(b)]);
                    p.col0 = 0;
                    p.ncols = static_cast<int>(nvs);
                    panels.push_back(p);
                    targets.push_back(Target{0, 0, b, -1, -1, false, cap});
                }
                range_batch(panels, targets, norm_tol);
            } else {
                VecLevel& v = sweep_part(*plan, nvec, out, c_.N_loc);
                for (int nth = nth_s; nth <= nth_e; ++nth) {
                    const int i = block_row(nth);
                    if (!owns_row(i)) continue;
                    const int ik = local_of(i, K.idx_r, K.inc_r, K.nr);
                    for (int jl = 0; jl < g.nc[lu]; ++jl) {
                        const int index_j = g.idx_c[lu] + jl * g.inc_c[lu];
                        const int jk = local_of(2 * index_j - 1, K.idx_c, K.inc_c, K.nc);
                        if (ik < 0 || jk < 0) throw std::runtime_error("bf::gpu::Randomized: a kernel block this rank resolves is not local");
                        const auto b = pair_of(v, i, jl);
                        PanelItem p;
                        p.s0 = v.data + v.offset[b.first];
                        p.ld0 = v.rows[b.first];
                        p.m0 = p.ld0;
                        p.s1 = v.data + v.offset[b.second];
                        p.ld1 = v.rows[b.second];
                        p.m1 = p.ld1;
                        p.col0 = (nth - nth_s) * static_cast<int>(nvs);
                        p.ncols = static_cast<int>(nvs);
                        panels.push_back(p);
                        targets.push_back(Target{l, ik, jk, ik, jk + 1, true, -1});
                    }
                }
                range_batch(panels, targets, 0.0);
                free_level(v);
            }
            give_back(own);
        }
    }

    void right_level_dist(int l, int cap) {
        const int L = c_.L, lh = c_.level_half, lls = lh + 1;
        const int nsub = nint(std::pow(2.0, std::ceil(L / 2.0)) / std::pow(2.0, l - lls));
        const int ng = (1 << L) / nsub;
        const int nbind = std::min(c_.nbundle, nsub);
        double norm_tol = 0.0;
        if (l == L + 1) {
            const double norm = norm_of('N');
            norm_tol = std::sqrt(norm * norm / static_cast<double>(1 << L)) * c_.tol_rdetect / 1e2;
        }
        // the 'T' sweep from U to level l+1, and where the halves meet the
        // 'N' sweep to level_half with the midpoint redistribution
        SweepPlan* pt = l <= L ? &partial_plan('T', l + 1) : nullptr;
        SweepPlan* pn = l == lls ? &partial_plan('N', lh) : nullptr;
        const SweepRanges& g = c_.layout.kernel_c;  // the host's loop grid ('C' ordering)
        const size_t lu = static_cast<size_t>(l);
        const RLevel& K = lv_[static_cast<size_t>(l)];
        auto owns_col = [&](int index_j) { return g.nr[lu] > 0 && local_of(index_j, g.idx_c[lu], g.inc_c[lu], g.nc[lu]) >= 0; };
        // the pair of the 'T' sweep's blocks the loop grid's row il reads for
        // column index_j: 2 index_i - 1, (index_j+1)/2, and the next local row
        auto pair_of = [&](const VecLevel& v, int il, int index_j) {
            const int index_i = g.idx_r[lu] + il * g.inc_r[lu];
            const size_t b = static_cast<size_t>(v.row_local(2 * index_i - 1)) * v.nc + v.col_local((index_j + 1) / 2);
            return std::make_pair(b, b + static_cast<size_t>(v.nc));
        };
        // the left half's block (index_i, index_j) where the halves meet
        auto left_of = [&](const VecLevel& v, int il, int index_j) {
            const int index_i = g.idx_r[lu] + il * g.inc_r[lu];
            return static_cast<size_t>(v.row_local(index_i)) * v.nc + v.col_local(index_j);
        };
        for (int s = 0; s < nsub / nbind; ++s) {
            const int nth_s = s * nbind + 1, nth_e = (s + 1) * nbind;
            auto block_col = [&](int nth) { return l == L + 1 ? 1 : static_cast<int>(std::ceil(nth * ng * 2.0 / (1 << l))); };
            int64_t nvs = 0;
            if (l == L + 1) {
                nvs = cap + 5;
            } else {
                const VecLevel& vt = pt->vec[static_cast<size_t>(pt->result)];
                for (int nth = nth_s; nth <= nth_e; ++nth) {
                    const int j = block_col(nth);
                    if (!owns_col(j)) continue;
                    for (int il = 0; il < g.nr[lu]; ++il) {
                        const auto b = pair_of(vt, il, j);
                        nvs = std::max<int64_t>(nvs, vt.rows[b.first] + vt.rows[b.second]);
                        if (l == lls) nvs = std::max<int64_t>(nvs, pn->vec[static_cast<size_t>(pn->result)].rows[left_of(pn->vec[static_cast<size_t>(pn->result)], il, j)]);
                    }
                }
                nvs = all_max(nvs);
            }
            const int nvec = static_cast<int>(nvs) * (nth_e - nth_s + 1);
            std::vector<char*> own;
            DataType *in = nullptr, *out = nullptr;
            bundle_room(c_.N_loc * nvec, c_.M_loc * nvec, in, out, own);
            fill_bundles(in, c_.N_loc, nvec, c_.col_off, lv_[0].idx_c - 1, c_.col0, nth_s, nth_e, ng,
                         static_cast<int>(nvs));
            apply('N', nvec, in, out);
            std::vector<PanelItem> panels, a_panels;
            std::vector<Target> targets;
            if (l == L + 1) {
                for (int a = 0; a < K.nr; ++a) {
                    PanelItem p;
                    p.s0 = out + c_.row_off[static_cast<size_t>(a)];
                    p.ld0 = static_cast<int>(std::max<int64_t>(c_.M_loc, 1));
                    p.m0 = static_cast<int>(c_.row_off[static_cast<size_t>(a) + 1] - c_.row_off[static_cast<size_t>(a)]);
                    p.col0 = 0;
                    p.ncols = static_cast<int>(nvs);
                    panels.push_back(p);
                    targets.push_back(Target{L + 1, a, 0, -1, -1, false, cap});
                }
                range_batch(panels, targets, norm_tol);
                give_back(own);
                continue;
            }
            VecLevel& vt = sweep_part(*pt, nvec, out, c_.M_loc);
            VecLevel* vn = l == lls ? &sweep_part(*pn, nvec, in, c_.N_loc) : nullptr;
            for (int nth = nth_s; nth <= nth_e; ++nth) {
                const int j = block_col(nth);
                if (!owns_col(j)) continue;
                const int jk = local_of(j, K.idx_c, K.inc_c, K.nc);
                for (int il = 0; il < g.nr[lu]; ++il) {
                    const int index_i = g.idx_r[lu] + il * g.inc_r[lu];
                    const int ik = local_of(2 * index_i - 1, K.idx_r, K.inc_r, K.nr);
                    if (ik < 0 || jk < 0) throw std::runtime_error("bf::gpu::Randomized: a kernel block this rank resolves is not local");
                    const auto b = pair_of(vt, il, j);
                    PanelItem p;
                    p.s0 = vt.data + vt.offset[b.first];
                    p.ld0 = vt.rows[b.first];
                    p.m0 = p.ld0;
                    p.s1 = vt.data + vt.offset[b.second];
                    p.ld1 = vt.rows[b.second];
                    p.m1 = p.ld1;
                    p.col0 = (nth - nth_s) * static_cast<int>(nvs);
                    p.ncols = static_cast<int>(nvs);
                    panels.push_back(p);
                    targets.push_back(Target{l, ik, jk, ik + 1, jk, false, -1});
                    if (vn != nullptr) {
                        const size_t bn = left_of(*vn, il, j);
                        PanelItem a;
                        a.s0 = vn->data + vn->offset[bn];
                        a.ld0 = vn->rows[bn];
                        a.m0 = a.ld0;
                        a.col0 = p.col0;
                        a.ncols = p.ncols;
                        a_panels.push_back(a);
                    }
                }
            }
            if (vn != nullptr) {
                solve_batch(a_panels, panels, targets);
                free_level(*vn);
            } else {
                range_batch(panels, targets, 0.0);
            }
            free_level(vt);
            give_back(own);
        }
    }

    // BF_Test_Reconstruction_Error on several ranks: the new butterfly's
    // matvec is the distributed sweep's
    void test_error_dist() {
        const auto t0 = clock::now();
        std::vector<char*> own;
        DataType* x = take(c_.N_loc, own);
        DataType* yref = take(c_.M_loc, own);
        DataType* y = take(c_.M_loc, own);
        launch_fill_value(x, c_.N_loc, DataType(1), stream());
        apply('N', 1, x, yref);
        refresh_layout();
        SweepPlan plan;
        build_sweep_plan(c_.bridge, dev_, 'N', 1, plan);
        SweepScratch unused;
        sweep(dev_, plan, 1, x, static_cast<int>(std::max<int64_t>(c_.N_loc, 1)), y,
              static_cast<int>(std::max<int64_t>(c_.M_loc, 1)), DataType(1), DataType(0), unused, stream());
        std::vector<DataType> h(static_cast<size_t>(c_.M_loc)), href(h.size());
        fmm::gpu::check_cuda(cudaMemcpyAsync(h.data(), y, h.size() * sizeof(DataType), cudaMemcpyDeviceToHost,
                                             stream()), "bf randomized test");
        fmm::gpu::check_cuda(cudaMemcpyAsync(href.data(), yref, href.size() * sizeof(DataType),
                                             cudaMemcpyDeviceToHost, stream()), "bf randomized test");
        fmm::gpu::check_cuda(cudaStreamSynchronize(stream()), "bf randomized test");
        plan.release();
        give_back(own);
        double d = 0.0, r = 0.0;
        for (size_t i = 0; i < h.size(); ++i) {
            d += (h[i] - href[i]) * (h[i] - href[i]);
            r += href[i] * href[i];
        }
        d = all_sum(d);
        r = all_sum(r);
        constexpr double kSafeUnderflow = 2.2250738585072014e-308;
        prof.error = (d < kSafeUnderflow && r < kSafeUnderflow) ? 0.0 : std::sqrt(d) / std::sqrt(r);
        prof.test += seconds_since(t0);
    }

    // the heights of vector level t's block (i, j), from the factors built so far
    int vec_height_n(int t, int i, int j) const {
        return t == 1 ? lv_[0].at(0, j).cols : lv_[t - 1].at(i, 2 * j).rows;
    }
    int vec_height_t(int t, int i, int j) const {
        return t == c_.L + 1 ? lv_[t].at(i, 0).cols : lv_[t].at(2 * i, j).cols;
    }
    uint64_t next_seed() { return c_.seed * 0x9e3779b97f4a7c15ull + (++draws_); }

    RandomizedConfig c_;
    BlackBox bb_;
    std::vector<RLevel> lv_;
    // several ranks: the butterfly as the sweeps read it, and their plans
    BfDevice dev_;
    std::map<int, SweepPlan> plans_n_, plans_t_;
    uint64_t draws_ = 0;
    double norm_t_ = -1.0, norm_n_ = -1.0;
};

}  // namespace gpu
}  // namespace bf

#endif  // H2_HAVE_GPU
