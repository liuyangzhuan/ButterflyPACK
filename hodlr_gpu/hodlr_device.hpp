#pragma once
// GPU backend of the HODLR format (format 1 with LRlevel 0: every
// off-diagonal block low rank).  The Fortran HODLR code keeps the tree, the
// process groups and host copies of all blocks; this class holds their device
// copies and runs the batched operations on them.  It reuses the runtime of
// the H2 GPU backend (h2_parallel/color_gpu): the device and MAGMA queue of
// the rank, the device heap, and the variable-size MAGMA batches (VBatch).
//
// Rows are local: row 0 is the first row this rank owns.  A low-rank block
// U V^T (plain transpose, also for complex data) covers rows
// [row0, row0 + m) and columns [col0, col0 + n).  In the symmetric HODLR
// (option%sym=1) only A21 = U1 V0^T of a node is stored, flagged sym, and
// stands for A12 = A21^T as well.
//
// With several ranks, the nodes of the top levels belong to process groups
// (shared levels).  A rank then holds of each block of its node the rows of
// U it owns (m of them, or none) and the rows of V it owns (n, or none):
// the Fortran side moves V to the ranks that own the rows V multiplies.  A
// product sums the partial coefficients V^T x (U^T x) over the node's ranks
// (add_shared_level: its communicator) before applying U (V).

#ifdef H2_HAVE_GPU

#include "color_gpu/device_heap.hpp"
#include "color_gpu/device_scalar.hpp"
#include "color_gpu/gpu_runtime.hpp"
#include "hodlr_kernels.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace bpack {
namespace gpu {

using fmm::gpu::check_cuda;
using fmm::gpu::Context;
using fmm::gpu::DeviceBuffer;
using fmm::gpu::DeviceHeap;
using fmm::gpu::MetaBuilder;
using fmm::gpu::PinnedBuffer;
using fmm::gpu::VBatch;

// Host-to-device copies of many column-major blocks into one contiguous
// device range, through two alternating pinned slots (pinning costs about a
// second per GB, so the slots stay small and live for the run).
class StreamUploader {
public:
    static constexpr size_t kSlot = size_t{64} << 20;

    StreamUploader() {
        for (int i = 0; i < 2; ++i) {
            check_cuda(cudaEventCreateWithFlags(&done_[i], cudaEventDisableTiming), "cudaEventCreate");
        }
    }
    StreamUploader(const StreamUploader&) = delete;
    StreamUploader& operator=(const StreamUploader&) = delete;
    // (the events live until process exit, as the H2 backend's)

    // Start writing at device address dst.
    void begin(char* dst, cudaStream_t stream) {
        dst_ = dst;
        stream_ = stream;
        cur_ = 0;
        used_ = 0;
        slot_base_ = static_cast<char*>(slots_[0].reserve(kSlot));
    }
    // Append the m x n column-major block a (leading dimension lda).
    template<typename T>
    void write_block(const T* a, int64_t m, int64_t n, int64_t lda) {
        if (m == lda) {
            write(a, static_cast<size_t>(m * n) * sizeof(T));
            return;
        }
        for (int64_t j = 0; j < n; ++j) write(a + j * lda, static_cast<size_t>(m) * sizeof(T));
    }
    void write(const void* src, size_t bytes) {
        const char* s = static_cast<const char*>(src);
        while (bytes > 0) {
            const size_t take = std::min(bytes, kSlot - used_);
            std::memcpy(slot_base_ + used_, s, take);
            used_ += take;
            s += take;
            bytes -= take;
            if (used_ == kSlot) flush();
        }
    }
    // Leave the next bytes of the destination as they are (written otherwise).
    void skip(size_t bytes) {
        flush();
        dst_ += bytes;
    }
    // Copy what is left and wait for all copies.
    void finish() {
        flush();
        for (int i = 0; i < 2; ++i) check_cuda(cudaEventSynchronize(done_[i]), "upload");
    }

private:
    void flush() {
        if (used_ == 0) return;
        check_cuda(cudaMemcpyAsync(dst_, slot_base_, used_, cudaMemcpyHostToDevice, stream_), "HODLR upload");
        check_cuda(cudaEventRecord(done_[cur_], stream_), "cudaEventRecord");
        dst_ += used_;
        used_ = 0;
        cur_ ^= 1;
        check_cuda(cudaEventSynchronize(done_[cur_]), "upload");  // the slot is free again
        slot_base_ = static_cast<char*>(slots_[cur_].reserve(kSlot));
    }

    PinnedBuffer slots_[2];
    cudaEvent_t done_[2] = {nullptr, nullptr};
    char* dst_ = nullptr;
    char* slot_base_ = nullptr;
    cudaStream_t stream_ = nullptr;
    size_t used_ = 0;
    int cur_ = 0;
};

inline StreamUploader& stream_uploader() {
    static StreamUploader* uploader = new StreamUploader;  // never freed (CUDA shutdown order)
    return *uploader;
}

template<typename T>
class HodlrDevice {
public:
    struct Leaf {
        int64_t row0;
        int m;
        T* d;              // m x m, leading dimension m
        const T* host;     // until commit_forward
        int host_ld;
    };
    struct LowRank {
        int level;
        int64_t row0, col0;
        int m, n, k;       // m or n is 0 when this rank holds no rows of U or V (shared levels)
        bool sym;
        T* u;              // m x k, leading dimension m
        T* v;              // n x k, leading dimension n
        const T* host_u;   // until commit_forward
        const T* host_v;
        int host_ldu, host_ldv;
        const T* dev_v;    // or V on the device (n x k, leading dimension n), until commit_forward
    };

    HodlrDevice() = default;
    HodlrDevice(const HodlrDevice&) = delete;
    HodlrDevice& operator=(const HodlrDevice&) = delete;
    ~HodlrDevice() {
        release_forward();
        clear_mirrors();
    }

    // Forget the forward blocks: a new matrix of n_loc local rows follows.
    void reset(int64_t n_loc, int maxlevel) {
        release_forward();
        n_loc_ = n_loc;
        maxlevel_ = maxlevel;
        level_comm_.assign(static_cast<size_t>(maxlevel + 1), MPI_COMM_NULL);
    }

    // Level `level` is shared: its node on this rank belongs to the ranks of
    // comm (the coefficients of its blocks are summed over comm).
    void add_shared_level(int level, MPI_Comm comm) {
        if (level < 1 || level > maxlevel_) throw std::invalid_argument("HODLR GPU: shared level out of range");
        level_comm_[static_cast<size_t>(level)] = comm;
    }
    bool shared_level(int level) const { return level_comm_[static_cast<size_t>(level)] != MPI_COMM_NULL; }
    MPI_Comm level_comm(int level) const { return level_comm_[static_cast<size_t>(level)]; }
    bool distributed() const {
        for (MPI_Comm c : level_comm_) {
            if (c != MPI_COMM_NULL) return true;
        }
        return false;
    }

    // A device copy (bytes) of a host array that add_lowrank will record, from
    // the construction on the GPU: commit_forward copies it on the device
    // instead of uploading the host array.  The device memory handed over by
    // own_mirror_memory is freed after the commit (or with the object).
    // stale: the host array was not filled (the device copy is the only one;
    // fetch_block brings it back after the commit)
    void add_mirror(const void* host, const void* dev, size_t bytes, bool stale = false) {
        mirrors_[host] = Mirror{dev, bytes, stale};
    }
    void own_mirror_memory(void* dev) { mirror_memory_.push_back(dev); }
    void clear_mirrors() {
        if (!mirror_memory_.empty()) {
            Context::instance().activate();
            for (void* p : mirror_memory_) DeviceHeap::instance().free(p);
        }
        mirror_memory_.clear();
        mirrors_.clear();
    }

    // Record a block; its host array must live until commit_forward.
    void add_leaf(int64_t row0, int m, const T* d, int ldd) {
        check_rows(row0, m, "leaf");
        leaves_.push_back(Leaf{row0, m, nullptr, d, ldd});
    }
    // (dev_v: V already on the device, n x k with leading dimension n,
    // alive until commit_forward, instead of the host array v)
    void add_lowrank(int level, int64_t row0, int m, int64_t col0, int n, int k, const T* u, int ldu, const T* v,
                     int ldv, bool sym, const T* dev_v = nullptr) {
        check_rows(row0, m, "low-rank block rows");
        check_rows(col0, n, "low-rank block columns");
        if (level < 1 || level > maxlevel_) throw std::invalid_argument("HODLR GPU: block level out of range");
        if (k < 0) throw std::invalid_argument("HODLR GPU: negative rank");
        blocks_.push_back(LowRank{level, row0, col0, m, n, k, sym, nullptr, nullptr, u, v, ldu, ldv, dev_v});
    }
    // the device copy of host (bytes) given by add_mirror, or null
    const void* mirror_ptr(const void* host, size_t bytes) const {
        auto it = mirrors_.find(host);
        return it != mirrors_.end() && it->second.bytes == bytes ? it->second.dev : nullptr;
    }

    // Upload the recorded blocks into one resident device allocation.
    void commit_forward() {
        Context& ctx = Context::instance();
        ctx.activate();
        size_t bytes = 0;
        for (const Leaf& l : leaves_) bytes += bytes_of(static_cast<int64_t>(l.m) * l.m);
        for (const LowRank& b : blocks_) {  // U and V are placed (and aligned) one after the other
            bytes += bytes_of(static_cast<int64_t>(b.m) * b.k) + bytes_of(static_cast<int64_t>(b.n) * b.k);
        }
        DeviceHeap& heap = DeviceHeap::instance();
        if (forward_ != nullptr) heap.free(forward_);
        forward_ = bytes > 0 ? heap.alloc_resident<char>(bytes) : nullptr;
        forward_bytes_ = bytes;
        if (bytes > 0) {  // (the padding of the blocks copied on the device)
            check_cuda(cudaMemsetAsync(forward_, 0, bytes, ctx.stream()), "HODLR forward blocks");
        }
        mirrored_bytes_ = 0;
        // the device copy of a host array of count elements (leading dimension ld, rows m), or null
        auto mirror_of = [&](const T* host, int64_t count, int64_t ld, int64_t m) -> const T* {
            if (host == nullptr || count <= 0) return nullptr;
            auto it = mirrors_.find(host);
            if (it == mirrors_.end()) return nullptr;
            if (ld != m || it->second.bytes != static_cast<size_t>(count) * sizeof(T)) {
                if (it->second.stale) {  // (the host array holds no data to upload instead)
                    throw std::logic_error("HODLR GPU: a block's device copy from the construction does not match it");
                }
                return nullptr;
            }
            return static_cast<const T*>(it->second.dev);
        };

        StreamUploader& up = stream_uploader();
        up.begin(forward_, ctx.stream());
        char* at = forward_;
        auto place = [&](int64_t count) {
            T* p = reinterpret_cast<T*>(at);
            at += bytes_of(count);
            return p;
        };
        // (each placement is followed by its data, padded to the same 256 bytes)
        auto pad = [&](int64_t count) {
            const size_t extra = bytes_of(count) - static_cast<size_t>(count) * sizeof(T);
            static const char zeros[256] = {0};
            if (extra > 0) up.write(zeros, extra);
        };
        for (Leaf& l : leaves_) {
            const int64_t count = static_cast<int64_t>(l.m) * l.m;
            l.d = place(count);
            up.write_block(l.host, l.m, l.m, l.host_ld);
            pad(count);
        }
        for (LowRank& b : blocks_) {
            const int64_t cu = static_cast<int64_t>(b.m) * b.k, cv = static_cast<int64_t>(b.n) * b.k;
            b.u = place(cu);
            if (const T* d = mirror_of(b.host_u, cu, b.host_ldu, b.m)) {
                up.skip(bytes_of(cu));
                check_cuda(cudaMemcpyAsync(b.u, d, static_cast<size_t>(cu) * sizeof(T), cudaMemcpyDeviceToDevice, ctx.stream()),
                           "HODLR forward block (device copy)");
                mirrored_bytes_ += static_cast<size_t>(cu) * sizeof(T);
            } else {
                up.write_block(b.host_u, b.m, b.k, b.host_ldu);
                pad(cu);
            }
            b.v = place(cv);
            const T* dv = b.dev_v != nullptr ? b.dev_v : mirror_of(b.host_v, cv, b.host_ldv, b.n);
            if (cv > 0 && dv != nullptr) {
                const T* d = dv;
                up.skip(bytes_of(cv));
                check_cuda(cudaMemcpyAsync(b.v, d, static_cast<size_t>(cv) * sizeof(T), cudaMemcpyDeviceToDevice, ctx.stream()),
                           "HODLR forward block (device copy)");
                mirrored_bytes_ += static_cast<size_t>(cv) * sizeof(T);
            } else {
                up.write_block(b.host_v, b.n, b.k, b.host_ldv);
                pad(cv);
            }
        }
        up.finish();
        check_cuda(cudaStreamSynchronize(ctx.stream()), "HODLR forward blocks");
        clear_mirrors();
        if (at != forward_ + bytes) throw std::logic_error("HODLR GPU: forward blocks placed outside their allocation");
        committed_ = true;
        if (fmm::env::check("hodlr") || fmm::env::check("hodlr-transpose")) verify_upload();  // BPACK_CHECK
        for (Leaf& l : leaves_) l.host = nullptr;
        for (LowRank& b : blocks_) {
            b.host_u = b.host_v = nullptr;
            b.dev_v = nullptr;
        }
    }

    // Compare the device copies with the host arrays they came from (exactly).
    void verify_upload() const {
        size_t bad_blocks = 0, bad_leaves = 0;
        auto same = [](const T* dev, const T* host, int64_t m, int64_t n, int64_t ld, int64_t& first_col) {
            std::vector<T> tmp(static_cast<size_t>(m * n));
            check_cuda(cudaMemcpy(tmp.data(), dev, tmp.size() * sizeof(T), cudaMemcpyDeviceToHost), "verify upload");
            for (int64_t j = 0; j < n; ++j) {
                if (std::memcmp(tmp.data() + j * m, host + j * ld, static_cast<size_t>(m) * sizeof(T)) != 0) {
                    first_col = j;
                    return false;
                }
            }
            return true;
        };
        for (const Leaf& l : leaves_) {
            int64_t col = -1;
            if (!same(l.d, l.host, l.m, l.m, l.host_ld, col)) ++bad_leaves;
        }
        for (const LowRank& b : blocks_) {
            int64_t cu = -1, cv = -1;
            const bool okU = same(b.u, b.host_u, b.m, b.k, b.host_ldu, cu);
            const bool okV = same(b.v, b.host_v, b.n, b.k, b.host_ldv, cv);
            if (!okU || !okV) {
                if (bad_blocks < 5) {
                    std::printf(" HODLR GPU check (upload): block level %d rows %lld+%d cols %lld+%d rank %d: U from column %lld,"
                                " V from column %lld differ (host ld %d, %d)\n",
                                b.level, static_cast<long long>(b.row0), b.m, static_cast<long long>(b.col0), b.n, b.k,
                                static_cast<long long>(cu), static_cast<long long>(cv), b.host_ldu, b.host_ldv);
                }
                ++bad_blocks;
            }
        }
        std::printf(" HODLR GPU check (upload): %zu of %zu low-rank blocks and %zu of %zu leaves differ from the host\n",
                    bad_blocks, blocks_.size(), bad_leaves, leaves_.size());
        std::fflush(stdout);
    }

    bool forward_ready() const { return committed_; }
    double mult_flops() const { return mult_flops_; }  // of the last mult()
    size_t mirrored_bytes() const { return mirrored_bytes_; }
    size_t forward_bytes() const { return forward_bytes_; }
    int64_t n_loc() const { return n_loc_; }
    int maxlevel() const { return maxlevel_; }
    const char* forward_region() const { return forward_; }
    const std::vector<Leaf>& leaves() const { return leaves_; }
    const std::vector<LowRank>& blocks() const { return blocks_; }
    // Block idx (in the order of add_lowrank) after the commit: its shape, and
    // its U, V copied back to host arrays of leading dimensions ldu, ldv (null:
    // not copied)
    const LowRank& committed_block(int idx) const {
        if (!committed_) throw std::logic_error("HODLR GPU: the forward blocks are not on the device");
        if (idx < 0 || static_cast<size_t>(idx) >= blocks_.size()) throw std::out_of_range("HODLR GPU: block index");
        return blocks_[static_cast<size_t>(idx)];
    }
    // Rows of U (role 0) or V (role 1) of block idx after the commit to the
    // host: out (n x k, leading dimension n) = the rows rows[0..n), numbered
    // as this rank's rows from 0 (the block's own start at row0 or col0);
    // returns the block's rank k
    int gather_rows(int idx, int role, int n, const int* rows, T* out) const {
        const LowRank& b = committed_block(idx);
        const T* a = role == 0 ? b.u : b.v;
        const int ld = role == 0 ? b.m : b.n;
        const int64_t first = role == 0 ? b.row0 : b.col0;
        if (n <= 0 || b.k <= 0) return b.k;
        std::vector<int> loc(static_cast<size_t>(n));
        for (int t = 0; t < n; ++t) {
            const int64_t r = rows[t] - first;
            if (r < 0 || r >= ld) throw std::out_of_range("HODLR GPU: a requested row outside a block's rows on this rank");
            loc[static_cast<size_t>(t)] = static_cast<int>(r);
        }
        Context& ctx = Context::instance();
        ctx.activate();
        const cudaStream_t stream = ctx.stream();
        DeviceHeap& heap = DeviceHeap::instance();
        const size_t count = static_cast<size_t>(n) * b.k;
        int* rd = heap.alloc<int>(loc.size() * sizeof(int));
        T* od = heap.alloc<T>(count * sizeof(T));
        check_cuda(cudaMemcpyAsync(rd, loc.data(), loc.size() * sizeof(int), cudaMemcpyHostToDevice, stream),
                   "HODLR forward block rows (indices)");
        launch_gather_rows(a, ld, b.k, rd, n, od, stream);
        check_cuda(cudaMemcpyAsync(out, od, count * sizeof(T), cudaMemcpyDeviceToHost, stream),
                   "HODLR forward block rows (to the host)");
        check_cuda(cudaStreamSynchronize(stream), "HODLR forward block rows");
        heap.free(od);
        heap.free(rd);
        return b.k;
    }
    void fetch_block(int idx, T* hu, int ldu, T* hv, int ldv) const {
        const LowRank& b = committed_block(idx);
        Context::instance().activate();
        auto down = [](T* h, int ld, const T* d, int rows, int cols) {
            if (h == nullptr || rows <= 0 || cols <= 0) return;
            if (ld < rows) throw std::invalid_argument("HODLR GPU: leading dimension below the rows of a block");
            check_cuda(cudaMemcpy2D(h, static_cast<size_t>(ld) * sizeof(T), d, static_cast<size_t>(rows) * sizeof(T),
                                    static_cast<size_t>(rows) * sizeof(T), static_cast<size_t>(cols), cudaMemcpyDeviceToHost),
                       "HODLR forward block (to the host)");
        };
        down(hu, ldu, b.u, b.m, b.k);
        down(hv, ldv, b.v, b.n, b.k);
    }

    // y = op(A) x for the host n_loc x nrhs arrays x and y (leading dimension
    // n_loc); op is 'N' or 'T' (a symmetric HODLR ignores it).  mode 2 runs
    // the double GEMM batches on the FP64 tensor cores.
    void mult(char trans, int nrhs, const T* x, T* y, int mode) {
        if (!committed_) throw std::logic_error("HODLR GPU: multiply before the forward blocks are uploaded");
        if (trans != 'N' && trans != 'T') throw std::invalid_argument("HODLR GPU: mult takes 'N' or 'T'");
        {  // (this rank's flops, counted as the CPU counts a gemm: 8 m n k in complex)
            const double scale = fmm::gpu::is_complex_scalar<T> ? 4.0 : 1.0;
            double f = 0.0;
            for (const Leaf& l : leaves_) f += 2.0 * l.m * static_cast<double>(l.m) * nrhs;
            for (const LowRank& b : blocks_) f += (b.sym ? 2.0 : 1.0) * 2.0 * (static_cast<double>(b.m) + b.n) * b.k * nrhs;
            mult_flops_ = scale * f;
        }
        Context& ctx = Context::instance();
        ctx.activate();
        fmm::gpu::tensor_core_gemm() = (mode == 2);
        const cudaStream_t stream = ctx.stream();
        const int64_t n = n_loc_;
        const size_t vec_bytes = static_cast<size_t>(n) * nrhs * sizeof(T);
        DeviceHeap& heap = DeviceHeap::instance();

        // coefficient workspace: k x nrhs per product (two per symmetric
        // block), level after level (a shared level's range is summed over
        // its node's ranks)
        std::vector<int64_t> coef_off(blocks_.size());
        std::vector<int64_t> level_begin(static_cast<size_t>(maxlevel_ + 2), 0);
        int64_t coef_count = 0;
        for (int l = 1; l <= maxlevel_; ++l) {
            level_begin[static_cast<size_t>(l)] = coef_count;
            for (size_t i = 0; i < blocks_.size(); ++i) {
                if (blocks_[i].level != l) continue;
                coef_off[i] = coef_count;
                coef_count += static_cast<int64_t>(blocks_[i].k) * nrhs * (blocks_[i].sym ? 2 : 1);
            }
        }
        level_begin[static_cast<size_t>(maxlevel_ + 1)] = coef_count;
        const bool dist = distributed();
        T* dx = heap.alloc<T>(vec_bytes);
        T* dy = heap.alloc<T>(vec_bytes);
        const size_t coef_bytes = std::max<size_t>(sizeof(T), static_cast<size_t>(coef_count) * sizeof(T));
        DeviceHeap& arena = DeviceHeap::exchange_arena();
        T* dc = nullptr;
        if (dist && fmm::gpu::device_exchange_enabled()) dc = reinterpret_cast<T*>(arena.try_alloc(coef_bytes));
        const bool dc_in_arena = dc != nullptr;
        if (dc == nullptr) dc = heap.alloc<T>(coef_bytes);
        // (a rank without the rows of one side of a shared block contributes zeros)
        if (dist) check_cuda(cudaMemsetAsync(dc, 0, coef_bytes, stream), "HODLR mult coefficients");
        check_cuda(cudaMemcpyAsync(dx, x, vec_bytes, cudaMemcpyHostToDevice, stream), "HODLR mult input");

        const bool nt = trans == 'N';
        VBatch<T> leaf_batch, coef_batch;
        std::vector<VBatch<T>> apply_batch(static_cast<size_t>(maxlevel_ + 1));
        const int ldn = static_cast<int>(n);
        for (const Leaf& l : leaves_) {
            leaf_batch.entries.push_back({l.d, dx + l.row0, dy + l.row0, l.m, nrhs, l.m, l.m, ldn, ldn});
        }
        for (size_t i = 0; i < blocks_.size(); ++i) {
            const LowRank& b = blocks_[i];
            if (b.k == 0) continue;
            T* c0 = dc + coef_off[i];
            VBatch<T>& apply = apply_batch[static_cast<size_t>(b.level)];
            // (the parts this rank holds: U when m > 0, V when n > 0)
            const bool hu = b.m > 0, hv = b.n > 0;
            if (b.sym) {
                // A21 x0 = U1 (V0^T x0) into rows of child 1; A12 x1 = V0 (U1^T x1) into child 0
                T* c1 = c0 + static_cast<int64_t>(b.k) * nrhs;
                if (hv) coef_batch.entries.push_back({b.v, dx + b.col0, c0, b.k, nrhs, b.n, b.n, ldn, b.k});
                if (hu) coef_batch.entries.push_back({b.u, dx + b.row0, c1, b.k, nrhs, b.m, b.m, ldn, b.k});
                if (hu) apply.entries.push_back({b.u, c0, dy + b.row0, b.m, nrhs, b.k, b.m, b.k, ldn});
                if (hv) apply.entries.push_back({b.v, c1, dy + b.col0, b.n, nrhs, b.k, b.n, b.k, ldn});
            } else if (nt) {
                if (hv) coef_batch.entries.push_back({b.v, dx + b.col0, c0, b.k, nrhs, b.n, b.n, ldn, b.k});
                if (hu) apply.entries.push_back({b.u, c0, dy + b.row0, b.m, nrhs, b.k, b.m, b.k, ldn});
            } else {
                if (hu) coef_batch.entries.push_back({b.u, dx + b.row0, c0, b.k, nrhs, b.m, b.m, ldn, b.k});
                if (hv) apply.entries.push_back({b.v, c0, dy + b.col0, b.n, nrhs, b.k, b.n, b.k, ldn});
            }
        }
        meta_.clear();
        leaf_batch.stage(meta_);
        coef_batch.stage(meta_);
        for (VBatch<T>& a : apply_batch) a.stage(meta_);
        char* meta = meta_.upload(meta_device_, stream);

        const magma_queue_t queue = ctx.queue();
        const T one(1.0), zero(0.0);
        // the leaves cover every local row, so their products initialize y
        leaf_batch.gemm(meta, nt ? MagmaNoTrans : MagmaTrans, MagmaNoTrans, one, zero, queue);
        if (coef_batch.count() > 0) coef_batch.gemm(meta, MagmaTrans, MagmaNoTrans, one, zero, queue);
        // the shared levels: the partial coefficients summed over each node's ranks
        if (dist) {
            for (int l = 1; l <= maxlevel_; ++l) {
                if (!shared_level(l)) continue;
                const int64_t c0 = level_begin[static_cast<size_t>(l)], c1 = level_begin[static_cast<size_t>(l + 1)];
                sum_over(level_comm_[static_cast<size_t>(l)], dc + c0, c1 - c0, dc_in_arena, stream);
            }
        }
        // one level at a time: the blocks of a level write disjoint rows, those of different levels do not
        for (const VBatch<T>& a : apply_batch) {
            if (a.count() > 0) a.gemm(meta, MagmaNoTrans, MagmaNoTrans, one, one, queue);
        }
        check_cuda(cudaMemcpyAsync(y, dy, vec_bytes, cudaMemcpyDeviceToHost, stream), "HODLR mult output");
        check_cuda(cudaStreamSynchronize(stream), "HODLR mult");
        if (dc_in_arena) {
            arena.free(dc);
        } else {
            heap.free(dc);
        }
        heap.free(dy);
        heap.free(dx);
    }

    // a[0 .. count) (device) = its copy on rank 0 of comm (as sum_over)
    static void bcast_from_root(MPI_Comm comm, T* a, int64_t count, bool device_buffer, cudaStream_t stream) {
        if (count <= 0) return;
        constexpr int per = fmm::gpu::is_complex_scalar<T> ? 2 : 1;
        const int64_t n = count * per;
        if (n > std::numeric_limits<int>::max()) throw std::runtime_error("HODLR GPU: message too long for MPI");
        check_cuda(cudaStreamSynchronize(stream), "HODLR broadcast");
        if (device_buffer) {
            MPI_Bcast(a, static_cast<int>(n), MPI_DOUBLE, 0, comm);
            return;
        }
        static PinnedBuffer* staging = new PinnedBuffer;  // never freed (CUDA shutdown order)
        double* h = static_cast<double*>(staging->reserve(static_cast<size_t>(n) * sizeof(double)));
        int rank = 0;
        MPI_Comm_rank(comm, &rank);
        if (rank == 0) check_cuda(cudaMemcpy(h, a, static_cast<size_t>(count) * sizeof(T), cudaMemcpyDeviceToHost), "HODLR broadcast download");
        MPI_Bcast(h, static_cast<int>(n), MPI_DOUBLE, 0, comm);
        if (rank != 0) check_cuda(cudaMemcpy(a, h, static_cast<size_t>(count) * sizeof(T), cudaMemcpyHostToDevice), "HODLR broadcast upload");
    }

    // a[0 .. count) (device) = its sum over the ranks of comm: in place in
    // device memory when MPI is CUDA-aware and a is in the exchange arena,
    // else through a pinned host buffer.  Waits for the stream.
    static void sum_over(MPI_Comm comm, T* a, int64_t count, bool device_buffer, cudaStream_t stream) {
        if (count <= 0) return;
        // (complex values sum as pairs of doubles)
        constexpr int per = fmm::gpu::is_complex_scalar<T> ? 2 : 1;
        const int64_t n = count * per;
        if (n > std::numeric_limits<int>::max()) throw std::runtime_error("HODLR GPU: message too long for MPI");
        check_cuda(cudaStreamSynchronize(stream), "HODLR sum over ranks");
        if (device_buffer) {
            MPI_Allreduce(MPI_IN_PLACE, a, static_cast<int>(n), MPI_DOUBLE, MPI_SUM, comm);
            return;
        }
        static PinnedBuffer* staging = new PinnedBuffer;  // never freed (CUDA shutdown order)
        double* h = static_cast<double*>(staging->reserve(static_cast<size_t>(n) * sizeof(double)));
        check_cuda(cudaMemcpy(h, a, static_cast<size_t>(count) * sizeof(T), cudaMemcpyDeviceToHost), "HODLR sum download");
        MPI_Allreduce(MPI_IN_PLACE, h, static_cast<int>(n), MPI_DOUBLE, MPI_SUM, comm);
        check_cuda(cudaMemcpy(a, h, static_cast<size_t>(count) * sizeof(T), cudaMemcpyHostToDevice), "HODLR sum upload");
    }

private:
    static size_t bytes_of(int64_t count) { return fmm::gpu::align_up(static_cast<size_t>(count) * sizeof(T)); }

    void check_rows(int64_t row0, int64_t m, const char* what) const {
        if (row0 < 0 || m < 0 || row0 + m > n_loc_) {
            throw std::invalid_argument(std::string("HODLR GPU: ") + what + " outside the local rows");
        }
    }

    void release_forward() {
        if (forward_ != nullptr) {
            Context::instance().activate();
            DeviceHeap::instance().free(forward_);
            forward_ = nullptr;
        }
        forward_bytes_ = 0;
        leaves_.clear();
        blocks_.clear();
        committed_ = false;
    }

    int64_t n_loc_ = 0;
    int maxlevel_ = 0;
    std::vector<MPI_Comm> level_comm_;  // per level: the node's communicator if the level is shared
    std::vector<Leaf> leaves_;
    std::vector<LowRank> blocks_;
    char* forward_ = nullptr;
    size_t forward_bytes_ = 0;
    double mult_flops_ = 0.0;
    struct Mirror {
        const void* dev;
        size_t bytes;
        bool stale;  // the host array was not filled
    };
    std::unordered_map<const void*, Mirror> mirrors_;  // device copies of host arrays (add_mirror)
    std::vector<void*> mirror_memory_;                  // their device memory, freed after the commit
    size_t mirrored_bytes_ = 0;                         // of the last commit, copied on the device
    bool committed_ = false;
    MetaBuilder meta_;
    DeviceBuffer meta_device_;
};

}  // namespace gpu
}  // namespace bpack

#endif  // H2_HAVE_GPU
