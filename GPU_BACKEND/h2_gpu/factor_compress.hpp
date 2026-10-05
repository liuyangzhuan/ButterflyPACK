#pragma once

// =============================================================================
// Lossy compression of the retained X_NR, the bulk of the H2 factors.
//
// BPACK_TRACE=memory's H2_MEM_FACTOR line shows X_NR is 86-88% of the factor
// bytes at every rank count (X_RN is never allocated: the symmetric path uses
// X_NR^T).  It is also smooth: its rows are a kernel coupling between a box's
// redundant points and its one-hop neighbours' points, both Morton-ordered,
// and only the small r x r X_RR^{-1} mixes its columns.  Lossy Bitcomp
// (nvCOMP) compresses it about 23x at a tolerance of tol_comp * maxabs, at
// several hundred GB/s on the device.
//
// The compression is per wave, inside keep_solve_factors: a wave's X_NR is
// still in the wave's working memory there, and the store that accumulates
// over the level then holds compressed bytes from the first wave.  Compressing
// at level end instead would need the whole level's X_NR uncompressed at once,
// which is the peak this is meant to remove.
//
// delta = 2 * tol_comp * maxabs, so the error is at most delta/2 = the
// tolerance; Bitcomp rounds delta down to a power of two, so the achieved
// bound lands in [tol/2, tol].  maxabs is per box (a box whose entries are all
// small must not be compressed to a needlessly tight absolute tolerance).
//
// Without H2_HAVE_NVCOMP every entry point reports "unavailable" and the
// factors are stored uncompressed, exactly as before.
// =============================================================================

#include "gpu_runtime.hpp"
#include "device_heap.hpp"
#include "device_kernels.hpp"

#ifdef H2_HAVE_NVCOMP
#include <nvcomp/native/bitcomp.h>
#endif

#include <cstdint>
#include <string>
#include <vector>

namespace fmm {
namespace gpu {

// Whether the build has a compressor at all.
inline bool factor_compression_built() {
#ifdef H2_HAVE_NVCOMP
    return true;
#else
    return false;
#endif
}

// One box's compressed X_NR inside a kept block.
struct CompressedXnr {
    size_t offset = 0;      // of the compressed bytes in the kept block
    size_t bytes = 0;       // their length
    int64_t rows = 0;       // the uncompressed shape, for the solve
    int64_t cols = 0;
    double delta = 0.0;     // the quantization step it was compressed with
};

#ifdef H2_HAVE_NVCOMP

#define H2_BK(x)                                                              \
    do {                                                                      \
        bitcompResult_t r_ = (x);                                             \
        if (r_ != BITCOMP_SUCCESS) {                                          \
            throw std::runtime_error("H2 GPU factor compression: bitcomp error " + \
                                     std::to_string(static_cast<int>(r_)));   \
        }                                                                     \
    } while (0)

// A wave's X_NR blocks, compressed in one batch.
//
// `srcs`/`bytes` describe the uncompressed blocks (device pointers, contiguous
// nt*r elements each).  On return `out` holds, per block, where its compressed
// bytes live in `staging` and how many there are; the caller copies them into
// the kept block it then sizes from them.  `staging` is transient: it is only
// as large as a wave, which is already resident while the wave runs.
//
// `tolerance` is the relative tolerance (tol_comp); the delta of each block is
// 2 * tolerance * max|block|, measured on the device.
template<typename S>
class WaveXnrCompressor {
public:
    explicit WaveXnrCompressor(DeviceHeap& heap) : heap_(heap) {}
    ~WaveXnrCompressor() { release(); }
    WaveXnrCompressor(const WaveXnrCompressor&) = delete;
    WaveXnrCompressor& operator=(const WaveXnrCompressor&) = delete;

    // Returns false when the scratch does not fit, so the caller stores the
    // wave uncompressed (as it would without a compressor at all).
    bool compress(const std::vector<const S*>& srcs, const std::vector<size_t>& elems,
                  double tolerance, std::vector<CompressedXnr>& out, cudaStream_t stream) {
        release();
        stream_ = stream;
        const size_t n = srcs.size();
        out.assign(n, CompressedXnr{});
        if (n == 0) return false;

        // The doubles of each block: a complex element is two of them, which
        // is the component-wise bound the delta applies to.
        const size_t per = sizeof(S) / sizeof(double);
        std::vector<size_t> nbytes(n), doubles(n);
        size_t stage_bytes = 0, meta_bytes = 0;
        for (size_t i = 0; i < n; ++i) {
            nbytes[i] = elems[i] * sizeof(S);
            doubles[i] = elems[i] * per;
            stage_bytes += align_up(bitcompMaxBuflen(nbytes[i]));
        }
        meta_bytes = align_up(n * sizeof(AbsMaxItem)) + align_up(n * sizeof(double)) +
                     3 * align_up(n * sizeof(void*)) + 2 * align_up(n * sizeof(size_t)) +
                     align_up(n * sizeof(double));
        staging_ = heap_.try_alloc(stage_bytes);
        if (staging_ == nullptr) return false;
        meta_ = heap_.try_alloc(meta_bytes);
        if (meta_ == nullptr) { release(); return false; }

        // ---- per-box max|.| on the device, for the deltas
        char* p = meta_;
        auto* d_absmax_items = reinterpret_cast<AbsMaxItem*>(p); p += align_up(n * sizeof(AbsMaxItem));
        auto* d_maxabs = reinterpret_cast<double*>(p);            p += align_up(n * sizeof(double));
        auto* d_in = reinterpret_cast<const double**>(p);         p += align_up(n * sizeof(void*));
        auto* d_out = reinterpret_cast<void**>(p);                p += align_up(n * sizeof(void*));
        auto* d_nb = reinterpret_cast<size_t*>(p);                p += align_up(n * sizeof(size_t));
        auto* d_csz = reinterpret_cast<size_t*>(p);               p += align_up(n * sizeof(size_t));
        auto* d_delta = reinterpret_cast<double*>(p);

        std::vector<AbsMaxItem> absmax(n);
        std::vector<const double*> h_in(n);
        std::vector<void*> h_out(n);
        size_t at = 0;
        for (size_t i = 0; i < n; ++i) {
            absmax[i].a = reinterpret_cast<const double*>(srcs[i]);
            absmax[i].n = static_cast<int64_t>(doubles[i]);
            absmax[i].out = d_maxabs + i;
            h_in[i] = reinterpret_cast<const double*>(srcs[i]);
            h_out[i] = staging_ + at;
            out[i].offset = at;  // provisional: the staging offset
            at += align_up(bitcompMaxBuflen(nbytes[i]));
        }
        check_cuda(cudaMemcpyAsync(d_absmax_items, absmax.data(), n * sizeof(AbsMaxItem),
                                   cudaMemcpyHostToDevice, stream_), "xnr absmax items");
        launch_abs_max(d_absmax_items, static_cast<int>(n), stream_);

        // maxabs comes back so the deltas (and the solve's dequantization)
        // are known on the host; a wave is a few hundred doubles.
        std::vector<double> maxabs(n, 0.0);
        check_cuda(cudaMemcpyAsync(maxabs.data(), d_maxabs, n * sizeof(double),
                                   cudaMemcpyDeviceToHost, stream_), "xnr maxabs");
        check_cuda(cudaStreamSynchronize(stream_), "xnr maxabs");

        std::vector<double> delta(n);
        for (size_t i = 0; i < n; ++i) {
            // A block of zeros has no scale: any positive delta is exact.
            delta[i] = maxabs[i] > 0.0 ? 2.0 * tolerance * maxabs[i] : 1.0;
            out[i].delta = delta[i];
        }

        check_cuda(cudaMemcpyAsync(d_in, h_in.data(), n * sizeof(void*), cudaMemcpyHostToDevice, stream_), "xnr in");
        check_cuda(cudaMemcpyAsync(d_out, h_out.data(), n * sizeof(void*), cudaMemcpyHostToDevice, stream_), "xnr out");
        check_cuda(cudaMemcpyAsync(d_nb, nbytes.data(), n * sizeof(size_t), cudaMemcpyHostToDevice, stream_), "xnr nb");
        check_cuda(cudaMemcpyAsync(d_delta, delta.data(), n * sizeof(double), cudaMemcpyHostToDevice, stream_), "xnr delta");

        H2_BK(bitcompCreateBatchPlan(&plan_, n, BITCOMP_FP64_DATA, BITCOMP_LOSSY_FP_TO_SIGNED,
                                     BITCOMP_DEFAULT_ALGO));
        have_plan_ = true;
        H2_BK(bitcompSetStream(plan_, stream_));
        H2_BK(bitcompBatchCompressLossy_fp64(plan_, d_in, d_out, d_nb, d_csz, d_delta));

        std::vector<size_t> csizes(n, 0);
        check_cuda(cudaMemcpyAsync(csizes.data(), d_csz, n * sizeof(size_t),
                                   cudaMemcpyDeviceToHost, stream_), "xnr sizes");
        check_cuda(cudaStreamSynchronize(stream_), "xnr sizes");
        for (size_t i = 0; i < n; ++i) out[i].bytes = csizes[i];
        return true;
    }

    // Where a block's compressed bytes sit, for the copy into the kept block.
    const char* staging() const { return staging_; }

    void release() {
        if (have_plan_) { bitcompDestroyPlan(plan_); have_plan_ = false; }
        if (meta_ != nullptr) { heap_.free(meta_); meta_ = nullptr; }
        if (staging_ != nullptr) { heap_.free(staging_); staging_ = nullptr; }
    }

private:
    DeviceHeap& heap_;
    cudaStream_t stream_ = nullptr;
    char* staging_ = nullptr;
    char* meta_ = nullptr;
    bitcompHandle_t plan_{};
    bool have_plan_ = false;
};

// Decompress one block, for the solve.  `dst` must hold rows*cols elements.
template<typename S>
inline void decompress_xnr(const char* src, S* dst, cudaStream_t stream) {
    bitcompHandle_t plan{};
    H2_BK(bitcompCreatePlanFromCompressedData(&plan, src));
    H2_BK(bitcompSetStream(plan, stream));
    H2_BK(bitcompUncompress(plan, src, dst));
    bitcompDestroyPlan(plan);
}

#undef H2_BK

#endif  // H2_HAVE_NVCOMP

}  // namespace gpu
}  // namespace fmm
