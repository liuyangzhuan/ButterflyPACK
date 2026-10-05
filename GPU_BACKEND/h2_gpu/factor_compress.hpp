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

#include <algorithm>
#include <cmath>
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
// |v| for the device element types (double, and the backend's dcomplex).
inline double xnr_magnitude(double v) { return std::abs(v); }
template<typename T>
inline double xnr_magnitude(const T& v) { return std::sqrt(v.re * v.re + v.im * v.im); }

// The compression's tolerance as a fraction of tol_comp.  At tol_comp itself
// the quantization error (about 6e-4 relative) is comparable to the H2
// approximation error, and measured at 96^3 it multiplied acc_mvp by 3.6 and
// acc_forward by 4.0.  A fifth of it costs little ratio (X_NR compresses about
// 23x at tol_comp, 12-15x here) and keeps the error well below the
// approximation's.
inline constexpr double kXnrToleranceFraction = 0.2;

// Lossy Bitcomp decompresses a batch in parallel across its blocks, so the
// rate follows the block count, not the bytes: on an A100 a wave of 37 boxes
// (436 MB, one block each) decompresses at 63 GB/s, the same bytes cut into
// 0.5 MB blocks at 548 GB/s.  So each box's X_NR is compressed in chunks of
// this many bytes - contiguous ranges of its column-major elements, which
// reconstruct the same buffer.  A chunk costs its 256-byte padding (about
// 0.9% of a compressed chunk) and the measured ratio loss is 0.04x of 19.0.
constexpr size_t kXnrChunkBytes = 512 * 1024;

// One box's compressed X_NR: the chunks it was cut into, and the shape and
// quantization the solve needs to put it back.
struct CompressedXnr {
    size_t bytes = 0;       // compressed bytes in all its chunks, as stored
    int64_t rows = 0;       // the uncompressed shape, for the solve
    int64_t cols = 0;
    double delta = 0.0;     // the quantization step it was compressed with
    size_t chunk0 = 0;      // its chunks in the compressor's chunk list
    size_t nchunks = 0;
};

// The compressed X_NR of one solve wave, decompressed into a shared scratch
// before the wave's sweep kernels run and discarded after, so only one wave's
// X_NR is uncompressed at a time.  Both sweeps read X_NR, so a wave is
// decompressed twice per solve; at a wave-sized batch that is a few hundred
// GB/s, against a solve of a quarter second.
//
// The batch plan is built once (from the compressed data itself) and reused by
// every sweep of every solve.
struct SolveXnrWave {
    // per box of the wave: where it goes in the level's SolveBox array, and
    // where its X_NR lands in the scratch (its first chunk's destination)
    std::vector<int> boxes;
    std::vector<void*> h_box_dsts;
    int* d_box_idx = nullptr;          // for launch_patch_xnr
    void** d_box_dsts = nullptr;
    // per chunk of the wave, in box order: the batch the decompression runs
    std::vector<const void*> h_srcs;   // compressed bytes, in the kept block
    std::vector<void*> h_dsts;         // into the level's scratch
    const void** d_srcs = nullptr;
    void** d_dsts = nullptr;
    size_t raw_bytes = 0;              // uncompressed bytes the wave decompresses to
#ifdef H2_HAVE_NVCOMP
    bitcompHandle_t plan{};
    bool have_plan = false;
#endif
    bool empty() const { return boxes.empty(); }
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
// nt*r elements each).  Each block is cut into chunks of kXnrChunkBytes, so
// the batch has enough blocks to fill the GPU; on return `out` holds, per
// block, its chunks in `chunks()` and the bytes they come to, and the caller
// copies each chunk into the kept block it then sizes from them.  `staging` is
// transient: it is only as large as a wave, which is already resident while
// the wave runs.
//
// `tolerance` is the relative tolerance (tol_comp); the delta of each block is
// 2 * tolerance * max|block|, measured on the device.
template<typename S>
class WaveXnrCompressor {
public:
    // One chunk in the staging buffer: where it is, how long it came out, and
    // the padded length the kept block gives it.
    struct StageChunk {
        size_t offset = 0;
        size_t bytes = 0;
        size_t stride = 0;
    };

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
        const size_t chunk_elems = std::max<size_t>(1, kXnrChunkBytes / sizeof(S));

        // ---- the chunks of the whole wave, in box order.  A chunk is a
        // contiguous range of a box's elements, so the box's buffer comes back
        // as its chunks laid end to end.
        chunks_.clear();
        std::vector<const double*> h_in;
        std::vector<size_t> nbytes;
        std::vector<size_t> owner;  // per chunk, its box
        for (size_t i = 0; i < n; ++i) {
            out[i].chunk0 = chunks_.size();
            for (size_t at = 0; at < elems[i]; at += chunk_elems) {
                const size_t take = std::min(chunk_elems, elems[i] - at);
                h_in.push_back(reinterpret_cast<const double*>(srcs[i] + at));
                nbytes.push_back(take * sizeof(S));
                owner.push_back(i);
                chunks_.push_back(StageChunk{});
            }
            out[i].nchunks = chunks_.size() - out[i].chunk0;
        }
        const size_t nc = chunks_.size();
        if (nc == 0) return false;

        size_t stage_bytes = 0;
        for (size_t c = 0; c < nc; ++c) {
            chunks_[c].offset = stage_bytes;
            stage_bytes += align_up(bitcompMaxBuflen(nbytes[c]));
        }
        const size_t meta_bytes = align_up(n * sizeof(AbsMaxItem)) + align_up(n * sizeof(double)) +
                                  3 * align_up(nc * sizeof(void*)) + 2 * align_up(nc * sizeof(size_t)) +
                                  align_up(nc * sizeof(double));
        staging_ = heap_.try_alloc(stage_bytes);
        if (staging_ == nullptr) return false;
        meta_ = heap_.try_alloc(meta_bytes);
        if (meta_ == nullptr) { release(); return false; }

        // ---- per-box max|.| on the device, for the deltas.  The scale stays
        // per box, not per chunk, so the tolerance a box is held to does not
        // depend on how it was cut.
        char* p = meta_;
        auto* d_absmax_items = reinterpret_cast<AbsMaxItem*>(p); p += align_up(n * sizeof(AbsMaxItem));
        auto* d_maxabs = reinterpret_cast<double*>(p);            p += align_up(n * sizeof(double));
        auto* d_in = reinterpret_cast<const double**>(p);         p += align_up(nc * sizeof(void*));
        auto* d_out = reinterpret_cast<void**>(p);                p += align_up(nc * sizeof(void*));
        auto* d_nb = reinterpret_cast<size_t*>(p);                p += align_up(nc * sizeof(size_t));
        auto* d_csz = reinterpret_cast<size_t*>(p);               p += align_up(nc * sizeof(size_t));
        auto* d_delta = reinterpret_cast<double*>(p);

        std::vector<AbsMaxItem> absmax(n);
        for (size_t i = 0; i < n; ++i) {
            absmax[i].a = reinterpret_cast<const double*>(srcs[i]);
            absmax[i].n = static_cast<int64_t>(elems[i] * per);
            absmax[i].out = d_maxabs + i;
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

        for (size_t i = 0; i < n; ++i) {
            // A block of zeros has no scale: any positive delta is exact.
            out[i].delta = maxabs[i] > 0.0 ? 2.0 * tolerance * maxabs[i] : 1.0;
        }
        std::vector<void*> h_out(nc);
        std::vector<double> delta(nc);
        for (size_t c = 0; c < nc; ++c) {
            h_out[c] = staging_ + chunks_[c].offset;
            delta[c] = out[owner[c]].delta;
        }

        check_cuda(cudaMemcpyAsync(d_in, h_in.data(), nc * sizeof(void*), cudaMemcpyHostToDevice, stream_), "xnr in");
        check_cuda(cudaMemcpyAsync(d_out, h_out.data(), nc * sizeof(void*), cudaMemcpyHostToDevice, stream_), "xnr out");
        check_cuda(cudaMemcpyAsync(d_nb, nbytes.data(), nc * sizeof(size_t), cudaMemcpyHostToDevice, stream_), "xnr nb");
        check_cuda(cudaMemcpyAsync(d_delta, delta.data(), nc * sizeof(double), cudaMemcpyHostToDevice, stream_), "xnr delta");

        H2_BK(bitcompCreateBatchPlan(&plan_, nc, BITCOMP_FP64_DATA, BITCOMP_LOSSY_FP_TO_SIGNED,
                                     BITCOMP_DEFAULT_ALGO));
        have_plan_ = true;
        H2_BK(bitcompSetStream(plan_, stream_));
        H2_BK(bitcompBatchCompressLossy_fp64(plan_, d_in, d_out, d_nb, d_csz, d_delta));

        std::vector<size_t> csizes(nc, 0);
        check_cuda(cudaMemcpyAsync(csizes.data(), d_csz, nc * sizeof(size_t),
                                   cudaMemcpyDeviceToHost, stream_), "xnr sizes");
        check_cuda(cudaStreamSynchronize(stream_), "xnr sizes");
        // The kept block holds a box's chunks end to end, each padded to 256
        // bytes: the gather copies 4-byte words, but bitcomp reads a chunk's
        // compressed bytes with wide loads and faults on a weaker start.  That
        // padded length is what the solve walks to find a chunk (bitcomp reads
        // the real one from the data itself).
        for (size_t c = 0; c < nc; ++c) {
            chunks_[c].bytes = csizes[c];
            chunks_[c].stride = align_up(csizes[c]);
            out[owner[c]].bytes += chunks_[c].stride;
        }
        return true;
    }

    // The wave's chunks, in box order (CompressedXnr::chunk0 indexes this).
    const std::vector<StageChunk>& chunks() const { return chunks_; }

    // Where a block's compressed bytes sit, for the copy into the kept block.
    const char* staging() const { return staging_; }

    void release() {
        if (have_plan_) { bitcompDestroyPlan(plan_); have_plan_ = false; }
        if (meta_ != nullptr) { heap_.free(meta_); meta_ = nullptr; }
        if (staging_ != nullptr) { heap_.free(staging_); staging_ = nullptr; }
        chunks_.clear();
    }

private:
    DeviceHeap& heap_;
    cudaStream_t stream_ = nullptr;
    char* staging_ = nullptr;
    char* meta_ = nullptr;
    std::vector<StageChunk> chunks_;
    bitcompHandle_t plan_{};
    bool have_plan_ = false;
};


#ifdef H2_HAVE_NVCOMP

// Build a wave's plan and upload its pointer arrays.  `temp` is scratch the
// plan keeps for the compressed headers.
inline bool build_solve_xnr_wave(SolveXnrWave& w, DeviceHeap& heap, cudaStream_t stream) {
    if (w.boxes.empty()) return true;
    const size_t nb = w.boxes.size();    // boxes, for the pointer patch
    const size_t nc = w.h_srcs.size();   // chunks, the decompression batch
    const size_t chunk_span = align_up(nc * sizeof(void*));
    const size_t box_span = align_up(nb * sizeof(void*));
    char* ptrs = heap.alloc_resident(chunk_span * 2 + box_span + align_up(nb * sizeof(int)));
    w.d_srcs = reinterpret_cast<const void**>(ptrs);
    w.d_dsts = reinterpret_cast<void**>(ptrs + chunk_span);
    w.d_box_dsts = reinterpret_cast<void**>(ptrs + chunk_span * 2);
    w.d_box_idx = reinterpret_cast<int*>(ptrs + chunk_span * 2 + box_span);
    check_cuda(cudaMemcpyAsync(w.d_srcs, w.h_srcs.data(), nc * sizeof(void*),
                               cudaMemcpyHostToDevice, stream), "xnr wave srcs");
    check_cuda(cudaMemcpyAsync(w.d_dsts, w.h_dsts.data(), nc * sizeof(void*),
                               cudaMemcpyHostToDevice, stream), "xnr wave dsts");
    check_cuda(cudaMemcpyAsync(w.d_box_dsts, w.h_box_dsts.data(), nb * sizeof(void*),
                               cudaMemcpyHostToDevice, stream), "xnr wave box dsts");
    check_cuda(cudaMemcpyAsync(w.d_box_idx, w.boxes.data(), nb * sizeof(int),
                               cudaMemcpyHostToDevice, stream), "xnr wave box indices");
    check_cuda(cudaStreamSynchronize(stream), "xnr wave pointers");
    if (bitcompCreateBatchPlanFromCompressedData(&w.plan, w.d_srcs, nc, nullptr, stream) != BITCOMP_SUCCESS) {
        return false;
    }
    w.have_plan = true;
    bitcompSetStream(w.plan, stream);
    return true;
}

// Decompress a wave's X_NR into the scratch its pointers were built for.
inline void decompress_solve_xnr_wave(const SolveXnrWave& w, cudaStream_t stream) {
    if (w.boxes.empty() || !w.have_plan) return;
    bitcompSetStream(w.plan, stream);
    if (bitcompBatchUncompress(w.plan, w.d_srcs, w.d_dsts) != BITCOMP_SUCCESS) {
        throw std::runtime_error("H2 GPU solve: decompressing a wave's X_NR failed");
    }
}

inline void release_solve_xnr_wave(SolveXnrWave& w) {
    if (w.have_plan) { bitcompDestroyPlan(w.plan); w.have_plan = false; }
}

#endif  // H2_HAVE_NVCOMP

// Decompress one block on the host (bitcompHost*), for the paths that read
// the factors there: the host solve fallback and a kept level solved on the
// host.  The format is the same as the device's.
template<typename S>
inline void decompress_xnr_host(const void* src, S* dst, const uint32_t* strides, size_t nchunks) {
    const size_t chunk_elems = std::max<size_t>(1, kXnrChunkBytes / sizeof(S));
    const char* at = static_cast<const char*>(src);
    for (size_t c = 0; c < nchunks; ++c) {
        bitcompHandle_t plan{};
        H2_BK(bitcompCreatePlanFromCompressedData(&plan, at));
        H2_BK(bitcompHostUncompress(plan, at, dst + c * chunk_elems));
        bitcompDestroyPlan(plan);
        at += strides[c];
    }
}

// Decompress one block, for the solve.  `dst` must hold rows*cols elements.
template<typename S>
inline void decompress_xnr(const char* src, S* dst, const uint32_t* strides, size_t nchunks,
                           cudaStream_t stream) {
    const size_t chunk_elems = std::max<size_t>(1, kXnrChunkBytes / sizeof(S));
    const char* at = src;
    for (size_t c = 0; c < nchunks; ++c) {
        bitcompHandle_t plan{};
        H2_BK(bitcompCreatePlanFromCompressedData(&plan, at));
        H2_BK(bitcompSetStream(plan, stream));
        H2_BK(bitcompUncompress(plan, at, dst + c * chunk_elems));
        bitcompDestroyPlan(plan);
        at += strides[c];
    }
}

#undef H2_BK

#endif  // H2_HAVE_NVCOMP

}  // namespace gpu
}  // namespace fmm
