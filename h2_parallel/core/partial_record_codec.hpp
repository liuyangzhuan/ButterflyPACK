#pragma once
// =============================================================================
// H2_compress_factors under H2_partial: the F22 records (PartialRecord,
// core/tree.hpp) stored compressed.
//
// The records (X_RR, X_RS, X_SR, X_RN and the neighbours' X_NR of every
// sparsified box) are what the CB entries and products read, and nearly all
// of a partial factorization's host memory (32 of 40 GB at 64^3).  They are
// kept as lossy Bitcomp chunks in the format of the GPU backend's kept
// factors (GPU_BACKEND/h2_gpu/factor_compress.hpp): contiguous ranges of
// kRecordChunkBytes of the column-major entries, each chunk's bytes padded to
// kRecordChunkAlign, quantized at delta = 2 * fraction * tol * max|record|
// (Bitcomp rounds delta down to a power of two, so every entry is within
// fraction * tol * max|record|; a complex entry's components each).  The
// device compresses a chunk of boxes' records before their copy to the host
// (LevelEliminator::f22_submit_records), the host path a box's as it forms
// them, and the CB entries and products decompress what they read.  T, which
// the partial solves read, stays as it is.
//
// Without H2_HAVE_NVCOMP the records stay uncompressed.
// =============================================================================

#include "tree.hpp"

#ifdef H2_HAVE_NVCOMP
#include <nvcomp/native/bitcomp.h>
#endif

#include <algorithm>
#include <memory>
#include <mutex>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace fmm {

// A chunk's uncompressed bytes and a stored chunk's alignment, as the kept
// factors' (kXnrChunkBytes; bitcomp's device reads are wide loads).
inline constexpr size_t kRecordChunkBytes = 512 * 1024;
inline constexpr size_t kRecordChunkAlign = 256;
// The quantization step as a fraction of the tolerance (of each record's
// largest entry).  The CB entries and products are what a parent front
// assembles, so their error is held near the uncompressed one: 32^3 (F11 =
// 4096) against dense LAPACK, schur_mult / schur_extract 2.00e-6 / 1.54e-7
// uncompressed; 0.2: 2.67e-6 / 6.82e-7; 0.05: 2.04e-6 / 2.32e-7; 0.02: 2.01e-6
// / 1.68e-7.  64^3 (F11 = 16384) records 31.6 GB uncompressed: 4.2 / 5.2 /
// 5.9 GB.
inline constexpr double kRecordToleranceFraction = 0.02;

inline bool partial_record_codec_built() {
#ifdef H2_HAVE_NVCOMP
    return true;
#else
    return false;
#endif
}

#ifdef H2_HAVE_NVCOMP
inline void partial_record_check(bitcompResult_t r, const char* what) {
    if (r != BITCOMP_SUCCESS) {
        throw std::runtime_error(std::string("H2_partial record compression: ") + what + ": bitcomp error " +
                                 std::to_string(static_cast<int>(r)));
    }
}
#endif

// The record compressed in place on the host (its entries freed); `tolerance`
// the relative tolerance, the step 2 * kRecordToleranceFraction * tolerance *
// max|record|.  Nothing without a compressor, or for an empty record.
template<typename DataType>
void compress_record_host(PartialRecord<DataType>& r, double tolerance) {
#ifdef H2_HAVE_NVCOMP
    const size_t n = r.data.size();
    if (n == 0 || r.compressed()) return;
    const double* d = reinterpret_cast<const double*>(r.data.data());
    const size_t nd = n * (sizeof(DataType) / sizeof(double));
    double maxabs = 0.0;
    for (size_t i = 0; i < nd; ++i) maxabs = std::max(maxabs, std::abs(d[i]));
    const double delta = maxabs > 0.0 ? 2.0 * kRecordToleranceFraction * tolerance * maxabs : 1.0;
    const size_t chunk_elems = std::max<size_t>(1, kRecordChunkBytes / sizeof(DataType));
    const size_t nchunks = (n + chunk_elems - 1) / chunk_elems;
    std::vector<uint32_t> strides(nchunks);
    std::vector<double> out((bitcompMaxBuflen(chunk_elems * sizeof(DataType)) + sizeof(double) - 1) / sizeof(double));
    std::vector<char, DefaultInitAllocator<char>> packed;
    size_t at = 0;
    for (size_t c = 0; c < nchunks; ++c) {
        const size_t e0 = c * chunk_elems, len = std::min(chunk_elems, n - e0);
        bitcompHandle_t plan{};
        partial_record_check(bitcompCreatePlan(&plan, len * sizeof(DataType), BITCOMP_FP64_DATA, BITCOMP_LOSSY_FP_TO_SIGNED,
                                               BITCOMP_DEFAULT_ALGO),
                             "plan");
        const bitcompResult_t res = bitcompHostCompressLossy_fp64(plan, d + e0 * (sizeof(DataType) / sizeof(double)),
                                                                  out.data(), delta);
        bitcompDestroyPlan(plan);
        partial_record_check(res, "compress");
        size_t bytes = 0;
        partial_record_check(bitcompGetCompressedSize(out.data(), &bytes), "size");
        const size_t stride = (bytes + kRecordChunkAlign - 1) / kRecordChunkAlign * kRecordChunkAlign;
        packed.resize(at + stride);
        std::memcpy(packed.data() + at, out.data(), bytes);
        std::memset(packed.data() + at + bytes, 0, stride - bytes);
        strides[c] = static_cast<uint32_t>(stride);
        at += stride;
    }
    r.packed = std::move(packed);
    r.strides = std::move(strides);
    decltype(r.data)().swap(r.data);
#else
    (void)r;
    (void)tolerance;
#endif
}

// The record's entries (column-major, lda = rows): in place, or decompressed
// into `scratch` (its size rows * cols).
template<typename DataType>
const DataType* record_values(const PartialRecord<DataType>& r, std::vector<DataType>& scratch) {
    if (!r.compressed()) return r.data.data();
#ifdef H2_HAVE_NVCOMP
    const size_t n = static_cast<size_t>(r.rows) * static_cast<size_t>(r.cols);
    scratch.resize(n);
    const size_t chunk_elems = std::max<size_t>(1, kRecordChunkBytes / sizeof(DataType));
    const char* at = r.packed.data();
    for (size_t c = 0; c < r.strides.size(); ++c) {
        bitcompHandle_t plan{};
        partial_record_check(bitcompCreatePlanFromCompressedData(&plan, at), "plan");
        const bitcompResult_t res = bitcompHostUncompress(plan, at, scratch.data() + c * chunk_elems);
        bitcompDestroyPlan(plan);
        partial_record_check(res, "decompress");
        at += r.strides[c];
    }
    return scratch.data();
#else
    throw std::logic_error("H2_partial: a compressed record without a compressor");
#endif
}

// The records' entries the CB entries gather (a few rows and columns of a
// box's records for every block pair), read by many threads: each 512 KB
// chunk decompressed once, by the first reader, and kept until release(group)
// frees the chunks of the records a group's readers read (a group: the row
// blocks of one box, the only readers of its records).  Not kept across
// calls: a later factorization's records may sit at the same addresses.
template<typename DataType>
class RecordChunkCache {
public:
    /// a reader's last chunk (one a task)
    struct Reader {
        const PartialRecord<DataType>* record = nullptr;
        size_t chunk = 0;
        const DataType* values = nullptr;
    };
    DataType at(const PartialRecord<DataType>& r, size_t i, int64_t group, Reader& reader) {
        if (!r.compressed()) return r.data[i];
        const size_t c = i / kChunkElems;
        if (&r != reader.record || c != reader.chunk) {
            reader.values = chunk(r, c, group);
            reader.record = &r;
            reader.chunk = c;
        }
        return reader.values[i - c * kChunkElems];
    }
    void release(int64_t group) {
        std::lock_guard<std::mutex> lock(mutex_);
        for (auto it = records_.begin(); it != records_.end();) it = it->second->group == group ? records_.erase(it) : std::next(it);
    }

private:
    static constexpr size_t kChunkElems = std::max<size_t>(1, kRecordChunkBytes / sizeof(DataType));
    struct Chunk {
        std::once_flag once;
        std::vector<DataType> values;
    };
    struct Record {
        int64_t group = 0;
        std::vector<size_t> offsets;  // each chunk's bytes in packed
        std::vector<std::unique_ptr<Chunk>> chunks;
    };
    const DataType* chunk(const PartialRecord<DataType>& r, size_t c, int64_t group) {
        Chunk* ch = nullptr;
        const char* src = nullptr;
        {
            std::lock_guard<std::mutex> lock(mutex_);
            auto& e = records_[&r];
            if (!e) {
                e = std::make_unique<Record>();
                e->group = group;
                e->offsets.assign(r.strides.size() + 1, 0);
                for (size_t k = 0; k < r.strides.size(); ++k) e->offsets[k + 1] = e->offsets[k] + r.strides[k];
                e->chunks.resize(r.strides.size());
            }
            if (!e->chunks[c]) e->chunks[c] = std::make_unique<Chunk>();
            ch = e->chunks[c].get();
            src = r.packed.data() + e->offsets[c];
        }
        std::call_once(ch->once, [&] {
            const size_t n = static_cast<size_t>(r.rows) * static_cast<size_t>(r.cols);
            ch->values.resize(std::min(kChunkElems, n - c * kChunkElems));
            decompress_chunk(src, ch->values.data());
        });
        return ch->values.data();
    }
    static void decompress_chunk(const char* at, DataType* dst) {
#ifdef H2_HAVE_NVCOMP
        bitcompHandle_t plan{};
        partial_record_check(bitcompCreatePlanFromCompressedData(&plan, at), "plan");
        const bitcompResult_t res = bitcompHostUncompress(plan, at, dst);
        bitcompDestroyPlan(plan);
        partial_record_check(res, "decompress");
#else
        (void)at;
        (void)dst;
        throw std::logic_error("H2_partial: a compressed record without a compressor");
#endif
    }
    std::mutex mutex_;
    std::unordered_map<const PartialRecord<DataType>*, std::unique_ptr<Record>> records_;
};

}  // namespace fmm
