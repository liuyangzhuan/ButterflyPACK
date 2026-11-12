#pragma once
// =============================================================================
// H2_compress_factors for a compression-only H2 (GPU compression,
// GPU_BACKEND/h2_gpu/compression_gpu.hpp): the coupling blocks K(S_t, S_s)
// and the leaf near blocks K(B_t, B_s) are not kept as dense blocks.  With a
// GPU entry evaluator they are not stored at all - the device matvec
// evaluates each block where it reads it - and otherwise they are stored as
// lossy Bitcomp chunks (H2Block::packed, the GPU backend's kept-factor
// format), on the host and on the device.  The host readers - the host
// matvec, entry extraction - take a block's entries through h2_block_values.
// =============================================================================

#include "tree.hpp"
#ifdef H2_HAVE_NVCOMP
#include <nvcomp/native/bitcomp.h>
#endif
#include <stdexcept>
#include <string>

namespace fmm {

// The quantization step of a stored block as a fraction of the tolerance (of
// the block's largest entry): small, since the blocks are the operator itself
// rather than a factor's intermediate (the F22 records' fraction).
inline constexpr double kH2BlockToleranceFraction = 0.02;

// A compression-only H2 block's entries (column-major, lda = rows): the block
// itself, or decompressed into `scratch`.  `site` names the reader for the
// error a block that is not stored raises.
template<typename DataType>
const MatrixStorage<DataType>& h2_block_values(const H2Block<DataType>& block, MatrixStorage<DataType>& scratch,
                                               const char* site) {
    if (!block.compressed()) {
        if (block.matrix.data.empty() && block.matrix.rows * block.matrix.cols > 0) {
            throw std::runtime_error(std::string(site) +
                                     ": this compression-only H2 does not store its blocks (H2_compress_factors with "
                                     "a GPU entry evaluator: the device matvec evaluates them); the host paths "
                                     "cannot read them");
        }
        return block.matrix;
    }
#ifdef H2_HAVE_NVCOMP
    const int64_t rows = block.matrix.rows, cols = block.matrix.cols;
    scratch.rows = rows;
    scratch.cols = cols;
    scratch.lda = rows;
    scratch.format = MatrixStorage<DataType>::FULL;
    scratch.data.resize(static_cast<size_t>(rows * cols));
    char* dst = reinterpret_cast<char*>(scratch.data.data());
    const char* at = block.packed.data();
    for (const uint32_t stride : block.strides) {
        bitcompHandle_t plan{};
        size_t bytes = 0;
        if (bitcompCreatePlanFromCompressedData(&plan, at) != BITCOMP_SUCCESS ||
            bitcompGetUncompressedSizeFromHandle(plan, &bytes) != BITCOMP_SUCCESS ||
            bitcompHostUncompress(plan, at, dst) != BITCOMP_SUCCESS) {
            bitcompDestroyPlan(plan);
            throw std::runtime_error(std::string(site) + ": decompressing a compression-only H2 block failed");
        }
        bitcompDestroyPlan(plan);
        dst += bytes;
        at += stride;
    }
    if (dst != reinterpret_cast<char*>(scratch.data.data() + scratch.data.size())) {
        throw std::runtime_error(std::string(site) + ": a compressed compression-only H2 block has the wrong size");
    }
    return scratch;
#else
    throw std::logic_error(std::string(site) + ": a compressed H2 block without a compressor");
#endif
}

}  // namespace fmm
