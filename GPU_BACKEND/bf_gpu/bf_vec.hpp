#pragma once

#ifdef H2_HAVE_GPU

// The intermediate vectors of a butterfly sweep, on the device
// (GPU_BACKEND/bf_gpu).
//
// The CPU sweep keeps BFvec%vec(level): an nr x nc grid of blocks, each
// (rows, nrhs) with `rows` the row count of the kernel block that produced
// it, so the heights vary across a level.  Here one level is one contiguous
// device buffer with an offset table, which is what a batched GEMM wants.
//
// The indexing follows the CPU exactly (Bplus_utilities.f90, the 'N' sweep):
// a level carries idx_r/inc_r and idx_c/inc_c, and a global (index_i,
// index_j) maps to the local pair
//     i_loc = (index_i - idx_r)/inc_r,   j_loc = (index_j - idx_c)/inc_c
// (zero based here, one based there).  An interior level reads its input at
//     index_ii = (index_i + 1)/2,        index_jj = 2*index_j - 1
// and gathers index_jj and index_jj + 1 into one output block - the two
// GEMMs, beta 0 then beta 1, that make the transform a butterfly rather than
// a product of block-diagonal factors.

#include <cstdint>
#include <stdexcept>
#include <vector>

#include "bf_device.hpp"

namespace bf {
namespace gpu {

// One level's vector blocks.  `rows[i*nc + j]` is a block's height and
// `row_offset[i*nc + j]` where it starts in `data`, in ROWS - the element
// offset is that times nrhs.  Keeping it in rows is what lets an exchange
// plan be built once and reused for any nrhs; `offset` is the same thing in
// elements, for the host-side GEMM setup.  A height of 0 marks a block this
// rank does not hold.
struct VecLevel {
    int nr = 0, nc = 0;
    int idx_r = 1, idx_c = 1, inc_r = 1, inc_c = 1;
    std::vector<int> rows;
    std::vector<int64_t> row_offset;
    std::vector<int64_t> offset;
    DataType* data = nullptr;     // transient, from the arena
    int64_t elements = 0;
    int64_t row_total = 0;        // elements / nrhs
    int nrhs = 0;

    int64_t index(int i_loc, int j_loc) const {
        if (i_loc < 0 || i_loc >= nr || j_loc < 0 || j_loc >= nc) {
            throw std::out_of_range("bf::gpu::VecLevel: block out of range");
        }
        return static_cast<int64_t>(i_loc) * nc + j_loc;
    }
    // zero based local indices from the global ones, as the CPU computes them
    int row_local(int index_i) const { return (index_i - idx_r) / inc_r; }
    int col_local(int index_j) const { return (index_j - idx_c) / inc_c; }

    // Pack the offsets once the heights are set; returns the element count.
    int64_t pack(int nrhs_in) {
        nrhs = nrhs_in;
        row_offset.assign(rows.size(), 0);
        offset.assign(rows.size(), 0);
        int64_t at = 0;
        for (size_t b = 0; b < rows.size(); ++b) {
            row_offset[b] = at;
            offset[b] = at * nrhs;
            at += rows[b];
        }
        row_total = at;
        elements = at * nrhs;
        return elements;
    }
};

}  // namespace gpu
}  // namespace bf

#endif  // H2_HAVE_GPU
