#pragma once
// Blocks of the kind-3 kernel (the EMSURF EFIE entry, emsurf_kernel.cuh)
// by triangle pairs (emsurf_blocks.cu).  An RWG edge lives on two triangles,
// so an entry is a sum over 4 triangle pairs, and a triangle pair serves the
// up to 9 edge pairs on it.  For a block of row edges by column edges the
// quadrature of every pair of their triangles runs once, into 9 sums per pair
// that do not depend on the edges; the entries then combine the sums of
// their 4 pairs with the edges' geometry.

#ifdef H2_HAVE_GPU

#include "device_kernels.hpp"

#include <cstdint>

namespace fmm {
namespace gpu {

constexpr int kEfieSums = 9;  // per triangle pair (emsurf_blocks.cu)

// One block: nrow x ncol entries out(r, c) = out[r * rs + c * cs] of the
// row and column edges.  Int lists at byte offsets of the metadata image:
// the edges (ids), each edge's two triangles as indices into the block's
// distinct triangles (-1: none), and those triangles (ids).  M: scratch of
// ntr * ntc * kEfieSums values.
struct EfieBlockItem {
    int nrow = 0, ncol = 0, ntr = 0, ntc = 0;
    int64_t row_edges = 0, col_edges = 0, row_tri = 0, col_tri = 0, tr = 0, tc = 0;
    dcomplex* M = nullptr;
    dcomplex* out = nullptr;
    int64_t rs = 0, cs = 0;
};

// The sums of every triangle pair of each block, then its entries.
void launch_efie_blocks(const EfieBlockItem* items, int count, int64_t max_pairs, int64_t max_entries,
                        const char* meta, KernelSpec spec, cudaStream_t stream);

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
