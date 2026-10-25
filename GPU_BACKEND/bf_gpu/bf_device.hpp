#pragma once

#ifdef H2_HAVE_GPU
// Device-resident butterfly block (GPU_BACKEND/bf_gpu).
//
// A butterfly block on the host is a tree of small matrices: ButterflyV
// (rightmost factor, nblk_loc blocks), ButterflyKerl(level) for the interior
// levels (an nr x nc grid of blocks each), and ButterflyU (leftmost factor).
// Each block is a column-major DT matrix reached through a pointer, so the
// host layout is scattered over as many allocations as there are blocks.
//
// Here the same factors live in one device allocation per butterfly, with a
// descriptor per block giving its offset, rows, columns and leading
// dimension.  Two reasons for the flat form: a batched GEMM wants arrays of
// pointers it can build without touching the data, and the whole butterfly is
// uploaded and released in one pair of calls, which matters when a
// factorization holds five operators per front over a hundred-odd fronts.
//
// What this header does NOT do: it keeps no opinion about the matvec or the
// extraction.  It owns the layout, the upload and the release only, so that a
// round trip (upload then download) is testable on its own before any kernel
// consumes it.

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "../h2_gpu/device_heap.hpp"

namespace bf {
namespace gpu {

// double or std::complex<double>, as the library's precision copy fixes it
using DataType = double;

// One small matrix of a butterfly, as the device sees it.
struct BlockDesc {
    int64_t offset = 0;  // in elements, from the butterfly's base pointer
    int rows = 0;
    int cols = 0;
    int lda = 0;         // == rows: the upload packs each block tightly
};

// One level's grid of blocks: ButterflyKerl(level) has nr x nc of them,
// stored row-major in `blocks` (index_i_loc * nc + index_j_loc), matching the
// host's own (index_i_loc, index_j_loc) addressing.  The outer factors U and
// V use nr == 1 or nc == 1 and the same indexing.
struct LevelDesc {
    int nr = 0, nc = 0;
    int idx_r = 0, idx_c = 0;  // global row/column number of the first local block
    int inc_r = 0, inc_c = 0;  // stride between local blocks
    int64_t first = 0;         // index into BfDevice::blocks of this level's first block
};

// A (idx, inc, n) grid range per level, as GetLocalBlockRange gives it: the
// global index of local entry k is k * inc + idx, zero based here.  Used for
// both the kernel loop bounds of a sweep level and the vector grids.
struct SweepRanges {
    std::vector<int> idx_r, inc_r, nr, idx_c, inc_c, nc;

    void resize(int levels) {
        const size_t n = static_cast<size_t>(std::max(levels, 0));
        idx_r.assign(n, 1); inc_r.assign(n, 1); nr.assign(n, 0);
        idx_c.assign(n, 1); inc_c.assign(n, 1); nc.assign(n, 0);
    }
    int levels() const { return static_cast<int>(nr.size()); }
};

// The levels of a butterfly, in the order the matvec sweeps them for 'N':
// V (level 0), Kerl(1) .. Kerl(level_butterfly), U (level_butterfly + 1).
// A low-rank block (level_butterfly == 0) has V and U only.
struct BfDevice {
    int level_butterfly = 0;
    int level_half = 0;
    int M = 0, N = 0, M_loc = 0, N_loc = 0;
    std::vector<LevelDesc> levels;    // level_butterfly + 2 of them
    std::vector<BlockDesc> blocks;    // every level's blocks, concatenated
    // The grids the sweep runs over, fetched from the Fortran at upload: the
    // kernel loop bounds in both orderings (levels 0 .. level_butterfly + 1)
    // and the vector grids of each orientation (vec levels 0 ..
    // level_butterfly + 2).  Properties of the block, so they are read once
    // rather than per matvec.
    SweepRanges kernel_r, kernel_c;    // GetLocalBlockRange, dir 'R' and 'C'
    SweepRanges vec_n, vec_t;          // BFvec grids for chara 'N' and 'T'
    int nproc = 1, myid = 0;
    // the factors themselves, one resident block of the shared arena
    // (fmm::gpu::DeviceHeap), so a butterfly draws from the same pool as the
    // H2s of F11 and F22 rather than a private allocation
    DataType* values = nullptr;
    int64_t elements = 0;

    bool resident() const { return values != nullptr; }
    const LevelDesc& level(int l) const {
        if (l < 0 || l >= static_cast<int>(levels.size())) {
            throw std::out_of_range("bf::gpu::BfDevice: level " + std::to_string(l) +
                                    " of " + std::to_string(levels.size()));
        }
        return levels[static_cast<size_t>(l)];
    }
    // The descriptor of local block (i, j) at a level, by the host's own
    // local indices.
    const BlockDesc& block(int l, int i, int j) const {
        const LevelDesc& d = level(l);
        if (i < 0 || i >= d.nr || j < 0 || j >= d.nc) {
            throw std::out_of_range("bf::gpu::BfDevice: block out of range at level " +
                                    std::to_string(l));
        }
        return blocks[static_cast<size_t>(d.first + static_cast<int64_t>(i) * d.nc + j)];
    }
};

}  // namespace gpu
}  // namespace bf

#endif  // H2_HAVE_GPU
