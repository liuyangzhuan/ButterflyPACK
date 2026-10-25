#pragma once

#ifdef H2_HAVE_GPU
// Laying a host butterfly's factors flat on the device, and taking them back
// (GPU_BACKEND/bf_gpu; the layout is bf_device.hpp).

#include <algorithm>
#include <cstring>
#include <numeric>
#include <sstream>
#include <vector>

#include "bf_device.hpp"
#include "../h2_gpu/gpu_runtime.hpp"

namespace bf {
namespace gpu {

inline void release(BfDevice& bf) {
    if (bf.values != nullptr) {
        fmm::gpu::DeviceHeap::instance().free(bf.values);
        bf.values = nullptr;
    }
    bf.elements = 0;
    bf.blocks.clear();
    bf.levels.clear();
}

inline bool upload(const std::vector<LevelDesc>& level_rows,
                   const std::vector<std::vector<const DataType*>>& host_ptrs,
                   const std::vector<std::vector<std::pair<int, int>>>& host_shapes,
                   int level_butterfly, int level_half,
                   int M, int N, int M_loc, int N_loc, BfDevice& out) {
    const size_t nlevels = level_rows.size();
    if (host_ptrs.size() != nlevels || host_shapes.size() != nlevels) {
        throw std::invalid_argument("bf::gpu::upload: a level is missing its blocks");
    }

    // ---- lay the blocks out, tightly packed, level by level
    BfDevice next;
    next.level_butterfly = level_butterfly;
    next.level_half = level_half;
    next.M = M; next.N = N; next.M_loc = M_loc; next.N_loc = N_loc;
    next.levels.reserve(nlevels);
    int64_t at = 0;
    for (size_t l = 0; l < nlevels; ++l) {
        LevelDesc d = level_rows[l];
        const size_t want = static_cast<size_t>(d.nr) * static_cast<size_t>(d.nc);
        if (host_ptrs[l].size() != want || host_shapes[l].size() != want) {
            std::ostringstream os;
            os << "bf::gpu::upload: level " << l << " has " << host_ptrs[l].size()
               << " blocks, expected " << d.nr << " x " << d.nc;
            throw std::invalid_argument(os.str());
        }
        d.first = static_cast<int64_t>(next.blocks.size());
        for (size_t b = 0; b < want; ++b) {
            BlockDesc bd;
            // an absent local block keeps a zero-sized descriptor, so the
            // matvec skips it on its size rather than on a null test
            if (host_ptrs[l][b] != nullptr) {
                bd.rows = host_shapes[l][b].first;
                bd.cols = host_shapes[l][b].second;
                bd.lda = bd.rows;
                bd.offset = at;
                at += static_cast<int64_t>(bd.rows) * bd.cols;
            }
            next.blocks.push_back(bd);
        }
        next.levels.push_back(d);
    }
    next.elements = at;
    if (at == 0) { out = std::move(next); return true; }

    // ---- one resident block of the shared arena, tagged as matvec data
    auto& heap = fmm::gpu::DeviceHeap::instance();
    heap.ensure_initialized();
    DataType* dev = reinterpret_cast<DataType*>(heap.try_alloc_resident(
        static_cast<size_t>(at) * sizeof(DataType), fmm::gpu::ResidentKind::Matvec));
    if (dev == nullptr) return false;  // no room: the caller stays on the host

    // ---- stage into one host buffer, then a single copy: a copy per block
    // would cost more in launches than the factors are large (F12 and F21 are
    // about 0.1 GB a rank, spread over thousands of small blocks)
    std::vector<DataType> stage(static_cast<size_t>(at));
    for (size_t l = 0; l < nlevels; ++l) {
        const LevelDesc& d = next.levels[l];
        const size_t want = static_cast<size_t>(d.nr) * static_cast<size_t>(d.nc);
        for (size_t b = 0; b < want; ++b) {
            const BlockDesc& bd = next.blocks[static_cast<size_t>(d.first) + b];
            if (bd.rows == 0 || bd.cols == 0) continue;
            const DataType* src = host_ptrs[l][b];
            // the host block is column major with lda == rows (its matrix is
            // allocated exactly), so it copies as one run
            std::memcpy(stage.data() + bd.offset, src,
                        static_cast<size_t>(bd.rows) * bd.cols * sizeof(DataType));
        }
    }
    fmm::gpu::check_cuda(cudaMemcpy(dev, stage.data(),
                                    static_cast<size_t>(at) * sizeof(DataType),
                                    cudaMemcpyHostToDevice),
                         "bf::gpu::upload");
    next.values = dev;
    release(out);          // a previous residency of the same operator
    out = std::move(next);
    return true;
}

inline void download(const BfDevice& bf,
                     const std::vector<std::vector<DataType*>>& host_ptrs) {
    if (!bf.resident() || bf.elements == 0) return;
    std::vector<DataType> stage(static_cast<size_t>(bf.elements));
    fmm::gpu::check_cuda(cudaMemcpy(stage.data(), bf.values,
                                    static_cast<size_t>(bf.elements) * sizeof(DataType),
                                    cudaMemcpyDeviceToHost),
                         "bf::gpu::download");
    for (size_t l = 0; l < bf.levels.size(); ++l) {
        const LevelDesc& d = bf.levels[l];
        const size_t want = static_cast<size_t>(d.nr) * static_cast<size_t>(d.nc);
        if (host_ptrs[l].size() != want) {
            throw std::invalid_argument("bf::gpu::download: a level's block count changed");
        }
        for (size_t b = 0; b < want; ++b) {
            const BlockDesc& bd = bf.blocks[static_cast<size_t>(d.first) + b];
            if (bd.rows == 0 || bd.cols == 0 || host_ptrs[l][b] == nullptr) continue;
            std::memcpy(host_ptrs[l][b], stage.data() + bd.offset,
                        static_cast<size_t>(bd.rows) * bd.cols * sizeof(DataType));
        }
    }
}


}  // namespace gpu
}  // namespace bf

#endif  // H2_HAVE_GPU
