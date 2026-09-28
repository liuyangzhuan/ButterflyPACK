#pragma once
// Device data of the compression-only H2 matvec (h2_matvec.hpp): per level,
// the box and block tables and the device copies of the coupling and near
// blocks the compression kept.  The matvec runs on the device only when
// every rank kept all of its levels (decided collectively at the end of the
// compression); otherwise the host blocks serve it.

#ifdef H2_HAVE_GPU

#include "device_heap.hpp"
#include "h2_matvec_kernels.hpp"

#include <cstdint>
#include <cstdlib>
#include <unordered_map>
#include <vector>

namespace fmm {
namespace gpu {

// H2_GPU_MATVEC=0: keep no blocks on the device (host matvec).
inline bool device_matvec_enabled() {
    static const bool enabled = [] {
        const char* v = std::getenv("H2_GPU_MATVEC");
        return v == nullptr || std::atoi(v) != 0;
    }();
    return enabled;
}
// The blocks are kept while the heap stays under this fraction of its
// capacity (H2_GPU_MATVEC_KEEP_FRACTION, default 0.9).
inline double device_matvec_keep_fraction() {
    static const double fraction = [] {
        const char* v = std::getenv("H2_GPU_MATVEC_KEEP_FRACTION");
        const double f = v != nullptr ? std::atof(v) : 0.9;
        return f > 0.0 && f <= 1.0 ? f : 0.9;
    }();
    return fraction;
}
// H2_GPU_MATVEC_DIRECT=1: messages of the device matvec straight from
// device memory when MPI is CUDA-aware (default: through pinned host
// buffers; the messages are small, and registering device memory with MPI
// costs seconds on the first use).
inline bool device_matvec_direct() {
    static const bool direct = [] {
        const char* v = std::getenv("H2_GPU_MATVEC_DIRECT");
        return v != nullptr && std::atoi(v) != 0;
    }();
    return direct;
}
// H2_GPU_MATVEC_SYMMETRIC=0: keep both near blocks of a pair of local boxes
// (default: the device kernels are symmetric, so one block serves both).
inline bool device_matvec_symmetric() {
    static const bool symmetric = [] {
        const char* v = std::getenv("H2_GPU_MATVEC_SYMMETRIC");
        return v == nullptr || std::atoi(v) != 0;
    }();
    return symmetric;
}
// H2_GPU_MATVEC_CHECK=1: also run the host matvec and print the difference.
inline bool device_matvec_check() {
    static const bool check = [] {
        const char* v = std::getenv("H2_GPU_MATVEC_CHECK");
        return v != nullptr && std::atoi(v) != 0;
    }();
    return check;
}

// Messages of one step with the neighbor ranks, fixed by the compression:
// per peer, a send and a receive segment (points; a message is
// size * nrhs values at base * nrhs), and the spans that pack the sends.
struct H2Exchange {
    std::vector<int> peers;
    std::vector<int64_t> send_base, send_size, recv_base, recv_size;
    int64_t send_points = 0, recv_points = 0;
    const H2Span* d_spans = nullptr;
    int nspans = 0, max_len = 0;
};

struct H2MatvecLevel {
    bool active = false;
    std::vector<H2Box> boxes;             // local boxes, in order
    int64_t points = 0;                   // level vector length (per right-hand side)
    int64_t q_points = 0;                 // compact skeleton vectors of the local boxes
    int max_n = 0, max_k = 0, max_r = 0, max_q_cols = 0, max_near_cols = 0;
    // the blocks of the boxes (coupling blocks of every box in order, then
    // near blocks): data, columns, source box, and for a near block of two
    // local boxes its partial vector (offset, or -1)
    std::vector<const void*> block_K;
    std::vector<int> block_cols;
    std::vector<int64_t> block_source;
    std::vector<int64_t> block_partial;
    int ncoupling = 0;
    // pairs of local boxes, one kept near block each (symmetric kernels)
    std::vector<H2Pair> pairs;
    int64_t partial_points = 0;
    int max_pair_cols = 0;
    const H2Pair* d_pairs = nullptr;
    // the plan (commit_device_matvec)
    H2Exchange coupling, near;            // remote skeleton vectors; remote leaf vectors
    bool handoff_device = false;          // to the parent level, on the device
    const H2Handoff* d_handoff = nullptr;
    int nhandoff = 0;
    std::vector<char*> allocations;       // blocks, factors, tables (heap)
    const H2Box* d_boxes = nullptr;
    const H2BlockRef* d_blocks = nullptr;
    double bytes = 0.0;
};

struct DeviceMatvecStore {
    const void* tree = nullptr;
    bool building = false;   // the compression keeps blocks
    bool failed = false;     // this rank could not keep everything
    bool usable = false;     // all ranks kept everything
    std::vector<H2MatvecLevel> levels;
    double bytes = 0.0;

    bool empty() const { return levels.empty(); }
    void release() {
        DeviceHeap& heap = DeviceHeap::instance();
        for (auto& lv : levels) {
            for (char* p : lv.allocations) heap.free(p);
        }
        levels.clear();
        tree = nullptr;
        building = false;
        failed = false;
        usable = false;
        bytes = 0.0;
    }
};

inline DeviceMatvecStore& device_matvec_store() {
    static DeviceMatvecStore store;
    return store;
}

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
