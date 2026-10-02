// Device construction of the sparse-sign sketch lists of the H2 Color GPU
// backend (see SketchListsItem / SourceListsItem in device_kernels.hpp).
//
// The draws reproduce the host's: std::mt19937_64 seeded with morton + 1,
// sk draws per sketched row in row order, destination block_start[b] +
// rand % block_range[b], positive when the top bit is set.  One warp runs a
// box's generator: the twist runs in its three dependency phases and the 312
// outputs of a state are tempered in parallel.  The destination lists come
// from a stable counting sort (each lane owns a contiguous range of the
// input), so entries of a destination keep the host's accumulation order.

#include "device_kernels.hpp"
#include "evaluator_device.hpp"

#include <stdexcept>
#include <string>

namespace fmm {
namespace gpu {

namespace {

constexpr int kStateWords = 312;  // std::mt19937_64: n
constexpr int kShift = 156;       // m
constexpr uint64_t kMatrixA = 0xB5026F5AA96619E9ULL;
constexpr uint64_t kUpper = ~((uint64_t{1} << 31) - 1);
constexpr uint64_t kLower = (uint64_t{1} << 31) - 1;

__device__ __forceinline__ uint64_t mt_temper(uint64_t y) {
    y ^= (y >> 29) & 0x5555555555555555ULL;
    y ^= (y << 17) & 0x71D67FFFEDA60000ULL;
    y ^= (y << 37) & 0xFFF7EEE000000000ULL;
    y ^= y >> 43;
    return y;
}

__device__ __forceinline__ uint64_t mt_mix(uint64_t upper, uint64_t lower, uint64_t shifted) {
    const uint64_t y = (upper & kUpper) | (lower & kLower);
    return shifted ^ (y >> 1) ^ ((y & 1ULL) ? kMatrixA : 0ULL);
}

// Regenerate the state in place (std::mersenne_twister_engine::_M_gen_rand):
// words [0, n - m) read only old words, [n - m, n - 1) read words of the
// first phase, and word n - 1 reads words 0 and m - 1 of the first phase.
__device__ void mt_twist(uint64_t* mt, int lane) {
    uint64_t next[(kStateWords - kShift + 31) / 32];
    int q = 0;
    for (int i = lane; i < kStateWords - kShift; i += 32, ++q) next[q] = mt_mix(mt[i], mt[i + 1], mt[i + kShift]);
    __syncwarp();
    q = 0;
    for (int i = lane; i < kStateWords - kShift; i += 32, ++q) mt[i] = next[q];
    __syncwarp();
    q = 0;
    for (int i = kStateWords - kShift + lane; i < kStateWords - 1; i += 32, ++q) {
        next[q] = mt_mix(mt[i], mt[i + 1], mt[i - (kStateWords - kShift)]);
    }
    __syncwarp();
    q = 0;
    for (int i = kStateWords - kShift + lane; i < kStateWords - 1; i += 32, ++q) mt[i] = next[q];
    __syncwarp();
    if (lane == 0) mt[kStateWords - 1] = mt_mix(mt[kStateWords - 1], mt[0], mt[kShift - 1]);
    __syncwarp();
}

// draws[q] of the sketched row q / sk, draw q % sk: destination | sign bit
__global__ void sketch_draws_kernel(const SketchListsItem* items) {
    __shared__ uint64_t mt[kStateWords];
    const SketchListsItem item = items[blockIdx.x];
    const int lane = threadIdx.x;
    if (lane == 0) {
        mt[0] = item.seed;
        for (int i = 1; i < kStateWords; ++i) {
            mt[i] = 6364136223846793005ULL * (mt[i - 1] ^ (mt[i - 1] >> 62)) + static_cast<uint64_t>(i);
        }
    }
    __syncwarp();
    const int sk = item.sk;
    const int block_size = item.d / sk;
    const int64_t total = static_cast<int64_t>(item.rows) * sk;
    for (int64_t base = 0; base < total; base += kStateWords) {
        mt_twist(mt, lane);
        for (int i = lane; i < kStateWords && base + i < total; i += 32) {
            const int64_t q = base + i;
            const uint64_t y = mt_temper(mt[i]);
            const int b = static_cast<int>(q % sk);
            const int start = b * block_size;
            const int range = (b == sk - 1 ? item.d : (b + 1) * block_size) - start;
            const int dest = start + static_cast<int>(y % static_cast<uint64_t>(range));
            const bool positive = ((y >> 63) & 1ULL) != 0;
            item.draws[q] = positive ? dest : static_cast<int>(static_cast<unsigned>(dest) | 0x80000000u);
        }
        __syncwarp();
    }
}

// rows[row_base + i] = slot_base + (list ? list[i] : i) for every row block
__global__ void sketch_row_slots_kernel(const SketchListsItem* items, const char* meta) {
    const SketchListsItem item = items[blockIdx.x];
    const RowBlockDesc* blocks = reinterpret_cast<const RowBlockDesc*>(meta + item.blocks_offset);
    for (int j = blockIdx.y; j < item.num_blocks; j += gridDim.y) {
        const RowBlockDesc blk = blocks[j];
        for (int i = threadIdx.x; i < blk.count; i += blockDim.x) {
            item.rows_out[blk.row_base + i] = blk.slot_base + (blk.list ? blk.list[i] : i);
        }
    }
}

// Stable counting sort of `total` inputs into d destination lists, by one
// warp.  Input f yields (dest, value) through `get`; lane l owns inputs
// [l * chunk, (l + 1) * chunk).  hist is scratch of 32 * d ints.
template<typename Get>
__device__ void warp_counting_sort(int64_t total, int d, int* hist, int* ptr, int* entries, Get get) {
    const int lane = threadIdx.x & 31;
    const int64_t chunk = (total + 31) / 32;
    const int64_t f0 = lane * chunk, f1 = f0 + chunk < total ? f0 + chunk : total;
    int* mine = hist + static_cast<int64_t>(lane) * d;
    for (int b = 0; b < d; ++b) mine[b] = 0;
    for (int64_t f = f0; f < f1; ++f) {
        int dest, value;
        get(f, dest, value);
        ++mine[dest];
    }
    __syncwarp();
    // per destination: exclusive offsets of the lanes; totals into ptr[dest + 1]
    for (int b = lane; b < d; b += 32) {
        int running = 0;
        for (int l = 0; l < 32; ++l) {
            const int c = hist[static_cast<int64_t>(l) * d + b];
            hist[static_cast<int64_t>(l) * d + b] = running;
            running += c;
        }
        ptr[b + 1] = running;
    }
    __syncwarp();
    if (lane == 0) {
        ptr[0] = 0;
        for (int b = 0; b < d; ++b) ptr[b + 1] += ptr[b];
    }
    __syncwarp();
    for (int b = 0; b < d; ++b) mine[b] += ptr[b];
    for (int64_t f = f0; f < f1; ++f) {
        int dest, value;
        get(f, dest, value);
        entries[mine[dest]++] = value;
    }
}

__device__ __forceinline__ int pack_entry(int row, int dest, int draw) {
    return static_cast<int>(static_cast<unsigned>(row) | (static_cast<unsigned>(dest) << kEntryRowBits) |
                            (static_cast<unsigned>(draw) & 0x80000000u));
}

// Kernel-row lists: the box's draws in (row, draw) order, partitioned by
// owner warp (destination mod kSketchOwners).
__global__ void sketch_row_lists_kernel(const SketchListsItem* items) {
    const SketchListsItem item = items[blockIdx.x];
    const int sk = item.sk;
    const int* draws = item.draws;
    warp_counting_sort(static_cast<int64_t>(item.rows) * sk, kSketchOwners, item.hist, item.ptr, item.entries,
                       [&](int64_t q, int& owner, int& value) {
                           const int draw = draws[q];
                           const int dest = draw & 0x7fffffff;
                           owner = dest % kSketchOwners;
                           value = pack_entry(static_cast<int>(q / sk), dest, draw);
                       });
}

// Stored-row lists of one fill source: its sketched rows (the runs, in
// order) numbered 0, 1, ...; row_index[row] = the source row it reads,
// src_offset + (list ? list[i] : i); draws partitioned by owner warp.
__global__ void sketch_source_lists_kernel(const SourceListsItem* items, const SketchListsItem* boxes,
                                           const char* meta) {
    __shared__ int64_t run_start[kMaxSourceRuns + 1];  // in draws
    __shared__ int run_row[kMaxSourceRuns + 1];        // in rows
    const SourceListsItem item = items[blockIdx.x];
    const SketchListsItem& box = boxes[item.box];
    const int sk = box.sk;
    const RunDesc* runs = reinterpret_cast<const RunDesc*>(meta + item.runs_offset);
    if (threadIdx.x == 0) {
        run_start[0] = 0;
        run_row[0] = 0;
        for (int r = 0; r < item.num_runs; ++r) {
            run_start[r + 1] = run_start[r] + static_cast<int64_t>(runs[r].count) * sk;
            run_row[r + 1] = run_row[r] + runs[r].count;
        }
    }
    __syncwarp();
    for (int r = 0; r < item.num_runs; ++r) {
        const RunDesc rd = runs[r];
        for (int i = threadIdx.x; i < rd.count; i += 32) {
            item.row_index[run_row[r] + i] = rd.src_offset + (rd.list ? rd.list[i] : i);
        }
    }
    const int* draws = box.draws;
    const int num_runs = item.num_runs;
    // inputs are visited in order within a lane, so the run search restarts rarely
    int run = 0;
    warp_counting_sort(run_start[num_runs], kSketchOwners, item.hist, item.ptr, item.entries,
                       [&](int64_t f, int& owner, int& value) {
                           if (f < run_start[run]) run = 0;
                           while (f >= run_start[run + 1]) ++run;
                           const RunDesc rd = runs[run];
                           const int64_t local = f - run_start[run];
                           const int i = static_cast<int>(local / sk), b = static_cast<int>(local % sk);
                           const int draw = draws[static_cast<int64_t>(rd.row_base + i) * sk + b];
                           const int dest = draw & 0x7fffffff;
                           owner = dest % kSketchOwners;
                           value = pack_entry(run_row[run] + i, dest, draw);
                       });
}

void check_launch(const char* what) {
    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA launch failed in ") + what + ": " + cudaGetErrorString(status));
    }
}

}  // namespace

// The sketch of stored rows (the kernel itself, shared with the sketch of
// kernel rows: bpack_gpu_kernels.cuh).
template<typename T>
void launch_ordered_sketch_stored(const OrderedSketchItemT<T>* items, int count, int max_d, int max_cols,
                                  cudaStream_t stream) {
    if (count <= 0 || max_d <= 0 || max_cols <= 0) return;
    const SketchGeometry g = sketch_geometry(count, max_d, max_cols, sizeof(T));
    auto kernel = ordered_sketch_kernel<T, false, NoEntry<T>>;
    cudaFuncSetAttribute(kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, static_cast<int>(g.shared));
    kernel<<<g.grid, g.block, g.shared, stream>>>(items, NoEntry<T>{}, PointTable{}, g.dest_block);
    check_launch("ordered_sketch_kernel");
}
template void launch_ordered_sketch_stored<double>(const OrderedSketchItemT<double>*, int, int, int, cudaStream_t);
template void launch_ordered_sketch_stored<dcomplex>(const OrderedSketchItemT<dcomplex>*, int, int, int, cudaStream_t);

void launch_sketch_lists(const SketchListsItem* boxes, int num_boxes, int max_blocks,
                         const SourceListsItem* sources, int num_sources, const char* meta, cudaStream_t stream) {
    if (num_boxes <= 0) return;
    sketch_draws_kernel<<<static_cast<unsigned>(num_boxes), 32, 0, stream>>>(boxes);
    check_launch("sketch_draws_kernel");
    const dim3 grid(static_cast<unsigned>(num_boxes), static_cast<unsigned>(max_blocks < 1 ? 1 : (max_blocks < 64 ? max_blocks : 64)));
    sketch_row_slots_kernel<<<grid, 128, 0, stream>>>(boxes, meta);
    check_launch("sketch_row_slots_kernel");
    sketch_row_lists_kernel<<<static_cast<unsigned>(num_boxes), 32, 0, stream>>>(boxes);
    check_launch("sketch_row_lists_kernel");
    if (num_sources > 0) {
        sketch_source_lists_kernel<<<static_cast<unsigned>(num_sources), 32, 0, stream>>>(sources, boxes, meta);
        check_launch("sketch_source_lists_kernel");
    }
}

}  // namespace gpu
}  // namespace fmm
