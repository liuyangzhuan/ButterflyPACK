// Device construction of the sparse-sign sketch lists of the H2 Color GPU
// backend (see SketchListsItem / SourceListsItem in device_kernels.hpp).
//
// The draws reproduce the host's: std::mt19937_64 seeded with morton + 1,
// sk draws per sketched row in row order, destination block_start[b] +
// rand % block_range[b], positive when the top bit is set.  One block runs
// a box's generator: the twist runs in its three dependency phases, a word
// per thread, and the 312 outputs of a state are tempered in parallel.  The
// destination lists come from a block-wide stable partition (each warp owns
// a contiguous range of the input, read in order), so entries of a
// destination keep the host's accumulation order.

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

// One block of kListThreads per box or fill source: the generator's state
// regenerates with one word per thread, and the destination lists come
// from a block-wide stable partition.  (One warp per box, with its lanes'
// histograms in global memory, took 1-6 ms per wave at a 4096-box leaf: a
// serial chain of thousands of latency-bound steps per box.)
constexpr int kListThreads = 512;
constexpr int kListWarps = kListThreads / 32;
static_assert(kListThreads >= kStateWords, "one state word per thread");

// draws[q] of the sketched row q / sk, draw q % sk: destination | sign bit.
// The state regenerates in its three dependency phases
// (std::mersenne_twister_engine::_M_gen_rand): words [0, n - m) read only
// old words, [n - m, n - 1) read old words and the new words of the first
// phase, and word n - 1 reads the old word and new words 0 and m - 1; the
// 312 outputs of a state are tempered and drawn in parallel.
__global__ void __launch_bounds__(kListThreads) sketch_draws_kernel(const SketchListsItem* items) {
    __shared__ uint64_t state[2][kStateWords];
    const SketchListsItem item = items[blockIdx.x];
    const int t = threadIdx.x;
    if (t == 0) {
        uint64_t* mt = state[0];
        mt[0] = item.seed;
        for (int i = 1; i < kStateWords; ++i) {
            mt[i] = 6364136223846793005ULL * (mt[i - 1] ^ (mt[i - 1] >> 62)) + static_cast<uint64_t>(i);
        }
    }
    __syncthreads();
    const int sk = item.sk;
    const int block_size = item.d / sk;
    const int64_t total = static_cast<int64_t>(item.rows) * sk;
    int cur = 0;
    for (int64_t base = 0; base < total; base += kStateWords) {
        const uint64_t* mt = state[cur];
        uint64_t* nxt = state[cur ^ 1];
        if (t < kStateWords - kShift) nxt[t] = mt_mix(mt[t], mt[t + 1], mt[t + kShift]);
        __syncthreads();
        if (t < kShift - 1) {
            const int i = kStateWords - kShift + t;
            nxt[i] = mt_mix(mt[i], mt[i + 1], nxt[i - (kStateWords - kShift)]);
        }
        __syncthreads();
        if (t == 0) nxt[kStateWords - 1] = mt_mix(mt[kStateWords - 1], nxt[0], nxt[kShift - 1]);
        __syncthreads();
        // (the next twist reads these new words and writes the old buffer,
        // whose last reads ended before the sync above)
        if (t < kStateWords && base + t < total) {
            const int64_t q = base + t;
            const uint64_t y = mt_temper(nxt[t]);
            const int b = static_cast<int>(q % sk);
            const int start = b * block_size;
            const int range = (b == sk - 1 ? item.d : (b + 1) * block_size) - start;
            const int dest = start + static_cast<int>(y % static_cast<uint64_t>(range));
            const bool positive = ((y >> 63) & 1ULL) != 0;
            item.draws[q] = positive ? dest : static_cast<int>(static_cast<unsigned>(dest) | 0x80000000u);
        }
        cur ^= 1;
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

// Stable partition of `total` inputs into the kSketchOwners destination
// lists of one box or fill source, by one block.  Input f yields (owner,
// value) through `get`; warp w takes the inputs [w * chunk, (w + 1) * chunk),
// 32 at a time in order, so every list keeps the input order (the host's
// accumulation order).  Two passes over the inputs: the counts, then the
// scatter (each warp's run of a list starts after the lower warps' runs).
struct PartitionScratch {
    int counts[kListWarps][kSketchOwners];
    int base[kListWarps][kSketchOwners];
    int total_of[kSketchOwners];
    int start[kSketchOwners + 1];
};

template<typename Get>
__device__ void block_partition(PartitionScratch& s, int total, int* ptr, int* entries, Get get) {
    const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5;
    const unsigned below = (1u << lane) - 1u;
    const int chunk = ((total + kListWarps - 1) / kListWarps + 31) / 32 * 32;
    const int f0 = min(warp * chunk, total), f1 = min(f0 + chunk, total);
    int count[kSketchOwners];
#pragma unroll
    for (int b = 0; b < kSketchOwners; ++b) count[b] = 0;
    for (int f = f0; f < f1; f += 32) {
        int owner = -1, value = 0;
        if (f + lane < f1) get(f + lane, owner, value);
#pragma unroll
        for (int b = 0; b < kSketchOwners; ++b) count[b] += __popc(__ballot_sync(0xffffffffu, owner == b));
    }
    if (lane == 0) {
#pragma unroll
        for (int b = 0; b < kSketchOwners; ++b) s.counts[warp][b] = count[b];
    }
    __syncthreads();
    if (threadIdx.x < kSketchOwners) {
        const int b = threadIdx.x;
        int running = 0;
        for (int w = 0; w < kListWarps; ++w) {
            s.base[w][b] = running;
            running += s.counts[w][b];
        }
        s.total_of[b] = running;
    }
    __syncthreads();
    if (threadIdx.x == 0) {
        int running = 0;
        for (int b = 0; b < kSketchOwners; ++b) {
            s.start[b] = running;
            ptr[b] = running;
            running += s.total_of[b];
        }
        s.start[kSketchOwners] = running;
        ptr[kSketchOwners] = running;
    }
    __syncthreads();
    int at[kSketchOwners];
#pragma unroll
    for (int b = 0; b < kSketchOwners; ++b) at[b] = s.start[b] + s.base[warp][b];
    for (int f = f0; f < f1; f += 32) {
        int owner = -1, value = 0;
        if (f + lane < f1) get(f + lane, owner, value);
#pragma unroll
        for (int b = 0; b < kSketchOwners; ++b) {
            const unsigned mask = __ballot_sync(0xffffffffu, owner == b);
            if (owner == b) entries[at[b] + __popc(mask & below)] = value;
            at[b] += __popc(mask);
        }
    }
}

__device__ __forceinline__ int pack_entry(int row, int dest, int draw) {
    return static_cast<int>(static_cast<unsigned>(row) | (static_cast<unsigned>(dest) << kEntryRowBits) |
                            (static_cast<unsigned>(draw) & 0x80000000u));
}

// Kernel-row lists: the box's draws in (row, draw) order, partitioned by
// owner warp (destination mod kSketchOwners).  (rows * sk fits an int: the
// host bounds the rows by kEntryRowMask.)
__global__ void __launch_bounds__(kListThreads) sketch_row_lists_kernel(const SketchListsItem* items) {
    __shared__ PartitionScratch scratch;
    const SketchListsItem item = items[blockIdx.x];
    const int sk = item.sk;
    const int* draws = item.draws;
    block_partition(scratch, item.rows * sk, item.ptr, item.entries, [&](int q, int& owner, int& value) {
        const int draw = draws[q];
        const int dest = draw & 0x7fffffff;
        owner = dest % kSketchOwners;
        value = pack_entry(q / sk, dest, draw);
    });
}

// Stored-row lists of one fill source: its sketched rows (the runs, in
// order) numbered 0, 1, ...; row_index[row] = the source row it reads,
// src_offset + (list ? list[i] : i); draws partitioned by owner warp.
__global__ void __launch_bounds__(kListThreads) sketch_source_lists_kernel(const SourceListsItem* items,
                                                                           const SketchListsItem* boxes,
                                                                           const char* meta) {
    __shared__ PartitionScratch scratch;
    __shared__ RunDesc run[kMaxSourceRuns];
    __shared__ int run_start[kMaxSourceRuns + 1];  // in draws
    __shared__ int run_row[kMaxSourceRuns + 1];    // in rows
    const SourceListsItem item = items[blockIdx.x];
    const int sk = boxes[item.box].sk;
    const int* draws = boxes[item.box].draws;
    const RunDesc* runs = reinterpret_cast<const RunDesc*>(meta + item.runs_offset);
    const int num_runs = item.num_runs;
    for (int r = threadIdx.x; r < num_runs; r += blockDim.x) run[r] = runs[r];
    if (threadIdx.x == 0) {
        run_start[0] = 0;
        run_row[0] = 0;
        for (int r = 0; r < num_runs; ++r) {
            run_start[r + 1] = run_start[r] + runs[r].count * sk;
            run_row[r + 1] = run_row[r] + runs[r].count;
        }
    }
    __syncthreads();
    for (int r = 0; r < num_runs; ++r) {
        const RunDesc rd = run[r];
        for (int i = threadIdx.x; i < rd.count; i += blockDim.x) {
            item.row_index[run_row[r] + i] = rd.src_offset + (rd.list ? rd.list[i] : i);
        }
    }
    // a thread's inputs come in order, so the run search restarts rarely
    int cur = 0;
    block_partition(scratch, run_start[num_runs], item.ptr, item.entries, [&](int f, int& owner, int& value) {
        if (f < run_start[cur]) cur = 0;
        while (f >= run_start[cur + 1]) ++cur;
        const int local = f - run_start[cur];
        const int i = local / sk, b = local - i * sk;
        const int draw = draws[(run[cur].row_base + i) * sk + b];
        const int dest = draw & 0x7fffffff;
        owner = dest % kSketchOwners;
        value = pack_entry(run_row[cur] + i, dest, draw);
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
    sketch_draws_kernel<<<static_cast<unsigned>(num_boxes), kListThreads, 0, stream>>>(boxes);
    check_launch("sketch_draws_kernel");
    const dim3 grid(static_cast<unsigned>(num_boxes), static_cast<unsigned>(max_blocks < 1 ? 1 : (max_blocks < 64 ? max_blocks : 64)));
    sketch_row_slots_kernel<<<grid, 128, 0, stream>>>(boxes, meta);
    check_launch("sketch_row_slots_kernel");
    sketch_row_lists_kernel<<<static_cast<unsigned>(num_boxes), kListThreads, 0, stream>>>(boxes);
    check_launch("sketch_row_lists_kernel");
    if (num_sources > 0) {
        sketch_source_lists_kernel<<<static_cast<unsigned>(num_sources), kListThreads, 0, stream>>>(sources, boxes, meta);
        check_launch("sketch_source_lists_kernel");
    }
}

}  // namespace gpu
}  // namespace fmm
