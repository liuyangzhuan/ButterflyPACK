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
#include "kernel_eval.cuh"

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

// out(i, c) = sum over the draws (row, b) with destination i, in (row, b)
// order, of sign * value(row, c).  A block owns 32 columns and the
// destinations [i0, i0 + dest_block); the output tile lives in shared memory
// and 32 rows of values at a time are formed once (kernel entries, or source
// rows read once) and applied by the warps, warp w owning the destinations
// i = w mod kSketchOwners.  Every destination sums in the draw order, so the
// result matches the per-destination lists exactly.
constexpr int kTileCols = 32;
constexpr int kTileRows = 64;
constexpr int kPitch = kTileCols + 1;  // padded rows of the shared tiles
constexpr int kRowsPerWarp = kTileRows / kSketchOwners;

// Values of this thread's rows of the chunk starting at r0 (rows
// r0 + warp + q * kSketchOwners, its column c): the row indices first, then
// the values, so the loads of a thread are independent.
template<typename T, bool kKernelRows, int Kind>
__device__ __forceinline__ void chunk_values(const OrderedSketchItemT<T>& item, const KernelSpec& spec,
                                             const PointTable& points, const double col_xyz[3], int64_t col_id,
                                             bool column, int c, int warp, int r0, T (&v)[kRowsPerWarp]) {
    int src_row[kRowsPerWarp];
#pragma unroll
    for (int q = 0; q < kRowsPerWarp; ++q) {
        const int row = r0 + warp + q * kSketchOwners;
        src_row[q] = row < item.rows ? (kKernelRows ? item.row_slots[row] : (item.row_index ? item.row_index[row] : row))
                                     : -1;
    }
#pragma unroll
    for (int q = 0; q < kRowsPerWarp; ++q) {
        v[q] = T(0.0);
        if (src_row[q] >= 0 && column) {
            if (kKernelRows) {
                v[q] = kernel_value<T, Kind>(spec, col_xyz, col_id, points.xyz + 3 * static_cast<int64_t>(src_row[q]),
                                       points.ids[src_row[q]]);
            } else {
                v[q] = item.src[c + static_cast<int64_t>(src_row[q]) * item.row_stride];
            }
        }
    }
}

// y + s v for a real s: y + s v componentwise for complex values (the host's
// axpy with a real coefficient)
__device__ __forceinline__ double sketch_fma(double s, double v, double y) { return fma(s, v, y); }
__device__ __forceinline__ dcomplex sketch_fma(double s, dcomplex v, dcomplex y) {
    return dcomplex(fma(s, v.re, y.re), fma(s, v.im, y.im));
}

template<typename T, bool kKernelRows, int Kind>
__global__ void __launch_bounds__(32 * kSketchOwners)
ordered_sketch_kernel(const OrderedSketchItemT<T>* items, KernelSpec spec, PointTable points, int dest_block) {
    extern __shared__ __align__(16) unsigned char shared_bytes[];
    T* shared = reinterpret_cast<T*>(shared_bytes);
    const OrderedSketchItemT<T> item = items[blockIdx.x];
    const int c0 = blockIdx.y * kTileCols;
    const int i0 = blockIdx.z * dest_block;
    if (c0 >= item.ncols || i0 >= item.d) return;
    const int i1 = min(item.d, i0 + dest_block);
    T* tile = shared;                                               // (i1 - i0) x kPitch
    T* values = shared + static_cast<size_t>(dest_block) * kPitch;  // kTileRows x kPitch
    const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5;
    for (int e = threadIdx.x; e < (i1 - i0) * kPitch; e += blockDim.x) tile[e] = T(0.0);
    const int c = c0 + lane;
    const bool column = c < item.ncols;
    double col_xyz[3] = {0.0, 0.0, 0.0};
    int64_t col_id = 0;
    if (kKernelRows && column) {
        const int sc = item.col_base + c;
        for (int k = 0; k < 3; ++k) col_xyz[k] = points.xyz[3 * static_cast<int64_t>(sc) + k];
        col_id = points.ids[sc];
    }
    // this warp's entries, 32 at a time: one coalesced load, then broadcasts
    const int* list = item.entries + item.ptr[warp];
    const int* end = item.entries + item.ptr[warp + 1];
    // +-scale by the entry's sign bit (the same fma as with a branch)
    const unsigned long long scale_bits = static_cast<unsigned long long>(__double_as_longlong(item.scale));
    auto signed_scale = [&](int e) {
        return __longlong_as_double(static_cast<long long>(
            scale_bits ^ (static_cast<unsigned long long>(static_cast<unsigned>(e) & 0x80000000u) << 32)));
    };
    const bool all_dests = i0 == 0 && i1 == item.d;  // no range test
    int batch = 0, batch_count = 0, batch_pos = 0;
    T next[kRowsPerWarp];
    chunk_values<T, kKernelRows, Kind>(item, spec, points, col_xyz, col_id, column, c, warp, 0, next);
    for (int r0 = 0; r0 < item.rows; r0 += kTileRows) {
        __syncthreads();  // the value tile is free
#pragma unroll
        for (int q = 0; q < kRowsPerWarp; ++q) values[(warp + q * kSketchOwners) * kPitch + lane] = next[q];
        __syncthreads();
        // the next chunk's loads are in flight while this one is applied
        if (r0 + kTileRows < item.rows) {
            chunk_values<T, kKernelRows, Kind>(item, spec, points, col_xyz, col_id, column, c, warp, r0 + kTileRows, next);
        }
        const int r1 = r0 + kTileRows;
        const T* vrow = values + lane;
        T* ycol = tile + lane;
        while (true) {
            if (batch_pos == batch_count) {
                if (list >= end) break;
                batch_count = end - list < 32 ? static_cast<int>(end - list) : 32;
                batch = lane < batch_count ? list[lane] : 0;
                list += batch_count;
                batch_pos = 0;
            }
            // entries are in row order: those of this chunk are a run from batch_pos
            const unsigned here = __ballot_sync(0xffffffffu, lane >= batch_pos && lane < batch_count &&
                                                                 (batch & kEntryRowMask) < r1);
            const int stop = batch_pos + __popc(here);
            int j = batch_pos;
            if (all_dests) {
                // two entries at a time: different destinations update independently
                for (; j + 1 < stop; j += 2) {
                    const int e1 = __shfl_sync(0xffffffffu, batch, j);
                    const int e2 = __shfl_sync(0xffffffffu, batch, j + 1);
                    const int a1 = ((e1 >> kEntryRowBits) & kEntryDestMask) * kPitch;
                    const int a2 = ((e2 >> kEntryRowBits) & kEntryDestMask) * kPitch;
                    const T v1 = vrow[((e1 & kEntryRowMask) - r0) * kPitch];
                    const T v2 = vrow[((e2 & kEntryRowMask) - r0) * kPitch];
                    if (a1 != a2) {
                        const T y1 = ycol[a1], y2 = ycol[a2];
                        ycol[a1] = sketch_fma(signed_scale(e1), v1, y1);
                        ycol[a2] = sketch_fma(signed_scale(e2), v2, y2);
                    } else {
                        ycol[a1] = sketch_fma(signed_scale(e2), v2, sketch_fma(signed_scale(e1), v1, ycol[a1]));
                    }
                }
                if (j < stop) {
                    const int e = __shfl_sync(0xffffffffu, batch, j);
                    const int a = ((e >> kEntryRowBits) & kEntryDestMask) * kPitch;
                    ycol[a] = sketch_fma(signed_scale(e), vrow[((e & kEntryRowMask) - r0) * kPitch], ycol[a]);
                    ++j;
                }
            } else {
                for (; j < stop; ++j) {
                    const int e = __shfl_sync(0xffffffffu, batch, j);
                    const int i = (e >> kEntryRowBits) & kEntryDestMask;
                    if (i >= i0 && i < i1) {
                        T& y = ycol[(i - i0) * kPitch];
                        y = sketch_fma(signed_scale(e), vrow[((e & kEntryRowMask) - r0) * kPitch], y);
                    }
                }
            }
            batch_pos = stop;
            if (batch_pos < batch_count) break;  // the rest belongs to later chunks
        }
    }
    __syncthreads();
    for (int col = warp; col < kTileCols; col += kSketchOwners) {
        if (c0 + col >= item.ncols) break;
        for (int i = lane; i < i1 - i0; i += 32) {
            item.out[(i0 + i) + static_cast<int64_t>(c0 + col) * item.ldo] = tile[i * kPitch + col];
        }
    }
}

void check_launch(const char* what) {
    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA launch failed in ") + what + ": " + cudaGetErrorString(status));
    }
}

}  // namespace

template<typename T>
void launch_ordered_sketch(const OrderedSketchItemT<T>* items, int count, int max_d, int max_cols, bool kernel_rows,
                           KernelSpec spec, PointTable points, cudaStream_t stream) {
    if (count <= 0 || max_d <= 0 || max_cols <= 0) return;
    // All destinations in one block when the tile fits the shared memory of
    // a block: a split recomputes every value and rescans every list per part.
    static const int max_block = [] {
        int device = 0, optin = 0;
        cudaGetDevice(&device);
        cudaDeviceGetAttribute(&optin, cudaDevAttrMaxSharedMemoryPerBlockOptin, device);
        return optin / static_cast<int>(kPitch * sizeof(T)) - kTileRows;
    }();
    const int parts = (max_d + max_block - 1) / max_block;
    const int dest_block = (max_d + parts - 1) / parts;
    const size_t shared = static_cast<size_t>(dest_block + kTileRows) * kPitch * sizeof(T);
    const dim3 grid(static_cast<unsigned>(count), static_cast<unsigned>((max_cols + kTileCols - 1) / kTileCols),
                    static_cast<unsigned>((max_d + dest_block - 1) / dest_block));
    // (rows read from stored values do not evaluate the kernel)
    constexpr int kPlain = is_complex_scalar<T> ? 2 : 1;
    auto launch = [&](auto kernel) {
        cudaFuncSetAttribute(kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, static_cast<int>(shared));
        kernel<<<grid, 32 * kSketchOwners, shared, stream>>>(items, spec, points, dest_block);
    };
    if (!kernel_rows) {
        launch(ordered_sketch_kernel<T, false, kPlain>);
    } else if (kernel_kind_of<T>(spec, "launch_ordered_sketch") == 3) {
        if constexpr (is_complex_scalar<T>) launch(ordered_sketch_kernel<T, true, 3>);
    } else if (spec.kind == 4) {
        if constexpr (!is_complex_scalar<T>) launch(ordered_sketch_kernel<T, true, 4>);
    } else {
        launch(ordered_sketch_kernel<T, true, kPlain>);
    }
    check_launch("ordered_sketch_kernel");
}
template void launch_ordered_sketch<double>(const OrderedSketchItemT<double>*, int, int, int, bool, KernelSpec,
                                            PointTable, cudaStream_t);
template void launch_ordered_sketch<dcomplex>(const OrderedSketchItemT<dcomplex>*, int, int, int, bool, KernelSpec,
                                              PointTable, cudaStream_t);

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
