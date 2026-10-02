#pragma once
// Device side of the application's evaluator (evaluator.hpp, evaluator.cu):
// the helpers of the block and entry-list forms, the stored-rows sketch, and
// entry evaluators compiled with NVRTC.  CUDA only (no MPI): included by
// evaluator.cu.

#ifdef H2_HAVE_GPU

#include "device_kernels.hpp"

#include <cstdint>
#include <memory>
#include <string>

namespace fmm {
namespace gpu {

// ---- device helpers (evaluator.cu)

// A block's index lists (as an eval item's) and where its points' ids and
// coordinates go in the gathered arrays.
struct BlockGatherItem {
    IndexList rows;
    IndexList cols;
    int m;
    int n;
    int64_t row_at;
    int64_t col_at;
};
// ids[at + i] = points.ids[slot of i]; coords (if not null): dim values per point
void launch_gather_block_points(const BlockGatherItem* items, int count, int max_mn, const char* meta,
                                PointTable points, int64_t* ids, double* coords, cudaStream_t stream);

// A block of an entry list: entries at .. at + m n - 1, column by column.
struct FlatBlockItem {
    int m;
    int n;
    int ld;
    void* out;
    const int64_t* row_ids;
    const int64_t* col_ids;
    int64_t at;
};
void launch_flatten_entries(const FlatBlockItem* items, int count, int max_m, int max_n, int64_t* rows, int64_t* cols,
                            cudaStream_t stream);
template<typename T>
void launch_scatter_entries(const FlatBlockItem* items, int count, int max_m, int max_n, const T* values,
                            cudaStream_t stream);

// The ordered sketch of stored rows (device_sketch.cu).
template<typename T>
void launch_ordered_sketch_stored(const OrderedSketchItemT<T>* items, int count, int max_d, int max_cols,
                                  cudaStream_t stream);

// An entry evaluator compiled with NVRTC (evaluator.cu).
struct NvrtcEntryKernels;
std::shared_ptr<NvrtcEntryKernels> compile_entry_source(const std::string& source, bool complex_values);
void nvrtc_launch_eval(const NvrtcEntryKernels& k, const void* items, int count, int max_m, int max_n,
                       const char* meta, const PointTable& points, const double* params, cudaStream_t stream);
void nvrtc_launch_sketch(const NvrtcEntryKernels& k, const void* items, int count, int max_d, int max_cols,
                         const PointTable& points, const double* params, cudaStream_t stream);

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
