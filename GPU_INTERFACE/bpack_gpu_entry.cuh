#pragma once
// ButterflyPACK GPU backends: an entry evaluator from C++ (doc/gpu_kernels.md).
// Include this header in a .cu file compiled by nvcc, define the evaluator
// and register it:
//
//   struct MyEntry {
//       using value_type = double;          // or bpack::gpu::dcomplex
//       double c;                           // parameters, device pointers, ...
//       __device__ double operator()(const double* x, int64_t i, const double* y, int64_t j) const {
//           ...                             // K(i, j): 0-based global ids, coordinates x, y
//       }
//   };
//   bpack::gpu::set_entry_evaluator(&bmat, MyEntry{c}, BPACK_GPU_SYMMETRIC | BPACK_GPU_COORDINATES);
//
// The call instantiates the library's entry and sketch kernels for MyEntry
// here, in the application, and hands their launches to the library, which
// copies the evaluator and passes it to them by value (so it must be
// trivially copyable, and at most 4000 bytes).  The code compiled here
// depends on the library's version: rebuild it with the library (a
// mismatch is refused at registration).

#include "bpack_gpu.h"
#include "bpack_gpu_kernels.cuh"

#include <type_traits>

namespace bpack {
namespace gpu {

using dcomplex = fmm::gpu::dcomplex;

namespace detail {

template<typename T, typename Entry>
void launch_eval(const void* items, int count, int max_m, int max_n, const char* meta, const void* points,
                 const void* entry, void* stream) {
    if (count <= 0 || max_m <= 0 || max_n <= 0) return;
    fmm::gpu::eval_kernel<T, Entry>
        <<<fmm::gpu::eval_grid(count, max_m, max_n), fmm::gpu::eval_block(), 0, static_cast<cudaStream_t>(stream)>>>(
            static_cast<const fmm::gpu::EvalItemT<T>*>(items), meta, *static_cast<const Entry*>(entry),
            *static_cast<const fmm::gpu::PointTable*>(points));
}

template<typename T, typename Entry>
void launch_sketch(const void* items, int count, int max_d, int max_cols, const void* points, const void* entry,
                   void* stream) {
    if (count <= 0 || max_d <= 0 || max_cols <= 0) return;
    const fmm::gpu::SketchGeometry g = fmm::gpu::sketch_geometry(count, max_d, max_cols, sizeof(T));
    auto kernel = fmm::gpu::ordered_sketch_kernel<T, true, Entry>;
    cudaFuncSetAttribute(kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, static_cast<int>(g.shared));
    kernel<<<g.grid, g.block, g.shared, static_cast<cudaStream_t>(stream)>>>(
        static_cast<const fmm::gpu::OrderedSketchItemT<T>*>(items), *static_cast<const Entry*>(entry),
        *static_cast<const fmm::gpu::PointTable*>(points), g.dest_block);
}

}  // namespace detail

// Register `entry` as the evaluator of the matrix `bmat` (the F2Cptr of
// c_bpack_construct_init): value_type double for the d_ library, dcomplex
// for the z_ one.  `flags`: BPACK_GPU_SYMMETRIC (required by H2),
// BPACK_GPU_COORDINATES (the evaluator reads x and y).
template<typename Entry>
void set_entry_evaluator(void** bmat, const Entry& entry, int flags) {
    using T = typename Entry::value_type;
    static_assert(std::is_same<T, double>::value || std::is_same<T, dcomplex>::value,
                  "an entry evaluator's value_type is double or bpack::gpu::dcomplex");
    static_assert(std::is_trivially_copyable<Entry>::value,
                  "an entry evaluator is copied to the device by value: it must be trivially copyable");
    static_assert(sizeof(Entry) <= 4000, "an entry evaluator is a kernel argument: at most 4000 bytes");
    bpack_gpu_entry_launchers launchers;
    launchers.version = BPACK_GPU_INTERFACE_VERSION;
    launchers.scalar_bytes = static_cast<int>(sizeof(T));
    launchers.flags = flags;
    launchers.entry = &entry;
    launchers.entry_bytes = sizeof(Entry);
    launchers.eval = &detail::launch_eval<T, Entry>;
    launchers.sketch = &detail::launch_sketch<T, Entry>;
    if constexpr (std::is_same<T, double>::value) {
        d_c_bpack_set_gpu_entry_launchers(bmat, &launchers);
    } else {
        z_c_bpack_set_gpu_entry_launchers(bmat, &launchers);
    }
}

}  // namespace gpu
}  // namespace bpack
