#pragma once
// The device form (KernelSpec) of the kernel an application registered in
// an H2Kernel's gpu_spec: kind, parameters, and the tables of a kind that
// has them (kind 3: the mesh of emsurf_kernel.cuh).  The tables go to device
// memory of their own, outside the heap (whose resets free everything), once
// per registration.

#ifdef H2_HAVE_GPU

#include "device_kernels.hpp"
#include "gpu_runtime.hpp"

#include <cstdint>
#include <cstdlib>
#include <type_traits>

namespace fmm {
namespace gpu {

struct DeviceKernelTables {
    const void* owner = nullptr;  // the gpu_spec uploaded
    uint64_t version = 0;
    double* real = nullptr;
    int* ints = nullptr;
    size_t real_count = 0, int_count = 0;

    void release() {
        if (real != nullptr) cudaFree(real);
        if (ints != nullptr) cudaFree(ints);
        real = nullptr;
        ints = nullptr;
        owner = nullptr;
        version = 0;
        real_count = int_count = 0;
    }
};

inline DeviceKernelTables& device_kernel_tables() {
    static DeviceKernelTables tables;
    return tables;
}

// Whether gpu_spec holds a device kernel for DataType (kind 1 or 4 real; kind 2,
// or kind 3, 4 or 5 with its tables, complex); `why` says what is missing otherwise.
template<typename DataType, typename GpuSpecT>
bool device_kernel_registered(const GpuSpecT& g, const char** why) {
    if constexpr (std::is_same_v<DataType, double>) {
        if (g.kind == 1 || g.kind == 4) return true;
        *why = "no real device kernel registered (c_bpack_h2_set_gpu_kernel, kind 1 or 4)";
    } else {
        if (g.kind == 2 || ((g.kind == 3 || g.kind == 5) && !g.table_real.empty() && !g.table_int.empty()) ||
            (g.kind == 4 && !g.table_real.empty())) {
            return true;
        }
        *why = g.kind >= 3 && g.kind <= 5 ? "device kernel without its tables (c_bpack_set_gpu_kernel_tables)"
                                          : "no complex device kernel registered (c_bpack_set_gpu_kernel, kind 2 to 5)";
    }
    return false;
}

template<typename GpuSpecT>
KernelSpec device_kernel_spec(const GpuSpecT& g) {
    KernelSpec spec;
    spec.kind = g.kind;
    for (int i = 0; i < kKernelParams; ++i) spec.p[i] = g.params[i];
    if (g.table_real.empty() && g.table_int.empty()) return spec;
    DeviceKernelTables& t = device_kernel_tables();
    if (t.owner != &g || t.version != g.table_version || t.real_count != g.table_real.size() ||
        t.int_count != g.table_int.size()) {
        Context::instance().activate();
        t.release();
        const size_t rb = g.table_real.size() * sizeof(double), ib = g.table_int.size() * sizeof(int);
        if (rb > 0) {
            check_cuda(cudaMalloc(reinterpret_cast<void**>(&t.real), rb), "cudaMalloc (kernel tables)");
            check_cuda(cudaMemcpy(t.real, g.table_real.data(), rb, cudaMemcpyHostToDevice), "kernel tables");
        }
        if (ib > 0) {
            check_cuda(cudaMalloc(reinterpret_cast<void**>(&t.ints), ib), "cudaMalloc (kernel tables)");
            check_cuda(cudaMemcpy(t.ints, g.table_int.data(), ib, cudaMemcpyHostToDevice), "kernel tables");
        }
        t.owner = &g;
        t.version = g.table_version;
        t.real_count = g.table_real.size();
        t.int_count = g.table_int.size();
    }
    spec.treal = t.real;
    spec.tint = t.ints;
    return spec;
}

// The device sketch covers a level's ID targets (factorization and compression).  Static training rows
// (H2_ID_radius > 2, H2_ID_proxy 1: points anywhere in the tree) join the
// level's point table: a kernel of coordinates needs their global
// coordinates (kept for H2_ID_proxy 1), kind 3 only their ids.  Adaptive rows
// (H2_ID_proxy 2) are not streamed.
template<typename Tree>
bool device_sketch_supported(const Tree* tree, int kernel_kind) {
    if (tree->id_proxy_mode == 2) return false;
    if (tree->id_neighborhood_radius <= 2 && tree->id_proxy_mode != 1) return true;
    const size_t coordinates = static_cast<size_t>(tree->num_points) * static_cast<size_t>(tree->dimension);
    return kernel_kind == 3 || tree->id_source_point_coords.size() == coordinates;
}

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
