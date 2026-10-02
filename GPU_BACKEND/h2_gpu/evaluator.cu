// The application's evaluator of the matrix's entries (evaluator.hpp): the
// device helpers of the block and entry-list forms, and the entry evaluators
// given as source text, compiled with NVRTC.
#include "evaluator_device.hpp"

#include "bpack_gpu_kernels_source.h"  // generated: the text of GPU_INTERFACE/bpack_gpu_kernels.cuh

#include <cuda.h>
#include <nvrtc.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace fmm {
namespace gpu {

namespace {

void check_launch(const char* what) {
    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("CUDA launch failed in ") + what + ": " + cudaGetErrorString(status));
    }
}

__global__ void gather_block_points_kernel(const BlockGatherItem* items, const char* meta, PointTable points,
                                           int64_t* ids, double* coords) {
    const BlockGatherItem it = items[blockIdx.x];
    for (int i = blockIdx.y * blockDim.x + threadIdx.x; i < it.m + it.n; i += gridDim.y * blockDim.x) {
        const bool row = i < it.m;
        const int k = row ? i : i - it.m;
        const int slot = index_at(meta, row ? it.rows : it.cols, k);
        const int64_t at = (row ? it.row_at : it.col_at) + k;
        ids[at] = points.ids[slot];
        if (coords != nullptr) {
            for (int d = 0; d < points.dim; ++d) {
                coords[at * points.dim + d] = points.xyz[static_cast<int64_t>(slot) * points.dim + d];
            }
        }
    }
}

__global__ void flatten_entries_kernel(const FlatBlockItem* items, int64_t* rows, int64_t* cols) {
    const FlatBlockItem it = items[blockIdx.x];
    const int64_t total = static_cast<int64_t>(it.m) * it.n;
    for (int64_t e = static_cast<int64_t>(blockIdx.y) * blockDim.x + threadIdx.x; e < total;
         e += static_cast<int64_t>(gridDim.y) * blockDim.x) {
        const int r = static_cast<int>(e % it.m), c = static_cast<int>(e / it.m);
        rows[it.at + e] = it.row_ids[r];
        cols[it.at + e] = it.col_ids[c];
    }
}

template<typename T>
__global__ void scatter_entries_kernel(const FlatBlockItem* items, const T* values) {
    const FlatBlockItem it = items[blockIdx.x];
    T* out = static_cast<T*>(it.out);
    const int64_t total = static_cast<int64_t>(it.m) * it.n;
    for (int64_t e = static_cast<int64_t>(blockIdx.y) * blockDim.x + threadIdx.x; e < total;
         e += static_cast<int64_t>(gridDim.y) * blockDim.x) {
        const int r = static_cast<int>(e % it.m), c = static_cast<int>(e / it.m);
        out[r + static_cast<int64_t>(c) * it.ld] = values[it.at + e];
    }
}

constexpr int kHelperThreads = 256;

unsigned helper_blocks(int64_t work) {
    const int64_t b = (work + kHelperThreads - 1) / kHelperThreads;
    return static_cast<unsigned>(std::max<int64_t>(1, std::min<int64_t>(b, 1024)));
}

}  // namespace

void launch_gather_block_points(const BlockGatherItem* items, int count, int max_mn, const char* meta,
                                PointTable points, int64_t* ids, double* coords, cudaStream_t stream) {
    if (count <= 0) return;
    gather_block_points_kernel<<<dim3(static_cast<unsigned>(count), helper_blocks(2 * static_cast<int64_t>(max_mn))),
                                 kHelperThreads, 0, stream>>>(items, meta, points, ids, coords);
    check_launch("gather_block_points_kernel");
}

void launch_flatten_entries(const FlatBlockItem* items, int count, int max_m, int max_n, int64_t* rows, int64_t* cols,
                            cudaStream_t stream) {
    if (count <= 0) return;
    flatten_entries_kernel<<<dim3(static_cast<unsigned>(count), helper_blocks(static_cast<int64_t>(max_m) * max_n)),
                             kHelperThreads, 0, stream>>>(items, rows, cols);
    check_launch("flatten_entries_kernel");
}

template<typename T>
void launch_scatter_entries(const FlatBlockItem* items, int count, int max_m, int max_n, const T* values,
                            cudaStream_t stream) {
    if (count <= 0) return;
    scatter_entries_kernel<T><<<dim3(static_cast<unsigned>(count), helper_blocks(static_cast<int64_t>(max_m) * max_n)),
                                kHelperThreads, 0, stream>>>(items, values);
    check_launch("scatter_entries_kernel");
}
template void launch_scatter_entries<double>(const FlatBlockItem*, int, int, int, const double*, cudaStream_t);
template void launch_scatter_entries<dcomplex>(const FlatBlockItem*, int, int, int, const dcomplex*, cudaStream_t);

// ---------------------------------------------------------------------------
// Entry evaluators from source text: the library's kernels (the text of
// bpack_gpu_kernels.cuh) with the application's bpack_entry, compiled with
// NVRTC for the current device, loaded with the driver API.

struct NvrtcEntryKernels {
    CUmodule module = nullptr;
    CUfunction eval = nullptr;
    CUfunction sketch = nullptr;
    bool complex_values = false;
    ~NvrtcEntryKernels() {
        if (module != nullptr) cuModuleUnload(module);
    }
};

namespace {

// The evaluator the compiled kernels take by value (BpackNvrtcEntry below).
struct NvrtcEntryArgs {
    const double* params;
    int dim;
};

void check_nvrtc(nvrtcResult r, const char* what) {
    if (r != NVRTC_SUCCESS) throw std::runtime_error(std::string("NVRTC: ") + what + ": " + nvrtcGetErrorString(r));
}
void check_runtime(cudaError_t r, const char* what) {
    if (r != cudaSuccess) throw std::runtime_error(std::string("CUDA: ") + what + ": " + cudaGetErrorString(r));
}
void check_driver(CUresult r, const char* what) {
    if (r != CUDA_SUCCESS) {
        const char* name = nullptr;
        cuGetErrorName(r, &name);
        throw std::runtime_error(std::string("CUDA driver: ") + what + ": " + (name != nullptr ? name : "error"));
    }
}

}  // namespace

std::shared_ptr<NvrtcEntryKernels> compile_entry_source(const std::string& source, bool complex_values) {
    check_runtime(cudaFree(nullptr), "CUDA context");  // (the primary context, current for the driver API)
    const std::string type = complex_values ? "fmm::gpu::dcomplex" : "double";
    std::string text = kBpackGpuKernelsSource;
    text += "\nusing bpack_dcomplex = fmm::gpu::dcomplex;\n#line 1 \"bpack_entry\"\n";
    text += source;
    text += "\nstruct BpackNvrtcEntry {\n"
            "    const double* params;\n"
            "    int dim;\n"
            "    using value_type = " + type + ";\n"
            "    __device__ value_type operator()(const double* x, int64_t i, const double* y, int64_t j) const {\n"
            "        return bpack_entry(x, i, y, j, params, dim);\n"
            "    }\n"
            "};\n";
    const std::string eval_name = "fmm::gpu::eval_kernel<" + type + ", BpackNvrtcEntry>";
    const std::string sketch_name = "fmm::gpu::ordered_sketch_kernel<" + type + ", true, BpackNvrtcEntry>";

    nvrtcProgram program = nullptr;
    check_nvrtc(nvrtcCreateProgram(&program, text.c_str(), "bpack_gpu_entry.cu", 0, nullptr, nullptr),
                "create program");
    auto destroy = [&] { nvrtcDestroyProgram(&program); };
    try {
        check_nvrtc(nvrtcAddNameExpression(program, eval_name.c_str()), "name expression");
        check_nvrtc(nvrtcAddNameExpression(program, sketch_name.c_str()), "name expression");
        int device = 0;
        check_runtime(cudaGetDevice(&device), "cudaGetDevice");
        cudaDeviceProp prop;
        check_runtime(cudaGetDeviceProperties(&prop, device), "cudaGetDeviceProperties");
        const std::string arch = "--gpu-architecture=sm_" + std::to_string(prop.major) + std::to_string(prop.minor);
        const char* options[] = {arch.c_str(), "--std=c++17"};
        const nvrtcResult compiled = nvrtcCompileProgram(program, 2, options);
        if (compiled != NVRTC_SUCCESS) {
            size_t log_size = 0;
            nvrtcGetProgramLogSize(program, &log_size);
            std::string log(log_size, '\0');
            if (log_size > 0) nvrtcGetProgramLog(program, &log[0]);
            throw std::runtime_error("GPU entry evaluator: the source did not compile (NVRTC):\n" + log);
        }
        size_t size = 0;
        check_nvrtc(nvrtcGetCUBINSize(program, &size), "CUBIN size");
        std::vector<char> cubin(size);
        check_nvrtc(nvrtcGetCUBIN(program, cubin.data()), "CUBIN");
        auto k = std::make_shared<NvrtcEntryKernels>();
        k->complex_values = complex_values;
        check_driver(cuModuleLoadData(&k->module, cubin.data()), "cuModuleLoadData");
        const char* lowered = nullptr;
        check_nvrtc(nvrtcGetLoweredName(program, eval_name.c_str(), &lowered), "lowered name");
        check_driver(cuModuleGetFunction(&k->eval, k->module, lowered), "cuModuleGetFunction");
        check_nvrtc(nvrtcGetLoweredName(program, sketch_name.c_str(), &lowered), "lowered name");
        check_driver(cuModuleGetFunction(&k->sketch, k->module, lowered), "cuModuleGetFunction");
        destroy();
        return k;
    } catch (...) {
        destroy();
        throw;
    }
}

void nvrtc_launch_eval(const NvrtcEntryKernels& k, const void* items, int count, int max_m, int max_n,
                       const char* meta, const PointTable& points, const double* params, cudaStream_t stream) {
    if (count <= 0 || max_m <= 0 || max_n <= 0) return;
    NvrtcEntryArgs entry{params, points.dim};
    PointTable table = points;
    const void* items_arg = items;
    const char* meta_arg = meta;
    void* args[] = {&items_arg, &meta_arg, &entry, &table};
    const dim3 g = eval_grid(count, max_m, max_n), b = eval_block();
    check_driver(cuLaunchKernel(k.eval, g.x, g.y, g.z, b.x, b.y, b.z, 0, reinterpret_cast<CUstream>(stream), args,
                                nullptr),
                 "eval_kernel");
}

void nvrtc_launch_sketch(const NvrtcEntryKernels& k, const void* items, int count, int max_d, int max_cols,
                         const PointTable& points, const double* params, cudaStream_t stream) {
    if (count <= 0 || max_d <= 0 || max_cols <= 0) return;
    const SketchGeometry g = sketch_geometry(count, max_d, max_cols, k.complex_values ? sizeof(dcomplex) : sizeof(double));
    check_driver(cuFuncSetAttribute(k.sketch, CU_FUNC_ATTRIBUTE_MAX_DYNAMIC_SHARED_SIZE_BYTES, static_cast<int>(g.shared)),
                 "cuFuncSetAttribute");
    NvrtcEntryArgs entry{params, points.dim};
    PointTable table = points;
    const void* items_arg = items;
    int dest_block = g.dest_block;
    void* args[] = {&items_arg, &entry, &table, &dest_block};
    check_driver(cuLaunchKernel(k.sketch, g.grid.x, g.grid.y, g.grid.z, g.block.x, g.block.y, g.block.z,
                                static_cast<unsigned>(g.shared), reinterpret_cast<CUstream>(stream), args, nullptr),
                 "ordered_sketch_kernel");
}

}  // namespace gpu
}  // namespace fmm
