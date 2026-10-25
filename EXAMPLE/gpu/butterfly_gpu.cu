// GPU entry evaluators of Butterfly_Driver.cpp (gpu_evaluators.h): the device
// forms of its two kernels, for a butterfly block (bpack::gpu::set_bf_entry_evaluator).
// A butterfly's evaluator gets a row id i and a column id j, each the CPU
// callback's index minus one in its own mesh's order; these read their points
// from tables by those ids, so they need no coordinates from the library.
#include "gpu_evaluators.h"

#include "bpack_gpu_entry.cuh"

#include <cuda_runtime.h>

#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kPi = 3.1415926535897932384626433832795;

// K(i, j) = cos(2 pi s x_i xi_j)
struct DftEntry {
    using value_type = double;
    const double* x;
    const double* xi;
    double s;

    __device__ double operator()(const double*, int64_t i, const double*, int64_t j) const {
        return cos(2.0 * kPi * s * x[i] * xi[j]);
    }
};

// K(i, j) = 1 / (4 pi |p_i - q_j|)
struct LaplaceBfEntry {
    using value_type = double;
    const double* p;  // 3 per row
    const double* q;  // 3 per column

    __device__ double operator()(const double*, int64_t i, const double*, int64_t j) const {
        const double dx = p[3 * i] - q[3 * j];
        const double dy = p[3 * i + 1] - q[3 * j + 1];
        const double dz = p[3 * i + 2] - q[3 * j + 2];
        return 1.0 / (4.0 * kPi * sqrt(dx * dx + dy * dy + dz * dz));
    }
};

std::vector<void*>& tables() {
    static std::vector<void*> t;
    return t;
}

// The tables must be on the GPU the library uses: with several visible to a
// rank, the node-local rank's, as ButterflyPACK picks it (gpu_runtime.hpp).
void use_library_device() {
    int count = 0;
    if (cudaGetDeviceCount(&count) != cudaSuccess || count <= 1) return;
    for (const char* name : {"SLURM_LOCALID", "PALS_LOCAL_RANKID", "OMPI_COMM_WORLD_LOCAL_RANK", "MPI_LOCALRANKID",
                             "MV2_COMM_WORLD_LOCAL_RANK", "PMI_LOCAL_RANK"}) {
        if (const char* v = std::getenv(name)) {
            cudaSetDevice(std::atoi(v) % count);
            return;
        }
    }
}

const double* to_device(const double* host, int64_t n) {
    use_library_device();
    void* d = nullptr;
    const size_t bytes = static_cast<size_t>(n > 0 ? n : 1) * sizeof(double);
    if (cudaMalloc(&d, bytes) != cudaSuccess) throw std::runtime_error("butterfly_gpu: cudaMalloc failed");
    if (n > 0 && cudaMemcpy(d, host, static_cast<size_t>(n) * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) {
        throw std::runtime_error("butterfly_gpu: cudaMemcpy failed");
    }
    tables().push_back(d);
    return static_cast<const double*>(d);
}

}  // namespace

extern "C" void butterfly_gpu_register_dft(void** bf, const double* x, int64_t m, const double* xi, int64_t n,
                                           double s) {
    bpack::gpu::set_bf_entry_evaluator(bf, DftEntry{to_device(x, m), to_device(xi, n), s}, 0);
}

extern "C" void butterfly_gpu_register_laplace(void** bf, const double* rows_xyz, int64_t m, const double* cols_xyz,
                                               int64_t n) {
    bpack::gpu::set_bf_entry_evaluator(bf, LaplaceBfEntry{to_device(rows_xyz, 3 * m), to_device(cols_xyz, 3 * n)}, 0);
}

extern "C" void butterfly_gpu_release(void) {
    for (void* d : tables()) cudaFree(d);
    tables().clear();
}
