// GPU entry evaluator of Laplace3D_H2Test_Driver.cpp (gpu_evaluators.h):
// the device form of its entry callback.
#include "gpu_evaluators.h"

#include "bpack_gpu_entry.cuh"

namespace {

// K(i, j) = scale / |x - y|, and the self-cell integral `diagonal` when i == j
struct LaplaceEntry {
    using value_type = double;
    double scale;
    double diagonal;

    __device__ double operator()(const double* x, int64_t i, const double* y, int64_t j) const {
        if (i == j) return diagonal;
        const double dx = x[0] - y[0];
        const double dy = x[1] - y[1];
        const double dz = x[2] - y[2];
        return scale / sqrt(dx * dx + dy * dy + dz * dz);
    }
};

}  // namespace

extern "C" void laplace3d_gpu_register(void** bmat, double scale, double diagonal) {
    bpack::gpu::set_entry_evaluator(bmat, LaplaceEntry{scale, diagonal},
                                    BPACK_GPU_SYMMETRIC | BPACK_GPU_COORDINATES);
}
