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

// The driver's nonsymmetric forms (--skew, --column_weight): times
// 1 + skew (d . (x - y)) / r off the diagonal, and every entry times
// w(y) = 1 + weight sin(2 pi (y0 + 2 y1 + 3 y2))
struct LaplaceGeneralEntry {
    using value_type = double;
    double scale;
    double diagonal;
    double skew;
    double d0, d1, d2;
    double weight;

    __device__ double operator()(const double* x, int64_t i, const double* y, int64_t j) const {
        const double w =
            weight == 0.0 ? 1.0 : 1.0 + weight * sin(2.0 * 3.14159265358979323846 * (y[0] + 2.0 * y[1] + 3.0 * y[2]));
        if (i == j) return diagonal * w;
        const double dx = x[0] - y[0];
        const double dy = x[1] - y[1];
        const double dz = x[2] - y[2];
        const double r = sqrt(dx * dx + dy * dy + dz * dz);
        double value = scale / r;
        if (skew != 0.0) value *= 1.0 + skew * (d0 * dx + d1 * dy + d2 * dz) / r;
        return value * w;
    }
};

}  // namespace

extern "C" void laplace3d_gpu_register(void** bmat, double scale, double diagonal) {
    bpack::gpu::set_entry_evaluator(bmat, LaplaceEntry{scale, diagonal},
                                    BPACK_GPU_SYMMETRIC | BPACK_GPU_COORDINATES);
}

extern "C" void laplace3d_gpu_register_general(void** bmat, double scale, double diagonal, double skew,
                                               const double* direction, double weight) {
    if (skew == 0.0 && weight == 0.0) {
        laplace3d_gpu_register(bmat, scale, diagonal);
        return;
    }
    bpack::gpu::set_entry_evaluator(
        bmat, LaplaceGeneralEntry{scale, diagonal, skew, direction[0], direction[1], direction[2], weight},
        BPACK_GPU_COORDINATES);
}
