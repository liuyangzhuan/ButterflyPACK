// GPU entry evaluator of RBF3D_H2Test_Driver.cpp (gpu_evaluators.h):
// the device form of its entry callback.
#include "gpu_evaluators.h"

#include "bpack_gpu_entry.cuh"

namespace {

// K(i, j) = scale / |x - y| with length_scale <= 0, and the self-cell integral
// `diagonal` when i == j.  With length_scale > 0 the kernel is instead the
// squared exponential exp(-|x - y|^2 / (2 l^2)), which decays over a length the
// caller chooses rather than over the whole domain; `scale` is then its
// amplitude.  The short-ranged form is what reproduces the factorization's
// collect_lazy_far_sources failure, which the long-ranged Laplace kernel never
// reaches.
struct RBFEntry {
    using value_type = double;
    double scale;
    double diagonal;
    double length_scale;

    __device__ double operator()(const double* x, int64_t i, const double* y, int64_t j) const {
        if (i == j) return diagonal;
        const double dx = x[0] - y[0];
        const double dy = x[1] - y[1];
        const double dz = x[2] - y[2];
        const double r2 = dx * dx + dy * dy + dz * dz;
        if (length_scale > 0.0) {
            return scale * exp(-r2 / (2.0 * length_scale * length_scale));
        }
        return scale / sqrt(r2);
    }
};

}  // namespace

extern "C" void rbf3d_gpu_register(void** bmat, double scale, double diagonal,
                                       double length_scale) {
    bpack::gpu::set_entry_evaluator(bmat, RBFEntry{scale, diagonal, length_scale},
                                    BPACK_GPU_SYMMETRIC | BPACK_GPU_COORDINATES);
}
