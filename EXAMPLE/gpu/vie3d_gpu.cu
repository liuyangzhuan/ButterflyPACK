// GPU entry evaluators of VIE3D_H2Test_Driver.cpp (gpu_evaluators.h): the
// device forms of its scatterer-scatterer kernel, as
// assemble_fromD1D2Tau_s2s_with_coef.
#include "gpu_evaluators.h"

#include "bpack_gpu_entry.cuh"

#include <cstdio>
#include <cstdlib>

namespace {

using bpack::gpu::dcomplex;

// a e^{i k r} / (four_pi r), and `self` when i == j (scaleGreen 1), in the
// host's order of operations (phase k r, amplitude 1 / (four_pi r))
struct HelmholtzEntry {
    using value_type = dcomplex;
    double k;
    dcomplex a;
    dcomplex self;
    double four_pi;

    __device__ dcomplex operator()(const double* x, int64_t i, const double* y, int64_t j) const {
        if (i == j) return self;
        const double dx = x[0] - y[0];
        const double dy = x[1] - y[1];
        const double dz = x[2] - y[2];
        const double r = sqrt(dx * dx + dy * dy + dz * dz);
        const double amplitude = 1.0 / (four_pi * r);
        double sn, cs;
        sincos(k * r, &sn, &cs);
        return a * dcomplex(amplitude * cs, amplitude * sn);
    }
};

// coef[j] times HelmholtzEntry, plus `diagonal` when i == j (scaleGreen 0):
// a table on the device, read at the column point's global id
struct CoefficientEntry {
    using value_type = dcomplex;
    HelmholtzEntry helmholtz;
    dcomplex diagonal;
    const double* coef;  // device, one value per point

    __device__ dcomplex operator()(const double* x, int64_t i, const double* y, int64_t j) const {
        const double c = coef[j];
        if (i == j) return c * helmholtz.self + diagonal;
        return c * helmholtz(x, i, y, j);
    }
};

double* coef_device = nullptr;

void check(cudaError_t status, const char* what) {
    if (status != cudaSuccess) {
        std::fprintf(stderr, "vie3d_gpu: %s: %s\n", what, cudaGetErrorString(status));
        std::abort();
    }
}

}  // namespace

extern "C" void vie3d_gpu_register_helmholtz(void** bmat, double k, double a_re, double a_im, double self_re,
                                             double self_im, double four_pi) {
    const HelmholtzEntry entry{k, dcomplex(a_re, a_im), dcomplex(self_re, self_im), four_pi};
    bpack::gpu::set_entry_evaluator(bmat, entry, BPACK_GPU_SYMMETRIC | BPACK_GPU_COORDINATES);
}

extern "C" void vie3d_gpu_register_coefficient(void** bmat, double k, double a_re, double a_im, double self_re,
                                               double self_im, double four_pi, double diag_re, double diag_im,
                                               const double* coef, int64_t n) {
    vie3d_gpu_release();
    check(cudaMalloc(reinterpret_cast<void**>(&coef_device), static_cast<size_t>(n) * sizeof(double)),
          "cudaMalloc");
    check(cudaMemcpy(coef_device, coef, static_cast<size_t>(n) * sizeof(double), cudaMemcpyHostToDevice),
          "cudaMemcpy");
    const CoefficientEntry entry{HelmholtzEntry{k, dcomplex(a_re, a_im), dcomplex(self_re, self_im), four_pi},
                                 dcomplex(diag_re, diag_im), coef_device};
    bpack::gpu::set_entry_evaluator(bmat, entry, BPACK_GPU_COORDINATES);  // not symmetric
}

extern "C" void vie3d_gpu_release(void) {
    if (coef_device != nullptr) cudaFree(coef_device);
    coef_device = nullptr;
}
