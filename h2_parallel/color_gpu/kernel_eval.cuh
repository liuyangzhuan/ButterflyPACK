#pragma once
// Device evaluation of the registered kernel (KernelSpec, device_kernels.hpp),
// shared by the element and sketch kernels.

#include "device_kernels.hpp"

namespace fmm {
namespace gpu {

template<typename T>
__device__ __forceinline__ T kernel_value(const KernelSpec& spec, const double* x, int64_t x_id, const double* y,
                                          int64_t y_id);

// kind 1: inverse distance
template<>
__device__ __forceinline__ double kernel_value<double>(const KernelSpec& spec, const double* x, int64_t x_id,
                                                       const double* y, int64_t y_id) {
    if (x_id == y_id) return spec.p[1];
    const double dx = x[0] - y[0];
    const double dy = x[1] - y[1];
    const double dz = x[2] - y[2];
    const double r = sqrt(dx * dx + dy * dy + dz * dz);
    return spec.p[0] / r;
}

// kind 2: symmetric Helmholtz, a e^{i k r} / (4 pi r), in the host's order of
// operations (phase k r, amplitude 1 / ((4 pi) r))
template<>
__device__ __forceinline__ dcomplex kernel_value<dcomplex>(const KernelSpec& spec, const double* x, int64_t x_id,
                                                           const double* y, int64_t y_id) {
    if (x_id == y_id) return dcomplex(spec.p[3], spec.p[4]);
    const double dx = x[0] - y[0];
    const double dy = x[1] - y[1];
    const double dz = x[2] - y[2];
    const double r = sqrt(dx * dx + dy * dy + dz * dz);
    const double phase = spec.p[0] * r;
    const double amplitude = 1.0 / (spec.p[5] * r);
    double sn, cs;
    sincos(phase, &sn, &cs);
    return dcomplex(spec.p[1], spec.p[2]) * dcomplex(amplitude * cs, amplitude * sn);
}

}  // namespace gpu
}  // namespace fmm
