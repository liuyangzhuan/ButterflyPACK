#pragma once
// Device evaluation of the registered kernel (KernelSpec, device_kernels.hpp),
// shared by the element and sketch kernels.  The kind is a template
// parameter: a launch picks the instantiation of its spec's kind (with
// kernel_kind_of), so a heavy kind does not weigh on the others' kernels.

#include "device_kernels.hpp"
#include "emsurf_kernel.cuh"

#include <stdexcept>
#include <string>

namespace fmm {
namespace gpu {

template<typename T, int Kind>
__device__ __forceinline__ T kernel_value(const KernelSpec& spec, const double* x, int64_t x_id, const double* y,
                                          int64_t y_id);

// kind 1: inverse distance
template<>
__device__ __forceinline__ double kernel_value<double, 1>(const KernelSpec& spec, const double* x, int64_t x_id,
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
__device__ __forceinline__ dcomplex kernel_value<dcomplex, 2>(const KernelSpec& spec, const double* x, int64_t x_id,
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

// kind 3: EFIE entry of the RWG edges x_id, y_id (emsurf_kernel.cuh)
template<>
__device__ __forceinline__ dcomplex kernel_value<dcomplex, 3>(const KernelSpec& spec, const double*, int64_t x_id,
                                                              const double*, int64_t y_id) {
    return emsurf::efie_entry(spec, x_id, y_id);
}

// kind 4 with real data: Gaussian-process squared exponential kernel (a constant times an
// axis-aligned squared exponential) of the point coordinates: with t_d = p[2+d] (x_d - y_d)^2,
// q = t_0 + t_1 + t_2 and e = p[0] exp(-q/2), the mode p[5] selects
//   0: e, plus p[1] when the global ids match (the diagonal: noise and jitter)
//   1: e (derivative in the log of the constant p[0])
//   2, 3, 4: e t_d / 2, d = mode - 2 (derivative in the log metric of dimension d, the metric
//      being 1 / p[2+d] in the units of the coordinates)
//   5: e q / 2 (derivative in the log of an isotropic metric)
template<>
__device__ __forceinline__ double kernel_value<double, 4>(const KernelSpec& spec, const double* x, int64_t x_id,
                                                          const double* y, int64_t y_id) {
    double t[3];
    double q = 0.0;
#pragma unroll
    for (int d = 0; d < 3; ++d) {
        const double dd = x[d] - y[d];
        t[d] = spec.p[2 + d] * dd * dd;
        q += t[d];
    }
    const double e = spec.p[0] * exp(-0.5 * q);
    const int mode = static_cast<int>(spec.p[5]);
    switch (mode) {
        case 0: return x_id == y_id ? e + spec.p[1] : e;
        case 1: return e;
        case 2: case 3: case 4: return 0.5 * e * t[mode - 2];
        default: return 0.5 * e * q;
    }
}

// kind 4 with complex data: kind 2 scaled by a coefficient of the column point, c = treal
// indexed by y_id: c[y] (p[1], p[2]) e^{i p[0] r} / (p[5] r), and c[y] (p[3],
// p[4]) + (p[6], p[7]) on the diagonal (VIE3D's S2S kernel with
// scaleGreen = 0, assemble_fromD1D2Tau_s2s_with_coef)
template<>
__device__ __forceinline__ dcomplex kernel_value<dcomplex, 4>(const KernelSpec& spec, const double* x, int64_t x_id,
                                                              const double* y, int64_t y_id) {
    const double c = spec.treal[y_id];
    if (x_id == y_id) return c * dcomplex(spec.p[3], spec.p[4]) + dcomplex(spec.p[6], spec.p[7]);
    return c * kernel_value<dcomplex, 2>(spec, x, x_id, y, y_id);
}

// kind 5: CFIE entry of the RWG edges x_id, y_id (emsurf_kernel.cuh)
template<>
__device__ __forceinline__ dcomplex kernel_value<dcomplex, 5>(const KernelSpec& spec, const double*, int64_t x_id,
                                                              const double*, int64_t y_id) {
    return emsurf::cfie_entry(spec, x_id, y_id);
}

// The kind of a launch's spec, checked against the element type.
template<typename T>
inline int kernel_kind_of(const KernelSpec& spec, const char* what) {
    const bool ok = is_complex_scalar<T> ? (spec.kind >= 2 && spec.kind <= 5) : (spec.kind == 1 || spec.kind == 4);
    if (!ok) {
        throw std::runtime_error(std::string(what) + ": device kernel kind " + std::to_string(spec.kind) +
                                 " does not match the data type");
    }
    return spec.kind;
}

}  // namespace gpu
}  // namespace fmm
