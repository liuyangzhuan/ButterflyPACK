#pragma once
// Scalar types of the H2 Color GPU backend: the device kernels are templates
// on the element type, double or dcomplex (complex double laid out as
// std::complex<double> and magmaDoubleComplex).  The factorization is
// symmetric (not Hermitian), so no conjugates appear outside the pivoted QR.

#ifdef H2_HAVE_GPU

#include <cuda_runtime.h>
#include "../../GPU_INTERFACE/bpack_gpu_kernels.cuh"

#include <cmath>
#include <complex>
#include <type_traits>

namespace fmm {
namespace gpu {

// (dcomplex and its arithmetic: GPU_INTERFACE/bpack_gpu_kernels.cuh)
static_assert(sizeof(dcomplex) == sizeof(std::complex<double>), "dcomplex must match std::complex<double>");
__host__ __device__ inline double real_part(dcomplex a) { return a.re; }
__host__ __device__ inline double real_part(double a) { return a; }
__host__ __device__ inline double imag_part(dcomplex a) { return a.im; }
__host__ __device__ inline double imag_part(double) { return 0.0; }
__host__ __device__ inline double abs2(dcomplex a) { return a.re * a.re + a.im * a.im; }
__host__ __device__ inline double abs2(double a) { return a * a; }
__host__ __device__ inline double magnitude(dcomplex a) { return hypot(a.re, a.im); }
__host__ __device__ inline double magnitude(double a) { return fabs(a); }
__host__ __device__ inline bool is_finite(double a) {
#ifdef __CUDA_ARCH__
    return isfinite(a);
#else
    return std::isfinite(a);
#endif
}
__host__ __device__ inline bool is_finite(dcomplex a) { return is_finite(a.re) && is_finite(a.im); }

// re + i im as T (im is dropped for double)
template<typename T> __host__ __device__ inline T make_scalar(double re, double im);
template<> __host__ __device__ inline double make_scalar<double>(double re, double) { return re; }
template<> __host__ __device__ inline dcomplex make_scalar<dcomplex>(double re, double im) { return {re, im}; }

// Factorization data types the GPU backend runs; others (single
// precision) stay on the host.
template<typename DataType>
constexpr bool gpu_data_type = std::is_same_v<DataType, double> || std::is_same_v<DataType, std::complex<double>>;

// The device element type of a factorization's DataType (void: not run on
// the GPU).
template<typename DataType> struct DeviceScalar { using type = void; };
template<> struct DeviceScalar<double> { using type = double; };
template<> struct DeviceScalar<std::complex<double>> { using type = dcomplex; };

template<typename T> constexpr bool is_complex_scalar = false;
template<> constexpr bool is_complex_scalar<dcomplex> = true;

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
