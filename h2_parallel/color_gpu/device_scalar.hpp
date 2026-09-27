#pragma once
// Scalar types of the H2 Color GPU backend: the device kernels are templates
// on the element type, double or dcomplex (complex double laid out as
// std::complex<double> and magmaDoubleComplex).  The factorization is
// symmetric (not Hermitian), so no conjugates appear outside the pivoted QR.

#ifdef H2_HAVE_GPU

#include <cuda_runtime.h>

#include <cmath>
#include <complex>
#include <type_traits>

namespace fmm {
namespace gpu {

struct __align__(16) dcomplex {
    double re;
    double im;
    dcomplex() = default;
    __host__ __device__ constexpr dcomplex(double r, double i = 0.0) : re(r), im(i) {}
};
static_assert(sizeof(dcomplex) == sizeof(std::complex<double>), "dcomplex must match std::complex<double>");

__host__ __device__ inline dcomplex operator+(dcomplex a, dcomplex b) { return {a.re + b.re, a.im + b.im}; }
__host__ __device__ inline dcomplex operator-(dcomplex a, dcomplex b) { return {a.re - b.re, a.im - b.im}; }
__host__ __device__ inline dcomplex operator-(dcomplex a) { return {-a.re, -a.im}; }
__host__ __device__ inline dcomplex operator*(dcomplex a, dcomplex b) {
    return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
}
__host__ __device__ inline dcomplex operator*(double a, dcomplex b) { return {a * b.re, a * b.im}; }
__host__ __device__ inline dcomplex operator*(dcomplex a, double b) { return {a.re * b, a.im * b}; }
__host__ __device__ inline dcomplex operator/(dcomplex a, double b) { return {a.re / b, a.im / b}; }
// Smith's algorithm (no overflow in |b|^2)
__host__ __device__ inline dcomplex operator/(dcomplex a, dcomplex b) {
    if (fabs(b.re) >= fabs(b.im)) {
        const double r = b.im / b.re, d = b.re + b.im * r;
        return {(a.re + a.im * r) / d, (a.im - a.re * r) / d};
    }
    const double r = b.re / b.im, d = b.re * r + b.im;
    return {(a.re * r + a.im) / d, (a.im * r - a.re) / d};
}
__host__ __device__ inline dcomplex& operator+=(dcomplex& a, dcomplex b) { a.re += b.re; a.im += b.im; return a; }
__host__ __device__ inline dcomplex& operator-=(dcomplex& a, dcomplex b) { a.re -= b.re; a.im -= b.im; return a; }
__host__ __device__ inline dcomplex& operator*=(dcomplex& a, dcomplex b) { a = a * b; return a; }
__host__ __device__ inline dcomplex& operator/=(dcomplex& a, double b) { a.re /= b; a.im /= b; return a; }
__host__ __device__ inline bool operator==(dcomplex a, dcomplex b) { return a.re == b.re && a.im == b.im; }
__host__ __device__ inline bool operator!=(dcomplex a, dcomplex b) { return !(a == b); }

__host__ __device__ inline dcomplex conj(dcomplex a) { return {a.re, -a.im}; }
__host__ __device__ inline double conj(double a) { return a; }
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
