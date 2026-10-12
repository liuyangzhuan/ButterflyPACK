#pragma once
// SVD of one dense device matrix (cuSOLVER gesvdp, polar decomposition
// based): the recompression of the HODLR blocks of large rank and the root
// of the merges' TSQR.

#ifdef H2_HAVE_GPU

#include "h2_gpu/device_heap.hpp"
#include "gpu_common/device_scalar.hpp"
#include "h2_gpu/gpu_runtime.hpp"

#include <cusolverDn.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace bpack {
namespace gpu {

// An m x n identity (ones on the diagonal, leading dimension m) from the
// heap: the operand of the warm-ups
template<typename T>
T* heap_identity(int m, int n) {
    std::vector<T> host(static_cast<size_t>(m) * n, T(0.0));
    for (int i = 0; i < std::min(m, n); ++i) host[static_cast<size_t>(i) * m + i] = T(1.0);
    T* d = fmm::gpu::DeviceHeap::instance().alloc<T>(host.size() * sizeof(T));
    fmm::gpu::check_cuda(cudaMemcpy(d, host.data(), host.size() * sizeof(T), cudaMemcpyHostToDevice), "warm-up");
    return d;
}

// The Gram-Cholesky form of the truncated SVDs of the construction (the
// recompression of BACA, the TSQR of the merges) in place of their
// Householder QR factorizations: for a tall X (here U and V), R the upper
// Cholesky factor of X^H X + s I with s = kGramShift tr(X^H X), the SVD of the
// small factors as before, and the new factors X R^-1 (...), without forming
// Q (X stays as it is).  X R^-1 is a contraction when the shift s covers the
// rounding of the Gram matrix and of the Cholesky factorization, so the
// truncation error stays within the threshold, also for X of lower rank than
// its columns (BACA's last batch); singular values below about sqrt(s/tr) of
// the largest are not resolved (1.6e-6 for 32768 x 2300), hence the
// tolerances it is taken for.  Most of its work is GEMM, where the QR's
// panels are memory bound.  BPACK_CHECK=hodlr-exact keeps the QR, which
// matches the CPU at round-off.  (Measured, 3D Laplace 32^3: at tol 1e-4 the
// blocks agree with the CPU's to 1e-7 with the same ranks; at 1e-6 a shift of
// 1e-13 tr gave ranks up to 17 larger and differences of 1e-5.)
constexpr double kGramMinTol = 1e-4;
// s / tr(X^H X) for X of m rows and k columns
inline double gram_shift(int k, double m) { return 10.0 * 1.1102230246251565e-16 * (k + std::sqrt(std::max(m, 0.0))); }
inline bool gram_truncation(double tol) { return tol >= kGramMinTol && !fmm::env::check("hodlr-exact"); }

template<typename T>
class DeviceSvd {
public:
    DeviceSvd() = default;
    DeviceSvd(const DeviceSvd&) = delete;
    DeviceSvd& operator=(const DeviceSvd&) = delete;
    ~DeviceSvd() {
        if (params_ != nullptr) cusolverDnDestroyParams(params_);
        if (handle_ != nullptr) cusolverDnDestroy(handle_);
    }

    // a (m x n, leading dimension m) = u diag(s) v^H with u (m x mn), s (mn,
    // descending), v (n x mn), mn = min(m, n), all on the device (a
    // destroyed); returns the perturbation estimate of gesvdp, after the
    // stream's work
    double svd_device(int m, int n, T* da, double* ds, T* du, T* dv) {
        using fmm::gpu::check_cuda;
        fmm::gpu::Context& ctx = fmm::gpu::Context::instance();
        ctx.activate();
        const cudaStream_t stream = ctx.stream();
        if (handle_ == nullptr) {
            check(cusolverDnCreate(&handle_), "cusolverDnCreate");
            check(cusolverDnCreateParams(&params_), "cusolverDnCreateParams");
        }
        check(cusolverDnSetStream(handle_, stream), "cusolverDnSetStream");
        const cudaDataType type = fmm::gpu::is_complex_scalar<T> ? CUDA_C_64F : CUDA_R_64F;
        fmm::gpu::DeviceHeap& heap = fmm::gpu::DeviceHeap::instance();
        int* dinfo = heap.alloc<int>(sizeof(int));
        size_t dbytes = 0, hbytes = 0;
        check(cusolverDnXgesvdp_bufferSize(handle_, params_, CUSOLVER_EIG_MODE_VECTOR, 1, m, n, type, da, m,
                                           CUDA_R_64F, ds, type, du, m, type, dv, n, type, &dbytes, &hbytes),
              "cusolverDnXgesvdp_bufferSize");
        char* dwork = heap.alloc<char>(std::max<size_t>(dbytes, 1));
        std::vector<char> hwork(std::max<size_t>(hbytes, 1));
        double err_sigma = 0.0;
        check(cusolverDnXgesvdp(handle_, params_, CUSOLVER_EIG_MODE_VECTOR, 1, m, n, type, da, m, CUDA_R_64F, ds, type,
                                du, m, type, dv, n, type, dwork, dbytes, hwork.data(), hbytes, dinfo, &err_sigma),
              "cusolverDnXgesvdp");
        int info = 0;
        check_cuda(cudaMemcpyAsync(&info, dinfo, sizeof(int), cudaMemcpyDeviceToHost, stream), "SVD info");
        check_cuda(cudaStreamSynchronize(stream), "SVD");
        heap.free(dwork);
        heap.free(dinfo);
        if (info != 0) throw std::runtime_error("cusolverDnXgesvdp failed, info " + std::to_string(info));
        return err_sigma;
    }

    // The upper Cholesky factor of a (n x n, leading dimension lda, device) in
    // place (its strictly lower part as it was); returns the info of potrf
    // (0: done), after the stream's work
    int potrf_device(int n, T* da, int lda) {
        using fmm::gpu::check_cuda;
        fmm::gpu::Context& ctx = fmm::gpu::Context::instance();
        ctx.activate();
        const cudaStream_t stream = ctx.stream();
        if (handle_ == nullptr) {
            check(cusolverDnCreate(&handle_), "cusolverDnCreate");
            check(cusolverDnCreateParams(&params_), "cusolverDnCreateParams");
        }
        check(cusolverDnSetStream(handle_, stream), "cusolverDnSetStream");
        const cudaDataType type = fmm::gpu::is_complex_scalar<T> ? CUDA_C_64F : CUDA_R_64F;
        fmm::gpu::DeviceHeap& heap = fmm::gpu::DeviceHeap::instance();
        size_t dbytes = 0, hbytes = 0;
        check(cusolverDnXpotrf_bufferSize(handle_, params_, CUBLAS_FILL_MODE_UPPER, n, type, da, lda, type, &dbytes, &hbytes),
              "cusolverDnXpotrf_bufferSize");
        const size_t wbytes = fmm::gpu::align_up(std::max<size_t>(dbytes, 1));
        char* dwork = heap.alloc<char>(wbytes + sizeof(int));
        int* dinfo = reinterpret_cast<int*>(dwork + wbytes);
        std::vector<char> hwork(std::max<size_t>(hbytes, 1));
        check(cusolverDnXpotrf(handle_, params_, CUBLAS_FILL_MODE_UPPER, n, type, da, lda, type, dwork, dbytes, hwork.data(),
                               hbytes, dinfo),
              "cusolverDnXpotrf");
        int info = 0;
        check_cuda(cudaMemcpyAsync(&info, dinfo, sizeof(int), cudaMemcpyDeviceToHost, stream), "potrf info");
        check_cuda(cudaStreamSynchronize(stream), "potrf");
        heap.free(dwork);
        return info;
    }

    // Before the first SVD (GpuState::warm_up): the cuSOLVER handle and one
    // small SVD, so that no level pays their first use
    void warm_up() {
        constexpr int n = 8;
        fmm::gpu::DeviceHeap& heap = fmm::gpu::DeviceHeap::instance();
        T* a = heap_identity<T>(n, n);
        double* s = heap.alloc<double>(n * sizeof(double));
        T* u = heap.alloc<T>(n * n * sizeof(T));
        T* v = heap.alloc<T>(n * n * sizeof(T));
        svd_device(n, n, a, s, u, v);  // (synchronizes)
        for (void* p : {static_cast<void*>(v), static_cast<void*>(u), static_cast<void*>(s), static_cast<void*>(a)}) {
            heap.free(p);
        }
    }

private:
    static void check(cusolverStatus_t status, const char* what) {
        if (status != CUSOLVER_STATUS_SUCCESS) {
            throw std::runtime_error(std::string(what) + " failed with status " + std::to_string(static_cast<int>(status)));
        }
    }

    cusolverDnHandle_t handle_ = nullptr;
    cusolverDnParams_t params_ = nullptr;
};

}  // namespace gpu
}  // namespace bpack

#endif  // H2_HAVE_GPU
