#pragma once
// SVD of one dense device matrix (cuSOLVER gesvdp, polar decomposition
// based, checked against the truncation that follows, else one-sided Jacobi,
// gesvdj): the recompression of the HODLR blocks of large rank and the root
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
// That holds for the factors one at a time, not for U V^T: the singular
// values of R_U R_V^T are those of [U; s_U^1/2 I] [V; s_V^1/2 I]^T, U V^T
// plus blocks of energy s_U |V|_F^2 + s_V |U|_F^2 + k s_U s_V, which is small
// against (tol s_1)^2 unless U and V are ill conditioned and cancel (BACA
// seeded by the nearest neighbours, whose columns of A are nearly parallel:
// STRUMPACK's F12 at tol 1e-4 came out of rank 266 instead of 360, 1e-2
// off).  A batch with a block where it is not goes to the QR
// (HodlrConstruct::recompress_gram).
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
        if (jinfo_ != nullptr) cusolverDnDestroyGesvdjInfo(jinfo_);
        if (params_ != nullptr) cusolverDnDestroyParams(params_);
        if (handle_ != nullptr) cusolverDnDestroy(handle_);
    }

    // a (m x n, leading dimension m) = u diag(s) v^H with u (m x mn), s (mn,
    // descending), v (n x mn), mn = min(m, n), all on the device (a
    // destroyed); returns the perturbation estimate of gesvdp, after the
    // stream's work.  For a truncation at tol (relative to s(0), tol > 0) the
    // SVD is redone by one-sided Jacobi (gesvdj) when gesvdp's estimate
    // exceeds kSvdCheck tol s(0): gesvdp can lose the small singular values
    // of a dense, numerically rank deficient a (measured on HODLR blocks:
    // estimates up to 0.6 s(0), a block of rank about 120 at 1e-3 returned as
    // of full rank)
    double svd_device(int m, int n, T* da, double* ds, T* du, T* dv, double tol = 0.0) {
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
        // (a copy of a for the check's second SVD)
        const size_t abytes = static_cast<size_t>(m) * n * sizeof(T);
        T* acopy = nullptr;
        if (tol > 0.0) {
            acopy = heap.alloc<T>(std::max<size_t>(abytes, 1));
            check_cuda(cudaMemcpyAsync(acopy, da, abytes, cudaMemcpyDeviceToDevice, stream), "SVD copy");
        }
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
        if (acopy != nullptr) {
            double s0 = 0.0;
            check_cuda(cudaMemcpy(&s0, ds, sizeof(double), cudaMemcpyDeviceToHost), "SVD s(0)");
            if (err_sigma > kSvdCheck * tol * s0) jacobi_device(m, n, acopy, ds, du, dv);
            heap.free(acopy);
        }
        return err_sigma;
    }

    // a (m x n, leading dimension m) = u diag(s) v^H by one-sided Jacobi
    // (cuSOLVER gesvdj, economy size, its default tolerance and sweeps), as
    // svd_device
    void jacobi_device(int m, int n, T* da, double* ds, T* du, T* dv) {
        using fmm::gpu::check_cuda;
        fmm::gpu::Context& ctx = fmm::gpu::Context::instance();
        const cudaStream_t stream = ctx.stream();
        fmm::gpu::DeviceHeap& heap = fmm::gpu::DeviceHeap::instance();
        if (jinfo_ == nullptr) check(cusolverDnCreateGesvdjInfo(&jinfo_), "cusolverDnCreateGesvdjInfo");
        int lwork = 0;
        if constexpr (fmm::gpu::is_complex_scalar<T>) {
            auto* a = reinterpret_cast<cuDoubleComplex*>(da);
            auto* u = reinterpret_cast<cuDoubleComplex*>(du);
            auto* v = reinterpret_cast<cuDoubleComplex*>(dv);
            check(cusolverDnZgesvdj_bufferSize(handle_, CUSOLVER_EIG_MODE_VECTOR, 1, m, n, a, m, ds, u, m, v, n, &lwork, jinfo_),
                  "cusolverDnZgesvdj_bufferSize");
            auto* work = heap.alloc<cuDoubleComplex>(static_cast<size_t>(std::max(lwork, 1)) * sizeof(cuDoubleComplex));
            int* dinfo = heap.alloc<int>(sizeof(int));
            check(cusolverDnZgesvdj(handle_, CUSOLVER_EIG_MODE_VECTOR, 1, m, n, a, m, ds, u, m, v, n, work, lwork, dinfo, jinfo_),
                  "cusolverDnZgesvdj");
            jacobi_info(dinfo, m, n, stream);
            heap.free(dinfo);
            heap.free(work);
        } else {
            check(cusolverDnDgesvdj_bufferSize(handle_, CUSOLVER_EIG_MODE_VECTOR, 1, m, n, da, m, ds, du, m, dv, n, &lwork, jinfo_),
                  "cusolverDnDgesvdj_bufferSize");
            double* work = heap.alloc<double>(static_cast<size_t>(std::max(lwork, 1)) * sizeof(double));
            int* dinfo = heap.alloc<int>(sizeof(int));
            check(cusolverDnDgesvdj(handle_, CUSOLVER_EIG_MODE_VECTOR, 1, m, n, da, m, ds, du, m, dv, n, work, lwork, dinfo, jinfo_),
                  "cusolverDnDgesvdj");
            jacobi_info(dinfo, m, n, stream);
            heap.free(dinfo);
            heap.free(work);
        }
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

    // gesvdj's info after the stream's work: negative, an argument; min(m, n) + 1, not converged in its
    // sweeps (accepted: it stops at a tolerance far below the truncation's)
    static void jacobi_info(const int* dinfo, int m, int n, cudaStream_t stream) {
        int info = 0;
        fmm::gpu::check_cuda(cudaMemcpyAsync(&info, dinfo, sizeof(int), cudaMemcpyDeviceToHost, stream), "gesvdj info");
        fmm::gpu::check_cuda(cudaStreamSynchronize(stream), "gesvdj");
        if (info < 0) throw std::runtime_error("cusolver gesvdj failed, info " + std::to_string(info));
        if (info > 0 && info != std::min(m, n) + 1) throw std::runtime_error("cusolver gesvdj failed, info " + std::to_string(info));
    }

    // (gesvdp's perturbation estimate allowed, relative to s(0), per unit of the truncation's tolerance)
    static constexpr double kSvdCheck = 0.1;

    cusolverDnHandle_t handle_ = nullptr;
    cusolverDnParams_t params_ = nullptr;
    gesvdjInfo_t jinfo_ = nullptr;
};

}  // namespace gpu
}  // namespace bpack

#endif  // H2_HAVE_GPU
