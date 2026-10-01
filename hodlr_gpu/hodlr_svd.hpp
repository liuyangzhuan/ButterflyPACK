#pragma once
// SVD of one dense device matrix (cuSOLVER gesvdp, polar decomposition
// based): the recompression of the HODLR blocks of large rank and the root
// of the merges' TSQR.

#ifdef H2_HAVE_GPU

#include "color_gpu/device_heap.hpp"
#include "color_gpu/device_scalar.hpp"
#include "color_gpu/gpu_runtime.hpp"

#include <cusolverDn.h>

#include <algorithm>
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
