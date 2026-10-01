#pragma once
// GPU state of one ButterflyPACK matrix of a format other than H2
// (Bmatrix%gpu): the device kernel the application registered
// (c_bpack_set_gpu_kernel) and, for HODLR, the device copies of the blocks.
// H2 keeps its own registration in its H2Kernel (butterfly_types.hpp).

#ifdef H2_HAVE_GPU

#include "hodlr_construct.hpp"
#include "hodlr_device.hpp"
#include "hodlr_distqr.hpp"
#include "hodlr_svd.hpp"
#include "hodlr_sym.hpp"
#include "hodlr_unsym.hpp"

#include <cstdint>
#include <memory>
#include <vector>

namespace bpack {
namespace gpu {

// Same fields as H2Kernel::GpuSpec, so the H2 backend's helpers
// (fmm::gpu::device_kernel_spec, device_kernel_registered) take it as is.
struct KernelRegistration {
    int kind = 0;  // 0: none; 1 real, 2 and 3 complex (see fmm::gpu::KernelSpec)
    double params[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<double> table_real;
    std::vector<int> table_int;
    uint64_t table_version = 0;
};

template<typename T>
struct GpuState {
    KernelRegistration kernel;
    std::unique_ptr<HodlrDevice<T>> hodlr;
    std::unique_ptr<HodlrSymFactor<T>> sym;      // after a symmetric factorization
    std::unique_ptr<HodlrUnsymFactor<T>> unsym;  // after an unsymmetric factorization
    std::unique_ptr<HodlrConstruct<T>> construct;  // HODLR construction (points, BACA of a level)
    std::unique_ptr<DeviceSvd<T>> svd;               // dense SVDs of the recompression
    std::unique_ptr<DistQr<T>> distqr;               // TSQR of the merges of shared blocks

    HodlrDevice<T>& hodlr_device() {
        if (!hodlr) hodlr = std::make_unique<HodlrDevice<T>>();
        return *hodlr;
    }
    HodlrConstruct<T>& hodlr_construct() {
        if (!construct) construct = std::make_unique<HodlrConstruct<T>>();
        return *construct;
    }
    DistQr<T>& dist_qr() {
        if (!distqr) distqr = std::make_unique<DistQr<T>>();
        return *distqr;
    }
    DeviceSvd<T>& device_svd() {
        if (!svd) svd = std::make_unique<DeviceSvd<T>>();
        return *svd;
    }
    // Drop the factors (a new forward matrix or factorization follows).
    void reset_factors() {
        sym.reset();
        unsym.reset();
    }
};

}  // namespace gpu
}  // namespace bpack

#endif  // H2_HAVE_GPU
