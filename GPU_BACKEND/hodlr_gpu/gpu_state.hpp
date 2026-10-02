#pragma once
// GPU state of one ButterflyPACK matrix of a format other than H2
// (Bmatrix%gpu): the evaluator of the entries the application registered
// and, for HODLR, the device copies of the blocks.  H2 keeps its own
// registration in its H2Kernel (butterfly_types.hpp).

#ifdef H2_HAVE_GPU

#include "h2_gpu/evaluator.hpp"
#include "hodlr_construct.hpp"
#include "hodlr_device.hpp"
#include "hodlr_distqr.hpp"
#include "hodlr_svd.hpp"
#include "hodlr_sym.hpp"
#include "hodlr_unsym.hpp"

#include <chrono>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

namespace bpack {
namespace gpu {

// The batched LU, triangular solves and GEMM forms of the HODLR
// factorization and its solves, once per size range (MAGMA picks its kernels
// by size), before the process's first HODLR GPU work, so that no step pays
// the first-use loading of their kernels.  Once per process and value type;
// returns its seconds (0 when it already ran).
template<typename T>
double warm_up_hodlr_kernels() {
    static bool done = false;
    if (done) return 0.0;
    done = true;
    const auto t0 = std::chrono::steady_clock::now();
    fmm::gpu::Context& ctx = fmm::gpu::Context::instance();
    ctx.activate();
    const cudaStream_t stream = ctx.stream();
    const magma_queue_t queue = ctx.queue();
    fmm::gpu::DeviceHeap& heap = fmm::gpu::DeviceHeap::instance();
    fmm::gpu::MetaBuilder& meta = fmm::gpu::pinned_pool().meta;  // (no H2 level runs now)
    fmm::gpu::DeviceBuffer meta_device, work;
    for (int n : {8, 16, 32, 64, 128, 256, 512, 1024}) {
        T* a = heap_identity<T>(n, n);
        T* b = heap_identity<T>(n, n);
        T* c = heap_identity<T>(n, n);
        int* ipiv = heap.alloc<int>(static_cast<size_t>(n) * sizeof(int));
        magma_int_t* info = heap.alloc<magma_int_t>(sizeof(magma_int_t));
        GetrfBatch<T> lu;
        lu.add(a, n, n, ipiv);
        TrsmBatch<T> solve;
        solve.add(a, n, b, n, n, n);
        fmm::gpu::VBatch<T> g;
        g.entries.push_back({a, b, c, n, n, n, n, n, n});
        meta.clear();
        lu.stage(meta);
        solve.stage(meta);
        g.stage(meta);
        char* md = meta.upload(meta_device, stream);
        lu.factor(md, info, work, queue);
        solve.solve(md, MagmaLower, MagmaUnit, queue);
        solve.solve(md, MagmaUpper, MagmaNonUnit, queue);
        g.gemm(md, MagmaNoTrans, MagmaNoTrans, T(1.0), T(0.0), queue);
        g.gemm(md, MagmaNoTrans, MagmaTrans, T(1.0), T(1.0), queue);
        g.gemm(md, MagmaTrans, MagmaNoTrans, T(1.0), T(1.0), queue);
        if constexpr (fmm::gpu::is_complex_scalar<T>) g.gemm(md, MagmaConjTrans, MagmaNoTrans, T(1.0), T(1.0), queue);
        fmm::gpu::check_cuda(cudaStreamSynchronize(stream), "warm-up");
        for (void* p : {static_cast<void*>(info), static_cast<void*>(ipiv), static_cast<void*>(c), static_cast<void*>(b),
                        static_cast<void*>(a)}) {
            heap.free(p);
        }
    }
    meta.clear();
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

template<typename T>
struct GpuState {
    // the application's evaluator of the entries (c_bpack_set_gpu_*_evaluator,
    // GPU_INTERFACE/bpack_gpu.h); null: none (the entries come from the host)
    std::shared_ptr<fmm::gpu::Evaluator> evaluator;
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
    // Before the matrix's first GPU work (c_bpack_hodlr_gpu_warm_up, ahead of
    // the construction and not in its time): the kernels of the
    // factorization (once per process), with `construct` the routines and
    // cuSOLVER handles of the GPU construction, and with `evaluate` the
    // application's evaluator; the seconds of the first two together, then
    // of the evaluator.
    std::pair<double, double> warm_up(bool construct, bool evaluate) {
        const auto t0 = std::chrono::steady_clock::now();
        warm_up_hodlr_kernels<T>();
        if (construct) {
            hodlr_construct().warm_up();
            device_svd().warm_up();
            dist_qr().warm_up();
        }
        const double kernels = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        double evaluated = 0.0;
        if (evaluate && evaluator) evaluated = hodlr_construct().warm_up_evaluator(*evaluator);
        return {kernels, evaluated};
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
