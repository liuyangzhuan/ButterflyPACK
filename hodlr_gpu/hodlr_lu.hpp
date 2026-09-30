#pragma once
// Batched LU helpers of the HODLR GPU backend, shared by the symmetric and
// unsymmetric factorizations: MAGMA's LU with the failed entries reported,
// the diagonals and pivots on the host (log-determinants, pivot checks), the
// explicit inverse from an LU, and the modified LU of the CPU code
// (getrfmodf90: pivots below jitter * |A|_F are raised to that size).

#ifdef H2_HAVE_GPU

#include "hodlr_batch.hpp"
#include "hodlr_device.hpp"

#include <cmath>
#include <complex>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace bpack {
namespace gpu {

inline std::complex<double> to_std(double v) { return {v, 0.0}; }
inline std::complex<double> to_std(fmm::gpu::dcomplex v) { return {v.re, v.im}; }
inline double from_std(std::complex<double> v, double*) { return v.real(); }
inline fmm::gpu::dcomplex from_std(std::complex<double> v, fmm::gpu::dcomplex*) { return {v.real(), v.imag()}; }

// log|det| and phase of an LU factor from its diagonal and 1-based pivots,
// as getrf_slogdet (SRC/MISC_DenseLA.f90); false for an exact zero pivot.
template<typename T>
bool lu_slogdet(const T* diag, const int* ipiv, int n, std::complex<double>& phase, double& logabs) {
    phase = {1.0, 0.0};
    logabs = 0.0;
    int nswap = 0;
    for (int i = 0; i < n; ++i) {
        if (ipiv[i] != i + 1) ++nswap;
        const std::complex<double> u = to_std(diag[i]);
        const double a = std::abs(u);
        if (a == 0.0) {
            phase = {0.0, 0.0};
            logabs = -std::numeric_limits<double>::max();
            return false;
        }
        phase *= u / a;
        logabs += std::log(a);
    }
    if (nswap % 2 == 1) phase = -phase;
    return true;
}

// Run a staged LU batch; the indices (in add() order) of its failed entries.
template<typename T>
std::vector<int> run_getrf(const GetrfBatch<T>& lu, char* meta_device, DeviceBuffer& work) {
    Context& ctx = Context::instance();
    const cudaStream_t stream = ctx.stream();
    const int count = lu.count();
    if (count == 0) return {};
    DeviceHeap& heap = DeviceHeap::instance();
    magma_int_t* d_info = heap.alloc<magma_int_t>(sizeof(magma_int_t) * count);
    lu.factor(meta_device, d_info, work, ctx.queue());
    std::vector<magma_int_t> info(static_cast<size_t>(count));
    check_cuda(cudaMemcpyAsync(info.data(), d_info, sizeof(magma_int_t) * count, cudaMemcpyDeviceToHost, stream),
               "getrf info");
    check_cuda(cudaStreamSynchronize(stream), "getrf");
    heap.free(d_info);
    std::vector<int> failed;
    for (int i = 0; i < count; ++i) {
        if (info[static_cast<size_t>(i)] != 0) failed.push_back(i);
    }
    return failed;
}

// The diagonals and pivots of the LU factors mats[i] (order sizes[i]), one
// after the other; stops on a pivot outside LAPACK's range.
template<typename T>
void lu_diagonals(const std::vector<const T*>& mats, const std::vector<const int*>& pivs,
                  const std::vector<int>& sizes, MetaBuilder& meta, DeviceBuffer& meta_device, std::vector<T>& diag,
                  std::vector<int>& piv, const char* what) {
    const cudaStream_t stream = Context::instance().stream();
    size_t total = 0;
    for (int n : sizes) total += static_cast<size_t>(n);
    diag.assign(total, T(0.0));
    piv.assign(total, 0);
    if (total == 0) return;
    DeviceHeap& heap = DeviceHeap::instance();
    T* d_diag = heap.alloc<T>(sizeof(T) * total);
    ItemList<DiagCopyItem<T>> items;
    size_t off = 0;
    for (size_t i = 0; i < mats.size(); ++i) {
        items.add({mats[i], sizes[i], sizes[i], d_diag + off}, sizes[i], sizes[i]);
        off += static_cast<size_t>(sizes[i]);
    }
    meta.clear();
    items.stage(meta);
    char* md = meta.upload(meta_device, stream);
    launch_copy_diagonal(items.device(md), items.count(), items.max_n, stream);
    check_cuda(cudaMemcpyAsync(diag.data(), d_diag, sizeof(T) * total, cudaMemcpyDeviceToHost, stream), "diag");
    off = 0;
    for (size_t i = 0; i < mats.size(); ++i) {
        check_cuda(cudaMemcpyAsync(piv.data() + off, pivs[i], sizeof(int) * sizes[i], cudaMemcpyDeviceToHost, stream),
                   "pivots");
        off += static_cast<size_t>(sizes[i]);
    }
    check_cuda(cudaStreamSynchronize(stream), "LU diagonals");
    heap.free(d_diag);
    off = 0;
    for (size_t i = 0; i < mats.size(); ++i) {  // LAPACK pivots: row i swapped with a row in [i, n)
        for (int r = 0; r < sizes[i]; ++r) {
            const int p = piv[off + static_cast<size_t>(r)];
            if (p < r + 1 || p > sizes[i]) {
                std::ostringstream oss;
                oss << "HODLR GPU: pivot " << p << " of row " << r + 1 << " in a " << what << " of order " << sizes[i];
                throw std::runtime_error(oss.str());
            }
        }
        off += static_cast<size_t>(sizes[i]);
    }
}

// inv[i] = A_i^{-1} (n x n, leading dimension n) from the LU factor lu[i]
// and its pivots: inv = P I, then the unit lower and the upper solves.
template<typename T>
void lu_inverses(const std::vector<const T*>& lu, const std::vector<const int*>& piv, const std::vector<int>& sizes,
                 const std::vector<T*>& inv, MetaBuilder& meta, DeviceBuffer& meta_device) {
    Context& ctx = Context::instance();
    const cudaStream_t stream = ctx.stream();
    ItemList<DiagItem<T>> ones;
    ItemList<RowSwapItem<T>> swaps;
    TrsmBatch<T> lower, upper;
    for (size_t i = 0; i < lu.size(); ++i) {
        const int n = sizes[i];
        if (n == 0) continue;
        check_cuda(cudaMemsetAsync(inv[i], 0, sizeof(T) * static_cast<size_t>(n) * n, stream), "inverse init");
        ones.add({inv[i], n, n, T(1.0)}, n, n);
        swaps.add({inv[i], n, n, n, piv[i]}, n, n);
        lower.add(const_cast<T*>(lu[i]), n, inv[i], n, n, n);
        upper.add(const_cast<T*>(lu[i]), n, inv[i], n, n, n);
    }
    meta.clear();
    ones.stage(meta);
    swaps.stage(meta);
    lower.stage(meta);
    upper.stage(meta);
    char* md = meta.upload(meta_device, stream);
    launch_add_diagonal(ones.device(md), ones.count(), ones.max_n, stream);
    launch_row_swaps(swaps.device(md), swaps.count(), swaps.max_n, stream);
    lower.solve(md, MagmaLower, MagmaUnit, ctx.queue());
    upper.solve(md, MagmaUpper, MagmaNonUnit, ctx.queue());
}

// LU with partial pivoting of the n x n column-major a on the host, with the
// pivot rule of getrfmodf90 and pgetrfmodf90 (LAPACK_?getrf2mod.f,
// SCALAPACK_p?getf2mod.f): the pivot is the entry of largest |re| + |im|
// (i?amax), and a pivot smaller than thresh in magnitude is set to thresh
// (with its phase; thresh itself for an exact zero).  For the rare matrices
// whose GPU LU has such a pivot.  ipiv: 1-based.
template<typename T>
void host_lu_threshold(std::vector<T>& a, int n, double thresh, std::vector<int>& ipiv) {
    auto at = [&](int i, int j) -> std::complex<double> { return to_std(a[static_cast<size_t>(i + j * n)]); };
    auto put = [&](int i, int j, std::complex<double> v) {
        a[static_cast<size_t>(i + j * n)] = from_std(v, static_cast<T*>(nullptr));
    };
    ipiv.assign(static_cast<size_t>(n), 0);
    for (int j = 0; j < n; ++j) {
        int p = j;
        double best = -1.0;
        for (int i = j; i < n; ++i) {  // i?amax: |re| + |im| for complex data
            const std::complex<double> c = at(i, j);
            const double v = std::abs(c.real()) + std::abs(c.imag());
            if (v > best) {
                best = v;
                p = i;
            }
        }
        ipiv[static_cast<size_t>(j)] = p + 1;
        if (p != j)
            for (int c = 0; c < n; ++c) {
                const std::complex<double> t = at(j, c);
                put(j, c, at(p, c));
                put(p, c, t);
            }
        std::complex<double> piv = at(j, j);
        if (std::abs(piv) < thresh) {
            piv = std::abs(piv) == 0.0 ? std::complex<double>(thresh) : piv / std::abs(piv) * thresh;
            put(j, j, piv);
        }
        for (int i = j + 1; i < n; ++i) put(i, j, at(i, j) / piv);
        for (int c = j + 1; c < n; ++c) {
            const std::complex<double> u = at(j, c);
            if (u == 0.0) continue;
            for (int i = j + 1; i < n; ++i) put(i, c, at(i, c) - at(i, j) * u);
        }
    }
}

}  // namespace gpu
}  // namespace bpack

#endif  // H2_HAVE_GPU
