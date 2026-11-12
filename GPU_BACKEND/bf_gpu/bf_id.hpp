#pragma once

#ifdef H2_HAVE_GPU

// The rank-reveal of butterfly construction, batched on the device
// (GPU_BACKEND/bf_gpu).
//
// Both halves of BF_compress_NlogN reduce to one operation per block, and the
// four branches of BF_compress_NlogN_oneblock_{R,C}_rankreveal differ only in
// where the result is stored and how the skeleton indices are mapped.  The
// operation is an interpolative decomposition of a panel A (m x n) of sampled
// entries:
//
//     geqp3modf90(core, jpvt, tau, tol_comp, SafeUnderflow, rank)  ! truncated QRCP
//     un_or_mqrf90(core, tau, dat, 'L', 'C', m, n, rank)           ! dat := Q^H dat
//     trsmf90(core, dat, 'L', 'U', 'N', 'N', rank, n)              ! dat := R11^-1 dat
//
// leaving A ~ A(:, J) X with J = jpvt(1:rank) and X = dat(1:rank, :).
//
// It also serves forwardN15flag 3: that path replaces only the choice of
// sampled ROWS (BF_compress_NlogN_tree_R_sample fills the same submats
// contract), after which the caller runs this same rank-reveal.  So one batch
// covers both settings.
//
// The factorization itself is fmm::gpu::launch_qrcp (h2_gpu/device_id.cu),
// which the H2 backend already uses for exactly this: a batched column-pivoted
// Householder QR with the host's rank rule
// (|R_kk| > max(tol |R_11|, abs_tol), the same rule as LAPACK_*geqp3mod.f's
// early exit), leaving the skeleton first in jpvt and T = R11^-1 R12 in rows
// 0..K-1 of columns K..n-1.  It is real and complex, bit-deterministic across
// its three launch configurations, and has a cooperative multi-block path -
// which matters here, because a level holds 2^level_butterfly blocks globally
// and so a rank may see only tens of them.
//
// What is left for this file is the butterfly-specific glue, two small
// kernels:
//
//   * load the level's panels into device matrices, transposing the ones the
//     'C' half factorizes (it works on A^T and transposes the result back);
//   * turn (jpvt, rank, T) into the X the butterfly stores.  The CPU gets X by
//     applying Q^H and solving, which yields the identity on the skeleton
//     columns only to rounding; here X(:, jpvt[k]) = e_k exactly and
//     X(:, jpvt[K+j]) = T(:, j), which is the same decomposition written
//     directly.
//
// Pivot order is not bit-reproducible against LAPACK: the rule is the same but
// the reduction order differs, so a near-tie can pick the other column.
// Validation therefore compares the rank and the reconstruction, not the
// pivots.

#include <cstdint>
#include <vector>

#include "bf_device.hpp"
#include "../h2_gpu/device_kernels.hpp"

namespace bf {
namespace gpu {

// One block of a level: where its entries are, its shape, and whether it is
// factorized transposed.  `entry_offset` is into the level's single contiguous
// buffer - the Fortran fills one array for the whole level
// (`alldat_loc_in`) and points each block at an offset into it, so the upload
// is one copy.
struct IdPanel {
    int64_t entry_offset = 0;   // into the level's entry buffer, in elements
    int rows = 0, cols = 0;     // as the Fortran laid it out (column major, lda = rows)
    int transposed = 0;         // factorize A^T ('C' half)
    int64_t mat_offset = 0;     // into the factorization buffer (m * n)
    int64_t x_offset = 0;       // into the X buffer (m * n, X is its first `rank` rows)
    int64_t jpvt_offset = 0;    // into the pivot buffer (n)

#ifdef __CUDACC__
    __host__ __device__
#endif
    int m() const { return transposed ? cols : rows; }
#ifdef __CUDACC__
    __host__ __device__
#endif
    int n() const { return transposed ? rows : cols; }
};

// A level's batch.  `rank` and `jpvt` come back to the host: the next level's
// sampling and the skeleton exchange need them, and the skeleton index mapping
// is integer bookkeeping that stays in the Fortran.
template<typename T>
struct IdBatch {
    using DataType = T;
    std::vector<IdPanel> panels;
    IdPanel* d_panels = nullptr;
    DataType* d_mat = nullptr;    // the panels, destroyed by the factorization
    DataType* d_x = nullptr;      // X, rank x n per panel with lda = m
    int* d_jpvt = nullptr;
    int* d_rank = nullptr;
    double* d_norm = nullptr;
    int* d_flag = nullptr;
    fmm::gpu::QrcpItemT<DataType>* d_items = nullptr;
    void* d_work = nullptr;
    int64_t mat_elements = 0, jpvt_elements = 0;
    int max_n = 0;

    int count() const { return static_cast<int>(panels.size()); }
};

// Lay the panels out and take the device buffers from the shared arena,
// filling in each panel's offsets; the caller sets only shape and entry
// offset.  False when the arena has no room.
// (each of these defined for double and fmm::gpu::dcomplex in bf_id.cu)
template<typename T>
bool plan_id_batch(IdBatch<T>& batch);
template<typename T>
void release_id_batch(IdBatch<T>& batch);

// Load, factorize and build X for a whole level.  `entries` is device
// resident; `rtol` and `atol` are option%tol_comp and BPACK_SafeUnderflow.
template<typename T>
void run_id_batch(const T* entries, IdBatch<T>& batch, double rtol, double atol, cudaStream_t stream);

// The two glue kernels, separately for the harness.
template<typename T>
void launch_id_load(const T* entries, const IdPanel* panels, int count, T* mat, cudaStream_t stream);
template<typename T>
void launch_id_scatter(const IdPanel* panels, int count, const T* mat, const int* jpvt, const int* rank, T* x,
                       cudaStream_t stream);

}  // namespace gpu
}  // namespace bf

#endif  // H2_HAVE_GPU
