#pragma once
// Device kernels of the adaptive ID training rows (H2_ID_proxy 2): what
// adaptive_rows.hpp launches besides the element kernels and the pivoted QR of
// device_kernels.hpp.  One thread block per item; every result depends only on
// its item (sums in a fixed order), so a box gives the same bits in any batch.

#ifdef H2_HAVE_GPU

#include "device_scalar.hpp"

#include <cuda_runtime.h>

#include <cstdint>

namespace fmm {
namespace gpu {

// An interpolative decomposition as launch_qrcp leaves it: the rank K (read
// on the device from rank_ptr when set, else `rank`), the column order jpvt
// (skeleton first) and T (K x (n - K)) in rows 0..K-1 of columns K..n-1 of f.
template<typename T>
struct IdRefT {
    const T* f;
    int ld;
    const int* jpvt;
    const int* rank_ptr;
    int rank;
};

// r (rows x n) := its residual under the IDs ids[0..nids), applied in order
// as subtract_id_approximation does (color_CA/factorization.hpp):
//   r(:, redundant j) -= r(:, skeleton) T(:, j), then r(:, skeleton) = 0
// (an ID of rank n zeroes r, one of rank 0 leaves it).  norm, when set: the
// Frobenius norm of the result (of r itself for nids = 0).
template<typename T>
struct IdResidualItemT {
    T* r;
    int ld;
    int rows;
    int n;
    int64_t ids_offset;  // IdRefT array in the metadata
    int nids;
    double* norm;
};
template<typename T>
void launch_id_residual(const IdResidualItemT<T>* items, int count, const char* meta, cudaStream_t stream);

// Rows of src (s x n) chosen by a row ID (launch_qrcp of the transposed
// residual: *rank rows, in the order jpvt), written to dst from row0: all s
// rows in their order when the rank is s (compute_id_complex's full-rank
// convention), else rows jpvt[0..rank).
template<typename T>
struct AppendRowsItemT {
    T* dst;
    int ldd;
    int row0;
    const T* src;
    int lds;
    int s;
    int n;
    const int* rank;
    const int* jpvt;
};
template<typename T>
void launch_append_rows(const AppendRowsItemT<T>* items, int count, cudaStream_t stream);

}  // namespace gpu
}  // namespace fmm

#endif  // H2_HAVE_GPU
