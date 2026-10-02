#pragma once

// H2 (format-7) integration: log-determinant accumulation over the factored
// hierarchy (Bunch-Kaufman and LU blocks), shared by both backends.

#include <cmath>
#include <complex>
#include <cstdint>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "tree.hpp"

namespace butterfly {
using namespace fmm;

// To Do: need to find better place for this function
// thresh: block-determinant magnitudes at or below this are treated as singular
// (throws instead of feeding 0 into std::log and producing -inf). Default 0.0
// guards only the exactly-zero case; pass a small positive value to also catch
// near-singular pivots.
template<typename DataType>
inline void accumulate_logdet_bunch_kaufman(const DataType* A, int64_t r, int64_t lda,
                                            const std::vector<int>& pivots,
                                            double& logabs, double& arg,
                                            double thresh = 0.0) {
    int64_t k = 0;
    while (k < r) {
        if (pivots[static_cast<size_t>(k)] > 0) {          // 1×1 pivot
            DataType d = A[k + k * lda];
            double ad = std::abs(d);
            if (!(ad > thresh) || !std::isfinite(ad)) {
                std::ostringstream oss;
                oss << "accumulate_logdet_bunch_kaufman: singular/near-singular "
                    << "1x1 pivot at k=" << k
                    << " (|d|=" << std::scientific << std::setprecision(17) << ad
                    << ", threshold=" << thresh << ", d=" << d << ')';
                throw std::runtime_error(oss.str());
            }
            logabs += std::log(ad);
            arg    += std::arg(d);
            k += 1;
        } else {                                            // 2×2 pivot (pivots[k]==pivots[k+1]<0)
            if (k + 1 >= r) {
                throw std::runtime_error(
                    "accumulate_logdet_bunch_kaufman: 2x2 pivot at last index k=" +
                    std::to_string(k) + " (r=" + std::to_string(r) + "), malformed pivots");
            }
            DataType a = A[k       + k       * lda];
            DataType c = A[(k + 1) + (k + 1) * lda];
            DataType b = A[(k + 1) + k       * lda];        // sub-diagonal
            DataType det2 = a * c - b * b;
            double ad2 = std::abs(det2);
            if (!(ad2 > thresh) || !std::isfinite(ad2)) {
                std::ostringstream oss;
                oss << "accumulate_logdet_bunch_kaufman: singular/near-singular "
                    << "2x2 pivot at k=" << k
                    << " (|det2|=" << std::scientific << std::setprecision(17) << ad2
                    << ", threshold=" << thresh
                    << ", a=" << a << ", b=" << b << ", c=" << c << ')';
                throw std::runtime_error(oss.str());
            }
            logabs += std::log(ad2);
            arg    += std::arg(det2);
            k += 2;
        }
    }
}

template<typename DataType>
inline void accumulate_logdet_lu(const DataType* A, int64_t r, int64_t lda,
                                 const std::vector<int>& pivots,
                                 double& logabs, double& arg,
                                 double thresh = 1e-14) {
    if (pivots.size() != static_cast<size_t>(r)) {
        throw std::runtime_error(
            "accumulate_logdet_lu: pivot count (" + std::to_string(pivots.size()) +
            ") does not match matrix size (" + std::to_string(r) + ")");
    }

    for (int64_t k = 0; k < r; ++k) {
        const DataType diagonal = A[k + k * lda];
        const double magnitude = std::abs(diagonal);
        if (magnitude <= thresh) {
            throw std::runtime_error(
                "accumulate_logdet_lu: singular/near-singular pivot (|u_kk| = " +
                std::to_string(magnitude) + ")");
        }

        const int pivot = pivots[static_cast<size_t>(k)];
        if (pivot < 1 || pivot > r) {
            throw std::runtime_error(
                "accumulate_logdet_lu: invalid LAPACK pivot " +
                std::to_string(pivot) + " at index " + std::to_string(k));
        }

        logabs += std::log(magnitude);
        arg += std::arg(diagonal);
        if (pivot != k + 1) {
            arg += std::acos(-1.0);
        }
    }
}


template<typename CoordType, typename DataType>
void hierarchical_logdet_parallel(fmm::ParallelTree<CoordType, DataType>* tree,
                                  double* logabsdet, DataType* phase) {

    int leaf_level = tree->num_levels - 1;

    double logabs_local = 0.0;   // Σ log|det|
    double arg_local    = 0.0;   // Σ arg(det)

    // iterate through all the levels
    for (int level = leaf_level; level >= 2; level--) {

        auto& tree_level = tree->levels[level];
        if (!tree_level.is_process_active) continue;

        std::exception_ptr diagonal_exception;
        std::mutex diagonal_exception_mutex;
        std::atomic<bool> diagonal_failed{false};

        #pragma omp parallel default(shared) if (tree_level.num_boxes_local > 1)
        {
            double t_logabs = 0.0, t_arg = 0.0;   // thread-local

            // iterate through all the boxes (in parallel)
            #pragma omp for schedule(static)
            for (int64_t box_idx = 0; box_idx < tree_level.num_boxes_local; ++box_idx) {
                if (diagonal_failed.load(std::memory_order_relaxed)) {
                    continue;
                }

                try {
                    auto& box = tree_level.local_boxes[static_cast<size_t>(box_idx)];

                    if (box.redundant_indices.empty()) {
                        continue;
                    }

                    int64_t r = static_cast<int64_t>(box.redundant_indices.size());
                    int64_t lda = box.X_RR.rows;
                    if (r != lda) {
                        throw std::runtime_error(
                            "logdet: redundant_indices.size() (" + std::to_string(r) +
                            ") != X_RR.rows (" + std::to_string(lda) + ")");
                    }
                    const DataType* A = box.X_RR.data.data();

                    // compute log|det| and phase for each box
                    if (box.X_RR.format == MatrixStorage<DataType>::CHOLESKY_L) {
                        throw std::runtime_error(
                            "logdet: CHOLESKY_L is not implemented yet");
                    } else if (box.X_RR.format == MatrixStorage<DataType>::LU_FACTORED) {
                        accumulate_logdet_lu(
                            A, r, lda, box.X_RR_pivots, t_logabs, t_arg);
                    } else if (box.X_RR.format == MatrixStorage<DataType>::BUNCH_KAUFMAN) {
                        if (box.X_RR_pivots.size() != static_cast<size_t>(r)) {
                            throw std::runtime_error(
                                "logdet: pivots is incorrect (X_RR_pivots.size() = " +
                                std::to_string(box.X_RR_pivots.size()) +
                                ", expected " + std::to_string(r) + ")");
                        }
                        accumulate_logdet_bunch_kaufman(A, r, lda, box.X_RR_pivots, t_logabs, t_arg);
                    } else {
                        throw std::runtime_error("Diagonal multiply: unsupported X_RR format");
                    }

                } catch (...) {
                    if (!diagonal_failed.exchange(true, std::memory_order_relaxed)) {
                        std::lock_guard<std::mutex> lock(diagonal_exception_mutex);
                        diagonal_exception = std::current_exception();
                    }
                }
            }

            #pragma omp atomic
            logabs_local += t_logabs;
            #pragma omp atomic
            arg_local    += t_arg;
        }
        if (diagonal_exception) {
            std::rethrow_exception(diagonal_exception);
        }
    }

    // handle root level
    auto& root_level = tree->levels[0];
    if (root_level.is_process_active && !root_level.local_boxes.empty()) {
        auto& rb = root_level.local_boxes[0];
        if (rb.X_RR.format == MatrixStorage<DataType>::BUNCH_KAUFMAN) {
            accumulate_logdet_bunch_kaufman(rb.X_RR.data.data(), rb.X_RR.rows, rb.X_RR.rows,
                                            rb.X_RR_pivots, logabs_local, arg_local);
        } else if (rb.X_RR.format == MatrixStorage<DataType>::LU_FACTORED) {
            accumulate_logdet_lu(rb.X_RR.data.data(), rb.X_RR.rows, rb.X_RR.rows,
                                 rb.X_RR_pivots, logabs_local, arg_local);
        } else {
            throw std::runtime_error(
                "logdet (root): unsupported X_RR factorization format");
        }
    }

    // accumulate logabsdet and phase
    double buf[2] = { logabs_local, arg_local };
    MPI_Allreduce(MPI_IN_PLACE, buf, 2, MPI_DOUBLE, MPI_SUM, tree->comm);
    *logabsdet = buf[0];
    if constexpr (std::is_same_v<DataType, std::complex<double>>) {
        *phase = DataType(std::cos(buf[1]), std::sin(buf[1]));   // e^{iθ}, any angle
    } else {
        *phase = std::cos(buf[1]);   // ±1: for real dets, θ is a multiple of π, so cos θ = ±1, sin θ ≈ 0
    }
}

} // namespace butterfly
