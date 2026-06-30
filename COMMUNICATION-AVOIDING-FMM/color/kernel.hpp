#ifndef KERNEL_HPP
#define KERNEL_HPP

#include <cmath>
#include <vector>
#include <stdexcept>
#include <cstdint>
#include <complex>
#include <fftw3.h>
#include <algorithm>
#include <immintrin.h>
#include <cassert>

namespace kernel {

/**
 * @brief 2D Laplace kernel implementation
 * 
 * Implements the 2D Laplace kernel: K(x,y) = -1/(2π) * log(||x-y||)
 */
template<typename T>
struct LaplaceKernel2D {
    private:
        int64_t N_total;  // Total number of points for scaling
        T precomputed_diag_val;  // Precomputed diagonal value

    public:
        static constexpr int dimension = 2;
        static constexpr T pi = 3.14159265358979323846;

        
        LaplaceKernel2D(int64_t N) : N_total(N) {
            precomputed_diag_val = static_cast<T>(evaluate_diagonal(N));
        }
        
        /**
        * @brief Evaluate 2D Laplace kernel between two points
        */
        static T evaluate(const T* x, const T* y) {
            T dx = x[0] - y[0];
            T dy = x[1] - y[1];
            T r = std::sqrt(dx * dx + dy * dy);
            
            if (r < 1e-14) {
                return 0.0;
            }
            
            return -1.0 / (2.0 * pi) * std::log(r);
        }
        
        /**
        * @brief Evaluate diagonal self-interaction using 5-point Gauss-Legendre quadrature
        * 
        * Computes: ∫∫_{[0,h/2]²} K(r,0) dr where h = 1/√N
        * 
        * @param N Total number of points (assumes uniform √N × √N grid)
        * @return Integrated diagonal value (negative)
        */
        static T evaluate_diagonal(int64_t N) {
            int64_t n = static_cast<int64_t>(std::sqrt(static_cast<double>(N)));
            T h = 1.0 / static_cast<T>(n);
            
            // Integrand: K(sqrt(x^2 + y^2), 0) = -1/(2π) * log(r)
            auto integrand = [](T x, T y) -> T {
                T r = std::sqrt(x * x + y * y);
                if (r < 1e-14) return 0.0;
                return -1.0 / (2.0 * pi) * std::log(r);
            };
            
            // 5-point Gauss-Legendre quadrature nodes on [-1,1]
            const std::vector<T> nodes = {
                -0.9061798459386640, -0.5384693101056831, 0.0,
                0.5384693101056831,  0.9061798459386640
            };
            
            // Corresponding weights
            const std::vector<T> weights = {
                0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
                0.4786286704993665, 0.2369268850561891
            };
            
            // Transform nodes from [-1,1] to [0, h/2]
            T a = 0.0, b = h / 2.0;
            T transform = (b - a) / 2.0;  // Scaling factor
            T shift = (b + a) / 2.0;      // Translation
            
            // Compute 2D integral using tensor product quadrature
            T integral = 0.0;
            for (size_t i = 0; i < nodes.size(); ++i) {
                T xi = transform * nodes[i] + shift;  // Map to [0, h/2]
                for (size_t j = 0; j < nodes.size(); ++j) {
                    T yj = transform * nodes[j] + shift;
                    integral += weights[i] * weights[j] * integrand(xi, yj);
                }
            }
            
            // Apply Jacobian for coordinate transformation (transform²)
            integral *= transform * transform;
            
            // Multiply by 4 to account for:
            // 1. Integration over 1/4 of cell (factor of 4)
            return 4.0 * integral;
        }
        


        // void evaluate_block_precomputed_diagonal(
        //     const T* __restrict__ x_coords, int64_t x_size,
        //     const T* __restrict__ y_coords, int64_t y_size,
        //     int64_t N,
        //     T* __restrict__ A, int64_t lda) const {
            
        //     const T scale_factor = -1.0 / (2.0 * pi * static_cast<T>(N));
        //     const bool is_diagonal_block = (x_coords == y_coords) && (x_size == y_size);
            
        //     #pragma omp simd
        //     for (int64_t j = 0; j < y_size; ++j) {
        //         const T yj_x = y_coords[j * 2];
        //         const T yj_y = y_coords[j * 2 + 1];
                
        //         #pragma omp simd
        //         for (int64_t i = 0; i < x_size; ++i) {
        //             const T xi_x = x_coords[i * 2];
        //             const T xi_y = x_coords[i * 2 + 1];
                    
        //             const T dx = xi_x - yj_x;
        //             const T dy = xi_y - yj_y;
        //             const T r_sq = std::fma(dx, dx, dy * dy);  // Fused multiply-add
                    
        //             const T r_sq_safe = std::max(r_sq, T(1e-28));
        //             A[i + j * lda] = scale_factor * T(0.5) * std::log(r_sq_safe);
        //         }
        //     }
            
        //     if (is_diagonal_block) {
        //         #pragma omp simd
        //         for (int64_t i = 0; i < x_size; ++i) {
        //             A[i + i * lda] = precomputed_diag_val;
        //         }
        //     }
        // }
        void evaluate_block_precomputed_diagonal(
            const T* x_coords, int64_t x_size,
            const T* y_coords, int64_t y_size,
            int64_t N,
            T* __restrict__ A, int64_t lda) const {

            const T scale_factor = -1.0 / (2.0 * pi * static_cast<T>(N));

            for (int64_t j = 0; j < y_size; ++j) {
                const T yj_x = y_coords[j * 2];
                const T yj_y = y_coords[j * 2 + 1];

                // #pragma omp simd
                for (int64_t i = 0; i < x_size; ++i) {
                    const T xi_x = x_coords[i * 2];
                    const T xi_y = x_coords[i * 2 + 1];

                    const T dx = xi_x - yj_x;
                    const T dy = xi_y - yj_y;
                    const T r_sq = std::fma(dx, dx, dy * dy);

                    const T r_sq_safe = std::max(r_sq, T(1e-28));
                    A[i + j * lda] = scale_factor * T(0.5) * std::log(r_sq_safe);
                }

                const bool same_layout_ptr = (x_coords == y_coords) && (x_size == y_size);
                for (int64_t i = 0; i < x_size; ++i) {
                    const bool pointer_match = same_layout_ptr && (i == j);
                    const bool coord_match =
                        (x_coords[i * 2] == yj_x) &&
                        (x_coords[i * 2 + 1] == yj_y);
                    T match = T(pointer_match || coord_match);
                    A[i + j * lda] = match * precomputed_diag_val + (T(1.0) - match) * A[i + j * lda];
                }
            }
        }

        /**
        * @brief Instance method wrapper for kernel wrapper compatibility
        */
        void evaluate_block(
            const T* x_coords, int64_t x_size,
            const T* y_coords, int64_t y_size,
            T* A, int64_t lda) const {
            
            // Delegate to static method
            // kernel::LaplaceKernel2D::evaluate_block(x_coords, x_size, y_coords, y_size, N_total, A, lda);
            evaluate_block_precomputed_diagonal(x_coords, x_size, y_coords, y_size, N_total, A, lda);
        }
};

/**
 * @brief FFT-based fast matrix-vector multiplication for 2D Laplace kernel
 * 
 * Uses real-to-complex FFT to exploit:
 * 1. Toeplitz structure on regular grids
 * 2. Real-valued kernel → Hermitian FFT symmetry
 * 
 * Complexity: O(N log N) vs O(N²) for direct evaluation
 * Memory: ~50% savings vs complex FFT
 */
template <typename T>
class LaplaceKernel2D_FFT {
private:
    int64_t n;           // Grid dimension (n × n = N)
    int64_t N;           // Total points
    int64_t padded_size; // 2n-1 (circulant embedding)
    
    std::vector<std::complex<T>> G;  // Precomputed FFT (Hermitian, half-size)
    std::vector<T> work_real;        // Real working buffer
    std::vector<std::complex<T>> work_complex;  // Complex working buffer
    
    fftw_plan plan_forward;
    fftw_plan plan_backward;

public:
    /**
     * @brief Constructor: precompute FFT of circulant embedding
     * 
     * @param grid_points Point coordinates (column-major: 2 × N)
     * @param grid_size Grid dimension n (must satisfy n² = N)
     * @param total_points Total number of points N
     */
    LaplaceKernel2D_FFT(const T* grid_points, int64_t grid_size, int64_t total_points)
        : n(grid_size), N(total_points), padded_size(2 * grid_size - 1) {
        
        if (n * n != N) {
            throw std::invalid_argument("Grid must be square: n² ≠ N");
        }
        
        // Allocate buffers
        // Real: (2n-1)²
        // Complex: (2n-1) × (n) due to Hermitian symmetry
        work_real.resize(padded_size * padded_size);
        work_complex.resize(padded_size * (padded_size / 2 + 1));
        G.resize(padded_size * (padded_size / 2 + 1));
        
        // Step 1: Compute first column of kernel matrix
        std::vector<T> first_col(N);
        T y_fixed[2] = {grid_points[0], grid_points[1]};
        
        for (int64_t i = 0; i < N; ++i) {
            T x[2] = {grid_points[i * 2], grid_points[i * 2 + 1]};
            first_col[i] = (i == 0) ? 
                LaplaceKernel2D<T>::evaluate_diagonal(N) :
                LaplaceKernel2D<T>::evaluate(x, y_fixed) / static_cast<T>(N);
        }
        
        // Step 2: Reshape to n × n
        std::vector<T> A(n * n);
        for (int64_t i = 0; i < n; ++i) {
            for (int64_t j = 0; j < n; ++j) {
                A[i + j * n] = first_col[i * n + j];
            }
        }
        
        // Step 3: Build circulant embedding B in work_real
        std::fill(work_real.begin(), work_real.end(), 0.0);
        
        // B(1:n, 1:n) = A
        for (int64_t i = 0; i < n; ++i) {
            for (int64_t j = 0; j < n; ++j) {
                work_real[i + j * padded_size] = A[i + j * n];
            }
        }
        
        // B(1:n, n+1:2n-1) = A(:, n-1:-1:2) (flipped columns)
        for (int64_t i = 0; i < n; ++i) {
            for (int64_t j = 1; j < n; ++j) {
                work_real[i + (padded_size - j) * padded_size] = A[i + j * n];
            }
        }
        
        // B(n+1:2n-1, 1:n) = A(n-1:-1:2, :) (flipped rows)
        for (int64_t i = 1; i < n; ++i) {
            for (int64_t j = 0; j < n; ++j) {
                work_real[(padded_size - i) + j * padded_size] = A[i + j * n];
            }
        }
        
        // B(n+1:2n-1, n+1:2n-1) = A(n-1:-1:2, n-1:-1:2) (flipped both)
        for (int64_t i = 1; i < n; ++i) {
            for (int64_t j = 1; j < n; ++j) {
                work_real[(padded_size - i) + (padded_size - j) * padded_size] = A[i + j * n];
            }
        }
        
        // Step 4: Compute G = rfft2(B)
        if constexpr (std::is_same_v<T, double>) {
            fftw_plan temp_plan = fftw_plan_dft_r2c_2d(
                padded_size, padded_size,
                work_real.data(),
                reinterpret_cast<fftw_complex*>(G.data()),
                FFTW_ESTIMATE
            );
            fftw_execute(temp_plan);
            fftw_destroy_plan(temp_plan);
            
            // Create reusable plans for matvec
            plan_forward = fftw_plan_dft_r2c_2d(
                padded_size, padded_size,
                work_real.data(),
                reinterpret_cast<fftw_complex*>(work_complex.data()),
                FFTW_ESTIMATE
            );
            
            plan_backward = fftw_plan_dft_c2r_2d(
                padded_size, padded_size,
                reinterpret_cast<fftw_complex*>(work_complex.data()),
                work_real.data(),
                FFTW_ESTIMATE
            );
        } else {
            throw std::runtime_error("Only double precision supported for FFT");
        }
    }
    
    ~LaplaceKernel2D_FFT() {
        if constexpr (std::is_same_v<T, double>) {
            fftw_destroy_plan(plan_forward);
            fftw_destroy_plan(plan_backward);
        }
    }
    
    /**
     * @brief Fast matrix-vector multiplication: y = A * x
     * 
     * @param x Input vector (size N)
     * @param y Output vector (size N, pre-allocated)
     */
    void matvec(const T* x, T* y) {
        // Step 1: Zero-pad x to (2n-1)²
        std::fill(work_real.begin(), work_real.end(), 0.0);
        
        for (int64_t i = 0; i < n; ++i) {
            for (int64_t j = 0; j < n; ++j) {
                work_real[i + j * padded_size] = x[i * n + j];
            }
        }
        
        // Step 2: Real-to-complex FFT
        if constexpr (std::is_same_v<T, double>) {
            fftw_execute(plan_forward);
        }
        
        // Step 3: Element-wise multiplication with G
        for (int64_t i = 0; i < padded_size * (padded_size / 2 + 1); ++i) {
            work_complex[i] *= G[i];
        }
        
        // Step 4: Complex-to-real inverse FFT
        if constexpr (std::is_same_v<T, double>) {
            fftw_execute(plan_backward);
        }
        
        // Step 5: Extract first n × n and normalize
        // FFTW's c2r doesn't normalize, so we divide by (2n-1)²
        T scale = 1.0 / static_cast<T>(padded_size * padded_size);
        for (int64_t i = 0; i < n; ++i) {
            for (int64_t j = 0; j < n; ++j) {
                y[i * n + j] = work_real[i + j * padded_size] * scale;
            }
        }
    }
    
    /**
     * @brief Get memory usage in bytes
     */
    size_t memory_usage() const {
        size_t real_mem = work_real.size() * sizeof(T);
        size_t complex_mem = (G.size() + work_complex.size()) * sizeof(std::complex<T>);
        return real_mem + complex_mem;
    }
};



/**
 * @brief 3D Laplace kernel implementation
 * 
 * Implements the 3D Laplace kernel: K(x,y) = 1/(4π||x-y||)
 */
template<typename T>
struct LaplaceKernel3D {
    private:
        int64_t N_total;  // Total number of points for scaling
        T precomputed_diag_val;  // Precomputed diagonal value

    public:
        static constexpr int dimension = 3;
        static constexpr T pi = 3.14159265358979323846;

        
        LaplaceKernel3D(int64_t N) : N_total(N) {
            precomputed_diag_val = static_cast<T>(evaluate_diagonal(N));
        }
        
        /**
        * @brief Evaluate 3D Laplace kernel between two points
        */
        static T evaluate(const T* x, const T* y) {
            T dx = x[0] - y[0];
            T dy = x[1] - y[1];
            T dz = x[2] - y[2];
            T r = std::sqrt(dx * dx + dy * dy + dz * dz);
            
            if (r < 1e-14) {
                return 0.0;
            }
            
            return 1.0 / (4.0 * pi * r);
        }
        
        /**
        * @brief Evaluate diagonal self-interaction using 5-point Gauss-Legendre quadrature
        * 
        * Computes: ∫∫∫_{[0,h/2]³} K(r,0) dr where h = 1/∛N
        * 
        * @param N Total number of points (assumes uniform ∛N × ∛N × ∛N grid)
        * @return Integrated diagonal value
        */
        static T evaluate_diagonal(int64_t N) {
            int64_t n = static_cast<int64_t>(std::cbrt(static_cast<double>(N)));
            T h = 1.0 / static_cast<T>(n);
            
            // Integrand: K(sqrt(x^2 + y^2 + z^2), 0) = 1/(4π * r)
            auto integrand = [](T x, T y, T z) -> T {
                T r = std::sqrt(x * x + y * y + z * z);
                if (r < 1e-14) return 0.0;
                return 1.0 / (4.0 * pi * r);
            };
            
            // 5-point Gauss-Legendre quadrature nodes on [-1,1]
            const std::vector<T> nodes = {
                -0.9061798459386640, -0.5384693101056831, 0.0,
                 0.5384693101056831,  0.9061798459386640
            };
            
            // Corresponding weights
            const std::vector<T> weights = {
                0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
                0.4786286704993665, 0.2369268850561891
            };
            
            // Transform nodes from [-1,1] to [0, h/2]
            T a = 0.0, b = h / 2.0;
            T transform = (b - a) / 2.0;  // Scaling factor
            T shift = (b + a) / 2.0;      // Translation
            
            // Compute 3D integral using tensor product quadrature
            T integral = 0.0;
            for (size_t i = 0; i < nodes.size(); ++i) {
                T xi = transform * nodes[i] + shift;
                for (size_t j = 0; j < nodes.size(); ++j) {
                    T yj = transform * nodes[j] + shift;
                    for (size_t k = 0; k < nodes.size(); ++k) {
                        T zk = transform * nodes[k] + shift;
                        integral += weights[i] * weights[j] * weights[k] * integrand(xi, yj, zk);
                    }
                }
            }
            
            // Apply Jacobian for coordinate transformation (transform³)
            integral *= transform * transform * transform;
            
            // Multiply by 8 to account for integration over 1/8 of cell
            return 8.0 * integral;
        }
        

        // /**
        // * @brief Optimized evaluation with precomputed diagonal
        // */
        // void evaluate_block_precomputed_diagonal(
        //     const T* __restrict__ x_coords, int64_t x_size,
        //     const T* __restrict__ y_coords, int64_t y_size,
        //     int64_t N,
        //     T* __restrict__ A, int64_t lda) const {
            
        //     const T scale_factor = 1.0 / (4.0 * pi * static_cast<T>(N));
        //     const bool is_diagonal_block = (x_coords == y_coords) && (x_size == y_size);
            
        //     #pragma omp simd
        //     for (int64_t j = 0; j < y_size; ++j) {
        //         const T yj_x = y_coords[j * 3];
        //         const T yj_y = y_coords[j * 3 + 1];
        //         const T yj_z = y_coords[j * 3 + 2];
                
        //         #pragma omp simd
        //         for (int64_t i = 0; i < x_size; ++i) {
        //             const T xi_x = x_coords[i * 3];
        //             const T xi_y = x_coords[i * 3 + 1];
        //             const T xi_z = x_coords[i * 3 + 2];
                    
        //             const T dx = xi_x - yj_x;
        //             const T dy = xi_y - yj_y;
        //             const T dz = xi_z - yj_z;
        //             const T r_sq = std::fma(dx, dx, std::fma(dy, dy, dz * dz));
                    
        //             const T r_safe = std::sqrt(std::max(r_sq, T(1e-28)));
        //             A[i + j * lda] = scale_factor / r_safe;
        //         }
        //     }
            
        //     if (is_diagonal_block) {
        //         #pragma omp simd
        //         for (int64_t i = 0; i < x_size; ++i) {
        //             A[i + i * lda] = precomputed_diag_val;
        //         }
        //     }
        // }

        void evaluate_block_precomputed_diagonal(
            const T* x_coords, int64_t x_size,
            const T* y_coords, int64_t y_size,
            int64_t N,
            T* __restrict__ A, int64_t lda) const {

            const T scale_factor = 1.0 / (4.0 * pi * static_cast<T>(N));

            for (int64_t j = 0; j < y_size; ++j) {
                const T yj_x = y_coords[j * 3];
                const T yj_y = y_coords[j * 3 + 1];
                const T yj_z = y_coords[j * 3 + 2];

                // #pragma omp simd
                for (int64_t i = 0; i < x_size; ++i) {
                    const T xi_x = x_coords[i * 3];
                    const T xi_y = x_coords[i * 3 + 1];
                    const T xi_z = x_coords[i * 3 + 2];

                    const T dx = xi_x - yj_x;
                    const T dy = xi_y - yj_y;
                    const T dz = xi_z - yj_z;
                    const T r_sq = std::fma(dx, dx, std::fma(dy, dy, dz * dz));

                    const T r_safe = std::sqrt(std::max(r_sq, T(1e-28)));
                    A[i + j * lda] = scale_factor / r_safe;
                }

                const bool same_layout_ptr = (x_coords == y_coords) && (x_size == y_size);
                for (int64_t i = 0; i < x_size; ++i) {
                    const bool pointer_match = same_layout_ptr && (i == j);
                    const bool coord_match =
                        (x_coords[i * 3]     == yj_x) &&
                        (x_coords[i * 3 + 1] == yj_y) &&
                        (x_coords[i * 3 + 2] == yj_z);
                    T match = T(pointer_match || coord_match);
                    A[i + j * lda] = match * precomputed_diag_val + (T(1.0) - match) * A[i + j * lda];
                }
            }
        }
        
        /**
        * @brief Instance method wrapper
        */
        void evaluate_block(
            const T* x_coords, int64_t x_size,
            const T* y_coords, int64_t y_size,
            T* A, int64_t lda) const {
            
            evaluate_block_precomputed_diagonal(x_coords, x_size, y_coords, y_size, N_total, A, lda);
        }
};

/**
 * @brief FFT-based fast matrix-vector multiplication for 3D Laplace kernel
 */
template <typename T>
class LaplaceKernel3D_FFT {
private:
    int64_t n;           // Grid dimension (n × n × n = N)
    int64_t N;           // Total points
    int64_t padded_size; // 2n-1 (circulant embedding)
    
    std::vector<std::complex<T>> G;  // Precomputed FFT (Hermitian, half-size)
    std::vector<T> work_real;
    std::vector<std::complex<T>> work_complex;
    
    fftw_plan plan_forward;
    fftw_plan plan_backward;

public:
    LaplaceKernel3D_FFT(const T* grid_points, int64_t grid_size, int64_t total_points)
        : n(grid_size), N(total_points), padded_size(2 * grid_size - 1) {
        
        if (n * n * n != N) {
            throw std::invalid_argument("Grid must be cubic: n³ ≠ N");
        }
        
        // Allocate buffers for 3D
        work_real.resize(padded_size * padded_size * padded_size);
        work_complex.resize(padded_size * padded_size * (padded_size / 2 + 1));
        G.resize(padded_size * padded_size * (padded_size / 2 + 1));
        
        // Step 1: Compute first column of kernel matrix
        std::vector<T> first_col(N);
        T y_fixed[3] = {grid_points[0], grid_points[1], grid_points[2]};
        
        for (int64_t i = 0; i < N; ++i) {
            T x[3] = {grid_points[i * 3], grid_points[i * 3 + 1], grid_points[i * 3 + 2]};
            first_col[i] = (i == 0) ? 
                LaplaceKernel3D<T>::evaluate_diagonal(N) :
                LaplaceKernel3D<T>::evaluate(x, y_fixed) / static_cast<T>(N);
        }
        
        // Step 2: Reshape to n × n × n
        std::vector<T> A(n * n * n);
        for (int64_t i = 0; i < n; ++i) {
            for (int64_t j = 0; j < n; ++j) {
                for (int64_t k = 0; k < n; ++k) {
                    A[i + j * n + k * n * n] = first_col[i + j * n + k * n * n];
                }
            }
        }
        
        // Step 3: Build 3D circulant embedding
        std::fill(work_real.begin(), work_real.end(), 0.0);
        
        for (int64_t i = 0; i < n; ++i) {
            for (int64_t j = 0; j < n; ++j) {
                for (int64_t k = 0; k < n; ++k) {
                    int64_t i_pad = (i == 0) ? 0 : (i <= n/2) ? i : padded_size - (n - i);
                    int64_t j_pad = (j == 0) ? 0 : (j <= n/2) ? j : padded_size - (n - j);
                    int64_t k_pad = (k == 0) ? 0 : (k <= n/2) ? k : padded_size - (n - k);
                    
                    if (i < n && j < n && k < n) {
                        i_pad = i;
                        j_pad = j;
                        k_pad = k;
                    }
                }
            }
        }
        
        // Simplified embedding for 3D
        for (int64_t k = 0; k < n; ++k) {
            for (int64_t j = 0; j < n; ++j) {
                for (int64_t i = 0; i < n; ++i) {
                    work_real[i + j * padded_size + k * padded_size * padded_size] = 
                        A[i + j * n + k * n * n];
                }
            }
        }
        
        // Mirror in all three dimensions
        for (int64_t k = 0; k < n; ++k) {
            for (int64_t j = 0; j < n; ++j) {
                for (int64_t i = 1; i < n; ++i) {
                    work_real[(padded_size - i) + j * padded_size + k * padded_size * padded_size] = 
                        work_real[i + j * padded_size + k * padded_size * padded_size];
                }
            }
        }
        
        for (int64_t k = 0; k < n; ++k) {
            for (int64_t j = 1; j < n; ++j) {
                for (int64_t i = 0; i < padded_size; ++i) {
                    work_real[i + (padded_size - j) * padded_size + k * padded_size * padded_size] = 
                        work_real[i + j * padded_size + k * padded_size * padded_size];
                }
            }
        }
        
        for (int64_t k = 1; k < n; ++k) {
            for (int64_t j = 0; j < padded_size; ++j) {
                for (int64_t i = 0; i < padded_size; ++i) {
                    work_real[i + j * padded_size + (padded_size - k) * padded_size * padded_size] = 
                        work_real[i + j * padded_size + k * padded_size * padded_size];
                }
            }
        }
        
        // Step 4: Compute G = rfft3(B)
        if constexpr (std::is_same_v<T, double>) {
            fftw_plan temp_plan = fftw_plan_dft_r2c_3d(
                padded_size, padded_size, padded_size,
                work_real.data(),
                reinterpret_cast<fftw_complex*>(G.data()),
                FFTW_ESTIMATE
            );
            fftw_execute(temp_plan);
            fftw_destroy_plan(temp_plan);
            
            plan_forward = fftw_plan_dft_r2c_3d(
                padded_size, padded_size, padded_size,
                work_real.data(),
                reinterpret_cast<fftw_complex*>(work_complex.data()),
                FFTW_ESTIMATE
            );
            
            plan_backward = fftw_plan_dft_c2r_3d(
                padded_size, padded_size, padded_size,
                reinterpret_cast<fftw_complex*>(work_complex.data()),
                work_real.data(),
                FFTW_ESTIMATE
            );
        } else {
            throw std::runtime_error("Only double precision supported for FFT");
        }
    }
    
    ~LaplaceKernel3D_FFT() {
        if constexpr (std::is_same_v<T, double>) {
            fftw_destroy_plan(plan_forward);
            fftw_destroy_plan(plan_backward);
        }
    }
    
    void matvec(const T* x, T* y) {
        // Zero-pad x to (2n-1)³
        std::fill(work_real.begin(), work_real.end(), 0.0);
        
        for (int64_t k = 0; k < n; ++k) {
            for (int64_t j = 0; j < n; ++j) {
                for (int64_t i = 0; i < n; ++i) {
                    work_real[i + j * padded_size + k * padded_size * padded_size] = 
                        x[i + j * n + k * n * n];
                }
            }
        }
        
        // Real-to-complex FFT
        if constexpr (std::is_same_v<T, double>) {
            fftw_execute(plan_forward);
        }
        
        // Element-wise multiplication
        for (int64_t i = 0; i < padded_size * padded_size * (padded_size / 2 + 1); ++i) {
            work_complex[i] *= G[i];
        }
        
        // Complex-to-real inverse FFT
        if constexpr (std::is_same_v<T, double>) {
            fftw_execute(plan_backward);
        }
        
        // Extract and normalize
        T scale = 1.0 / static_cast<T>(padded_size * padded_size * padded_size);
        for (int64_t k = 0; k < n; ++k) {
            for (int64_t j = 0; j < n; ++j) {
                for (int64_t i = 0; i < n; ++i) {
                    y[i + j * n + k * n * n] = 
                        work_real[i + j * padded_size + k * padded_size * padded_size] * scale;
                }
            }
        }
    }
    
    size_t memory_usage() const {
        size_t real_mem = work_real.size() * sizeof(T);
        size_t complex_mem = (G.size() + work_complex.size()) * sizeof(std::complex<T>);
        return real_mem + complex_mem;
    }
};


/**
 * @brief Template-based kernel wrapper for block evaluation
 * 
 * Requirements for KernelType:
 * - Must have a static or member function: evaluate_block(x_coords, x_size, y_coords, y_size, A, lda, ...)
 * - Additional kernel-specific parameters can be stored as member variables
 * 
 * Example usage:
 *   KernelWrapper<LaplaceKernel2D<double>> wrapper(kernel_instance);
 *   wrapper.evaluate_block(x, nx, y, ny, A, lda);
 */
template<typename KernelType, typename CoordType, typename DataType>
class KernelWrapper {
private:
    KernelType* kernel_instance;  // Non-owning pointer to actual kernel
    int dimension;
    
public:
    /**
     * @brief Constructor with kernel instance
     * @param kernel Pointer to kernel object (must outlive this wrapper)
     * @param dim Spatial dimension (2 or 3)
     */
    KernelWrapper(KernelType* kernel, int dim) 
        : kernel_instance(kernel), dimension(dim) {
        if (kernel == nullptr) {
            throw std::invalid_argument("Kernel instance cannot be null");
        }
        if (dim != 2 && dim != 3) {
            throw std::invalid_argument("Dimension must be 2 or 3");
        }
    }
    
    /**
     * @brief Evaluate kernel matrix block: A = K(X, Y)
     * 
     * @param x_coords Source point coordinates (column-major: dim × x_size)
     * @param x_size Number of source points
     * @param y_coords Target point coordinates (column-major: dim × y_size)
     * @param y_size Number of target points
     * @param A Output kernel matrix (column-major: x_size × y_size)
     * @param lda Leading dimension of A
     */
    void evaluate_block(
        const CoordType* x_coords, int64_t x_size,
        const CoordType* y_coords, int64_t y_size,
        DataType* A, int64_t lda) const {
        
        // Delegate to the actual kernel's evaluate_block
        kernel_instance->evaluate_block(x_coords, x_size, y_coords, y_size, A, lda);
    }
    
    /**
     * @brief Get spatial dimension
     */
    int get_dimension() const { return dimension; }
    
    /**
     * @brief Get pointer to underlying kernel (for specialized operations)
     */
    KernelType* get_kernel() { return kernel_instance; }
    const KernelType* get_kernel() const { return kernel_instance; }
};


/**
 * @brief 2D Helmholtz kernel implementation
 *
 * Implements the 2D Helmholtz fundamental solution:
 *   K(x,y) = (i/4) * H_0^(1)(k * ||x - y||)
 *
 * With optional Gaussian potential scaling:
 *   A_ij = V(x_i) * K(x_i, y_j) / N * V(y_j)
 *
 * And identity shift on diagonal:
 *   A_ii += 1
 *
 * The wavenumber k is derived from the grid as:
 *   k = 2π * √N / wave_number_divisor
 *
 * @tparam T Floating-point type (float or double)
 * @tparam CoordType Coordinate type (double or float)
 */
template<typename CoordType, typename T>
struct HelmholtzKernel2D {
private:
    int64_t N_total;
    T wavenumber;
    T wave_number_divisor;
    std::complex<T> precomputed_diag_val;
    bool use_potential;

    /**
     * @brief Wrapper around POSIX j0 (Bessel J_0)
     */
    static T bessel_j0(T x) {
        return static_cast<T>(::j0(static_cast<double>(x)));
    }

    /**
     * @brief Wrapper around POSIX y0 (Bessel Y_0)
     */
    static T bessel_y0(T x) {
        return static_cast<T>(::y0(static_cast<double>(x)));
    }

public:
    static constexpr int dimension = 2;
    static constexpr T pi = T(3.14159265358979323846);

    /**
     * @brief Construct Helmholtz kernel
     *
     * @param N Total number of points
     * @param divisor Wavenumber divisor (k = 2π√N / divisor). Default 32.
     * @param apply_potential Whether to apply Gaussian potential V(x). Default true.
     */
    HelmholtzKernel2D(int64_t N, T divisor = T(32.0), bool apply_potential = true)
        : N_total(N),
          wave_number_divisor(divisor),
          use_potential(apply_potential)
    {
        wavenumber = T(2.0) * pi * std::sqrt(static_cast<T>(N)) / wave_number_divisor;
        precomputed_diag_val = evaluate_diagonal(N);
    }

    /**
     * @brief Get the wavenumber k
     */
    T get_wavenumber() const { return wavenumber; }

    /**
     * @brief Evaluate the Helmholtz Green's function between two points
     *
     * K(x,y) = (i/4) * H_0^(1)(k * ||x - y||)
     */
    std::complex<T> evaluate(const CoordType* x, const CoordType* y) const {
        CoordType dx = x[0] - y[0];
        CoordType dy = x[1] - y[1];
        T r = std::sqrt(dx * dx + dy * dy);

        if (std::abs(r) < 1e-14) {
            return std::complex<T>(0, 0);
        }

        T kr = wavenumber * r;
        T j0_val = bessel_j0(kr);
        T y0_val = bessel_y0(kr);
        // (i/4) * H_0^(1)(kr) = (i/4) * (J0 + iY0) = -Y0/4 + iJ0/4
        return std::complex<T>(T(-0.25) * y0_val, T(0.25) * j0_val);
    }

    /**
     * @brief Evaluate the multiplicative weight applied on each side of the kernel
     *
     * V(x) = k * exp(-16 * ((x₁ - 0.5)² + (x₂ - 0.5)²))
     *
     * The dense and FFT operators both apply this same factor on the left and
     * right of the Green's function, i.e. A = I + diag(V) * C * diag(V).
     */
    T evaluate_potential(const CoordType* x) const {
        T dx = x[0] - T(0.5);
        T dy = x[1] - T(0.5);
        return wavenumber * std::exp(T(-16.0) * (dx * dx + dy * dy));
    }

    /**
     * @brief Evaluate diagonal self-interaction using 5-point Gauss-Legendre quadrature
     *
     * Computes: ∫∫_{[0,h/2]²} (i/4) * H_0^(1)(k * r) dA
     * Multiplied by 4 for symmetry over the full cell [-h/2, h/2]².
     *
     * Real and imaginary parts are integrated separately.
     *
     * @param N Total number of points
     * @return Complex diagonal value (not including potential or identity)
     */
    std::complex<T> evaluate_diagonal(int64_t N) const {
        int64_t n = static_cast<int64_t>(std::sqrt(static_cast<double>(N)));
        T h = T(1.0) / static_cast<T>(n);

        T k = wavenumber;

        // Re((i/4) * H_0^(1)(kr)) = -Y_0(kr)/4
        auto integrand_real = [k](T x, T y) -> T {
            T r = std::sqrt(x * x + y * y);
            if (std::abs(r) < 1e-14) return T(0.0);
            return T(-0.25) * bessel_y0(k * r);
        };

        // Im((i/4) * H_0^(1)(kr)) = J_0(kr)/4
        auto integrand_imag = [k](T x, T y) -> T {
            T r = std::sqrt(x * x + y * y);
            if (std::abs(r) < 1e-14) return T(0.0);
            return T(0.25) * bessel_j0(k * r);
        };

        // 5-point Gauss-Legendre quadrature nodes on [-1,1]
        const T nodes[] = {
            T(-0.9061798459386640), T(-0.5384693101056831), T(0.0),
            T(0.5384693101056831),  T(0.9061798459386640)
        };

        const T weights[] = {
            T(0.2369268850561891), T(0.4786286704993665), T(0.5688888888888889),
            T(0.4786286704993665), T(0.2369268850561891)
        };

        constexpr int npts = 5;

        // Transform from [-1,1] to [0, h/2]
        T b = h / T(2.0);
        T transform = b / T(2.0);
        T shift = b / T(2.0);

        T integral_real = T(0.0);
        T integral_imag = T(0.0);

        for (int i = 0; i < npts; ++i) {
            T xi = transform * nodes[i] + shift;
            for (int j = 0; j < npts; ++j) {
                T yj = transform * nodes[j] + shift;
                T w = weights[i] * weights[j];
                integral_real += w * integrand_real(xi, yj);
                integral_imag += w * integrand_imag(xi, yj);
            }
        }

        // Jacobian for coordinate transformation
        integral_real *= transform * transform;
        integral_imag *= transform * transform;

        // Factor of 4 for symmetry (integrated over 1/4 of cell)
        integral_real *= T(4.0);
        integral_imag *= T(4.0);

        return std::complex<T>(integral_real, integral_imag);
    }

    /**
     * @brief Evaluate a block of the kernel matrix (complex-valued)
     *
     * Fills A with:
     *   A[i,j] = V(x_i) * (i/4) * H_0^(1)(k * ||x_i - y_j||) / N * V(y_j)
     *
     * Diagonal entries (when x_coords == y_coords) are replaced with the
     * precomputed quadrature value, and the identity is added.
     *
     * @param x_coords  Points, stored as [x0, y0, x1, y1, ...] (interleaved 2D)
     * @param x_size    Number of x points
     * @param y_coords  Points, same layout
     * @param y_size    Number of y points
     * @param A         Output matrix (column-major), complex-valued
     * @param lda       Leading dimension of A
     */
    
    void evaluate_block(
        const CoordType *x_coords, int64_t x_size,
        const CoordType *y_coords, int64_t y_size,
        std::complex<T>* __restrict__ A, int64_t lda) const
    {
        const T inv_N = T(1.0) / static_cast<T>(N_total);
        const bool same_layout_ptr = (x_coords == y_coords) && (x_size == y_size);

        // 1. Deinterleave coordinates for better cache access patterns
        //    Instead of [x0,y0, x1,y1, ...], use separate x[] and y[] arrays
        std::vector<T> xi_x(x_size), xi_y(x_size);
        std::vector<T> yj_x(y_size), yj_y(y_size);

        for (int64_t i = 0; i < x_size; ++i) {
            xi_x[i] = static_cast<T>(x_coords[i * 2]);
            xi_y[i] = static_cast<T>(x_coords[i * 2 + 1]);
        }
        for (int64_t j = 0; j < y_size; ++j) {
            yj_x[j] = static_cast<T>(y_coords[j * 2]);
            yj_y[j] = static_cast<T>(y_coords[j * 2 + 1]);
        }

        // 2. Precompute potential values if needed
        std::vector<T> Vx, Vy;
        if (use_potential) {
            Vx.resize(x_size);
            Vy.resize(y_size);
            for (int64_t i = 0; i < x_size; ++i) {
                Vx[i] = evaluate_potential(&x_coords[i * 2]);
            }
            for (int64_t j = 0; j < y_size; ++j) {
                Vy[j] = evaluate_potential(&y_coords[j * 2]);
            }
        }

        // 3. Batch compute: separate distance computation from Bessel evaluation
        //    This allows the distance computation to vectorize even if Bessel can't
        // #pragma omp parallel for schedule(static) if(y_size > 16)
        for (int64_t j = 0; j < y_size; ++j) {
            const T yjx = yj_x[j];
            const T yjy = yj_y[j];
            const T vy = use_potential ? Vy[j] : T(1.0);

            // Compute all kr values for this column — this part can vectorize
            std::vector<T> kr_vals(x_size);

            // #pragma omp simd
            for (int64_t i = 0; i < x_size; ++i) {
                const T dx = xi_x[i] - yjx;
                const T dy = xi_y[i] - yjy;
                kr_vals[i] = wavenumber * std::sqrt(dx * dx + dy * dy);
            }

            // Evaluate Bessel functions in a tight scalar loop.
            for (int64_t i = 0; i < x_size; ++i) {
                const T kr = kr_vals[i];
                T re, im;

                if (kr < T(1e-14)) {
                    re = T(0.0);
                    im = T(0.0);
                } else {
                    // Direct calls avoid the static member function overhead
                    double kr_d = static_cast<double>(kr);
                    re = static_cast<T>(T(-0.25) * ::y0(kr_d));
                    im = static_cast<T>(T(0.25)  * ::j0(kr_d));
                }

                if (use_potential) {
                    const T pair_weight = Vx[i] * vy;
                    const T scale = inv_N * pair_weight;
                    A[i + j * lda] = std::complex<T>(re * scale, im * scale);
                } else {
                    A[i + j * lda] = std::complex<T>(re * inv_N, im * inv_N);
                }
            }

            for (int64_t i = 0; i < x_size; ++i) {
                const bool pointer_match = same_layout_ptr && (i == j);
                const bool coord_match =
                    (x_coords[i * 2] == yjx) &&
                    (x_coords[i * 2 + 1] == yjy);
                if (pointer_match || coord_match) {
                    std::complex<T> diag = precomputed_diag_val;
                    if (use_potential) {
                        const T pair_weight = Vx[i] * vy;
                        diag *= pair_weight;
                    }
                    A[i + j * lda] = diag + std::complex<T>(T(1.0), T(0.0));
                }
            }
        }
    }

    /**
     * @brief Evaluate block WITHOUT potential scaling or identity shift
     *
     * Pure Green's function evaluation:
     *   A[i,j] = (i/4) * H_0^(1)(k * ||x_i - y_j||) / N
     *
     * Useful when potential and identity are handled externally.
     */

    void evaluate_block_pure(
        const CoordType* x_coords, int64_t x_size,
        const CoordType* y_coords, int64_t y_size,
        std::complex<T>* __restrict__ A, int64_t lda) const
    {
        const T inv_N = T(1.0) / static_cast<T>(N_total);
        const bool same_layout_ptr = (x_coords == y_coords) && (x_size == y_size);

        for (int64_t j = 0; j < y_size; ++j) {
            const T yj_x = y_coords[j * 2];
            const T yj_y = y_coords[j * 2 + 1];

            for (int64_t i = 0; i < x_size; ++i) {
                const T xi_x = x_coords[i * 2];
                const T xi_y = x_coords[i * 2 + 1];

                const T dx = xi_x - yj_x;
                const T dy = xi_y - yj_y;
                const T r = std::sqrt(dx * dx + dy * dy);
                const T kr = wavenumber * r;

                std::complex<T> val;
                if (kr < T(1e-14)) {
                    val = std::complex<T>(0, 0);
                } else {
                    T j0_val = bessel_j0(kr);
                    T y0_val = bessel_y0(kr);
                    val = std::complex<T>(T(-0.25) * y0_val, T(0.25) * j0_val);
                }

                A[i + j * lda] = val * inv_N;
            }
        }

        for (int64_t j = 0; j < y_size; ++j) {
            const T yj_x = y_coords[j * 2];
            const T yj_y = y_coords[j * 2 + 1];
            for (int64_t i = 0; i < x_size; ++i) {
                const bool pointer_match = same_layout_ptr && (i == j);
                const bool coord_match =
                    (x_coords[i * 2] == yj_x) &&
                    (x_coords[i * 2 + 1] == yj_y);
                if (pointer_match || coord_match) {
                    A[i + j * lda] = precomputed_diag_val;
                }
            }
        }
    }
};

/**
 * @brief FFT-based fast matrix-vector multiplication for 2D Helmholtz kernel
 *
 * Uses complex-to-complex FFT to exploit Toeplitz structure on regular grids.
 * The Helmholtz kernel is complex-valued, so Hermitian FFT optimization
 * cannot be used (unlike Laplace).
 *
 * The matvec computes:
 *   y = x + V .* (C * (V .* x))
 *
 * where C is the Toeplitz kernel matrix (embedded as circulant) and
 * V(x) = k * exp(-16 * ((x₁ - 0.5)² + (x₂ - 0.5)²)).
 *
 * This matches the dense block assembly:
 *   A_ij = V(x_i) * K(x_i, y_j) / N * V(y_j) + δ_ij
 *
 * Complexity: O(N log N) vs O(N²) for direct evaluation
 *
 * @tparam T Floating-point type (double)
 */
template <typename CoordType, typename T>
class HelmholtzKernel2D_FFT {
private:
    int64_t n;           // Grid dimension (n × n = N)
    int64_t N;           // Total points
    int64_t padded_size; // 2n-1 (circulant embedding)

    T wavenumber;
    bool use_potential;

    std::vector<std::complex<T>> G;            // Precomputed FFT of circulant embedding
    mutable std::vector<std::complex<T>> work; // Scratch buffer for FFT (mutable for const matvec)
    std::vector<T> potential_weights;          // V(x) factor applied on both sides of C

    fftw_plan plan_forward;
    fftw_plan plan_backward;

public:
    static constexpr T pi = T(3.14159265358979323846);

    /**
     * @brief Constructor: precompute FFT of circulant embedding
     *
     * @param grid_points Point coordinates (interleaved: [x0,y0, x1,y1, ...])
     * @param grid_size Grid dimension n (must satisfy n² = N)
     * @param total_points Total number of points N
     * @param divisor Wavenumber divisor (k = 2π√N / divisor). Default 32.
     * @param apply_potential Whether to apply Gaussian potential. Default true.
     */
    HelmholtzKernel2D_FFT(const CoordType* grid_points, int64_t grid_size, int64_t total_points,
                           T divisor = T(32.0), bool apply_potential = true)
        : n(grid_size), N(total_points), padded_size(2 * grid_size - 1),
          use_potential(apply_potential),
          plan_forward(nullptr), plan_backward(nullptr)
    {
        if (n * n != N) {
            throw std::invalid_argument("Grid must be square: n² ≠ N");
        }

        wavenumber = T(2.0) * pi * std::sqrt(static_cast<T>(N)) / divisor;

        // Allocate complex buffers: full (2n-1)² since kernel is complex
        G.resize(padded_size * padded_size);
        work.resize(padded_size * padded_size);

        // Precompute the multiplicative weight used on both sides of the kernel.
        if (use_potential) {
            potential_weights.resize(N);
            for (int64_t i = 0; i < N; ++i) {
                T pt[2] = {grid_points[i * 2], grid_points[i * 2 + 1]};
                T dx = pt[0] - T(0.5);
                T dy = pt[1] - T(0.5);
                potential_weights[i] = wavenumber * std::exp(T(-16.0) * (dx * dx + dy * dy));
            }
        }

        // Step 1: Compute first column of kernel matrix
        // K(x_i, x_0) / N for all i, with diagonal replaced by quadrature value
        std::vector<std::complex<T>> first_col(N);
        T y_fixed[2] = {grid_points[0], grid_points[1]};

        std::complex<T> diag = HelmholtzKernel2D<CoordType, T>(N, divisor, false).evaluate_diagonal(N);

        for (int64_t i = 0; i < N; ++i) {
            if (i == 0) {
                first_col[i] = diag;
            } else {
                T x[2] = {grid_points[i * 2], grid_points[i * 2 + 1]};
                T dx = x[0] - y_fixed[0];
                T dy = x[1] - y_fixed[1];
                T r = std::sqrt(dx * dx + dy * dy);
                T kr = wavenumber * r;

                // (i/4) * H_0^(1)(kr) / N
                T j0_val = static_cast<T>(::j0(static_cast<double>(kr)));
                T y0_val = static_cast<T>(::y0(static_cast<double>(kr)));
                first_col[i] = std::complex<T>(T(-0.25) * y0_val, T(0.25) * j0_val)
                               / static_cast<T>(N);
            }
        }

        // Step 2: Reshape first_col to n × n (column-major)
        std::vector<std::complex<T>> A(n * n);
        for (int64_t i = 0; i < n; ++i) {
            for (int64_t j = 0; j < n; ++j) {
                A[i + j * n] = first_col[i * n + j];
            }
        }

        // Step 3: Build circulant embedding B in G buffer
        // Mirrors Julia:
        //   B[1:n, 1:n] = a
        //   B[1:n, n+1:end] = a[:, 2:n]
        //   B[n+1:end, 1:n] = a[2:n, :]
        //   B[n+1:end, n+1:end] = a[2:n, 2:n]
        //   B[n+1:end, :] = reverse(B[n+1:end, :], dims=1)
        //   B[:, n+1:end] = reverse(B[:, n+1:end], dims=2)
        std::fill(G.begin(), G.end(), std::complex<T>(0, 0));

        // B(1:n, 1:n) = A
        for (int64_t j = 0; j < n; ++j) {
            for (int64_t i = 0; i < n; ++i) {
                G[i + j * padded_size] = A[i + j * n];
            }
        }

        // B(1:n, n+1:2n-1) = A(:, 2:n) reversed → column k maps to 2n-1-k
        for (int64_t j = 1; j < n; ++j) {
            for (int64_t i = 0; i < n; ++i) {
                G[i + (padded_size - j) * padded_size] = A[i + j * n];
            }
        }

        // B(n+1:2n-1, 1:n) = A(2:n, :) reversed → row k maps to 2n-1-k
        for (int64_t i = 1; i < n; ++i) {
            for (int64_t j = 0; j < n; ++j) {
                G[(padded_size - i) + j * padded_size] = A[i + j * n];
            }
        }

        // B(n+1:2n-1, n+1:2n-1) = A(2:n, 2:n) reversed both
        for (int64_t i = 1; i < n; ++i) {
            for (int64_t j = 1; j < n; ++j) {
                G[(padded_size - i) + (padded_size - j) * padded_size] = A[i + j * n];
            }
        }

        // Step 4: Compute G = fft2(B) in-place
        if constexpr (std::is_same_v<T, double>) {
            fftw_plan temp_plan = fftw_plan_dft_2d(
                padded_size, padded_size,
                reinterpret_cast<fftw_complex*>(G.data()),
                reinterpret_cast<fftw_complex*>(G.data()),
                FFTW_FORWARD, FFTW_ESTIMATE
            );
            fftw_execute(temp_plan);
            fftw_destroy_plan(temp_plan);

            // Create reusable plans for matvec (operate on work buffer)
            plan_forward = fftw_plan_dft_2d(
                padded_size, padded_size,
                reinterpret_cast<fftw_complex*>(work.data()),
                reinterpret_cast<fftw_complex*>(work.data()),
                FFTW_FORWARD, FFTW_ESTIMATE
            );

            plan_backward = fftw_plan_dft_2d(
                padded_size, padded_size,
                reinterpret_cast<fftw_complex*>(work.data()),
                reinterpret_cast<fftw_complex*>(work.data()),
                FFTW_BACKWARD, FFTW_ESTIMATE
            );
        } else {
            throw std::runtime_error("Only double precision currently supported for FFT");
        }
    }

    ~HelmholtzKernel2D_FFT() {
        if constexpr (std::is_same_v<T, double>) {
            if (plan_forward) fftw_destroy_plan(plan_forward);
            if (plan_backward) fftw_destroy_plan(plan_backward);
        }
    }

    // Non-copyable due to FFTW plans
    HelmholtzKernel2D_FFT(const HelmholtzKernel2D_FFT&) = delete;
    HelmholtzKernel2D_FFT& operator=(const HelmholtzKernel2D_FFT&) = delete;

    // Movable
    HelmholtzKernel2D_FFT(HelmholtzKernel2D_FFT&& other) noexcept
        : n(other.n), N(other.N), padded_size(other.padded_size),
          wavenumber(other.wavenumber), use_potential(other.use_potential),
          G(std::move(other.G)), work(std::move(other.work)),
          potential_weights(std::move(other.potential_weights)),
          plan_forward(other.plan_forward), plan_backward(other.plan_backward)
    {
        other.plan_forward = nullptr;
        other.plan_backward = nullptr;
    }

    HelmholtzKernel2D_FFT& operator=(HelmholtzKernel2D_FFT&& other) noexcept {
        if (this != &other) {
            if constexpr (std::is_same_v<T, double>) {
                if (plan_forward) fftw_destroy_plan(plan_forward);
                if (plan_backward) fftw_destroy_plan(plan_backward);
            }
            n = other.n;
            N = other.N;
            padded_size = other.padded_size;
            wavenumber = other.wavenumber;
            use_potential = other.use_potential;
            G = std::move(other.G);
            work = std::move(other.work);
            potential_weights = std::move(other.potential_weights);
            plan_forward = other.plan_forward;
            plan_backward = other.plan_backward;
            other.plan_forward = nullptr;
            other.plan_backward = nullptr;
        }
        return *this;
    }

    /**
     * @brief Fast matrix-vector multiplication: y = x + V .* (C * (V .* x))
     *
     * Mirrors the Julia fft_mv_:
     *   larger[1:n, 1:n] = reshape(weight .* x, n, n)
     *   fft!(larger)
     *   larger .= F .* larger
     *   ifft!(larger)
     *   y = weight .* reshape(larger[1:n, 1:n], N) + x
     *
     * @param x Input vector (size N), complex-valued
     * @param y Output vector (size N, pre-allocated), complex-valued
     */
    void matvec(const std::complex<T>* x, std::complex<T>* y) const {
        // Step 1: Apply the left/right weight and zero-pad to (2n-1)².
        std::fill(work.begin(), work.end(), std::complex<T>(0, 0));

        if (use_potential) {
            for (int64_t i = 0; i < n; ++i) {
                for (int64_t j = 0; j < n; ++j) {
                    int64_t idx = i * n + j;
                    work[i + j * padded_size] = potential_weights[idx] * x[idx];
                }
            }
        } else {
            for (int64_t i = 0; i < n; ++i) {
                for (int64_t j = 0; j < n; ++j) {
                    int64_t idx = i * n + j;
                    work[i + j * padded_size] = x[idx];
                }
            }
        }

        // Step 2: Forward FFT (use new-array execute for safety)
        if constexpr (std::is_same_v<T, double>) {
            fftw_execute_dft(plan_forward,
                reinterpret_cast<fftw_complex*>(work.data()),
                reinterpret_cast<fftw_complex*>(work.data()));
        }

        // Step 3: Element-wise multiplication with precomputed G
        int64_t total_padded = padded_size * padded_size;
        for (int64_t i = 0; i < total_padded; ++i) {
            work[i] *= G[i];
        }

        // Step 4: Inverse FFT
        if constexpr (std::is_same_v<T, double>) {
            fftw_execute_dft(plan_backward,
                reinterpret_cast<fftw_complex*>(work.data()),
                reinterpret_cast<fftw_complex*>(work.data()));
        }

        // Step 5: Extract first n × n, normalize, apply the weight again, add identity.
        // FFTW does not normalize inverse FFT, so divide by (2n-1)²
        T scale = T(1.0) / static_cast<T>(total_padded);

        if (use_potential) {
            for (int64_t i = 0; i < n; ++i) {
                for (int64_t j = 0; j < n; ++j) {
                    int64_t idx = i * n + j;
                    y[idx] = potential_weights[idx] * work[i + j * padded_size] * scale + x[idx];
                }
            }
        } else {
            for (int64_t i = 0; i < n; ++i) {
                for (int64_t j = 0; j < n; ++j) {
                    int64_t idx = i * n + j;
                    y[idx] = work[i + j * padded_size] * scale + x[idx];
                }
            }
        }
    }

    /**
     * @brief Convenience overload for real-valued input
     *
     * @param x_real Input vector (size N), real-valued
     * @param y Output vector (size N, pre-allocated), complex-valued
     */
    void matvec(const T* x_real, std::complex<T>* y) const {
        std::vector<std::complex<T>> x_complex(N);
        for (int64_t i = 0; i < N; ++i) {
            x_complex[i] = std::complex<T>(x_real[i], T(0));
        }
        matvec(x_complex.data(), y);
    }

    /**
     * @brief Get the wavenumber k
     */
    T get_wavenumber() const { return wavenumber; }

    /**
     * @brief Get memory usage in bytes
     */
    size_t memory_usage() const {
        size_t complex_mem = (G.size() + work.size()) * sizeof(std::complex<T>);
        size_t potential_mem = potential_weights.size() * sizeof(T);
        return complex_mem + potential_mem;
    }
};


/**
 * @brief 3D Helmholtz kernel implementation
 *
 * Implements the 3D Helmholtz fundamental solution:
 *   K(x,y) = exp(i * k * ||x - y||) / (4π * ||x - y||)
 *
 * With optional Gaussian potential scaling:
 *   A_ij = V(x_i) * K(x_i, y_j) / N * V(y_j)
 *
 * And identity shift on diagonal:
 *   A_ii += 1
 *
 * The wavenumber k is derived from the grid as:
 *   k = 2π * N^(1/3) / wave_number_divisor
 *
 * @tparam T Floating-point type (float or double)
 * @tparam CoordType Coordinate type (double or float)
 */
template<typename CoordType, typename T>
struct HelmholtzKernel3D {
private:
    int64_t N_total;
    T wavenumber;
    T wave_number_divisor;
    std::complex<T> precomputed_diag_val;
    // s1 / inverse_coef storage:
    //   * If s1 is constant across all points, only the scalar fields are
    //     populated and s1_values / inverse_coef_values stay empty.
    //   * If s1 varies, scalar fields are unused and the vectors hold the
    //     per-point values (same layout as before).
    bool s1_is_constant;
    T s1_constant;
    T inverse_coef_constant;
    std::vector<T> s1_values;
    std::vector<T> inverse_coef_values;
    bool use_potential;

    int64_t grid_points_per_dim() const {
        return static_cast<int64_t>(
            std::llround(std::cbrt(static_cast<long double>(N_total))));
    }

    int64_t global_index_from_coords(T x, T y, T z) const {
        const int64_t n = grid_points_per_dim();
        const auto clamp_index = [n](int64_t idx) {
            return std::max<int64_t>(0, std::min<int64_t>(n - 1, idx));
        };

        // Leaf points live at cell centers x = (ix + 0.5) / n on the unit cube.
        const int64_t ix = clamp_index(static_cast<int64_t>(std::llround(x * n - T(0.5))));
        const int64_t iy = clamp_index(static_cast<int64_t>(std::llround(y * n - T(0.5))));
        const int64_t iz = clamp_index(static_cast<int64_t>(std::llround(z * n - T(0.5))));

        // Match the tree's global ordering: (ix * n + iy) * n + iz.
        return (ix * n + iy) * n + iz;
    }

public:
    static constexpr int dimension = 3;
    static constexpr T pi = T(3.14159265358979323846);

    /**
     * @brief Construct 3D Helmholtz kernel
     *
     * @param N Total number of points
     * @param divisor Wavenumber divisor (k = 2π * N^(1/3) / divisor). Default 32.
     * @param apply_potential Whether to apply Gaussian potential V(x). Default true.
     * @param s1_input Length-N vector of slowness values. If empty, defaults to 3 everywhere.
     */
    HelmholtzKernel3D(int64_t N,
                    T divisor = T(32.0),
                    bool apply_potential = true,
                    const std::vector<T>& s1_input = {})
        : N_total(N),
        wave_number_divisor(divisor),
        s1_is_constant(true),
        s1_constant(T(3)),
        inverse_coef_constant(T(0)),
        use_potential(apply_potential)
    {
        int64_t n = std::cbrt(static_cast<T>(N));
        assert(n > 0 && static_cast<long double>(n)*n*n == static_cast<long double>(N) && "N must be a perfect cube (N = n^3).");
        wavenumber = T(2.0) * pi * n / wave_number_divisor;

        if (s1_input.empty()) {
            s1_is_constant = true;
            s1_constant = T(3);
        } else {
            if (static_cast<int64_t>(s1_input.size()) != N_total) {
                throw std::invalid_argument("s1_input must have length N");
            }
            s1_is_constant = true;
            s1_constant = s1_input[0];
            for (int64_t i = 1; i < N_total; ++i) {
                if (s1_input[i] != s1_constant) { s1_is_constant = false; break; }
            }
            if (!s1_is_constant) s1_values = s1_input;
        }

        precomputed_diag_val = evaluate_diagonal(N);

        const T s0 = T(2.0);
        const T k0 = wavenumber;
        auto coef_of = [&](T s1v) {
            return (k0 * k0) * ((s1v / s0) * (s1v / s0) - T(1));
        };
        if (s1_is_constant) {
            inverse_coef_constant = T(1) / coef_of(s1_constant);
        } else {
            inverse_coef_values.resize(N_total);
            for (int64_t i = 0; i < N_total; ++i) {
                inverse_coef_values[i] = T(1) / coef_of(s1_values[i]);
            }
        }
    }

    /**
     * @brief Get the wavenumber k
     */
    T get_wavenumber() const { return wavenumber; }

    /**
     * @brief Evaluate the 3D Helmholtz Green's function between two points
     *
     * K(x,y) = exp(i * k * r) / (4π * r),  r = ||x - y||
     */
    std::complex<T> evaluate(const CoordType* x, const CoordType* y) const {
        CoordType dx = x[0] - y[0];
        CoordType dy = x[1] - y[1];
        CoordType dz = x[2] - y[2];
        T r = std::sqrt(static_cast<T>(dx * dx + dy * dy + dz * dz));

        if (r < T(1e-14)) {
            return std::complex<T>(0, 0);
        }

        T kr = wavenumber * r;
        T inv_4pi_r = T(1.0) / (T(4.0) * pi * r);
        // exp(ikr) / (4πr) = (cos(kr) + i*sin(kr)) / (4πr)
        return std::complex<T>(std::cos(kr) * inv_4pi_r, std::sin(kr) * inv_4pi_r);
    }

    /**
     * @brief Evaluate Gaussian potential at a point in 3D
     *
     * V(x) = k * exp(-16 * ((x₁ - 0.5)² + (x₂ - 0.5)² + (x₃ - 0.5)²))
     */
    T evaluate_potential(const T* x) const {
        T dx = x[0] - T(0.5);
        T dy = x[1] - T(0.5);
        T dz = x[2] - T(0.5);
        return wavenumber * std::exp(T(-16.0) * (dx * dx + dy * dy + dz * dz));
    }

    /**
     * @brief Evaluate diagonal self-interaction using 5-point Gauss-Legendre quadrature
     *
     * Computes: ∫∫∫_{[0,h/2]³} exp(ikr) / (4πr) dV
     * Multiplied by 8 for symmetry over the full cell [-h/2, h/2]³.
     *
     * Real and imaginary parts are integrated separately.
     *
     * @param N Total number of points
     * @return Complex diagonal value (not including potential or identity)
     */
    std::complex<T> evaluate_diagonal(int64_t N) {
        // Grid spacing from your original convention:
        // N ~ total grid points, n ~ points per dimension, h = 1/n
        T s0 = T(2.0);
        int64_t n = static_cast<int64_t>(std::cbrt(static_cast<double>(N)));
        T h = T(1) / static_cast<T>(n);

        // If your member 'wavenumber' is already k = omega*s0, then:
        T k = wavenumber;

        // If instead wavenumber is omega and you want k = omega*s0, use:
        // T k = wavenumber * s0;

        T tt = h * k;                 // dimensionless
        constexpr T pi = T(3.1415926535897932384626433832795L);
        const std::complex<T> I(0, 1);

        // 7-pt Gauss nodes/weights already mapped to [0, pi/4] and [0, 1]
        const T nodes_phi[7] = {
            T(0.765412887308718), T(0.683897697334573), T(0.552074099956508),
            T(0.392699081698724), T(0.233324063440941), T(0.101500466062876),
            T(0.019985276088730)
        };
        const T weights_phi[7] = {
            T(0.050848627308305), T(0.109840050384021), T(0.149944310198338),
            T(0.164132187616120), T(0.149944310198338), T(0.109840050384021),
            T(0.050848627308305)
        };

        const T nodes_rhop[7] = {
            T(0.974553956171379), T(0.870765592799697), T(0.702922575688699),
            T(0.500000000000000), T(0.297077424311301), T(0.129234407200303),
            T(0.025446043828621)
        };
        const T weights_rhop[7] = {
            T(0.064742483084435), T(0.139852695744638), T(0.190915025252560),
            T(0.208979591836735), T(0.190915025252560), T(0.139852695744638),
            T(0.064742483084435)
        };

        std::complex<T> sum = 0;

        for (int i = 0; i < 7; ++i) {
            T rhop = nodes_rhop[i];
            for (int j = 0; j < 7; ++j) {
                T phi = nodes_phi[j];

                T c = std::cos(phi);
                T rho = rhop / c;

                T s = std::sqrt(T(1) + rho * rho);
                T r1 = (tt / T(2)) * s;

                // fun = ( (e^{i r1}(1 - i r1) - 1) * rho/(1+rho^2)^{3/2} ) / cos(phi)
                std::complex<T> e = std::exp(I * r1);
                std::complex<T> term = e * (T(1) - I * r1) - T(1);

                T geom = (rho / std::pow(T(1) + rho * rho, T(1.5))) / c;

                std::complex<T> fun = term * geom;

                sum += weights_rhop[i] * weights_phi[j] * fun;
            }
        }

        // Their scaling: *val *= -48 i; then *= i/(4π tt^2 h)
        // Combine to +48/(4π tt^2 h) = 12/(π tt^2 h)
        std::complex<T> val = sum * (T(12) / (pi * tt * tt * h));
        std::complex<T> Aii = val * (h*h*h);

        return Aii;
    }

    void evaluate_block(
        const CoordType *x_coords, int64_t x_size,
        const CoordType *y_coords, int64_t y_size,
        std::complex<T>* __restrict__ A, int64_t lda) const
    {
        const T inv_N = T(1.0) / static_cast<T>(N_total);
        const bool same_layout_ptr = (x_coords == y_coords) && (x_size == y_size);
        const T inv_4pi = T(1.0) / (T(4.0) * pi);

        // Deinterleave coordinates for contiguous SIMD access
        std::vector<T> xi_x(x_size), xi_y(x_size), xi_z(x_size);
        std::vector<T> yj_x(y_size), yj_y(y_size), yj_z(y_size);

        for (int64_t i = 0; i < x_size; ++i) {
            xi_x[i] = static_cast<T>(x_coords[i * 3]);
            xi_y[i] = static_cast<T>(x_coords[i * 3 + 1]);
            xi_z[i] = static_cast<T>(x_coords[i * 3 + 2]);
        }
        for (int64_t j = 0; j < y_size; ++j) {
            yj_x[j] = static_cast<T>(y_coords[j * 3]);
            yj_y[j] = static_cast<T>(y_coords[j * 3 + 1]);
            yj_z[j] = static_cast<T>(y_coords[j * 3 + 2]);
        }

        // Precompute potential values
        std::vector<T> Vx, Vy;
        if (use_potential) {
            Vx.resize(x_size);
            Vy.resize(y_size);
            for (int64_t i = 0; i < x_size; ++i) {
                T pt[3] = {xi_x[i], xi_y[i], xi_z[i]};
                Vx[i] = evaluate_potential(pt);
            }
            for (int64_t j = 0; j < y_size; ++j) {
                T pt[3] = {yj_x[j], yj_y[j], yj_z[j]};
                Vy[j] = evaluate_potential(pt);
            }
        }

        // Temporary contiguous buffers for real and imaginary parts
        // (avoids interleaved complex stores which break SIMD)
        std::vector<T> col_real(x_size);
        std::vector<T> col_imag(x_size);

        for (int64_t j = 0; j < y_size; ++j) {
            const T yjx = yj_x[j];
            const T yjy = yj_y[j];
            const T yjz = yj_z[j];

            // Phase 1: Compute kr values — vectorizes cleanly
            std::vector<T> kr(x_size);

            // #pragma omp simd
            for (int64_t i = 0; i < x_size; ++i) {
                T dx = xi_x[i] - yjx;
                T dy = xi_y[i] - yjy;
                T dz = xi_z[i] - yjz;
                kr[i] = wavenumber * std::sqrt(dx * dx + dy * dy + dz * dz);
            }

            // Phase 2: Compute sin/cos — vectorizable with appropriate flags
            std::vector<T> ck(x_size), sk(x_size);

            // #pragma omp simd
            for (int64_t i = 0; i < x_size; ++i) {
                ck[i] = std::cos(kr[i]);
                sk[i] = std::sin(kr[i]);
            }

            // Phase 3: Assemble kernel values — vectorizes cleanly
            const T vy = use_potential ? Vy[j] : T(1.0);
            if (use_potential) {
                // #pragma omp simd
                for (int64_t i = 0; i < x_size; ++i) {
                    T r = kr[i] / wavenumber;
                    const T pair_weight = Vx[i] * vy;
                    const T common_scale =
                        (r < T(1e-14)) ? T(0.0) : inv_4pi / r * inv_N;
                    const T scale = common_scale * pair_weight;
                    col_real[i] = ck[i] * scale;
                    col_imag[i] = sk[i] * scale;
                }
            } else {
                // #pragma omp simd
                for (int64_t i = 0; i < x_size; ++i) {
                    T r = kr[i] / wavenumber;
                    T scale = (r < T(1e-14)) ? T(0.0) : inv_4pi / r * inv_N;
                    col_real[i] = ck[i] * scale;
                    col_imag[i] = sk[i] * scale;
                }
            }

            // Phase 4: Write to complex output (interleaved store, not SIMD-friendly
            // but only one pass over the data)
            for (int64_t i = 0; i < x_size; ++i) {
                A[i + j * lda] = std::complex<T>(col_real[i], col_imag[i]);
            }

            for (int64_t i = 0; i < x_size; ++i) {
                const bool pointer_match = same_layout_ptr && (i == j);
                const bool coord_match =
                    (x_coords[i * 3] == yjx) &&
                    (x_coords[i * 3 + 1] == yjy) &&
                    (x_coords[i * 3 + 2] == yjz);
                if (pointer_match || coord_match) {
                    std::complex<T> diag = precomputed_diag_val;
                    if (use_potential) {
                        const T pair_weight = Vx[i] * vy;
                        diag *= pair_weight;
                    }
                    const T ic = s1_is_constant
                        ? inverse_coef_constant
                        : inverse_coef_values[
                              global_index_from_coords(xi_x[i], xi_y[i], xi_z[i])];
                    A[i + j * lda] = diag - std::complex<T>(ic, T(0));
                }
            }
        }
    }

    void evaluate_block_no_diag(
        const CoordType* x_coords, int64_t x_size,
        const CoordType* y_coords, int64_t y_size,
        std::complex<T>* __restrict__ A, int64_t lda) const
    {
        const T inv_N = T(1.0) / static_cast<T>(N_total);

        for (int64_t j = 0; j < y_size; ++j) {
            const T yj_x = static_cast<T>(y_coords[j * 3]);
            const T yj_y = static_cast<T>(y_coords[j * 3 + 1]);
            const T yj_z = static_cast<T>(y_coords[j * 3 + 2]);

            for (int64_t i = 0; i < x_size; ++i) {
                const T xi_x = static_cast<T>(x_coords[i * 3]);
                const T xi_y = static_cast<T>(x_coords[i * 3 + 1]);
                const T xi_z = static_cast<T>(x_coords[i * 3 + 2]);

                const T dx = xi_x - yj_x;
                const T dy = xi_y - yj_y;
                const T dz = xi_z - yj_z;
                const T r = std::sqrt(dx * dx + dy * dy + dz * dz);

                std::complex<T> val;
                if (r < T(1e-14)) {
                    val = std::complex<T>(0, 0);
                } else {
                    T kr = wavenumber * r;
                    T inv_4pi_r = T(1.0) / (T(4.0) * pi * r);
                    val = std::complex<T>(std::cos(kr) * inv_4pi_r,
                                          std::sin(kr) * inv_4pi_r);
                }

                A[i + j * lda] = val * inv_N;
            }
        }
    }

    const std::vector<T>& get_s1_values() const { return s1_values; }
    const std::vector<T>& get_inverse_coef_values() const { return inverse_coef_values; }
};


template <typename CoordType, typename T>
class HelmholtzKernel3D_FFT {
private:
    int64_t n;
    int64_t N;
    int64_t padded;
    int64_t padded3;

    T wavenumber;
    bool use_potential;

    // Adds the extra diagonal shift:  - inverse_coef * x
    // (set from coef = k0^2 * ((s1/s0)^2 - 1), now varying pointwise with s1[i])
    // When s1 is constant across all points, only the scalar fields are used
    // and the length-N vectors stay empty.
    bool s1_is_constant;
    T s1_constant;
    T inverse_coef_constant;
    std::vector<T> s1_values;
    std::vector<T> inverse_coef_values;

    std::vector<std::complex<T>> G;
    mutable std::vector<std::complex<T>> work;
    std::vector<T> sqrt_potential;

    fftw_plan plan_forward;
    fftw_plan plan_backward;

public:
    static constexpr T pi = T(3.14159265358979323846);

    HelmholtzKernel3D_FFT(const CoordType* grid_points,
                        int64_t grid_size,
                        int64_t total_points,
                        T divisor = T(32.0),
                        bool apply_potential = true,
                        const std::vector<T>& s1_input = {})
        : n(grid_size),
        N(total_points),
        padded(2 * grid_size - 1),
        padded3(0),
        wavenumber(T(0)),
        use_potential(apply_potential),
        s1_is_constant(true),
        s1_constant(T(3)),
        inverse_coef_constant(T(0)),
        plan_forward(nullptr),
        plan_backward(nullptr)
    {
        if (n * n * n != N) {
            throw std::invalid_argument("Grid must be cubic: n³ ≠ N");
        }

        padded3 = padded * padded * padded;

        // Here wavenumber is k0 (spatial wavenumber)
        wavenumber = T(2.0) * pi * std::cbrt(static_cast<T>(N)) / divisor;

        if (s1_input.empty()) {
            s1_is_constant = true;
            s1_constant = T(3);
        } else {
            if (static_cast<int64_t>(s1_input.size()) != N) {
                throw std::invalid_argument("s1_input must have length N");
            }
            s1_is_constant = true;
            s1_constant = s1_input[0];
            for (int64_t i = 1; i < N; ++i) {
                if (s1_input[i] != s1_constant) { s1_is_constant = false; break; }
            }
            if (!s1_is_constant) s1_values = s1_input;
        }

        // ---- compute pointwise inverse_coef (CHANGE s0,s1 HERE as needed) ----
        const T s0 = T(2);
        auto coef_of = [&](T s1v) {
            return (wavenumber * wavenumber) * ((s1v / s0) * (s1v / s0) - T(1));
        };
        if (s1_is_constant) {
            T c = coef_of(s1_constant);
            if (std::abs(c) == T(0))
                throw std::runtime_error("coef is zero; cannot form inverse_coef");
            inverse_coef_constant = T(1) / c;
        } else {
            inverse_coef_values.resize(N);
            for (int64_t i = 0; i < N; ++i) {
                T c = coef_of(s1_values[i]);
                if (std::abs(c) == T(0))
                    throw std::runtime_error("coef is zero; cannot form inverse_coef");
                inverse_coef_values[i] = T(1) / c;
            }
        }
        // -------------------------------------------------------------------

        G.resize(padded3);
        work.resize(padded3);

        // Precompute potential values (sqrtV)
        if (use_potential) {
            sqrt_potential.resize(N);
            for (int64_t i = 0; i < N; ++i) {
                T dx = static_cast<T>(grid_points[i * 3])     - T(0.5);
                T dy = static_cast<T>(grid_points[i * 3 + 1]) - T(0.5);
                T dz = static_cast<T>(grid_points[i * 3 + 2]) - T(0.5);
                sqrt_potential[i] = wavenumber *
                    std::exp(T(-16.0) * (dx * dx + dy * dy + dz * dz));
            }
        }

        // Step 1: Compute first column of kernel matrix (WITHOUT potential)
        std::vector<std::complex<T>> first_col(N);

        // IMPORTANT: evaluate_diagonal(N) must return ONLY the kernel self-cell integral,
        // not including identity or -1/coef. The -1/coef is applied in matvec() below.
        std::complex<T> diag =
            HelmholtzKernel3D<CoordType, T>(N, divisor, false, s1_values).evaluate_diagonal(N);

        T y0[3] = {
            static_cast<T>(grid_points[0]),
            static_cast<T>(grid_points[1]),
            static_cast<T>(grid_points[2])
        };

        for (int64_t i = 0; i < N; ++i) {
            if (i == 0) {
                first_col[i] = diag;
            } else {
                T x[3] = {
                    static_cast<T>(grid_points[i * 3]),
                    static_cast<T>(grid_points[i * 3 + 1]),
                    static_cast<T>(grid_points[i * 3 + 2])
                };
                T dx = x[0] - y0[0];
                T dy = x[1] - y0[1];
                T dz = x[2] - y0[2];
                T r = std::sqrt(dx * dx + dy * dy + dz * dz);

                T kr = wavenumber * r;
                T inv_4pi_r = T(1.0) / (T(4.0) * pi * r);

                // Multiply by h^3 = 1/N
                first_col[i] = std::complex<T>(std::cos(kr) * inv_4pi_r,
                                            std::sin(kr) * inv_4pi_r)
                            / static_cast<T>(N);
            }
        }

        // Step 2: Place first_col into top-left n×n×n corner of padded buffer
        std::fill(G.begin(), G.end(), std::complex<T>(0, 0));
        for (int64_t k = 0; k < n; ++k) {
            for (int64_t j = 0; j < n; ++j) {
                for (int64_t i = 0; i < n; ++i) {
                    G[i + j * padded + k * padded * padded] =
                        first_col[i + j * n + k * n * n];
                }
            }
        }

        // Step 3: Build circulant embedding by sequential mirroring
        for (int64_t k = 0; k < n; ++k) {
            for (int64_t j = 0; j < n; ++j) {
                for (int64_t i = 1; i < n; ++i) {
                    G[(padded - i) + j * padded + k * padded * padded] =
                        G[i + j * padded + k * padded * padded];
                }
            }
        }

        for (int64_t k = 0; k < n; ++k) {
            for (int64_t j = 1; j < n; ++j) {
                for (int64_t i = 0; i < padded; ++i) {
                    G[i + (padded - j) * padded + k * padded * padded] =
                        G[i + j * padded + k * padded * padded];
                }
            }
        }

        for (int64_t k = 1; k < n; ++k) {
            for (int64_t j = 0; j < padded; ++j) {
                for (int64_t i = 0; i < padded; ++i) {
                    G[i + j * padded + (padded - k) * padded * padded] =
                        G[i + j * padded + k * padded * padded];
                }
            }
        }

        // Step 4: FFT of circulant embedding
        if constexpr (std::is_same_v<T, double>) {
            fftw_plan temp_plan = fftw_plan_dft_3d(
                static_cast<int>(padded), static_cast<int>(padded), static_cast<int>(padded),
                reinterpret_cast<fftw_complex*>(G.data()),
                reinterpret_cast<fftw_complex*>(G.data()),
                FFTW_FORWARD, FFTW_ESTIMATE
            );
            fftw_execute(temp_plan);
            fftw_destroy_plan(temp_plan);

            plan_forward = fftw_plan_dft_3d(
                static_cast<int>(padded), static_cast<int>(padded), static_cast<int>(padded),
                reinterpret_cast<fftw_complex*>(work.data()),
                reinterpret_cast<fftw_complex*>(work.data()),
                FFTW_FORWARD, FFTW_ESTIMATE
            );

            plan_backward = fftw_plan_dft_3d(
                static_cast<int>(padded), static_cast<int>(padded), static_cast<int>(padded),
                reinterpret_cast<fftw_complex*>(work.data()),
                reinterpret_cast<fftw_complex*>(work.data()),
                FFTW_BACKWARD, FFTW_ESTIMATE
            );
        } else {
            throw std::runtime_error("Only double precision currently supported for FFT");
        }
    }

    ~HelmholtzKernel3D_FFT() {
        if constexpr (std::is_same_v<T, double>) {
            if (plan_forward) fftw_destroy_plan(plan_forward);
            if (plan_backward) fftw_destroy_plan(plan_backward);
        }
    }

    HelmholtzKernel3D_FFT(const HelmholtzKernel3D_FFT&) = delete;
    HelmholtzKernel3D_FFT& operator=(const HelmholtzKernel3D_FFT&) = delete;

    HelmholtzKernel3D_FFT(HelmholtzKernel3D_FFT&& other) noexcept
        : n(other.n), N(other.N), padded(other.padded), padded3(other.padded3),
        wavenumber(other.wavenumber), use_potential(other.use_potential),
        s1_is_constant(other.s1_is_constant),
        s1_constant(other.s1_constant),
        inverse_coef_constant(other.inverse_coef_constant),
        s1_values(std::move(other.s1_values)),
        inverse_coef_values(std::move(other.inverse_coef_values)),
        G(std::move(other.G)), work(std::move(other.work)),
        sqrt_potential(std::move(other.sqrt_potential)),
        plan_forward(other.plan_forward), plan_backward(other.plan_backward)
    {
        other.plan_forward = nullptr;
        other.plan_backward = nullptr;
    }

    HelmholtzKernel3D_FFT& operator=(HelmholtzKernel3D_FFT&& other) noexcept {
        if (this != &other) {
            if constexpr (std::is_same_v<T, double>) {
                if (plan_forward) fftw_destroy_plan(plan_forward);
                if (plan_backward) fftw_destroy_plan(plan_backward);
            }
            n = other.n;
            N = other.N;
            padded = other.padded;
            padded3 = other.padded3;
            wavenumber = other.wavenumber;
            use_potential = other.use_potential;
            s1_is_constant = other.s1_is_constant;
            s1_constant = other.s1_constant;
            inverse_coef_constant = other.inverse_coef_constant;
            s1_values = std::move(other.s1_values);
            inverse_coef_values = std::move(other.inverse_coef_values);
            G = std::move(other.G);
            work = std::move(other.work);
            sqrt_potential = std::move(other.sqrt_potential);
            plan_forward = other.plan_forward;
            plan_backward = other.plan_backward;
            other.plan_forward = nullptr;
            other.plan_backward = nullptr;
        }
        return *this;
    }

    /**
    * @brief Fast matvec:
    * y = sqrtV .* (C * (sqrtV .* x)) - inverse_coef * x
    */
    void matvec(const std::complex<T>* x, std::complex<T>* y) const {
        // Step 1: Apply potential and zero-pad
        std::fill(work.begin(), work.end(), std::complex<T>(0, 0));

        if (use_potential) {
            for (int64_t k = 0; k < n; ++k) {
                for (int64_t j = 0; j < n; ++j) {
                    for (int64_t i = 0; i < n; ++i) {
                        int64_t src = i + j * n + k * n * n;
                        work[i + j * padded + k * padded * padded] =
                            sqrt_potential[src] * x[src];
                    }
                }
            }
        } else {
            for (int64_t k = 0; k < n; ++k) {
                for (int64_t j = 0; j < n; ++j) {
                    for (int64_t i = 0; i < n; ++i) {
                        int64_t src = i + j * n + k * n * n;
                        work[i + j * padded + k * padded * padded] = x[src];
                    }
                }
            }
        }

        // Step 2: Forward FFT
        if constexpr (std::is_same_v<T, double>) {
            fftw_execute_dft(plan_forward,
                reinterpret_cast<fftw_complex*>(work.data()),
                reinterpret_cast<fftw_complex*>(work.data()));
        }

        // Step 3: Hadamard product
        for (int64_t i = 0; i < padded3; ++i) {
            work[i] *= G[i];
        }

        // Step 4: Inverse FFT
        if constexpr (std::is_same_v<T, double>) {
            fftw_execute_dft(plan_backward,
                reinterpret_cast<fftw_complex*>(work.data()),
                reinterpret_cast<fftw_complex*>(work.data()));
        }

        // Step 5: Extract, normalize, apply potential, add identity, subtract inverse_coef*x
        T scale = T(1.0) / static_cast<T>(padded3);

        if (use_potential) {
            if (s1_is_constant) {
                const T ic = inverse_coef_constant;
                for (int64_t k = 0; k < n; ++k) {
                    for (int64_t j = 0; j < n; ++j) {
                        for (int64_t i = 0; i < n; ++i) {
                            int64_t idx = i + j * n + k * n * n;
                            y[idx] = sqrt_potential[idx] *
                                    work[i + j * padded + k * padded * padded] * scale
                                    - ic * x[idx];
                        }
                    }
                }
            } else {
                for (int64_t k = 0; k < n; ++k) {
                    for (int64_t j = 0; j < n; ++j) {
                        for (int64_t i = 0; i < n; ++i) {
                            int64_t idx = i + j * n + k * n * n;
                            y[idx] = sqrt_potential[idx] *
                                    work[i + j * padded + k * padded * padded] * scale
                                    - inverse_coef_values[idx] * x[idx];
                        }
                    }
                }
            }
        } else {
            if (s1_is_constant) {
                const T ic = inverse_coef_constant;
                for (int64_t k = 0; k < n; ++k) {
                    for (int64_t j = 0; j < n; ++j) {
                        for (int64_t i = 0; i < n; ++i) {
                            int64_t idx = i + j * n + k * n * n;
                            y[idx] = work[i + j * padded + k * padded * padded] * scale
                                    - ic * x[idx];
                        }
                    }
                }
            } else {
                for (int64_t k = 0; k < n; ++k) {
                    for (int64_t j = 0; j < n; ++j) {
                        for (int64_t i = 0; i < n; ++i) {
                            int64_t idx = i + j * n + k * n * n;
                            y[idx] = work[i + j * padded + k * padded * padded] * scale
                                    - inverse_coef_values[idx] * x[idx];
                        }
                    }
                }
            }
        }
    }

    void matvec(const T* x_real, std::complex<T>* y) const {
        std::vector<std::complex<T>> x_complex(N);
        for (int64_t i = 0; i < N; ++i) {
            x_complex[i] = std::complex<T>(x_real[i], T(0));
        }
        matvec(x_complex.data(), y);
    }

    T get_wavenumber() const { return wavenumber; }
    const std::vector<T>& get_s1_values() const { return s1_values; }
    const std::vector<T>& get_inverse_coef_values() const { return inverse_coef_values; }

    size_t memory_usage() const {
        size_t complex_mem = (G.size() + work.size()) * sizeof(std::complex<T>);
        size_t potential_mem = sqrt_potential.size() * sizeof(T);
        return complex_mem + potential_mem;
    }
};

/**
* @brief Helper function to create kernel wrapper with type deduction
*/
template<typename KernelType, typename CoordType, typename DataType>
KernelWrapper<KernelType, CoordType, DataType> make_kernel_wrapper(
    KernelType* kernel, int dimension) {
    return KernelWrapper<KernelType, CoordType, DataType>(kernel, dimension);
}

// ============================================================================
// Matérn 5/2 Kernel
// ============================================================================
//
// Covariance kernel for Gaussian process regression:
//
//   k(r) = (1 + √5·r/ℓ + 5r²/(3ℓ²)) · exp(−√5·r/ℓ)
//
// where ℓ is the length scale and r = ‖x_i − x_j‖.
//
// Matrix convention (covariance matrix, NOT an integral equation):
//   A_ij = k(‖x_i − x_j‖)           i ≠ j
//   A_ii = k(0) + nugget = 1 + nugget
//
// Framework convention:
//   evaluate(x, y) returns N · k(r)   (caller divides by N)
//   evaluate_diagonal(N) returns 1 + nugget
//
// No proxy points should be used — the Matérn kernel is NOT a PDE Green's
// function, so proxy surfaces do not accurately capture the far field.
// ============================================================================


/**
 * @brief 2D Matérn 5/2 kernel for Gaussian process covariance matrices
 *
 * Constructs the covariance matrix on a √N × √N regular grid in [0,1]²:
 *
 *   A_ij = k(‖x_i − x_j‖)          i ≠ j
 *   A_ii = k(0) + nugget = 1 + nugget
 *
 * with k(r) = (1 + √5·r/ℓ + 5r²/(3ℓ²)) exp(−√5·r/ℓ).
 *
 * @tparam T Floating-point type (double recommended)
 */
template<typename T>
struct Matern52Kernel2D {
private:
    int64_t N_total;         // Total grid points (for evaluate() convention)
    T length_scale;          // ℓ: characteristic correlation length
    T nugget;                // σ_n²: diagonal regularisation (noise variance)
    T precomputed_diag_val;  // k(0) + nugget = 1 + nugget

public:
    static constexpr int dimension = 2;
    static constexpr T sqrt5 = T(2.2360679774997896964091736687747632);

    /**
     * @brief Constructor
     * @param N  Total grid points (must satisfy n² = N for some integer n)
     * @param ls Length scale ℓ — controls how fast correlations decay with distance.
     *           Smaller ℓ → faster decay → higher off-diagonal block rank → harder problem.
     * @param nug Nugget σ_n² added to the diagonal for numerical stability.
     */
    Matern52Kernel2D(int64_t N, T ls = T(0.1), T nug = T(1e-6))
        : N_total(N), length_scale(ls), nugget(nug) {
        precomputed_diag_val = T(1) + nugget;
    }

    /**
     * @brief Pure Matérn 5/2 covariance at Euclidean distance r
     *
     * k(r) = (1 + √5·r/ℓ + 5r²/(3ℓ²)) · exp(−√5·r/ℓ)
     *
     * @param r Non-negative distance
     * @return Covariance value in (0, 1]
     */
    T kernel_value(T r) const {
        const T s = sqrt5 * r / length_scale;
        return (T(1) + s + s * s / T(3)) * std::exp(-s);
    }

    /**
     * @brief Framework-compatible point evaluation: returns N · k(‖x−y‖)
     *
     * The caller (verify_solution_direct, FFT constructor) divides by N,
     * yielding the true covariance k(r).
     */
    T evaluate(const T* x, const T* y) const {
        const T dx = x[0] - y[0];
        const T dy = x[1] - y[1];
        const T r  = std::sqrt(dx * dx + dy * dy);
        return static_cast<T>(N_total) * kernel_value(r);
    }

    /**
     * @brief Diagonal matrix entry: k(0) + nugget = 1 + nugget
     *
     * Unlike Green's-function kernels, the Matérn kernel is non-singular
     * at r = 0, so no quadrature integration is needed.
     */
    T evaluate_diagonal(int64_t /*N*/) const {
        return precomputed_diag_val;
    }

    /**
     * @brief Fill a dense submatrix block A_ij = k(x_i, y_j)
     *
     * Off-diagonal entries: k(‖x_i − y_j‖)
     * Diagonal entries (when x_i == y_j): k(0) + nugget
     *
     * @param x_coords Source coordinates, point-major: [x0_0, x0_1, x1_0, x1_1, ...]
     * @param x_size   Number of source points
     * @param y_coords Target coordinates (same layout)
     * @param y_size   Number of target points
     * @param N        Total grid points (unused — scaling is handled internally)
     * @param A        Output matrix, column-major (x_size × y_size)
     * @param lda      Leading dimension of A
     */
    void evaluate_block_precomputed_diagonal(
        const T* x_coords, int64_t x_size,
        const T* y_coords, int64_t y_size,
        int64_t /*N*/,
        T* __restrict__ A, int64_t lda) const {

        for (int64_t j = 0; j < y_size; ++j) {
            const T yj_x = y_coords[j * 2];
            const T yj_y = y_coords[j * 2 + 1];

            for (int64_t i = 0; i < x_size; ++i) {
                const T dx = x_coords[i * 2]     - yj_x;
                const T dy = x_coords[i * 2 + 1] - yj_y;
                const T r  = std::sqrt(dx * dx + dy * dy);
                A[i + j * lda] = kernel_value(r);
            }

            // Override diagonal entries where source == target
            const bool same_ptr = (x_coords == y_coords) && (x_size == y_size);
            for (int64_t i = 0; i < x_size; ++i) {
                const bool pointer_match = same_ptr && (i == j);
                const bool coord_match =
                    (x_coords[i * 2]     == yj_x) &&
                    (x_coords[i * 2 + 1] == yj_y);
                const T match = T(pointer_match || coord_match);
                A[i + j * lda] = match * precomputed_diag_val
                               + (T(1) - match) * A[i + j * lda];
            }
        }
    }

    /** @brief Instance-method wrapper for KernelWrapper compatibility */
    void evaluate_block(
        const T* x_coords, int64_t x_size,
        const T* y_coords, int64_t y_size,
        T* A, int64_t lda) const {
        evaluate_block_precomputed_diagonal(x_coords, x_size, y_coords, y_size, N_total, A, lda);
    }

    T get_length_scale() const { return length_scale; }
    T get_nugget() const { return nugget; }
};


/**
 * @brief FFT-based fast matvec for the 2D Matérn 5/2 covariance on a regular grid
 *
 * Exploits the BCCB (Block Circulant with Circulant Blocks) structure that
 * arises because the Matérn kernel is stationary (translation-invariant) on
 * a uniform n × n grid.
 *
 * Algorithm (circulant embedding):
 *   1. Compute the first column of the Toeplitz matrix (n² entries)
 *   2. Reshape to n × n and embed in a (2n−1) × (2n−1) circulant
 *   3. Precompute G = rfft2(circulant)
 *   4. Matvec: zero-pad x → rfft2 → pointwise multiply by G → irfft2 → extract
 *
 * Complexity: O(N log N) per matvec  vs  O(N²) direct
 * Memory:     ~3 × (2n−1)² reals  (≈ 12N for the real + complex buffers)
 *
 * @tparam T Floating-point type (double)
 */
template<typename T>
class Matern52Kernel2D_FFT {
private:
    int64_t n;           // Grid dimension per axis (n × n = N)
    int64_t N;           // Total points
    int64_t padded_size; // 2n − 1

    T length_scale;
    T nugget;

    std::vector<std::complex<T>> G;           // Precomputed FFT of circulant (Hermitian half)
    std::vector<T> work_real;                 // Real work buffer  [(2n-1)²]
    std::vector<std::complex<T>> work_complex; // Complex work buffer [(2n-1)(n)]

    fftw_plan plan_forward;
    fftw_plan plan_backward;

    static constexpr T sqrt5 = T(2.2360679774997896964091736687747632);

    /** @brief Pure Matérn 5/2 covariance k(r) */
    T kernel_value(T r) const {
        const T s = sqrt5 * r / length_scale;
        return (T(1) + s + s * s / T(3)) * std::exp(-s);
    }

public:
    /**
     * @brief Constructor — precomputes the FFT of the circulant embedding
     *
     * @param grid_points Point coordinates, point-major: [x0,y0, x1,y1, ...]
     * @param grid_size   Grid dimension n (must satisfy n² = N)
     * @param total_points Total number of points N
     * @param ls          Length scale ℓ
     * @param nug         Nugget σ_n²
     */
    Matern52Kernel2D_FFT(const T* grid_points, int64_t grid_size, int64_t total_points,
                          T ls = T(0.1), T nug = T(1e-6))
        : n(grid_size), N(total_points), padded_size(2 * grid_size - 1),
          length_scale(ls), nugget(nug) {

        if (n * n != N) {
            throw std::invalid_argument("Matern52Kernel2D_FFT: grid must be square (n² ≠ N)");
        }

        // --- Allocate buffers ---
        work_real.resize(padded_size * padded_size);
        work_complex.resize(padded_size * (padded_size / 2 + 1));
        G.resize(padded_size * (padded_size / 2 + 1));

        // --- Step 1: First column of the covariance matrix ---
        std::vector<T> first_col(N);
        const T y_fixed[2] = {grid_points[0], grid_points[1]};

        for (int64_t i = 0; i < N; ++i) {
            if (i == 0) {
                first_col[i] = T(1) + nugget;
            } else {
                const T dx = grid_points[i * 2]     - y_fixed[0];
                const T dy = grid_points[i * 2 + 1] - y_fixed[1];
                const T r  = std::sqrt(dx * dx + dy * dy);
                first_col[i] = kernel_value(r);
            }
        }

        // --- Step 2: Reshape first_col to n × n (column-major) ---
        std::vector<T> A(n * n);
        for (int64_t i = 0; i < n; ++i) {
            for (int64_t j = 0; j < n; ++j) {
                A[i + j * n] = first_col[i * n + j];
            }
        }

        // --- Step 3: Build circulant embedding B in work_real ---
        std::fill(work_real.begin(), work_real.end(), T(0));

        for (int64_t j = 0; j < n; ++j)
            for (int64_t i = 0; i < n; ++i)
                work_real[i + j * padded_size] = A[i + j * n];

        for (int64_t j = 1; j < n; ++j)
            for (int64_t i = 0; i < n; ++i)
                work_real[i + (padded_size - j) * padded_size] = A[i + j * n];

        for (int64_t i = 1; i < n; ++i)
            for (int64_t j = 0; j < n; ++j)
                work_real[(padded_size - i) + j * padded_size] = A[i + j * n];

        for (int64_t i = 1; i < n; ++i)
            for (int64_t j = 1; j < n; ++j)
                work_real[(padded_size - i) + (padded_size - j) * padded_size] = A[i + j * n];

        // --- Step 4: G = rfft2(B) ---
        if constexpr (std::is_same_v<T, double>) {
            fftw_plan temp_plan = fftw_plan_dft_r2c_2d(
                padded_size, padded_size,
                work_real.data(),
                reinterpret_cast<fftw_complex*>(G.data()),
                FFTW_ESTIMATE);
            fftw_execute(temp_plan);
            fftw_destroy_plan(temp_plan);

            plan_forward = fftw_plan_dft_r2c_2d(
                padded_size, padded_size,
                work_real.data(),
                reinterpret_cast<fftw_complex*>(work_complex.data()),
                FFTW_ESTIMATE);
            plan_backward = fftw_plan_dft_c2r_2d(
                padded_size, padded_size,
                reinterpret_cast<fftw_complex*>(work_complex.data()),
                work_real.data(),
                FFTW_ESTIMATE);
        } else {
            throw std::runtime_error("Matern52Kernel2D_FFT: only double precision supported");
        }
    }

    ~Matern52Kernel2D_FFT() {
        if constexpr (std::is_same_v<T, double>) {
            fftw_destroy_plan(plan_forward);
            fftw_destroy_plan(plan_backward);
        }
    }

    /**
     * @brief Fast matvec: y = A · x   where A is the Matérn 5/2 covariance matrix
     */
    void matvec(const T* x, T* y) {
        std::fill(work_real.begin(), work_real.end(), T(0));
        for (int64_t i = 0; i < n; ++i)
            for (int64_t j = 0; j < n; ++j)
                work_real[i + j * padded_size] = x[i * n + j];

        if constexpr (std::is_same_v<T, double>) {
            fftw_execute(plan_forward);
        }

        for (int64_t i = 0; i < padded_size * (padded_size / 2 + 1); ++i) {
            work_complex[i] *= G[i];
        }

        if constexpr (std::is_same_v<T, double>) {
            fftw_execute(plan_backward);
        }

        const T scale = T(1) / static_cast<T>(padded_size * padded_size);
        for (int64_t i = 0; i < n; ++i)
            for (int64_t j = 0; j < n; ++j)
                y[i * n + j] = work_real[i + j * padded_size] * scale;
    }

    size_t memory_usage() const {
        const size_t real_mem    = work_real.size() * sizeof(T);
        const size_t complex_mem = (G.size() + work_complex.size()) * sizeof(std::complex<T>);
        return real_mem + complex_mem;
    }
};


/**
 * @brief 3D Matérn 5/2 kernel for Gaussian process covariance matrices
 *
 * Constructs the covariance matrix on an n × n × n regular grid in [0,1]³:
 *
 *   A_ij = k(‖x_i − x_j‖)          i ≠ j
 *   A_ii = k(0) + nugget = 1 + nugget
 *
 * with k(r) = (1 + √5·r/ℓ + 5r²/(3ℓ²)) exp(−√5·r/ℓ).
 *
 * @tparam T Floating-point type (double recommended)
 */
template<typename T>
struct Matern52Kernel3D {
private:
    int64_t N_total;
    T length_scale;
    T nugget;
    T precomputed_diag_val;

public:
    static constexpr int dimension = 3;
    static constexpr T sqrt5 = T(2.2360679774997896964091736687747632);

    Matern52Kernel3D(int64_t N, T ls = T(0.1), T nug = T(1e-6))
        : N_total(N), length_scale(ls), nugget(nug) {
        precomputed_diag_val = T(1) + nugget;
    }

    T kernel_value(T r) const {
        const T s = sqrt5 * r / length_scale;
        return (T(1) + s + s * s / T(3)) * std::exp(-s);
    }

    T evaluate(const T* x, const T* y) const {
        const T dx = x[0] - y[0];
        const T dy = x[1] - y[1];
        const T dz = x[2] - y[2];
        const T r  = std::sqrt(dx * dx + dy * dy + dz * dz);
        return static_cast<T>(N_total) * kernel_value(r);
    }

    T evaluate_diagonal(int64_t /*N*/) const {
        return precomputed_diag_val;
    }

    void evaluate_block_precomputed_diagonal(
        const T* x_coords, int64_t x_size,
        const T* y_coords, int64_t y_size,
        int64_t /*N*/,
        T* __restrict__ A, int64_t lda) const {

        for (int64_t j = 0; j < y_size; ++j) {
            const T yj_x = y_coords[j * 3];
            const T yj_y = y_coords[j * 3 + 1];
            const T yj_z = y_coords[j * 3 + 2];

            for (int64_t i = 0; i < x_size; ++i) {
                const T dx = x_coords[i * 3]     - yj_x;
                const T dy = x_coords[i * 3 + 1] - yj_y;
                const T dz = x_coords[i * 3 + 2] - yj_z;
                const T r  = std::sqrt(dx * dx + dy * dy + dz * dz);
                A[i + j * lda] = kernel_value(r);
            }

            const bool same_ptr = (x_coords == y_coords) && (x_size == y_size);
            for (int64_t i = 0; i < x_size; ++i) {
                const bool pointer_match = same_ptr && (i == j);
                const bool coord_match =
                    (x_coords[i * 3]     == yj_x) &&
                    (x_coords[i * 3 + 1] == yj_y) &&
                    (x_coords[i * 3 + 2] == yj_z);
                const T match = T(pointer_match || coord_match);
                A[i + j * lda] = match * precomputed_diag_val
                               + (T(1) - match) * A[i + j * lda];
            }
        }
    }

    void evaluate_block(
        const T* x_coords, int64_t x_size,
        const T* y_coords, int64_t y_size,
        T* A, int64_t lda) const {
        evaluate_block_precomputed_diagonal(x_coords, x_size, y_coords, y_size, N_total, A, lda);
    }

    T get_length_scale() const { return length_scale; }
    T get_nugget() const { return nugget; }
};


/**
 * @brief FFT-based fast matvec for the 3D Matérn 5/2 covariance on a regular grid
 *
 * Exploits the BTTB (Block Toeplitz with Toeplitz Blocks) structure on an
 * n × n × n uniform grid.  Uses circulant embedding and real-to-complex 3D FFT.
 *
 * @tparam T Floating-point type (double)
 */
template<typename T>
class Matern52Kernel3D_FFT {
private:
    int64_t n;
    int64_t N;
    int64_t padded_size;

    T length_scale;
    T nugget;

    std::vector<std::complex<T>> G;
    std::vector<T> work_real;
    std::vector<std::complex<T>> work_complex;

    fftw_plan plan_forward;
    fftw_plan plan_backward;

    static constexpr T sqrt5 = T(2.2360679774997896964091736687747632);

    T kernel_value(T r) const {
        const T s = sqrt5 * r / length_scale;
        return (T(1) + s + s * s / T(3)) * std::exp(-s);
    }

public:
    Matern52Kernel3D_FFT(const T* grid_points, int64_t grid_size, int64_t total_points,
                          T ls = T(0.1), T nug = T(1e-6))
        : n(grid_size), N(total_points), padded_size(2 * grid_size - 1),
          length_scale(ls), nugget(nug) {

        if (n * n * n != N) {
            throw std::invalid_argument("Matern52Kernel3D_FFT: grid must be cubic (n³ ≠ N)");
        }

        work_real.resize(padded_size * padded_size * padded_size);
        work_complex.resize(padded_size * padded_size * (padded_size / 2 + 1));
        G.resize(padded_size * padded_size * (padded_size / 2 + 1));

        std::vector<T> first_col(N);
        const T y_fixed[3] = {grid_points[0], grid_points[1], grid_points[2]};

        for (int64_t i = 0; i < N; ++i) {
            if (i == 0) {
                first_col[i] = T(1) + nugget;
            } else {
                const T dx = grid_points[i * 3]     - y_fixed[0];
                const T dy = grid_points[i * 3 + 1] - y_fixed[1];
                const T dz = grid_points[i * 3 + 2] - y_fixed[2];
                const T r  = std::sqrt(dx * dx + dy * dy + dz * dz);
                first_col[i] = kernel_value(r);
            }
        }

        std::vector<T> A(N);
        for (int64_t i = 0; i < n; ++i)
            for (int64_t j = 0; j < n; ++j)
                for (int64_t k = 0; k < n; ++k)
                    A[i + j * n + k * n * n] = first_col[i + j * n + k * n * n];

        std::fill(work_real.begin(), work_real.end(), T(0));

        for (int64_t k = 0; k < n; ++k)
            for (int64_t j = 0; j < n; ++j)
                for (int64_t i = 0; i < n; ++i)
                    work_real[i + j * padded_size + k * padded_size * padded_size] =
                        A[i + j * n + k * n * n];

        for (int64_t k = 0; k < n; ++k)
            for (int64_t j = 0; j < n; ++j)
                for (int64_t i = 1; i < n; ++i)
                    work_real[(padded_size - i) + j * padded_size + k * padded_size * padded_size] =
                        work_real[i + j * padded_size + k * padded_size * padded_size];

        for (int64_t k = 0; k < n; ++k)
            for (int64_t j = 1; j < n; ++j)
                for (int64_t i = 0; i < padded_size; ++i)
                    work_real[i + (padded_size - j) * padded_size + k * padded_size * padded_size] =
                        work_real[i + j * padded_size + k * padded_size * padded_size];

        for (int64_t k = 1; k < n; ++k)
            for (int64_t j = 0; j < padded_size; ++j)
                for (int64_t i = 0; i < padded_size; ++i)
                    work_real[i + j * padded_size + (padded_size - k) * padded_size * padded_size] =
                        work_real[i + j * padded_size + k * padded_size * padded_size];

        if constexpr (std::is_same_v<T, double>) {
            fftw_plan temp_plan = fftw_plan_dft_r2c_3d(
                padded_size, padded_size, padded_size,
                work_real.data(),
                reinterpret_cast<fftw_complex*>(G.data()),
                FFTW_ESTIMATE);
            fftw_execute(temp_plan);
            fftw_destroy_plan(temp_plan);

            plan_forward = fftw_plan_dft_r2c_3d(
                padded_size, padded_size, padded_size,
                work_real.data(),
                reinterpret_cast<fftw_complex*>(work_complex.data()),
                FFTW_ESTIMATE);
            plan_backward = fftw_plan_dft_c2r_3d(
                padded_size, padded_size, padded_size,
                reinterpret_cast<fftw_complex*>(work_complex.data()),
                work_real.data(),
                FFTW_ESTIMATE);
        } else {
            throw std::runtime_error("Matern52Kernel3D_FFT: only double precision supported");
        }
    }

    ~Matern52Kernel3D_FFT() {
        if constexpr (std::is_same_v<T, double>) {
            fftw_destroy_plan(plan_forward);
            fftw_destroy_plan(plan_backward);
        }
    }

    void matvec(const T* x, T* y) {
        std::fill(work_real.begin(), work_real.end(), T(0));
        for (int64_t k = 0; k < n; ++k)
            for (int64_t j = 0; j < n; ++j)
                for (int64_t i = 0; i < n; ++i)
                    work_real[i + j * padded_size + k * padded_size * padded_size] =
                        x[i + j * n + k * n * n];

        if constexpr (std::is_same_v<T, double>) {
            fftw_execute(plan_forward);
        }

        for (int64_t i = 0; i < padded_size * padded_size * (padded_size / 2 + 1); ++i) {
            work_complex[i] *= G[i];
        }

        if constexpr (std::is_same_v<T, double>) {
            fftw_execute(plan_backward);
        }

        const T scale = T(1) / static_cast<T>(padded_size * padded_size * padded_size);
        for (int64_t k = 0; k < n; ++k)
            for (int64_t j = 0; j < n; ++j)
                for (int64_t i = 0; i < n; ++i)
                    y[i + j * n + k * n * n] =
                        work_real[i + j * padded_size + k * padded_size * padded_size] * scale;
    }

    size_t memory_usage() const {
        const size_t real_mem    = work_real.size() * sizeof(T);
        const size_t complex_mem = (G.size() + work_complex.size()) * sizeof(std::complex<T>);
        return real_mem + complex_mem;
    }
};


// ============================================================================
// Yukawa (Screened Coulomb) Kernel
// ============================================================================
//
// Green's function for the screened Poisson equation:
//
//   (∇² − κ²) u(x) = −δ(x)
//
// where κ > 0 is the screening parameter (inverse Debye length).
//
// 3D:  K(r) = exp(−κr) / (4πr)
// 2D:  K(r) = K₀(κr) / (2π)
//
// where K₀ is the modified Bessel function of the second kind of order 0.
//
// Behaviour with screening parameter κ:
//   κ → 0:     recovers the Laplace Green's function
//   κ large:   rapid exponential decay → diagonally dominant, easy to compress
//   κ ≈ 1/h:   transition region at grid scale → interesting benchmark
//
// Singularity at r = 0:
//   3D:  1/(4πr)  — same as Laplace 3D
//   2D:  −ln(κr/2)/(2π) + const — log singularity like Laplace 2D
// Both are handled by Gauss-Legendre quadrature on the diagonal.
//
// Applications: Debye-Hückel electrostatics (ions in solution), DLVO colloidal
// interactions, plasma physics, nuclear physics (Yukawa's original meson exchange).
//
// Framework convention (same as Laplace):
//   evaluate(x, y)        returns the raw kernel value K(r)
//   evaluate_diagonal(N)  returns the quadrature-integrated self-interaction
//   Off-diagonal matrix entry:  A_ij = K(‖x_i − x_j‖) / N
//   Diagonal matrix entry:      A_ii = evaluate_diagonal(N)
// ============================================================================


/**
 * @brief Modified Bessel function of the second kind, K₀(x)
 *
 * Polynomial approximations from Abramowitz & Stegun §9.8.
 * Absolute error |ε| < 1.9 × 10⁻⁷, sufficient for kernel evaluation
 * where hierarchical compression tolerance dominates the error budget.
 *
 * @param x Argument (must be > 0)
 * @return K₀(x)
 */
template<typename T>
T bessel_k0(T x) {
    if (x <= T(2)) {
        // A&S 9.8.1: I₀(x) for |x/3.75| < 1
        const T t2 = (x / T(3.75)) * (x / T(3.75));
        const T i0 = T(1) + t2 * (T(3.5156229) + t2 * (T(3.0899424)
                   + t2 * (T(1.2067492) + t2 * (T(0.2659732)
                   + t2 * (T(0.0360768) + t2 * T(0.0045813))))));

        // A&S 9.8.5: K₀(x) = −ln(x/2)·I₀(x) + P((x/2)²)
        const T t = x * x / T(4);
        const T p = T(-0.57721566) + t * (T(0.42278420) + t * (T(0.23069756)
                  + t * (T(0.03488590) + t * (T(0.00262698)
                  + t * (T(0.00010750) + t * T(0.00000740))))));
        return -std::log(x / T(2)) * i0 + p;
    } else {
        // A&S 9.8.6: K₀(x) = exp(−x)/√x · Q(2/x)
        const T t = T(2) / x;
        const T q = T(1.25331414) + t * (T(-0.07832358) + t * (T(0.02189568)
                  + t * (T(-0.01062446) + t * (T(0.00587872)
                  + t * (T(-0.00251540) + t * T(0.00053208))))));
        return std::exp(-x) / std::sqrt(x) * q;
    }
}


/**
 * @brief 2D Yukawa (screened Coulomb) kernel
 *
 * Implements the 2D screened Poisson Green's function:
 *
 *   K(r) = K₀(κr) / (2π)
 *
 * where K₀ is the modified Bessel function of the second kind.
 * As κ → 0 this approaches the 2D Laplace kernel −ln(r)/(2π).
 *
 * @tparam T Floating-point type (double recommended)
 */
template<typename T>
struct YukawaKernel2D {
private:
    int64_t N_total;
    T kappa;                 // Screening parameter κ (inverse Debye length)
    T precomputed_diag_val;

public:
    static constexpr int dimension = 2;
    static constexpr T pi = T(3.14159265358979323846);

    /**
     * @brief Constructor
     * @param N     Total grid points (n² for an n × n grid)
     * @param k     Screening parameter κ. Larger κ → faster decay → easier problem.
     *              Debye length = 1/κ. On [0,1]², κ = 10 gives screening at 10% of domain.
     */
    YukawaKernel2D(int64_t N, T k = T(10))
        : N_total(N), kappa(k) {
        precomputed_diag_val = evaluate_diagonal(N);
    }

    /**
     * @brief Evaluate 2D Yukawa kernel: K₀(κr) / (2π)
     */
    T evaluate(const T* x, const T* y) const {
        const T dx = x[0] - y[0];
        const T dy = x[1] - y[1];
        const T r  = std::sqrt(dx * dx + dy * dy);

        if (r < T(1e-14)) {
            return T(0);  // Singular — diagonal handled by quadrature
        }

        return bessel_k0(kappa * r) / (T(2) * pi);
    }

    /**
     * @brief Diagonal self-interaction using 5-point Gauss-Legendre quadrature
     *
     * Computes: ∫∫_{[0,h/2]²} K₀(κ√(x²+y²)) / (2π) dA   where h = 1/√N
     * then multiplies by 4 to account for integration over the full cell.
     */
    T evaluate_diagonal(int64_t N) const {
        const int64_t n = static_cast<int64_t>(std::sqrt(static_cast<double>(N)));
        const T h = T(1) / static_cast<T>(n);

        auto integrand = [this](T x, T y) -> T {
            const T r = std::sqrt(x * x + y * y);
            if (r < T(1e-14)) return T(0);
            return bessel_k0(kappa * r) / (T(2) * pi);
        };

        // 5-point Gauss-Legendre on [-1,1]
        const T nodes[]   = {T(-0.9061798459386640), T(-0.5384693101056831), T(0),
                             T( 0.5384693101056831), T( 0.9061798459386640)};
        const T weights[] = {T(0.2369268850561891), T(0.4786286704993665), T(0.5688888888888889),
                             T(0.4786286704993665), T(0.2369268850561891)};

        // Map from [-1,1] to [0, h/2]
        const T a = T(0), b = h / T(2);
        const T transform = (b - a) / T(2);
        const T shift     = (b + a) / T(2);

        T integral = T(0);
        for (int i = 0; i < 5; ++i) {
            const T xi = transform * nodes[i] + shift;
            for (int j = 0; j < 5; ++j) {
                const T yj = transform * nodes[j] + shift;
                integral += weights[i] * weights[j] * integrand(xi, yj);
            }
        }

        integral *= transform * transform;  // Jacobian
        return T(4) * integral;             // 4 quadrants
    }

    /**
     * @brief Fill dense submatrix block
     *
     * Off-diagonal: A_ij = K₀(κ‖x_i − y_j‖) / (2πN)
     * Diagonal:     A_ii = precomputed quadrature value
     */
    void evaluate_block_precomputed_diagonal(
        const T* x_coords, int64_t x_size,
        const T* y_coords, int64_t y_size,
        int64_t N,
        T* __restrict__ A, int64_t lda) const {

        const T scale = T(1) / (T(2) * pi * static_cast<T>(N));

        for (int64_t j = 0; j < y_size; ++j) {
            const T yj_x = y_coords[j * 2];
            const T yj_y = y_coords[j * 2 + 1];

            for (int64_t i = 0; i < x_size; ++i) {
                const T dx = x_coords[i * 2]     - yj_x;
                const T dy = x_coords[i * 2 + 1] - yj_y;
                const T r  = std::sqrt(dx * dx + dy * dy);

                if (r < T(1e-14)) {
                    A[i + j * lda] = T(0);  // Will be overridden below
                } else {
                    A[i + j * lda] = scale * bessel_k0(kappa * r);
                }
            }

            // Override diagonal entries
            const bool same_ptr = (x_coords == y_coords) && (x_size == y_size);
            for (int64_t i = 0; i < x_size; ++i) {
                const bool pointer_match = same_ptr && (i == j);
                const bool coord_match =
                    (x_coords[i * 2]     == yj_x) &&
                    (x_coords[i * 2 + 1] == yj_y);
                const T match = T(pointer_match || coord_match);
                A[i + j * lda] = match * precomputed_diag_val
                               + (T(1) - match) * A[i + j * lda];
            }
        }
    }

    void evaluate_block(
        const T* x_coords, int64_t x_size,
        const T* y_coords, int64_t y_size,
        T* A, int64_t lda) const {
        evaluate_block_precomputed_diagonal(x_coords, x_size, y_coords, y_size, N_total, A, lda);
    }

    T get_kappa() const { return kappa; }
};


/**
 * @brief 3D Yukawa (screened Coulomb) kernel
 *
 * Implements the 3D screened Poisson Green's function:
 *
 *   K(r) = exp(−κr) / (4πr)
 *
 * As κ → 0 this recovers the 3D Laplace kernel 1/(4πr).
 *
 * @tparam T Floating-point type (double recommended)
 */
template<typename T>
struct YukawaKernel3D {
private:
    int64_t N_total;
    T kappa;
    T precomputed_diag_val;

public:
    static constexpr int dimension = 3;
    static constexpr T pi = T(3.14159265358979323846);

    YukawaKernel3D(int64_t N, T k = T(10))
        : N_total(N), kappa(k) {
        precomputed_diag_val = evaluate_diagonal(N);
    }

    /**
     * @brief Evaluate 3D Yukawa kernel: exp(−κr) / (4πr)
     */
    T evaluate(const T* x, const T* y) const {
        const T dx = x[0] - y[0];
        const T dy = x[1] - y[1];
        const T dz = x[2] - y[2];
        const T r  = std::sqrt(dx * dx + dy * dy + dz * dz);

        if (r < T(1e-14)) {
            return T(0);
        }

        return std::exp(-kappa * r) / (T(4) * pi * r);
    }

    /**
     * @brief Diagonal self-interaction using 5-point Gauss-Legendre quadrature
     */
    T evaluate_diagonal(int64_t N) const {
        const int64_t n = static_cast<int64_t>(std::cbrt(static_cast<double>(N)));
        const T h = T(1) / static_cast<T>(n);

        auto integrand = [this](T x, T y, T z) -> T {
            const T r = std::sqrt(x * x + y * y + z * z);
            if (r < T(1e-14)) return T(0);
            return std::exp(-kappa * r) / (T(4) * pi * r);
        };

        const T nodes[]   = {T(-0.9061798459386640), T(-0.5384693101056831), T(0),
                             T( 0.5384693101056831), T( 0.9061798459386640)};
        const T weights[] = {T(0.2369268850561891), T(0.4786286704993665), T(0.5688888888888889),
                             T(0.4786286704993665), T(0.2369268850561891)};

        const T a = T(0), b = h / T(2);
        const T transform = (b - a) / T(2);
        const T shift     = (b + a) / T(2);

        T integral = T(0);
        for (int i = 0; i < 5; ++i) {
            const T xi = transform * nodes[i] + shift;
            for (int j = 0; j < 5; ++j) {
                const T yj = transform * nodes[j] + shift;
                for (int k = 0; k < 5; ++k) {
                    const T zk = transform * nodes[k] + shift;
                    integral += weights[i] * weights[j] * weights[k] * integrand(xi, yj, zk);
                }
            }
        }

        integral *= transform * transform * transform;  // Jacobian
        return T(8) * integral;                          // 8 octants
    }

    /**
     * @brief Fill dense submatrix block
     *
     * Off-diagonal: A_ij = exp(−κr) / (4πrN)
     * Diagonal:     A_ii = precomputed quadrature value
     */
    void evaluate_block_precomputed_diagonal(
        const T* x_coords, int64_t x_size,
        const T* y_coords, int64_t y_size,
        int64_t N,
        T* __restrict__ A, int64_t lda) const {

        const T scale = T(1) / (T(4) * pi * static_cast<T>(N));

        for (int64_t j = 0; j < y_size; ++j) {
            const T yj_x = y_coords[j * 3];
            const T yj_y = y_coords[j * 3 + 1];
            const T yj_z = y_coords[j * 3 + 2];

            for (int64_t i = 0; i < x_size; ++i) {
                const T dx = x_coords[i * 3]     - yj_x;
                const T dy = x_coords[i * 3 + 1] - yj_y;
                const T dz = x_coords[i * 3 + 2] - yj_z;
                const T r_sq = dx * dx + dy * dy + dz * dz;
                const T r = std::sqrt(std::max(r_sq, T(1e-28)));
                A[i + j * lda] = scale * std::exp(-kappa * r) / r;
            }

            const bool same_ptr = (x_coords == y_coords) && (x_size == y_size);
            for (int64_t i = 0; i < x_size; ++i) {
                const bool pointer_match = same_ptr && (i == j);
                const bool coord_match =
                    (x_coords[i * 3]     == yj_x) &&
                    (x_coords[i * 3 + 1] == yj_y) &&
                    (x_coords[i * 3 + 2] == yj_z);
                const T match = T(pointer_match || coord_match);
                A[i + j * lda] = match * precomputed_diag_val
                               + (T(1) - match) * A[i + j * lda];
            }
        }
    }

    void evaluate_block(
        const T* x_coords, int64_t x_size,
        const T* y_coords, int64_t y_size,
        T* A, int64_t lda) const {
        evaluate_block_precomputed_diagonal(x_coords, x_size, y_coords, y_size, N_total, A, lda);
    }

    T get_kappa() const { return kappa; }
};


/**
 * @brief FFT-based fast matvec for the 2D Yukawa kernel on a regular grid
 *
 * Uses circulant embedding and real-to-complex FFT.
 * The Yukawa kernel is real, symmetric, and translation-invariant,
 * so the same BCCB approach as Laplace applies.
 *
 * @tparam T Floating-point type (double)
 */
template<typename T>
class YukawaKernel2D_FFT {
private:
    int64_t n;
    int64_t N;
    int64_t padded_size;

    T kappa;

    std::vector<std::complex<T>> G;
    std::vector<T> work_real;
    std::vector<std::complex<T>> work_complex;

    fftw_plan plan_forward;
    fftw_plan plan_backward;

    static constexpr T pi = T(3.14159265358979323846);

public:
    YukawaKernel2D_FFT(const T* grid_points, int64_t grid_size, int64_t total_points,
                        T k = T(10))
        : n(grid_size), N(total_points), padded_size(2 * grid_size - 1), kappa(k) {

        if (n * n != N) {
            throw std::invalid_argument("YukawaKernel2D_FFT: grid must be square (n² ≠ N)");
        }

        work_real.resize(padded_size * padded_size);
        work_complex.resize(padded_size * (padded_size / 2 + 1));
        G.resize(padded_size * (padded_size / 2 + 1));

        std::vector<T> first_col(N);
        const T y_fixed[2] = {grid_points[0], grid_points[1]};

        YukawaKernel2D<T> temp_kernel(N, kappa);
        const T diag_val = temp_kernel.evaluate_diagonal(N);

        for (int64_t i = 0; i < N; ++i) {
            if (i == 0) {
                first_col[i] = diag_val;
            } else {
                const T dx = grid_points[i * 2]     - y_fixed[0];
                const T dy = grid_points[i * 2 + 1] - y_fixed[1];
                const T r  = std::sqrt(dx * dx + dy * dy);
                first_col[i] = bessel_k0(kappa * r) / (T(2) * pi * static_cast<T>(N));
            }
        }

        std::vector<T> A(n * n);
        for (int64_t i = 0; i < n; ++i)
            for (int64_t j = 0; j < n; ++j)
                A[i + j * n] = first_col[i * n + j];

        std::fill(work_real.begin(), work_real.end(), T(0));

        for (int64_t j = 0; j < n; ++j)
            for (int64_t i = 0; i < n; ++i)
                work_real[i + j * padded_size] = A[i + j * n];

        for (int64_t j = 1; j < n; ++j)
            for (int64_t i = 0; i < n; ++i)
                work_real[i + (padded_size - j) * padded_size] = A[i + j * n];

        for (int64_t i = 1; i < n; ++i)
            for (int64_t j = 0; j < n; ++j)
                work_real[(padded_size - i) + j * padded_size] = A[i + j * n];

        for (int64_t i = 1; i < n; ++i)
            for (int64_t j = 1; j < n; ++j)
                work_real[(padded_size - i) + (padded_size - j) * padded_size] = A[i + j * n];

        if constexpr (std::is_same_v<T, double>) {
            fftw_plan temp_plan = fftw_plan_dft_r2c_2d(
                padded_size, padded_size, work_real.data(),
                reinterpret_cast<fftw_complex*>(G.data()), FFTW_ESTIMATE);
            fftw_execute(temp_plan);
            fftw_destroy_plan(temp_plan);

            plan_forward = fftw_plan_dft_r2c_2d(
                padded_size, padded_size, work_real.data(),
                reinterpret_cast<fftw_complex*>(work_complex.data()), FFTW_ESTIMATE);
            plan_backward = fftw_plan_dft_c2r_2d(
                padded_size, padded_size,
                reinterpret_cast<fftw_complex*>(work_complex.data()),
                work_real.data(), FFTW_ESTIMATE);
        } else {
            throw std::runtime_error("YukawaKernel2D_FFT: only double precision supported");
        }
    }

    ~YukawaKernel2D_FFT() {
        if constexpr (std::is_same_v<T, double>) {
            fftw_destroy_plan(plan_forward);
            fftw_destroy_plan(plan_backward);
        }
    }

    void matvec(const T* x, T* y) {
        std::fill(work_real.begin(), work_real.end(), T(0));
        for (int64_t i = 0; i < n; ++i)
            for (int64_t j = 0; j < n; ++j)
                work_real[i + j * padded_size] = x[i * n + j];

        if constexpr (std::is_same_v<T, double>) { fftw_execute(plan_forward); }

        for (int64_t i = 0; i < padded_size * (padded_size / 2 + 1); ++i)
            work_complex[i] *= G[i];

        if constexpr (std::is_same_v<T, double>) { fftw_execute(plan_backward); }

        const T scale = T(1) / static_cast<T>(padded_size * padded_size);
        for (int64_t i = 0; i < n; ++i)
            for (int64_t j = 0; j < n; ++j)
                y[i * n + j] = work_real[i + j * padded_size] * scale;
    }

    size_t memory_usage() const {
        const size_t real_mem    = work_real.size() * sizeof(T);
        const size_t complex_mem = (G.size() + work_complex.size()) * sizeof(std::complex<T>);
        return real_mem + complex_mem;
    }
};


/**
 * @brief FFT-based fast matvec for the 3D Yukawa kernel on a regular grid
 *
 * Circulant embedding + real-to-complex 3D FFT, identical structure to
 * LaplaceKernel3D_FFT but using exp(−κr)/(4πr) as the kernel.
 *
 * @tparam T Floating-point type (double)
 */
template<typename T>
class YukawaKernel3D_FFT {
private:
    int64_t n;
    int64_t N;
    int64_t padded_size;

    T kappa;

    std::vector<std::complex<T>> G;
    std::vector<T> work_real;
    std::vector<std::complex<T>> work_complex;

    fftw_plan plan_forward;
    fftw_plan plan_backward;

    static constexpr T pi = T(3.14159265358979323846);

public:
    YukawaKernel3D_FFT(const T* grid_points, int64_t grid_size, int64_t total_points,
                        T k = T(10))
        : n(grid_size), N(total_points), padded_size(2 * grid_size - 1), kappa(k) {

        if (n * n * n != N) {
            throw std::invalid_argument("YukawaKernel3D_FFT: grid must be cubic (n³ ≠ N)");
        }

        work_real.resize(padded_size * padded_size * padded_size);
        work_complex.resize(padded_size * padded_size * (padded_size / 2 + 1));
        G.resize(padded_size * padded_size * (padded_size / 2 + 1));

        std::vector<T> first_col(N);
        const T y_fixed[3] = {grid_points[0], grid_points[1], grid_points[2]};

        YukawaKernel3D<T> temp_kernel(N, kappa);
        const T diag_val = temp_kernel.evaluate_diagonal(N);

        for (int64_t i = 0; i < N; ++i) {
            if (i == 0) {
                first_col[i] = diag_val;
            } else {
                const T dx = grid_points[i * 3]     - y_fixed[0];
                const T dy = grid_points[i * 3 + 1] - y_fixed[1];
                const T dz = grid_points[i * 3 + 2] - y_fixed[2];
                const T r  = std::sqrt(dx * dx + dy * dy + dz * dz);
                first_col[i] = std::exp(-kappa * r) / (T(4) * pi * r * static_cast<T>(N));
            }
        }

        std::vector<T> A(N);
        for (int64_t i = 0; i < n; ++i)
            for (int64_t j = 0; j < n; ++j)
                for (int64_t k = 0; k < n; ++k)
                    A[i + j * n + k * n * n] = first_col[i + j * n + k * n * n];

        std::fill(work_real.begin(), work_real.end(), T(0));

        for (int64_t k = 0; k < n; ++k)
            for (int64_t j = 0; j < n; ++j)
                for (int64_t i = 0; i < n; ++i)
                    work_real[i + j * padded_size + k * padded_size * padded_size] =
                        A[i + j * n + k * n * n];

        for (int64_t k = 0; k < n; ++k)
            for (int64_t j = 0; j < n; ++j)
                for (int64_t i = 1; i < n; ++i)
                    work_real[(padded_size - i) + j * padded_size + k * padded_size * padded_size] =
                        work_real[i + j * padded_size + k * padded_size * padded_size];

        for (int64_t k = 0; k < n; ++k)
            for (int64_t j = 1; j < n; ++j)
                for (int64_t i = 0; i < padded_size; ++i)
                    work_real[i + (padded_size - j) * padded_size + k * padded_size * padded_size] =
                        work_real[i + j * padded_size + k * padded_size * padded_size];

        for (int64_t k = 1; k < n; ++k)
            for (int64_t j = 0; j < padded_size; ++j)
                for (int64_t i = 0; i < padded_size; ++i)
                    work_real[i + j * padded_size + (padded_size - k) * padded_size * padded_size] =
                        work_real[i + j * padded_size + k * padded_size * padded_size];

        if constexpr (std::is_same_v<T, double>) {
            fftw_plan temp_plan = fftw_plan_dft_r2c_3d(
                padded_size, padded_size, padded_size, work_real.data(),
                reinterpret_cast<fftw_complex*>(G.data()), FFTW_ESTIMATE);
            fftw_execute(temp_plan);
            fftw_destroy_plan(temp_plan);

            plan_forward = fftw_plan_dft_r2c_3d(
                padded_size, padded_size, padded_size, work_real.data(),
                reinterpret_cast<fftw_complex*>(work_complex.data()), FFTW_ESTIMATE);
            plan_backward = fftw_plan_dft_c2r_3d(
                padded_size, padded_size, padded_size,
                reinterpret_cast<fftw_complex*>(work_complex.data()),
                work_real.data(), FFTW_ESTIMATE);
        } else {
            throw std::runtime_error("YukawaKernel3D_FFT: only double precision supported");
        }
    }

    ~YukawaKernel3D_FFT() {
        if constexpr (std::is_same_v<T, double>) {
            fftw_destroy_plan(plan_forward);
            fftw_destroy_plan(plan_backward);
        }
    }

    void matvec(const T* x, T* y) {
        std::fill(work_real.begin(), work_real.end(), T(0));
        for (int64_t k = 0; k < n; ++k)
            for (int64_t j = 0; j < n; ++j)
                for (int64_t i = 0; i < n; ++i)
                    work_real[i + j * padded_size + k * padded_size * padded_size] =
                        x[i + j * n + k * n * n];

        if constexpr (std::is_same_v<T, double>) { fftw_execute(plan_forward); }

        for (int64_t i = 0; i < padded_size * padded_size * (padded_size / 2 + 1); ++i)
            work_complex[i] *= G[i];

        if constexpr (std::is_same_v<T, double>) { fftw_execute(plan_backward); }

        const T scale = T(1) / static_cast<T>(padded_size * padded_size * padded_size);
        for (int64_t k = 0; k < n; ++k)
            for (int64_t j = 0; j < n; ++j)
                for (int64_t i = 0; i < n; ++i)
                    y[i + j * n + k * n * n] =
                        work_real[i + j * padded_size + k * padded_size * padded_size] * scale;
    }

    size_t memory_usage() const {
        const size_t real_mem    = work_real.size() * sizeof(T);
        const size_t complex_mem = (G.size() + work_complex.size()) * sizeof(std::complex<T>);
        return real_mem + complex_mem;
    }
};


} // namespace kernel

 #endif // KERNEL_HPP
