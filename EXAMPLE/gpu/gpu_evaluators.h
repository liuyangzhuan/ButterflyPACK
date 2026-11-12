/*
 * GPU evaluators of the examples' matrices: the device counterparts of their
 * CPU entry callbacks, for H2_use_gpu (format 7) and HODLR_use_gpu (format
 * 1).  Each registers its evaluator with the matrix after
 * c_bpack_construct_init and before the construction; doc/gpu_kernels.md is
 * the user guide.  Built with the examples when ButterflyPACK has its GPU
 * backends (enable_h2_gpu, which defines BPACK_EXAMPLE_GPU for the drivers).
 *
 *   laplace3d_gpu.cu  entry evaluator (C++): 3D Laplace
 *   vie3d_gpu.cu      entry evaluators (C++): Helmholtz, Helmholtz
 *                     scaled by a coefficient table on the device, and
 *                     Helmholtz plus a diagonal table on the device
 *   emsurf_gpu.cu     block evaluators: EFIE and CFIE of RWG edges by
 *                     triangle pairs
 *   butterfly_gpu.cu  entry evaluators (C++) of a butterfly block
 *
 * The Python example (EXAMPLE/user_block_funcs_gp_gpu.py) shows the two
 * routes from Python: entry source text compiled with NVRTC, and an
 * entry-list evaluator with CuPy.
 */
#ifndef BPACK_GPU_EVALUATORS_EXAMPLE_H
#define BPACK_GPU_EVALUATORS_EXAMPLE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Laplace3D_H2Test_Driver.cpp: scale / |x - y|, and `diagonal` when the ids
 * match (real, symmetric). */
void laplace3d_gpu_register(void** bmat, double scale, double diagonal);
/* The same with the driver's nonsymmetric forms: times 1 + skew (d . (x - y))
 * / r off the diagonal (direction: d, 3 values), and every entry times w(y) =
 * 1 + weight sin(2 pi (y0 + 2 y1 + 3 y2)); symmetric when both are 0. */
void laplace3d_gpu_register_general(void** bmat, double scale, double diagonal, double skew,
                                    const double* direction, double weight);

/* RBF3D_H2Test_Driver.cpp: scale * exp(-|x - y|^2 / (2 l^2)) with length_scale > 0,
 * and scale / |x - y| otherwise; `diagonal` when the ids match (real, symmetric). */
void rbf3d_gpu_register(void** bmat, double scale, double diagonal,
                        double length_scale);

/* VIE3D_H2Test_Driver.cpp with scaleGreen 1: a e^{i k r} / (four_pi r), and
 * `self` when the ids match (complex, symmetric). */
void vie3d_gpu_register_helmholtz(void** bmat, double k, double a_re, double a_im, double self_re, double self_im,
                                  double four_pi);

/* VIE3D_H2Test_Driver.cpp with scaleGreen 0: coef[j] times the Helmholtz
 * kernel above, plus (diag_re, diag_im) when the ids match; coef: one value
 * per point, by 0-based global id (complex, not symmetric: HODLR only).
 * The table stays on the device until vie3d_gpu_release. */
void vie3d_gpu_register_coefficient(void** bmat, double k, double a_re, double a_im, double self_re, double self_im,
                                    double four_pi, double diag_re, double diag_im, const double* coef, int64_t n);
/* VIE3D_tensor_Driver_scattereronly.cpp's diagonal blocks (bdiag_format 7):
 * the Helmholtz kernel above, plus diag[i] when the ids match; diag: one
 * value per point, by 0-based global id (complex, symmetric: H2 and HODLR).
 * Each call keeps its own table on the device until vie3d_gpu_release. */
void vie3d_gpu_register_diagonal(void** bmat, double k, double a_re, double a_im, double self_re, double self_im,
                                 double four_pi, const double* diag, int64_t n);
void vie3d_gpu_release(void);

/* Butterfly_Driver.cpp, for a butterfly block (bf: the handle of
 * c_bf_construct_init): cos(2 pi s x_i xi_j), and 1 / (4 pi |p_i - q_j|) with
 * 3 coordinates a point; the tables by row id (m values or points) and by
 * column id (n), copied to the device until butterfly_gpu_release (real, not
 * symmetric).  The _complex forms register with a z_ butterfly:
 * exp(2 pi i s x_i xi_j), and the Laplace kernel's real values. */
void butterfly_gpu_register_dft(void** bf, const double* x, int64_t m, const double* xi, int64_t n, double s);
void butterfly_gpu_register_laplace(void** bf, const double* rows_xyz, int64_t m, const double* cols_xyz, int64_t n);
void butterfly_gpu_register_dft_complex(void** bf, const double* x, int64_t m, const double* xi, int64_t n, double s);
void butterfly_gpu_register_laplace_complex(void** bf, const double* rows_xyz, int64_t m, const double* cols_xyz,
                                            int64_t n);
void butterfly_gpu_release(void);

/* EMSURF_Driver.cpp: the EFIE (cfie_alpha == 1: complex symmetric) or CFIE
 * (not symmetric) entry of RWG edges, a block evaluator either way, from the
 * mesh tables of emsurf_get_gpu_tables_c (EMSURF_C_Bindings.f90).  The tables stay on the device until
 * emsurf_gpu_release. */
void emsurf_gpu_register(void** bmat, const double* params, int nparams, const double* reals, int64_t nreals,
                         const int* ints, int64_t nints, double cfie_alpha);
void emsurf_gpu_release(void);

#ifdef __cplusplus
}
#endif

#endif /* BPACK_GPU_EVALUATORS_EXAMPLE_H */
