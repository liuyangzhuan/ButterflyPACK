/*
 * GPU evaluators of the examples' matrices: the device counterparts of their
 * CPU entry callbacks, for H2_use_gpu (format 7) and HODLR_use_gpu (format
 * 1).  Each registers its evaluator with the matrix after
 * c_bpack_construct_init and before the construction; doc/gpu_kernels.md is
 * the user guide.  Built with the examples when ButterflyPACK has its GPU
 * backends (enable_h2_gpu, which defines BPACK_EXAMPLE_GPU for the drivers).
 *
 *   laplace3d_gpu.cu  entry evaluator (C++): 3D Laplace
 *   vie3d_gpu.cu      entry evaluators (C++): Helmholtz, and Helmholtz
 *                     scaled by a coefficient table on the device
 *   emsurf_gpu.cu     block evaluator: EFIE of RWG edges by triangle pairs;
 *                     entry evaluator: CFIE
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
void vie3d_gpu_release(void);

/* EMSURF_Driver.cpp: the EFIE (cfie_alpha == 1: a block evaluator, complex
 * symmetric) or CFIE (an entry evaluator, not symmetric: HODLR only) entry
 * of RWG edges, from the mesh tables of emsurf_get_gpu_tables_c
 * (EMSURF_C_Bindings.f90).  The tables stay on the device until
 * emsurf_gpu_release. */
void emsurf_gpu_register(void** bmat, const double* params, int nparams, const double* reals, int64_t nreals,
                         const int* ints, int64_t nints, double cfie_alpha);
void emsurf_gpu_release(void);

#ifdef __cplusplus
}
#endif

#endif /* BPACK_GPU_EVALUATORS_EXAMPLE_H */
