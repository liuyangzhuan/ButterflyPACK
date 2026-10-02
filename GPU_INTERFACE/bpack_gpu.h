/*
 * ButterflyPACK GPU backends (H2 with the option H2_use_gpu, HODLR with
 * HODLR_use_gpu): the matrix's entries on the device come from an evaluator
 * the application registers, the GPU counterpart of its CPU entry callback.
 * doc/gpu_kernels.md is the user guide; EXAMPLE/gpu/ has one example of
 * each kind of evaluator.
 *
 * Three ways to give the entries, registered after c_bpack_construct_init
 * and before the construction (one per matrix; a later one replaces it):
 *
 *  - entry evaluator: a device function K(i, j) that the library's kernels
 *    inline.  From C++ (CUDA): bpack::gpu::set_entry_evaluator
 *    (bpack_gpu_entry.cuh), compiled with the application by nvcc.  From any
 *    language: its CUDA source text, compiled at run time with NVRTC
 *    (c_bpack_set_gpu_entry_source).  The fastest.
 *
 *  - block evaluator: the application's own launches that fill a batch of
 *    blocks (c_bpack_set_gpu_block_evaluator): for entries that are cheaper
 *    to compute block by block (shared work between entries).
 *
 *  - entry-list evaluator: the application's function that fills the values
 *    of a list of (row, column) pairs (c_bpack_set_gpu_list_evaluator): a
 *    block evaluator in flat form, e.g. vectorized CuPy from Python.
 *
 * Ids are the points' 0-based global ids in the application's ordering
 * (the 1-based ids of the CPU callbacks minus one).  Only the real double
 * (d_) and double complex (z_) libraries have a GPU backend; elsewhere the
 * registrations are ignored.
 */
#ifndef BPACK_GPU_H
#define BPACK_GPU_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Properties of an evaluator (the `flags` of the registrations, or-ed). */
#define BPACK_GPU_SYMMETRIC 1   /* K(i, j) == K(j, i): required by the H2 format */
#define BPACK_GPU_COORDINATES 2 /* reads the points' coordinates (entries: x and y; blocks: row_coords and
                                   col_coords); without it the library may pass none */
#define BPACK_GPU_HOST_IDS 4    /* a block evaluator reads row_ids_host and col_ids_host (the library then
                                   waits for their copy to the host before each batch) */

/* One block of a batch: out(r, c) = K(row_ids[r], col_ids[c]), column-major,
 * out[r + c * ld], r < m, c < n.  `out` is device memory of the matrix's
 * type (double or double complex, the layout of C's double _Complex). */
typedef struct bpack_gpu_block {
    int m, n, ld;
    void* out;
    const int64_t* row_ids;      /* device, m ids */
    const int64_t* col_ids;      /* device, n ids */
    const int64_t* row_ids_host; /* host copies of the same ids; BPACK_GPU_HOST_IDS */
    const int64_t* col_ids_host;
    const double* row_coords;    /* device, dim per point (point r at row_coords[r * dim]); BPACK_GPU_COORDINATES */
    const double* col_coords;
} bpack_gpu_block;

/* A batch of blocks for a block evaluator.  Launch on `stream` (a
 * cudaStream_t); the library orders its later work after it and does not
 * wait for it.  Device scratch should come from allocate / release: the
 * library holds most of the device memory, and memory released here may be
 * reused by its later launches on the same stream. */
typedef struct bpack_gpu_block_batch {
    int count;
    const bpack_gpu_block* blocks;        /* host array of `count` blocks */
    const bpack_gpu_block* blocks_device; /* the same array in device memory */
    int max_m, max_n;                     /* the largest m and n of the batch */
    int dim;                              /* coordinates per point */
    void* stream;
    void* (*allocate)(void* allocator, size_t bytes); /* NULL when out of memory */
    void (*release)(void* allocator, void* ptr);
    void* allocator;
} bpack_gpu_block_batch;

typedef void (*bpack_gpu_block_evaluator)(const bpack_gpu_block_batch* batch, void* user);

/* An entry-list evaluator: values[e] = K(rows[e], cols[e]) for e < count,
 * all three device arrays (values of the matrix's type), on `stream`. */
typedef void (*bpack_gpu_list_evaluator)(int64_t count, const int64_t* rows, const int64_t* cols, void* values,
                                         void* stream, void* user);

/* The registrations, per matrix type (d_: double, z_: double complex). */
void d_c_bpack_set_gpu_block_evaluator(void** bmat, bpack_gpu_block_evaluator evaluate, void* user, const int* flags);
void z_c_bpack_set_gpu_block_evaluator(void** bmat, bpack_gpu_block_evaluator evaluate, void* user, const int* flags);
void d_c_bpack_set_gpu_list_evaluator(void** bmat, bpack_gpu_list_evaluator evaluate, void* user, const int* flags);
void z_c_bpack_set_gpu_list_evaluator(void** bmat, bpack_gpu_list_evaluator evaluate, void* user, const int* flags);

/* An entry evaluator from CUDA source text, compiled with NVRTC at the first
 * use on each rank.  The source defines
 *
 *   __device__ double bpack_entry(const double* x, long long i, const double* y, long long j,
 *                                 const double* params, int dim);
 *
 * (z_: returning bpack_dcomplex, {re, im}), and may define helpers.  params:
 * a device copy of params[0 .. nparams - 1]. */
void d_c_bpack_set_gpu_entry_source(void** bmat, const char* source, const double* params, const int* nparams,
                                    const int* flags);
void z_c_bpack_set_gpu_entry_source(void** bmat, const char* source, const double* params, const int* nparams,
                                    const int* flags);

/* Optional: the registered evaluator's warm-up now (its first-use costs: an
 * NVRTC compile, the application's own), on the GPU the matrix will use,
 * rather than before the first level that calls it, where the library does
 * it otherwise (out of the reported times).  For H2 with H2_use_gpu; HODLR
 * warms up as its construction starts.  seconds[0]: its time on this rank,
 * seconds[1]: of which the NVRTC compile.  Collective over the matrix's ranks. */
void d_c_bpack_gpu_warm_up(void** bmat, double* seconds);
void z_c_bpack_gpu_warm_up(void** bmat, double* seconds);

/* Used by bpack::gpu::set_entry_evaluator (bpack_gpu_entry.cuh): the launches of the
 * library's kernels instantiated for the application's entry evaluator. */
typedef struct bpack_gpu_entry_launchers {
    int version;      /* BPACK_GPU_INTERFACE_VERSION of the headers they were compiled with */
    int scalar_bytes; /* sizeof the value type: 8 (double) or 16 (double complex) */
    int flags;
    const void* entry; /* the evaluator, copied by the library */
    size_t entry_bytes;
    void (*eval)(const void* items, int count, int max_m, int max_n, const char* meta, const void* points,
                 const void* entry, void* stream);
    void (*sketch)(const void* items, int count, int max_d, int max_cols, const void* points, const void* entry,
                   void* stream);
} bpack_gpu_entry_launchers;

void d_c_bpack_set_gpu_entry_launchers(void** bmat, const bpack_gpu_entry_launchers* launchers);
void z_c_bpack_set_gpu_entry_launchers(void** bmat, const bpack_gpu_entry_launchers* launchers);

#ifdef __cplusplus
}
#endif

#endif /* BPACK_GPU_H */
