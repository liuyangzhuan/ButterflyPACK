#ifndef BPACK_HMAT_PTG_BRIDGE_H
#define BPACK_HMAT_PTG_BRIDGE_H

#include <stdint.h>

#ifdef HAVE_PARSEC
#include <parsec/data_distribution.h>
#endif

#if defined(DAT)
#if DAT == 0
#define BPACK_HMAT_PTG_PREFIX z_
#elif DAT == 1
#define BPACK_HMAT_PTG_PREFIX d_
#elif DAT == 2
#define BPACK_HMAT_PTG_PREFIX cpx_
#elif DAT == 3
#define BPACK_HMAT_PTG_PREFIX s_
#endif
#endif

#ifdef BPACK_HMAT_PTG_PREFIX
#define BPACK_HMAT_PTG_CAT2(a, b) a ## b
#define BPACK_HMAT_PTG_CAT(a, b) BPACK_HMAT_PTG_CAT2(a, b)
#define bpack_hmat_ptg_create BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_create)
#define bpack_hmat_ptg_destroy BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_destroy)
#define bpack_hmat_ptg_owner_of BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_owner_of)
#define bpack_hmat_ptg_is_local BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_is_local)
#define bpack_hmat_ptg_a_token BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_a_token)
#define bpack_hmat_ptg_diag BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_diag)
#define bpack_hmat_ptg_upanel BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_upanel)
#define bpack_hmat_ptg_lpanel BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_lpanel)
#define bpack_hmat_ptg_update BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_update)
#define bpack_hmat_ptg_diag_bytes BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_diag_bytes)
#define bpack_hmat_ptg_u_bytes BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_u_bytes)
#define bpack_hmat_ptg_l_bytes BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_l_bytes)
#define bpack_hmat_ptg_diag_set_copy BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_diag_set_copy)
#define bpack_hmat_ptg_u_set_copy BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_u_set_copy)
#define bpack_hmat_ptg_l_set_copy BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_l_set_copy)
#define bpack_hmat_ptg_diag_release BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_diag_release)
#define bpack_hmat_ptg_u_release BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_u_release)
#define bpack_hmat_ptg_l_release BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_l_release)
#define bpack_hmat_ptg_desc_create BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_desc_create)
#define bpack_hmat_ptg_desc_destroy BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_desc_destroy)
#if !defined(BPACK_HMAT_PTG_NO_PTR_RENAME) && \
    !defined(PARSEC_zparsec_hmat_factorization_ptg_NB_TASK_CLASSES) && \
    !defined(PARSEC_dparsec_hmat_factorization_ptg_NB_TASK_CLASSES) && \
    !defined(PARSEC_cparsec_hmat_factorization_ptg_NB_TASK_CLASSES) && \
    !defined(PARSEC_sparsec_hmat_factorization_ptg_NB_TASK_CLASSES)
#define bpack_hmat_ptg_diag_ptr BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_diag_ptr)
#define bpack_hmat_ptg_u_ptr BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_u_ptr)
#define bpack_hmat_ptg_l_ptr BPACK_HMAT_PTG_CAT(BPACK_HMAT_PTG_PREFIX, bpack_hmat_ptg_l_ptr)
#endif
#endif

#if defined(DAT)
#if DAT == 0
#define BPACK_HMAT_PTG_TASKPOOL_NEW parsec_zparsec_hmat_factorization_ptg_new
#elif DAT == 1
#define BPACK_HMAT_PTG_TASKPOOL_NEW parsec_dparsec_hmat_factorization_ptg_new
#elif DAT == 2
#define BPACK_HMAT_PTG_TASKPOOL_NEW parsec_cparsec_hmat_factorization_ptg_new
#elif DAT == 3
#define BPACK_HMAT_PTG_TASKPOOL_NEW parsec_sparsec_hmat_factorization_ptg_new
#endif
#endif

#ifndef BPACK_HMAT_PTG_TASKPOOL_NEW
#define BPACK_HMAT_PTG_TASKPOOL_NEW parsec_parsec_hmat_factorization_ptg_new
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef void *F2Cptr;

typedef struct bpack_hmat_ptg_slot_s {
    void    *ptr;
    uint64_t bytes;
    int      ready;
#ifdef HAVE_PARSEC
    parsec_data_copy_t *copy;
#endif
} bpack_hmat_ptg_slot_t;

typedef struct bpack_hmat_ptg_s {
    F2Cptr hmat;
    F2Cptr option;
    F2Cptr stats;
    F2Cptr ptree;
    F2Cptr msh;

    int nb;
    int myrank;

    bpack_hmat_ptg_slot_t *diag;
    bpack_hmat_ptg_slot_t *u;
    bpack_hmat_ptg_slot_t *l;
} bpack_hmat_ptg_t;

#ifdef HAVE_PARSEC
typedef enum bpack_hmat_ptg_desc_kind_e {
    BPACK_HMAT_PTG_DESC_A,
    BPACK_HMAT_PTG_DESC_D,
    BPACK_HMAT_PTG_DESC_U,
    BPACK_HMAT_PTG_DESC_L
} bpack_hmat_ptg_desc_kind_t;

typedef struct bpack_hmat_ptg_desc_s {
    parsec_data_collection_t super;
    bpack_hmat_ptg_t *ctx;
    bpack_hmat_ptg_desc_kind_t kind;
    int ntile;
    parsec_data_t **tiles;
    uint8_t *tile_storage;
} bpack_hmat_ptg_desc_t;

typedef struct bpack_hmat_ptg_desc_bundle_s {
    bpack_hmat_ptg_desc_t *a;
    bpack_hmat_ptg_desc_t *d;
    bpack_hmat_ptg_desc_t *u;
    bpack_hmat_ptg_desc_t *l;
} bpack_hmat_ptg_desc_bundle_t;
#endif

/*
 * Packed payload lifetime:
 *
 * - Each D(k), U(k,j), and L(i,k) payload slot is produced exactly once.
 * - Producer tasks allocate/fill their slot and never overwrite it.
 * - PTG release tasks free sender-side packed buffers after their last
 *   consumers complete.  bpack_hmat_ptg_destroy releases anything left after
 *   parsec_context_wait, which also keeps error paths simple.
 */
bpack_hmat_ptg_t *bpack_hmat_ptg_create(F2Cptr hmat, F2Cptr option,
                                         F2Cptr stats, F2Cptr ptree,
                                         F2Cptr msh, int nb);
void bpack_hmat_ptg_destroy(bpack_hmat_ptg_t *ctx);

int   bpack_hmat_ptg_owner_of(bpack_hmat_ptg_t *ctx, int row, int col);
int   bpack_hmat_ptg_is_local(bpack_hmat_ptg_t *ctx, int row, int col);
void *bpack_hmat_ptg_a_token(bpack_hmat_ptg_t *ctx, int row, int col);

void bpack_hmat_ptg_diag(bpack_hmat_ptg_t *ctx, int k);
void bpack_hmat_ptg_upanel(bpack_hmat_ptg_t *ctx, int k, int j,
                           const void *diag_bytes);
void bpack_hmat_ptg_lpanel(bpack_hmat_ptg_t *ctx, int k, int i,
                           const void *diag_bytes);
void bpack_hmat_ptg_update(bpack_hmat_ptg_t *ctx, int k, int i, int j,
                           const void *l_bytes, const void *u_bytes);
void bpack_hmat_ptg_diag_release(bpack_hmat_ptg_t *ctx, int k);
void bpack_hmat_ptg_u_release(bpack_hmat_ptg_t *ctx, int k, int j);
void bpack_hmat_ptg_l_release(bpack_hmat_ptg_t *ctx, int i, int k);

uint64_t bpack_hmat_ptg_diag_bytes(bpack_hmat_ptg_t *ctx, int k);
uint64_t bpack_hmat_ptg_u_bytes(bpack_hmat_ptg_t *ctx, int k, int j);
uint64_t bpack_hmat_ptg_l_bytes(bpack_hmat_ptg_t *ctx, int i, int k);

void *bpack_hmat_ptg_diag_ptr(bpack_hmat_ptg_t *ctx, int k);
void *bpack_hmat_ptg_u_ptr(bpack_hmat_ptg_t *ctx, int k, int j);
void *bpack_hmat_ptg_l_ptr(bpack_hmat_ptg_t *ctx, int i, int k);

#ifdef HAVE_PARSEC
void bpack_hmat_ptg_diag_set_copy(bpack_hmat_ptg_t *ctx, int k,
                                  parsec_data_copy_t *copy);
void bpack_hmat_ptg_u_set_copy(bpack_hmat_ptg_t *ctx, int k, int j,
                               parsec_data_copy_t *copy);
void bpack_hmat_ptg_l_set_copy(bpack_hmat_ptg_t *ctx, int i, int k,
                               parsec_data_copy_t *copy);

bpack_hmat_ptg_desc_bundle_t *bpack_hmat_ptg_desc_create(bpack_hmat_ptg_t *ctx);
void bpack_hmat_ptg_desc_destroy(bpack_hmat_ptg_desc_bundle_t *bundle);
#endif

void z_c_bpack_hmat_factorization_ptg(F2Cptr hmat, F2Cptr option,
                                    F2Cptr stats, F2Cptr ptree,
                                    F2Cptr msh, int nb, int *ierr);

#ifdef __cplusplus
}
#endif

#endif
