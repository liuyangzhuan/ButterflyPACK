#include "zButterflyPACK_config.fi"
#ifdef assert
#undef assert
#endif
#include "zparsec_hmat_ptg_bridge.h"

#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>

#ifdef HAVE_PARSEC
#include <parsec.h>
#include <parsec/data.h>
#include <parsec/data_internal.h>
#if defined(DAT) && DAT == 0
#include "zparsec_hmat_factorization_ptg.h"
#elif defined(DAT) && DAT == 1
#include "dparsec_hmat_factorization_ptg.h"
#elif defined(DAT) && DAT == 2
#include "cparsec_hmat_factorization_ptg.h"
#elif defined(DAT) && DAT == 3
#include "sparsec_hmat_factorization_ptg.h"
#else
#include "parsec_hmat_factorization_ptg.h"
#endif

#if defined(DAT) && DAT == 0
#define BPACK_HMAT_PTG_TASKPOOL_T parsec_zparsec_hmat_factorization_ptg_taskpool_t
#define BPACK_HMAT_PTG_DEFAULT_ADT_IDX PARSEC_zparsec_hmat_factorization_ptg_DEFAULT_ADT_IDX
#define BPACK_HMAT_PTG_AUTO_ADT_IDX PARSEC_zparsec_hmat_factorization_ptg_AUTO_ADT_IDX
#elif defined(DAT) && DAT == 1
#define BPACK_HMAT_PTG_TASKPOOL_T parsec_dparsec_hmat_factorization_ptg_taskpool_t
#define BPACK_HMAT_PTG_DEFAULT_ADT_IDX PARSEC_dparsec_hmat_factorization_ptg_DEFAULT_ADT_IDX
#define BPACK_HMAT_PTG_AUTO_ADT_IDX PARSEC_dparsec_hmat_factorization_ptg_AUTO_ADT_IDX
#elif defined(DAT) && DAT == 2
#define BPACK_HMAT_PTG_TASKPOOL_T parsec_cparsec_hmat_factorization_ptg_taskpool_t
#define BPACK_HMAT_PTG_DEFAULT_ADT_IDX PARSEC_cparsec_hmat_factorization_ptg_DEFAULT_ADT_IDX
#define BPACK_HMAT_PTG_AUTO_ADT_IDX PARSEC_cparsec_hmat_factorization_ptg_AUTO_ADT_IDX
#elif defined(DAT) && DAT == 3
#define BPACK_HMAT_PTG_TASKPOOL_T parsec_sparsec_hmat_factorization_ptg_taskpool_t
#define BPACK_HMAT_PTG_DEFAULT_ADT_IDX PARSEC_sparsec_hmat_factorization_ptg_DEFAULT_ADT_IDX
#define BPACK_HMAT_PTG_AUTO_ADT_IDX PARSEC_sparsec_hmat_factorization_ptg_AUTO_ADT_IDX
#else
#define BPACK_HMAT_PTG_TASKPOOL_T parsec_parsec_hmat_factorization_ptg_taskpool_t
#define BPACK_HMAT_PTG_DEFAULT_ADT_IDX PARSEC_parsec_hmat_factorization_ptg_DEFAULT_ADT_IDX
#define BPACK_HMAT_PTG_AUTO_ADT_IDX PARSEC_parsec_hmat_factorization_ptg_AUTO_ADT_IDX
#endif
#endif

extern void c_bpack_hmat_ptg_rank(F2Cptr *ptree, int *rank);
extern void c_bpack_hmat_ptg_size(F2Cptr *ptree, int *size);
extern void c_bpack_hmat_ptg_owner_of(F2Cptr *hmat, F2Cptr *ptree,
                                       int *row, int *col, int *owner);

extern void c_bpack_hmat_ptg_diag(F2Cptr *hmat, F2Cptr *option,
                                  F2Cptr *stats, F2Cptr *ptree,
                                  F2Cptr *msh, int *k,
                                  void **out_ptr, uint64_t *out_bytes);
extern void c_bpack_hmat_ptg_upanel(F2Cptr *hmat, F2Cptr *option,
                                    F2Cptr *stats, F2Cptr *ptree,
                                    F2Cptr *msh, int *k, int *j,
                                    const void *diag_ptr,
                                    void **out_ptr, uint64_t *out_bytes);
extern void c_bpack_hmat_ptg_lpanel(F2Cptr *hmat, F2Cptr *option,
                                    F2Cptr *stats, F2Cptr *ptree,
                                    F2Cptr *msh, int *k, int *i,
                                    const void *diag_ptr,
                                    void **out_ptr, uint64_t *out_bytes);
extern void c_bpack_hmat_ptg_update(F2Cptr *hmat, F2Cptr *option,
                                    F2Cptr *stats, F2Cptr *ptree,
                                    F2Cptr *msh, int *k, int *i, int *j,
                                    const void *l_ptr, const void *u_ptr);
extern void c_bpack_hmat_ptg_release_buffer(void **ptr);

static int slot2d(const bpack_hmat_ptg_t *ctx, int row, int col)
{
    assert(ctx != NULL);
    assert(row >= 0 && row < ctx->nb);
    assert(col >= 0 && col < ctx->nb);
    return row * ctx->nb + col;
}

static void release_slot(bpack_hmat_ptg_slot_t *slot)
{
    if(NULL != slot->ptr) {
        c_bpack_hmat_ptg_release_buffer(&slot->ptr);
    }
#ifdef HAVE_PARSEC
    if(NULL != slot->copy) {
        slot->copy->device_private = NULL;
        slot->copy = NULL;
    }
#endif
    slot->ptr = NULL;
    slot->bytes = 0;
    slot->ready = 0;
}

bpack_hmat_ptg_t *bpack_hmat_ptg_create(F2Cptr hmat, F2Cptr option,
                                         F2Cptr stats, F2Cptr ptree,
                                         F2Cptr msh, int nb)
{
    bpack_hmat_ptg_t *ctx = calloc(1, sizeof(*ctx));
    assert(NULL != ctx);
    assert(nb > 0);

    ctx->hmat = hmat;
    ctx->option = option;
    ctx->stats = stats;
    ctx->ptree = ptree;
    ctx->msh = msh;
    ctx->nb = nb;

    ctx->diag = calloc((size_t)nb, sizeof(*ctx->diag));
    ctx->u = calloc((size_t)nb * (size_t)nb, sizeof(*ctx->u));
    ctx->l = calloc((size_t)nb * (size_t)nb, sizeof(*ctx->l));
    assert(NULL != ctx->diag);
    assert(NULL != ctx->u);
    assert(NULL != ctx->l);

    c_bpack_hmat_ptg_rank(&ctx->ptree, &ctx->myrank);
    return ctx;
}

void bpack_hmat_ptg_destroy(bpack_hmat_ptg_t *ctx)
{
    int idx;

    if(NULL == ctx) {
        return;
    }
    /* Release any packed payloads that were not reached by release tasks. */
    for(idx = 0; idx < ctx->nb; idx++) {
        release_slot(&ctx->diag[idx]);
    }
    for(idx = 0; idx < ctx->nb * ctx->nb; idx++) {
        release_slot(&ctx->u[idx]);
        release_slot(&ctx->l[idx]);
    }
    free(ctx->diag);
    free(ctx->u);
    free(ctx->l);
    free(ctx);
}

int bpack_hmat_ptg_owner_of(bpack_hmat_ptg_t *ctx, int row, int col)
{
    int owner = -1;
    assert(NULL != ctx);
    c_bpack_hmat_ptg_owner_of(&ctx->hmat, &ctx->ptree, &row, &col, &owner);
    return owner;
}

int bpack_hmat_ptg_is_local(bpack_hmat_ptg_t *ctx, int row, int col)
{
    return bpack_hmat_ptg_owner_of(ctx, row, col) == ctx->myrank;
}

void *bpack_hmat_ptg_a_token(bpack_hmat_ptg_t *ctx, int row, int col)
{
    (void)ctx;
    (void)row;
    (void)col;
    return NULL;
}

void bpack_hmat_ptg_diag(bpack_hmat_ptg_t *ctx, int k)
{
    bpack_hmat_ptg_slot_t *slot = &ctx->diag[k];
    assert(!slot->ready);
    c_bpack_hmat_ptg_diag(&ctx->hmat, &ctx->option, &ctx->stats,
                          &ctx->ptree, &ctx->msh, &k,
                          &slot->ptr, &slot->bytes);
    slot->ready = 1;
}

void bpack_hmat_ptg_upanel(bpack_hmat_ptg_t *ctx, int k, int j,
                           const void *diag_bytes)
{
    bpack_hmat_ptg_slot_t *slot = &ctx->u[slot2d(ctx, k, j)];
    assert(!slot->ready);
    c_bpack_hmat_ptg_upanel(&ctx->hmat, &ctx->option, &ctx->stats,
                            &ctx->ptree, &ctx->msh, &k, &j,
                            diag_bytes, &slot->ptr, &slot->bytes);
    slot->ready = 1;
}

void bpack_hmat_ptg_lpanel(bpack_hmat_ptg_t *ctx, int k, int i,
                           const void *diag_bytes)
{
    bpack_hmat_ptg_slot_t *slot = &ctx->l[slot2d(ctx, i, k)];
    assert(!slot->ready);
    c_bpack_hmat_ptg_lpanel(&ctx->hmat, &ctx->option, &ctx->stats,
                            &ctx->ptree, &ctx->msh, &k, &i,
                            diag_bytes, &slot->ptr, &slot->bytes);
    slot->ready = 1;
}

void bpack_hmat_ptg_update(bpack_hmat_ptg_t *ctx, int k, int i, int j,
                           const void *l_bytes, const void *u_bytes)
{
    c_bpack_hmat_ptg_update(&ctx->hmat, &ctx->option, &ctx->stats,
                            &ctx->ptree, &ctx->msh, &k, &i, &j,
                            l_bytes, u_bytes);
}

void bpack_hmat_ptg_diag_release(bpack_hmat_ptg_t *ctx, int k)
{
    assert(ctx != NULL);
    assert(k >= 0 && k < ctx->nb);
    release_slot(&ctx->diag[k]);
}

void bpack_hmat_ptg_u_release(bpack_hmat_ptg_t *ctx, int k, int j)
{
    assert(ctx != NULL);
    release_slot(&ctx->u[slot2d(ctx, k, j)]);
}

void bpack_hmat_ptg_l_release(bpack_hmat_ptg_t *ctx, int i, int k)
{
    assert(ctx != NULL);
    release_slot(&ctx->l[slot2d(ctx, i, k)]);
}

uint64_t bpack_hmat_ptg_diag_bytes(bpack_hmat_ptg_t *ctx, int k)
{
    assert(ctx->diag[k].ready);
    return ctx->diag[k].bytes;
}

uint64_t bpack_hmat_ptg_u_bytes(bpack_hmat_ptg_t *ctx, int k, int j)
{
    bpack_hmat_ptg_slot_t *slot = &ctx->u[slot2d(ctx, k, j)];
    assert(slot->ready);
    return slot->bytes;
}

uint64_t bpack_hmat_ptg_l_bytes(bpack_hmat_ptg_t *ctx, int i, int k)
{
    bpack_hmat_ptg_slot_t *slot = &ctx->l[slot2d(ctx, i, k)];
    assert(slot->ready);
    return slot->bytes;
}

void *bpack_hmat_ptg_diag_ptr(bpack_hmat_ptg_t *ctx, int k)
{
    assert(ctx->diag[k].ready);
    return ctx->diag[k].ptr;
}

void *bpack_hmat_ptg_u_ptr(bpack_hmat_ptg_t *ctx, int k, int j)
{
    bpack_hmat_ptg_slot_t *slot = &ctx->u[slot2d(ctx, k, j)];
    assert(slot->ready);
    return slot->ptr;
}

void *bpack_hmat_ptg_l_ptr(bpack_hmat_ptg_t *ctx, int i, int k)
{
    bpack_hmat_ptg_slot_t *slot = &ctx->l[slot2d(ctx, i, k)];
    assert(slot->ready);
    return slot->ptr;
}

#ifdef HAVE_PARSEC
static int bpack_hmat_ptg_taskpool_init_arenas(parsec_taskpool_t *taskpool)
{
    BPACK_HMAT_PTG_TASKPOOL_T *ptg = (BPACK_HMAT_PTG_TASKPOOL_T *)taskpool;
    int rc;

    assert(ptg != NULL);

    rc = parsec_arena_datatype_construct(
        &ptg->arenas_datatypes[BPACK_HMAT_PTG_DEFAULT_ADT_IDX],
        sizeof(uint8_t), PARSEC_ARENA_ALIGNMENT_SSE, parsec_datatype_int8_t);
    if(rc != PARSEC_SUCCESS) {
        return rc;
    }

    return parsec_arena_datatype_construct(
        &ptg->arenas_datatypes[BPACK_HMAT_PTG_AUTO_ADT_IDX],
        sizeof(uint8_t), PARSEC_ARENA_ALIGNMENT_SSE, parsec_datatype_int8_t);
}

static void hmat_ptg_set_copy_ptr(parsec_data_copy_t *copy, void *ptr)
{
    if(copy != NULL) {
        copy->device_private = ptr;
    }
}

void bpack_hmat_ptg_diag_set_copy(bpack_hmat_ptg_t *ctx, int k,
                                  parsec_data_copy_t *copy)
{
    assert(ctx->diag[k].ready);
    ctx->diag[k].copy = copy;
    hmat_ptg_set_copy_ptr(copy, ctx->diag[k].ptr);
}

void bpack_hmat_ptg_u_set_copy(bpack_hmat_ptg_t *ctx, int k, int j,
                               parsec_data_copy_t *copy)
{
    bpack_hmat_ptg_slot_t *slot = &ctx->u[slot2d(ctx, k, j)];
    assert(slot->ready);
    slot->copy = copy;
    hmat_ptg_set_copy_ptr(copy, slot->ptr);
}

void bpack_hmat_ptg_l_set_copy(bpack_hmat_ptg_t *ctx, int i, int k,
                               parsec_data_copy_t *copy)
{
    bpack_hmat_ptg_slot_t *slot = &ctx->l[slot2d(ctx, i, k)];
    assert(slot->ready);
    slot->copy = copy;
    hmat_ptg_set_copy_ptr(copy, slot->ptr);
}

static int hmat_ptg_data_key_from_args(const bpack_hmat_ptg_desc_t *desc,
                                       va_list ap, int *row, int *col);

static parsec_data_key_t hmat_ptg_data_key(parsec_data_collection_t *dc, ...)
{
    bpack_hmat_ptg_desc_t *desc = (bpack_hmat_ptg_desc_t *)dc;
    va_list ap;
    int row;
    int col;
    int index;

    va_start(ap, dc);
    index = hmat_ptg_data_key_from_args(desc, ap, &row, &col);
    va_end(ap);

    (void)row;
    (void)col;
    return (parsec_data_key_t)index;
}

static int hmat_ptg_data_key_from_args(const bpack_hmat_ptg_desc_t *desc,
                                       va_list ap, int *row, int *col)
{
    assert(desc != NULL);
    switch(desc->kind) {
    case BPACK_HMAT_PTG_DESC_A:
        *row = va_arg(ap, int);
        *col = va_arg(ap, int);
        break;
    case BPACK_HMAT_PTG_DESC_D:
        (void)va_arg(ap, bpack_hmat_ptg_t *);
        *row = va_arg(ap, int);
        *col = *row;
        assert(*row >= 0 && *row < desc->ctx->nb);
        return *row;
    case BPACK_HMAT_PTG_DESC_U:
        (void)va_arg(ap, bpack_hmat_ptg_t *);
        *row = va_arg(ap, int);
        *col = va_arg(ap, int);
        break;
    case BPACK_HMAT_PTG_DESC_L:
        (void)va_arg(ap, bpack_hmat_ptg_t *);
        *row = va_arg(ap, int);
        *col = va_arg(ap, int);
        break;
    default:
        assert(0);
    }

    return slot2d(desc->ctx, *row, *col);
}

static void hmat_ptg_key_to_ij(const bpack_hmat_ptg_desc_t *desc,
                               parsec_data_key_t key, int *row, int *col)
{
    assert(desc != NULL);
    assert(key < (parsec_data_key_t)desc->ntile);
    if(desc->kind == BPACK_HMAT_PTG_DESC_D) {
        *row = (int)key;
        *col = (int)key;
    } else {
        *row = (int)(key / (parsec_data_key_t)desc->ctx->nb);
        *col = (int)(key % (parsec_data_key_t)desc->ctx->nb);
    }
}

static uint32_t hmat_ptg_rank_of(parsec_data_collection_t *dc, ...)
{
    bpack_hmat_ptg_desc_t *desc = (bpack_hmat_ptg_desc_t *)dc;
    va_list ap;
    int row;
    int col;

    va_start(ap, dc);
    (void)hmat_ptg_data_key_from_args(desc, ap, &row, &col);
    va_end(ap);

    return (uint32_t)bpack_hmat_ptg_owner_of(desc->ctx, row, col);
}

static uint32_t hmat_ptg_rank_of_key(parsec_data_collection_t *dc,
                                     parsec_data_key_t key)
{
    bpack_hmat_ptg_desc_t *desc = (bpack_hmat_ptg_desc_t *)dc;
    int row;
    int col;

    hmat_ptg_key_to_ij(desc, key, &row, &col);
    return (uint32_t)bpack_hmat_ptg_owner_of(desc->ctx, row, col);
}

static int32_t hmat_ptg_vpid_of(parsec_data_collection_t *dc, ...)
{
    (void)dc;
    return 0;
}

static int32_t hmat_ptg_vpid_of_key(parsec_data_collection_t *dc,
                                    parsec_data_key_t key)
{
    (void)dc;
    (void)key;
    return 0;
}

static parsec_data_t *hmat_ptg_data_of(parsec_data_collection_t *dc, ...)
{
    bpack_hmat_ptg_desc_t *desc = (bpack_hmat_ptg_desc_t *)dc;
    va_list ap;
    int row;
    int col;
    int index;

    va_start(ap, dc);
    index = hmat_ptg_data_key_from_args(desc, ap, &row, &col);
    va_end(ap);

    (void)row;
    (void)col;
    return desc->tiles[index];
}

static parsec_data_t *hmat_ptg_data_of_key(parsec_data_collection_t *dc,
                                           parsec_data_key_t key)
{
    bpack_hmat_ptg_desc_t *desc = (bpack_hmat_ptg_desc_t *)dc;
    assert(key < (parsec_data_key_t)desc->ntile);
    return desc->tiles[key];
}

static int hmat_ptg_desc_ntile(const bpack_hmat_ptg_t *ctx,
                               bpack_hmat_ptg_desc_kind_t kind)
{
    return (kind == BPACK_HMAT_PTG_DESC_D) ? ctx->nb : ctx->nb * ctx->nb;
}

static bpack_hmat_ptg_desc_t *bpack_hmat_ptg_one_desc_create(
    bpack_hmat_ptg_t *ctx, bpack_hmat_ptg_desc_kind_t kind)
{
    bpack_hmat_ptg_desc_t *desc;
    int nodes = 1;
    int idx;

    assert(ctx != NULL);
    desc = calloc(1, sizeof(*desc));
    assert(desc != NULL);

    c_bpack_hmat_ptg_size(&ctx->ptree, &nodes);
    parsec_data_collection_init(&desc->super, nodes, ctx->myrank);
    desc->super.data_key = hmat_ptg_data_key;
    desc->super.rank_of = hmat_ptg_rank_of;
    desc->super.rank_of_key = hmat_ptg_rank_of_key;
    desc->super.vpid_of = hmat_ptg_vpid_of;
    desc->super.vpid_of_key = hmat_ptg_vpid_of_key;
    desc->super.data_of = hmat_ptg_data_of;
    desc->super.data_of_key = hmat_ptg_data_of_key;

    desc->ctx = ctx;
    desc->kind = kind;
    desc->ntile = hmat_ptg_desc_ntile(ctx, kind);
    desc->tiles = calloc((size_t)desc->ntile, sizeof(*desc->tiles));
    desc->tile_storage = calloc((size_t)desc->ntile, sizeof(*desc->tile_storage));
    assert(desc->tiles != NULL);
    assert(desc->tile_storage != NULL);

    for(idx = 0; idx < desc->ntile; idx++) {
        parsec_data_create(&desc->tiles[idx], &desc->super,
                           (parsec_data_key_t)idx,
                           &desc->tile_storage[idx],
                           sizeof(desc->tile_storage[idx]), 0);
    }

    return desc;
}

static void bpack_hmat_ptg_one_desc_destroy(bpack_hmat_ptg_desc_t *desc)
{
    int idx;

    if(desc == NULL) {
        return;
    }
    for(idx = 0; idx < desc->ntile; idx++) {
        if(desc->tiles[idx] != NULL) {
            parsec_data_destroy(desc->tiles[idx]);
        }
    }
    parsec_data_collection_destroy(&desc->super);
    free(desc->tiles);
    free(desc->tile_storage);
    free(desc);
}

bpack_hmat_ptg_desc_bundle_t *bpack_hmat_ptg_desc_create(bpack_hmat_ptg_t *ctx)
{
    bpack_hmat_ptg_desc_bundle_t *bundle = calloc(1, sizeof(*bundle));
    assert(bundle != NULL);

    bundle->a = bpack_hmat_ptg_one_desc_create(ctx, BPACK_HMAT_PTG_DESC_A);
    bundle->d = bpack_hmat_ptg_one_desc_create(ctx, BPACK_HMAT_PTG_DESC_D);
    bundle->u = bpack_hmat_ptg_one_desc_create(ctx, BPACK_HMAT_PTG_DESC_U);
    bundle->l = bpack_hmat_ptg_one_desc_create(ctx, BPACK_HMAT_PTG_DESC_L);

    return bundle;
}

void bpack_hmat_ptg_desc_destroy(bpack_hmat_ptg_desc_bundle_t *bundle)
{
    if(bundle == NULL) {
        return;
    }
    bpack_hmat_ptg_one_desc_destroy(bundle->a);
    bpack_hmat_ptg_one_desc_destroy(bundle->d);
    bpack_hmat_ptg_one_desc_destroy(bundle->u);
    bpack_hmat_ptg_one_desc_destroy(bundle->l);
    free(bundle);
}

void c_bpack_hmat_factorization_ptg(F2Cptr hmat, F2Cptr option,
                                    F2Cptr stats, F2Cptr ptree,
                                    F2Cptr msh, int nb, int *ierr)
{
    parsec_context_t *parsec = NULL;
    bpack_hmat_ptg_t *ctx = NULL;
    bpack_hmat_ptg_desc_bundle_t *desc = NULL;
    parsec_taskpool_t *taskpool = NULL;

    assert(ierr != NULL);
    *ierr = 0;

    /*
     * The HMAT factor kernels use process-global scratch/state and are not
     * safe to execute concurrently inside one MPI rank.  Keep one PaRSEC
     * worker per rank; distributed parallelism still comes from MPI ranks.
     */
    parsec = parsec_init(1, NULL, NULL);
    if(parsec == NULL) {
        *ierr = -1;
        return;
    }

    ctx = bpack_hmat_ptg_create(hmat, option, stats, ptree, msh, nb);
    desc = bpack_hmat_ptg_desc_create(ctx);
    taskpool = (parsec_taskpool_t *)BPACK_HMAT_PTG_TASKPOOL_NEW(&desc->l->super,
                                                               &desc->u->super,
                                                               &desc->d->super,
                                                               &desc->a->super,
                                                               ctx, nb);
    if(taskpool == NULL) {
        *ierr = -2;
        bpack_hmat_ptg_desc_destroy(desc);
        bpack_hmat_ptg_destroy(ctx);
        parsec_fini(&parsec);
        return;
    }

    *ierr = bpack_hmat_ptg_taskpool_init_arenas(taskpool);
    if(*ierr != PARSEC_SUCCESS) {
        parsec_taskpool_free(taskpool);
        bpack_hmat_ptg_desc_destroy(desc);
        bpack_hmat_ptg_destroy(ctx);
        parsec_fini(&parsec);
        return;
    }

    *ierr = parsec_context_add_taskpool(parsec, taskpool);
    if(*ierr == 0) {
        *ierr = parsec_context_start(parsec);
    }
    if(*ierr == 0) {
        *ierr = parsec_context_wait(parsec);
    }

    parsec_taskpool_free(taskpool);
    bpack_hmat_ptg_desc_destroy(desc);
    bpack_hmat_ptg_destroy(ctx);
    parsec_fini(&parsec);
}
#else
void c_bpack_hmat_factorization_ptg(F2Cptr hmat, F2Cptr option,
                                    F2Cptr stats, F2Cptr ptree,
                                    F2Cptr msh, int nb, int *ierr)
{
    (void)hmat;
    (void)option;
    (void)stats;
    (void)ptree;
    (void)msh;
    (void)nb;
    if(ierr != NULL) {
        *ierr = -1;
    }
}
#endif
