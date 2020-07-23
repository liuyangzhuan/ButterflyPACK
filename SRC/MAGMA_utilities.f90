!
!   -- MAGMA (version 2.5.3) --
!      Univ. of Tennessee, Knoxville
!      Univ. of California, Berkeley
!      Univ. of Colorado, Denver
!      @date March 2020
!
#include "ButterflyPACK_config.fi"
module magma_utilities
    use iso_c_binding
    implicit none

#ifdef HAVE_MAGMA


!! =====================================================================
!! Parameter constants
    real(c_float),             parameter :: sdummy = 0
    real(c_double),            parameter :: ddummy = 0
    complex(c_float_complex),  parameter :: cdummy = 0
    complex(c_double_complex), parameter :: zdummy = 0
    integer(c_int),            parameter :: idummy = 0
    type(c_ptr),               parameter :: ptr_dummy = c_null_ptr

    !! Intel ifort chokes on c_sizeof here, so use extension sizeof
    !! see https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/495001
    integer(c_size_t), parameter :: &
        sizeof_real      = sizeof(sdummy), &
        sizeof_double    = sizeof(ddummy), &
        sizeof_complex   = sizeof(cdummy), &
        sizeof_complex16 = sizeof(zdummy), &
        sizeof_int       = sizeof(idummy), &
        sizeof_ptr       = sizeof(ptr_dummy)

    !! =============================================================================
    !! Parameter constants from magma_types.h
    integer(c_int), parameter ::   &
        MagmaFalse         = 0,    &
        MagmaTrue          = 1,    &

        MagmaRowMajor      = 101,  &
        MagmaColMajor      = 102,  &

        MagmaNoTrans       = 111,  &
        MagmaTrans         = 112,  &
        MagmaConjTrans     = 113,  &

        MagmaUpper         = 121,  &
        MagmaLower         = 122,  &
        MagmaGeneral       = 123,  &
        MagmaFull          = 123,  &  !! deprecated, use MagmaGeneral

        MagmaNonUnit       = 131,  &
        MagmaUnit          = 132,  &

        MagmaLeft          = 141,  &
        MagmaRight         = 142,  &
        MagmaBothSides     = 143
    !! todo all the rest


    !! =============================================================================
    !! Fortran interfaces to C functions
    interface

    !! -------------------------------------------------------------------------
    !! MAGMA_malloc (GPU memory)
    integer(c_int) function MAGMA_malloc( ptr, bytes ) &
    bind(C, name="magma_malloc")
        use iso_c_binding
        type(c_ptr), target :: ptr  !! void**
        integer(c_size_t), value :: bytes
    end function MAGMA_malloc

    integer(c_int) function MAGMA_free_internal( ptr, func, file, line ) &
    bind(C, name="magma_free_internal")
        use iso_c_binding
        type(c_ptr), value :: ptr  !! void*
        character(c_char) :: func, file
        integer(c_int), value :: line
    end function MAGMA_free_internal

    !! -------------------------------------------------------------------------
    !! MAGMA_malloc_pinned (pinned CPU main memory)
    integer(c_int) function MAGMA_malloc_pinned( ptr, bytes ) &
    bind(C, name="magma_malloc_pinned")
        use iso_c_binding
        type(c_ptr), target :: ptr  !! void**
        integer(c_size_t), value :: bytes
    end function MAGMA_malloc_pinned

    integer(c_int) function MAGMA_free_pinned_internal( ptr, func, file, line ) &
    bind(C, name="magma_free_pinned_internal")
        use iso_c_binding
        type(c_ptr), value :: ptr  !! void*
        character(c_char), value :: func, file
        integer(c_int), value :: line
    end function MAGMA_free_pinned_internal

    !! -------------------------------------------------------------------------
    !! set/get
    subroutine MAGMA_setmatrix_internal( &
        m, n, elemsize, hA_src, lda, dB_dst, ldb, queue, func, file, line ) &
    bind(C, name="magma_setmatrix_internal")
        use iso_c_binding
        integer(c_int),    value  :: m, n, elemsize, lda, ldb
        type(c_ptr),       value  :: hA_src
        type(c_ptr),       value  :: dB_dst
        type(c_ptr),       value  :: queue
        character(c_char), value  :: func, file
        integer(c_int),    value  :: line
    end subroutine MAGMA_setmatrix_internal

    subroutine MAGMA_getmatrix_internal( &
        m, n, elemsize, dA_src, lda, hB_dst, ldb, queue, func, file, line ) &
    bind(C, name="magma_getmatrix_internal")
        use iso_c_binding
        integer(c_int),    value  :: m, n, elemsize, lda, ldb
        type(c_ptr),       value  :: dA_src
        type(c_ptr),       value  :: hB_dst
        type(c_ptr),       value  :: queue
        character(c_char), value  :: func, file
        integer(c_int),    value  :: line
    end subroutine MAGMA_getmatrix_internal

    subroutine MAGMA_setvector_internal( &
        n, elemsize, hx_src, incx, dy_dst, incy, queue, func, file, line ) &
    bind(C, name="magma_setvector_internal")
        use iso_c_binding
        integer(c_int),    value  :: n, elemsize, incx, incy
        type(c_ptr),       value  :: hx_src
        type(c_ptr),       value  :: dy_dst
        type(c_ptr),       value  :: queue
        character(c_char), value  :: func, file
        integer(c_int),    value  :: line
    end subroutine MAGMA_setvector_internal

    subroutine MAGMA_getvector_internal( &
        n, elemsize, dx_src, incx, hy_dst, incy, queue, func, file, line ) &
    bind(C, name="magma_getvector_internal")
        use iso_c_binding
        integer(c_int),    value  :: n, elemsize, incx, incy
        type(c_ptr),       value  :: dx_src
        type(c_ptr),       value  :: hy_dst
        type(c_ptr),       value  :: queue
        character(c_char), value  :: func, file
        integer(c_int),    value  :: line
    end subroutine MAGMA_getvector_internal



    !! -------------------------------------------------------------------------
    !! initialize
    subroutine MAGMA_init() &
    bind(C, name="magma_init")
        use iso_c_binding
    end subroutine MAGMA_init

    subroutine MAGMA_finalize() &
    bind(C, name="magma_finalize")
        use iso_c_binding
    end subroutine MAGMA_finalize

    !! -------------------------------------------------------------------------
    !! version
    subroutine MAGMA_version( major, minor, micro ) &
    bind(C, name="magma_version")
        use iso_c_binding
        integer(c_int), target :: major, minor, micro
    end subroutine MAGMA_version

    subroutine MAGMA_print_environment() &
    bind(C, name="magma_print_environment")
        use iso_c_binding
    end subroutine MAGMA_print_environment

    !! -------------------------------------------------------------------------
    !! timing
    real(c_double) function MAGMA_wtime() &
    bind(C, name="magma_wtime")
        use iso_c_binding
    end function MAGMA_wtime

    !! -------------------------------------------------------------------------
    !! device support
    integer(c_int) function MAGMA_num_gpus() &
    bind(C, name="magma_num_gpus")
        use iso_c_binding
    end function MAGMA_num_gpus

    integer(c_int) function MAGMA_get_device_arch() &
    bind(C, name="magma_getdevice_arch")
        use iso_c_binding
    end function MAGMA_get_device_arch

    subroutine MAGMA_get_device( dev ) &
    bind(C, name="magma_getdevice")
        use iso_c_binding
        integer(c_int), target :: dev
    end subroutine MAGMA_get_device

    subroutine MAGMA_set_device( dev ) &
    bind(C, name="magma_setdevice")
        use iso_c_binding
        integer(c_int), value :: dev
    end subroutine MAGMA_set_device

    integer(c_size_t) function MAGMA_mem_size( queue ) &
    bind(C, name="magma_mem_size")
        use iso_c_binding
        type(c_ptr), value :: queue
    end function MAGMA_mem_size

    !! -------------------------------------------------------------------------
    !! queue support
    subroutine MAGMA_queue_create_internal( dev, queue_ptr, func, file, line ) &
    bind(C, name="magma_queue_create_internal")
        use iso_c_binding
        integer(c_int), value :: dev
        type(c_ptr), target :: queue_ptr  !! queue_t*
        character(c_char) :: func, file
        integer(c_int), value :: line
    end subroutine MAGMA_queue_create_internal

    subroutine MAGMA_queue_destroy_internal( queue, func, file, line ) &
    bind(C, name="magma_queue_destroy_internal")
        use iso_c_binding
        type(c_ptr), value :: queue  !! queue_t
        character(c_char) :: func, file
        integer(c_int), value :: line
    end subroutine MAGMA_queue_destroy_internal

    subroutine MAGMA_queue_sync_internal( queue, func, file, line ) &
    bind(C, name="magma_queue_sync_internal")
        use iso_c_binding
        type(c_ptr), value :: queue  !! queue_t
        character(c_char) :: func, file
        integer(c_int), value :: line
    end subroutine MAGMA_queue_sync_internal

    integer(c_int) function MAGMA_queue_get_device( queue ) &
    bind(C, name="magma_queue_get_device")
        use iso_c_binding
        type(c_ptr), value :: queue  !! queue_t
    end function MAGMA_queue_get_device

    !! -------------------------------------------------------------------------
    !! offsets pointers -- 1D vectors with inc
    !! see MAGMA_offset.c
    type(c_ptr) function MAGMA_offset_1d( ptr, inc, i ) &
    bind(C, name="c_magma_offset_1d")
        use iso_c_binding
        type(c_ptr),    value :: ptr
        integer(c_int), value :: inc, i
    end function MAGMA_offset_1d

    !! -------------------------------------------------------------------------
    !! offsets pointers -- 2D matrices with lda
    !! see offset.c
    type(c_ptr) function MAGMA_offset_2d( ptr, lda, i, j ) &
    bind(C, name="c_magma_offset_2d")
        use iso_c_binding
        type(c_ptr),    value:: ptr
        integer(c_int), value :: lda, i, j
    end function MAGMA_offset_2d

    !! batched GPU gemm interfaces
    subroutine MAGMA_gemm_vbatched( &
        transA, transB, m, n, k, alpha, dA_array, &
        ldda, dB_array, lddb,beta,dC_array, lddc, batchCount, queue) &
    bind(C, name="magmablas_gemm_vbatched")
        use iso_c_binding
        integer(c_int),             value :: transA, transB
        type(c_ptr),    value  :: dA_array,dB_array,dC_array    !! double_real**
        type(c_ptr),    value  :: m,n,k  !! int*
        type(c_ptr),    value  :: ldda,lddb,lddc  !! int*
        CBIND_DT,  value :: alpha, beta
        integer(c_int), value  :: batchcount
        type(c_ptr),    value  :: queue
    end subroutine MAGMA_gemm_vbatched

end interface

!! =============================================================================
!! Fortran routines & functions
contains

    !! -------------------------------------------------------------------------
    !! MAGMA_free wrappers
    integer(c_int) function MAGMA_free( ptr )
        type(c_ptr) :: ptr

        MAGMA_free = MAGMA_free_internal( &
                ptr, &
                "magma_free" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end function MAGMA_free

    integer(c_int) function MAGMA_free_pinned( ptr )
        type(c_ptr) :: ptr

        MAGMA_free_pinned = MAGMA_free_internal( &
                ptr, &
                "magma_free_pinned" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end function MAGMA_free_pinned

    !! -------------------------------------------------------------------------
    !! set/get wrappers
    subroutine MAGMA_setmatrix( &
        m, n, elemsize, hA_src, lda, dB_dst, ldb, queue )
        use iso_c_binding
        integer(c_int),    value  :: m, n, elemsize, lda, ldb
        type(c_ptr),       value  :: hA_src
        type(c_ptr),       value  :: dB_dst
        type(c_ptr),       value  :: queue

        call MAGMA_setmatrix_internal( &
                m, n, elemsize, hA_src, lda, dB_dst, ldb, queue, &
                "magma_setmatrix" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine MAGMA_setmatrix

    subroutine MAGMA_getmatrix( &
        m, n, elemsize, dA_src, lda, hB_dst, ldb, queue )
        use iso_c_binding
        integer(c_int),    value  :: m, n, elemsize, lda, ldb
        type(c_ptr),       value  :: dA_src
        type(c_ptr),       value  :: hB_dst
        type(c_ptr),       value  :: queue

        call MAGMA_getmatrix_internal( &
                m, n, elemsize, dA_src, lda, hB_dst, ldb, queue, &
                "magma_getmatrix" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine MAGMA_getmatrix

    subroutine MAGMA_setvector( &
        n, elemsize, hx_src, incx, dy_dst, incy, queue )
        use iso_c_binding
        integer(c_int),    value  :: n, elemsize, incx, incy
        type(c_ptr),       value  :: hx_src
        type(c_ptr),       value  :: dy_dst
        type(c_ptr),       value  :: queue

        call MAGMA_setvector_internal( &
                n, elemsize, hx_src, incx, dy_dst, incy, queue, &
                "magma_setvector" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine MAGMA_setvector

    subroutine MAGMA_getvector( &
        n, elemsize, dx_src, incx, hy_dst, incy, queue )
        use iso_c_binding
        integer(c_int),    value  :: n, elemsize, incx, incy
        type(c_ptr),       value  :: dx_src
        type(c_ptr),       value  :: hy_dst
        type(c_ptr),       value  :: queue

        call MAGMA_getvector_internal( &
                n, elemsize, dx_src, incx, hy_dst, incy, queue, &
                "magma_getvector" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine MAGMA_getvector


    !! -------------------------------------------------------------------------
    !! queue support
    subroutine MAGMA_queue_create( dev, queue_ptr )
        use iso_c_binding
        integer(c_int), value :: dev
        type(c_ptr), target :: queue_ptr  !! queue_t*

        call MAGMA_queue_create_internal( &
                dev, queue_ptr, &
                "magma_queue_create" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine MAGMA_queue_create

    subroutine MAGMA_queue_destroy( queue )
        use iso_c_binding
        type(c_ptr), value :: queue  !! queue_t

        call MAGMA_queue_destroy_internal( &
                queue, &
                "magma_queue_destroy" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine MAGMA_queue_destroy

    subroutine MAGMA_queue_sync( queue )
        use iso_c_binding
        type(c_ptr), value :: queue  !! queue_t

        call MAGMA_queue_sync_internal( &
                queue, &
                "magma_queue_sync" // c_null_char, &
                __FILE__ // c_null_char, &
                __LINE__ )
    end subroutine MAGMA_queue_sync
#endif

end module magma_utilities
