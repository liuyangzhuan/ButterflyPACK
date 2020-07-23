#ifdef HAVE_MAGMA
#include "magma.h"

#include <stdio.h>

#if DAT==1

double* d_c_magma_offset_1d(
    double* x, magma_int_t inc,
    magma_int_t i )
{
    return x + (i-1)*inc;
}


double* d_c_magma_offset_2d(
    double* A, magma_int_t lda,
    magma_int_t i, magma_int_t j )
{
    return A + (i-1) + (j-1)*lda;
}

#else

magmaDoubleComplex* z_c_magma_offset_1d(
    magmaDoubleComplex* x, magma_int_t inc,
    magma_int_t i )
{
    return x + (i-1)*inc;
}

magmaDoubleComplex* z_c_magma_offset_2d(
    magmaDoubleComplex* A, magma_int_t lda,
    magma_int_t i, magma_int_t j )
{
    return A + (i-1) + (j-1)*lda;
}

#endif
#endif
