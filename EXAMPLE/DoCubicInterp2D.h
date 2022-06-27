/*! @file
 * @brief Header contains declaration of C functions for 2D interpolation 
*/

#include "math.h"
#include "stdlib.h"
#include "stdio.h" 

#if defined (__cplusplus)
extern "C" {
#endif

double kernelu(double s);

void  CubicInterp2D(double *x, double *y, double *u, int Nx, int Ny, double *xx1, double *yy1, double *v, int No);

#if defined (__cplusplus)
}
#endif
