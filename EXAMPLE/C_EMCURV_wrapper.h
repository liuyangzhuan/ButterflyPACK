#ifndef EMCURV_WRAP /* allow multiple inclusions */
#define EMCURV_WRAP

typedef void* F2Cptr;  // pointer passing fortran derived types to c
typedef void* C2Fptr;  // pointer passing c objects to fortran

//------------------------------------------------------------------------------
// Declartion of FORTRAN subroutines to C code
extern "C" {
void c_emcurv_init(int* Npo,double* Locations,F2Cptr* quant_emcurv, int* model2d, double* wavelength, MPI_Fint* MPIcomm);
void c_emcurv_sample(int* m,int* n,_Complex double *value, F2Cptr quant_emcurv);
}
// -----------------------------------------------------------------------------

#endif
