/* “ButterflyPACK” Copyright (c) 2018, The Regents of the University of California, through
  Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the
  U.S. Dept. of Energy). All rights reserved.

  If you have questions about your rights to use or distribute this software, please contact
  Berkeley Lab's Intellectual Property Office at  IPO@lbl.gov.

  NOTICE.  This Software was developed under funding from the U.S. Department of Energy and the
  U.S. Government consequently retains certain rights. As such, the U.S. Government has been
  granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable
  worldwide license in the Software to reproduce, distribute copies to the public, prepare
  derivative works, and perform publicly and display publicly, and to permit other to do so. 

  Developers: Yang Liu, Xiaoye S. Li.
             (Lawrence Berkeley National Lab, Computational Research Division).
*/


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
