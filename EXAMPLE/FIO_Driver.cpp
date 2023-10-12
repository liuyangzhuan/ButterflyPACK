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

  Developers: Yang Liu
             (Lawrence Berkeley National Lab, Computational Research Division).
*/
/*! @file
 * @brief This c++ example compresses 1D Fourier integral operators (FIO) using entry-evaluation-based APIs, then compresses composition of 1D FIOs using matvec-based APIs. The example works on the double-complex data type.
*/
//------------------------------------------------------------------------------
#include <iostream> // header that defines the standard input/output stream objects
#include <math.h> // <cmath> (math.h) header that declares a set of functions to compute common mathematical operations and transformations
#include <fstream> // header providing file stream classes
#include <time.h> // <ctime> header contains definitions of functions to get and manipulate date and time information
#include <stdlib.h>
/* <cstdlib> defines several general purpose functions, including dynamic memory management, random number generation,
   communication with the environment, integer arithmetics, searching, sorting and converting */
#include <sstream> // header provide string stream classes
#include <string> // header introducing string types, character traits and a set of converting functions
#include <iomanip> // header providing parametric manipulators
#include <memory> // header defines general utilities to manage dynamic memory

#include <pthread.h> //

#include <cmath>
#include <cassert> // defines one macro function that can be used as a standard debugging tool
#include <iostream>
#include <random> // introduce random number generation facilities
#include <vector>  // defines the vector container class
#include <atomic> // Atomic types are types that encapsulate a value whose access is guaranteed to not cause data races and can be used to synchronize memory accesses among different threads
#include <mpi.h>
#include <complex.h>

#include <sstream>
#include <cstring>
#include <getopt.h>
#include <unistd.h>



#include "zC_BPACK_wrapper.h"



//------------------------------------------------------------------------------
using namespace std;

const _Complex double Im={0.0,1.0};


#ifdef HAVE_MPI
extern "C" {
      ///////////////////////////////////////////////
      ////// BLACS //////////////////////////////////
      ///////////////////////////////////////////////
      void Cblacs_exit(int);
}
#else
    void Cblacs_exit(int){};
#endif



// The object handling kernel parameters and sampling function
class C_QuantApp {
public:
  vector<double> _data;
  int _d = 0;
  int _n = 0;
  int _ker=1; // 1 is IFFT, 2 is FIO

  int _rank_rand;
  int _m_rand;
  int _n_rand;
  int _ni1_rand;
  int _ni2_rand;
  std::vector<_Complex double> _MatU;
  std::vector<_Complex double> _MatV;

  std::vector<int> _Hperm;
  std::vector<int> _iHperm;
  int _nloc = 0;

  F2Cptr* bmat;  //hierarchical matrix returned by Fortran code
  F2Cptr* bf_a, *bf_b, *bf_c;  //BF returned by Fortran code
  F2Cptr* stats_a, *stats_b, *stats_c;      //statistics structure returned by Fortran code
  F2Cptr* msh_a, *msh_b, *msh_c;		   //mesh structure returned by Fortran code
  F2Cptr* ptree_a, *ptree_b, *ptree_c;      //process tree returned by Fortran code
  F2Cptr* option_a, *option_b, *option_c;      //option structure returned by Fortran code


  C_QuantApp() = default;
  C_QuantApp(int m_rand, int n_rand, int rank_rand, int ker, vector<_Complex double> MatU, vector<_Complex double> MatV)
    : _m_rand(m_rand), _n_rand(n_rand), _rank_rand(rank_rand), _ker(ker), _MatU(move(MatU)), _MatV(move(MatV)){
	// cout<<_n_rand<<_rank_rand<<_MatU.size()<<endl;
    assert(size_t(_m_rand * _rank_rand) == _MatU.size());
    assert(size_t(_n_rand * _rank_rand) == _MatV.size());
	}

  C_QuantApp(int m_rand, int n_rand, int ni1_rand, int ni2_rand, int ker)
    : _m_rand(m_rand), _n_rand(n_rand), _ni1_rand(ni1_rand), _ni2_rand(ni2_rand),_ker(ker){
	}

  inline void Sample(int m, int n, _Complex double* val){

	  	if(_ker==0){
		*val =0;
		for (int k = 0; k < _rank_rand; k++){
			*val += _MatU[k*_m_rand+m]*_MatV[k*_n_rand+n];
		}
		}else if(_ker==1){
			double xi = m/((double)_m_rand);
			double ki = n-_n_rand/2.0;
			double ci = (2+sin(2*M_PI*xi))/8.0;
			double phi = xi*ki+ci*abs(ki);
			// *val = cexp(2*M_PI*Im*phi);
			*val = cos(2*M_PI*phi)+Im*sin(2*M_PI*phi);

			// double xi = m/((double)_n_rand);
			// double ki = n;
			// double phi = xi*ki;
			// *val = cexp(2*M_PI*Im*phi);

		}else if(_ker==2){
			double xi = n/((double)_n_rand);
			double ki = m;
			double phi = -xi*ki;
			// *val = cexp(2*M_PI*Im*phi);
      *val = cos(2*M_PI*phi)+Im*sin(2*M_PI*phi);
		}else if(_ker==3){
			double xi = m/((double)_m_rand);
			double ki = n-_n_rand/2.0;
			double ci = (2+sin(2*M_PI*xi))/8.0;
			double phi = ci*abs(ki);
			// *val = cexp(2*M_PI*Im*phi);
      *val = cos(2*M_PI*phi)+Im*sin(2*M_PI*phi);

			// double xi = m/((double)_n_rand);
			// double ki = n;
			// double phi = xi*ki;
			// *val = cexp(2*M_PI*Im*phi);
		}else if(_ker==4){
      double a=-1.0;
      double b=1.0;
      double h=(b-a)/((double)_m_rand);
			double xj = a + h*m;
      double mu_k = 2.0*M_PI*(n-_n_rand/2)/(b-a);
      // double sj= 1 + 0.3*sin(M_PI*xj/8.0);
      double sj= 1 + 0.2*tanh(cos(M_PI*xj/8.0));

      *val = pow(abs(mu_k),2*sj) * (cos(mu_k*(xj))+Im*sin(mu_k*(xj)))/(double)_m_rand;
		}else if(_ker==5){
      double a=-1.0;
      double b=1.0;
      double h=(b-a)/((double)_n_rand);
			double xl = a + h*n;
      double mu_k = 2.0*M_PI*(m-_m_rand/2)/(b-a);
      *val = (cos(-mu_k*xl)+Im*sin(-mu_k*xl));
		}
  }
};


// The sampling function wrapper required by the Fortran HODLR code
inline void C_FuncZmn(int *m, int *n, _Complex double *val, C2Fptr quant) {

  C_QuantApp* Q = (C_QuantApp*) quant;
  Q->Sample(*m-1,*n-1,val);
}

// The sampling function wrapper required by the Fortran HODLR code
inline void C_FuncBZmn(int *m, int *n, _Complex double *val, C2Fptr quant) {

  C_QuantApp* Q = (C_QuantApp*) quant;
// here positve inidex means row, negative index means column
  int m1, n1;
  if(*m>0){
	  m1=*m;
	  n1=-*n;
  }else{
	  m1=*n;
	  n1=-*m;
  }
  Q->Sample(m1-1,n1-1,val);
}

// The distance function wrapper required by the Fortran HODLR code
inline void C_FuncDistmn_dummy(int *m, int *n, double *val, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;

}

// The compressibility function wrapper required by the Fortran HODLR code
inline void C_FuncNearFar_dummy(int *m, int *n, int *val, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;

}

// The extraction sampling function wrapper required by the Fortran HODLR code
inline void C_FuncZmnBlock(int* Ninter, int* Nallrows, int* Nallcols, int* Nalldat_loc, int* allrows, int* allcols, _Complex double* alldat_loc, int* rowidx,int* colidx, int* pgidx, int* Npmap, int* pmaps, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;
}

// The extraction sampling function wrapper required by the Fortran HODLR code
inline void C_FuncBZmnBlock(int* Ninter, int* Nallrows, int* Nallcols, int* Nalldat_loc, int* allrows, int* allcols, _Complex double* alldat_loc, int* rowidx,int* colidx, int* pgidx, int* Npmap, int* pmaps, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;

//   z_c_bf_extractelement(Q->bf,Q->option,Q->msh,Q->stats,Q->ptree,Ninter,Nallrows, Nallcols, Nalldat_loc, allrows,allcols,alldat_loc,rowidx,colidx,pgidx,Npmap,pmaps);

}



// The matvec sampling function wrapper required by the Fortran HODLR code
inline void C_FuncHMatVec(char const *trans, int *nin, int *nout, int *nvec, _Complex double const *xin, _Complex double *xout, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;
  int cnt = (*nvec)*(*nout);
  _Complex double* xbuf1 = new _Complex double[(*nvec)*(Q->_ni1_rand)];

  if(*trans=='N'){
	z_c_bf_mult(trans, xin, xbuf1, nin, &(Q->_ni1_rand), nvec, Q->bf_b,Q->option_b,Q->stats_b,Q->ptree_b);
	z_c_bf_mult(trans, xbuf1, xout, &(Q->_ni1_rand), nout, nvec, Q->bf_a,Q->option_a,Q->stats_a,Q->ptree_a);
  }else if(*trans=='T'){
	z_c_bf_mult(trans, xin, xbuf1, nin, &(Q->_ni1_rand), nvec, Q->bf_a,Q->option_a,Q->stats_a,Q->ptree_a);
	z_c_bf_mult(trans, xbuf1, xout, &(Q->_ni1_rand), nout, nvec, Q->bf_b,Q->option_b,Q->stats_b,Q->ptree_b);
  }


// for (int nth=0; nth<*nvec; nth++){
//   for(int i=0; i<*nout; i++)
//     xout[i+nth*(Q->_ni1_rand)]=xin[i+nth*(Q->_ni1_rand)];
//   xout[*nout-1+nth*(Q->_ni1_rand)]=xin[*nout-1+nth*(Q->_ni1_rand)]*1e-10;
// }

// for (int nth=0; nth<*nvec; nth++){
//   for(int i=0; i<*nout; i++)
//     xout[i+nth*(Q->_ni1_rand)]=xin[i+nth*(Q->_ni1_rand)]*(i+1);
// }

  delete[] xbuf1;
}


// The matvec sampling function wrapper required by the Fortran HODLR code
inline void C_FuncBMatVec(char const *trans, int *nin, int *nout, int *nvec, _Complex double const *xin, _Complex double *xout, C2Fptr quant, _Complex double *a, _Complex double *b) {
  C_QuantApp* Q = (C_QuantApp*) quant;
  int cnt = (*nvec)*(*nout);
  _Complex double* xout1 = new _Complex double[cnt];
  _Complex double* xbuf1 = new _Complex double[(*nvec)*(Q->_ni1_rand)];
  _Complex double* xbuf2 = new _Complex double[(*nvec)*(Q->_ni2_rand)];


  if(*trans=='N'){
	z_c_bf_mult(trans, xin, xbuf2, nin, &(Q->_ni2_rand), nvec, Q->bf_c,Q->option_c,Q->stats_c,Q->ptree_c);
	z_c_bf_mult(trans, xbuf2, xbuf1, &(Q->_ni2_rand), &(Q->_ni1_rand), nvec, Q->bf_b,Q->option_b,Q->stats_b,Q->ptree_b);
	z_c_bf_mult(trans, xbuf1, xout1, &(Q->_ni1_rand), nout, nvec, Q->bf_a,Q->option_a,Q->stats_a,Q->ptree_a);
  }else if(*trans=='T'){
	z_c_bf_mult(trans, xin, xbuf1, nin, &(Q->_ni1_rand), nvec, Q->bf_a,Q->option_a,Q->stats_a,Q->ptree_a);
	z_c_bf_mult(trans, xbuf1, xbuf2, &(Q->_ni1_rand), &(Q->_ni2_rand), nvec, Q->bf_b,Q->option_b,Q->stats_b,Q->ptree_b);
	z_c_bf_mult(trans, xbuf2, xout1, &(Q->_ni2_rand), nout, nvec, Q->bf_c,Q->option_c,Q->stats_c,Q->ptree_c);
  }

// z_c_bf_mult(trans, xin, xout1, nin, nout, nvec, Q->bf_b,Q->option_b,Q->stats_b,Q->ptree_b);


  for (int ii=0; ii<cnt; ii++){
	xout[ii] = *b*xout[ii] + *a*xout1[ii];
  }
  delete[] xout1;
  delete[] xbuf1;
  delete[] xbuf2;
}


// The command line parser for the example related parameters
void set_option_from_command_line(int argc, const char* const* cargv,F2Cptr option0) {
    double opt_d;
    int opt_i;
    std::vector<std::unique_ptr<char[]>> argv_data(argc);
    std::vector<char*> argv(argc);
    for (int i=0; i<argc; i++) {
      argv_data[i].reset(new char[strlen(cargv[i])+1]);
      argv[i] = argv_data[i].get();
      strcpy(argv[i], cargv[i]);
    }
    option long_options[] =
      {{"nmin_leaf",                     required_argument, 0, 1},
       {"tol_comp",                   required_argument, 0, 2},
       {"tol_rand",                   required_argument, 0, 3},
       {"tol_Rdetect",             required_argument, 0, 4},
       {"tol_itersol",             required_argument, 0, 5},
       {"n_iter",          required_argument, 0, 6},
       {"level_check",         required_argument, 0, 7},
       {"precon",                  required_argument, 0, 8},
       {"xyzsort",      required_argument, 0, 9},
       {"lrlevel",     required_argument, 0, 10},
       {"errfillfull",       required_argument, 0, 11},
       {"baca_batch",      required_argument, 0, 12},
       {"reclr_leaf",      required_argument, 0, 13},
       {"nogeo",     required_argument, 0, 14},
       {"less_adapt",            required_argument, 0, 15},
       {"errsol",           required_argument, 0, 16},
       {"lr_blk_num",                  required_argument, 0, 17},
       {"rank0",  required_argument, 0, 18},
       {"rankrate", required_argument, 0, 19},
       {"itermax",               required_argument, 0, 20},
       {"powiter",  required_argument, 0, 21},
       {"ilu", required_argument, 0, 22},
       {"nbundle",     required_argument, 0, 23},
       {"near_para",  required_argument, 0, 24},
       {"format",  required_argument, 0, 25},
       {"verbosity", required_argument, 0, 26},
       {"rmax", required_argument, 0, 27},
       {"sample_para", required_argument, 0, 28},
       {"pat_comp",    required_argument, 0, 29},
       {"knn",         required_argument, 0, 30},
       {"forwardN15flag",         required_argument, 0, 31},
       {NULL, 0, NULL, 0}
      };
    int c, option_index = 0;
    // bool unrecognized_options = false;
    opterr = optind = 0;
    while ((c = getopt_long_only
            (argc, argv.data(), "",
             long_options, &option_index)) != -1) {
      switch (c) {
      case 1: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "Nmin_leaf", opt_i);
      } break;
      case 2: {
        std::istringstream iss(optarg);
        iss >> opt_d;
        z_c_bpack_set_D_option(&option0, "tol_comp", opt_d);
        z_c_bpack_set_D_option(&option0, "tol_rand", opt_d);
        z_c_bpack_set_D_option(&option0, "tol_Rdetect", opt_d*0.1);
      } break;
      case 3: {
        std::istringstream iss(optarg);
        iss >> opt_d;
        z_c_bpack_set_D_option(&option0, "tol_rand", opt_d);
      } break;
      case 4: {
        std::istringstream iss(optarg);
        iss >> opt_d;
        z_c_bpack_set_D_option(&option0, "tol_Rdetect", opt_d);
      } break;
      case 5: {
        std::istringstream iss(optarg);
        iss >> opt_d;
        z_c_bpack_set_D_option(&option0, "tol_itersol", opt_d);
      } break;
      case 6: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "n_iter", opt_i);
      } break;
      case 7: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "level_check", opt_i);
      } break;
      case 8: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "precon", opt_i);
      } break;
      case 9: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "xyzsort", opt_i);
      } break;
      case 10: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "LRlevel", opt_i);
      } break;
      case 11: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "ErrFillFull", opt_i);
      } break;
      case 12: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "BACA_Batch", opt_i);
      } break;
      case 13: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "RecLR_leaf", opt_i);
      } break;
      case 14: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "nogeo", opt_i);
      } break;
      case 15: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "less_adapt", opt_i);
      } break;
      case 16: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "ErrSol", opt_i);
      } break;
      case 17: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "LR_BLK_NUM", opt_i);
      } break;
      case 18: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "rank0", opt_i);
      } break;
      case 19: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "rankrate", opt_i);
      } break;
      case 20: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "itermax", opt_i);
      } break;
      case 21: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "powiter", opt_i);
      } break;
      case 22: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "ILU", opt_i);
      } break;
      case 23: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "Nbundle", opt_i);
      } break;
      case 24: {
        std::istringstream iss(optarg);
        iss >> opt_d;
        z_c_bpack_set_D_option(&option0, "near_para", opt_d);
      } break;
      case 25: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "format", opt_i);
      } break;
      case 26: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "verbosity", opt_i);
      } break;
      case 27: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "rmax", opt_i);
      } break;
      case 28: {
        std::istringstream iss(optarg);
        iss >> opt_d;
        z_c_bpack_set_D_option(&option0, "sample_para", opt_d);
      } break;
      case 29: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "pat_comp", opt_i);
      } break;
      case 30: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "knn", opt_i);
      } break;
      case 31: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        z_c_bpack_set_I_option(&option0, "forwardN15flag", opt_i);
      } break;
      default: break;
      }
    }
  }

// This example uses entry evaluation to compute three butterflies and use matvec to compress their products. The three butterflies are of size MxN, NxK, KxL, the fourth butterfly is of size MxL
////////////////////////////////////////////////////////////////////////////////
// --------------------------- Main Code Starts Here ------------------------ //

int main(int argc, char* argv[])
{

    int myrank, size;                     // Store values of processor rank and total no of procs requestedss
    int master_rank = 0;
	MPI_Init(&argc, &argv); 	                            // Initialize MPI, called only once
    MPI_Comm_size(MPI_COMM_WORLD, &size); 	                // Get no of procs
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 	                // Get no of procs
    MPI_Op op;
	double h=0.1; //kernel parameter
	double lambda=10.0 ; //kernel parameter
	int Npo=5000;   // matrix size
	int Ndim=1; //data dimension
	double starttime, endtime;
	double* dat_ptr_m, *dat_ptr_n, *dat_ptr_k, *dat_ptr_l;
	int* nns_ptr_m, *nns_ptr_n, *nns_ptr_k, *nns_ptr_l;
	int nogeo;  // 1: no geometrical information passed to hodlr, dat_ptr and Ndim are dummy
	int ker=0 ; // kernel choice

	int Nmin=4; //finest leafsize
	double tol=1e-4; //compression tolerance
	double sample_para=2.0; //oversampling factor in entry evaluation
	int com_opt=5; //1:SVD 2:RRQR 3:ACA 4:BACA 5:BACA_improved 6:Pseudo-skeleton
	int sort_opt=1; //0:natural order 1:kd-tree 2:cobble-like ordering 3:gram distance-based cobble-like ordering
	int checkerr = 0; //1: check compression quality
	int batch = 100; //batch size for BACA
	int bnum = 1; //sqrt of #of subblocks in H-BACA
	int knn=0; //k nearest neighbours stored per point
	C_QuantApp *quant_ptr_a, *quant_ptr_b, *quant_ptr_c, *quant_ptr_dummy;
	int v_major,v_minor,v_bugfix; //version numbers

  int tst = 1;
	int lrlevel=100;

if(myrank==master_rank){
	z_c_bpack_getversionnumber(&v_major,&v_minor,&v_bugfix);
	std::cout<<"ButterflyPACK Version: "<<v_major<<"."<<v_minor<<"."<<v_bugfix<<std::endl;
}



	/*****************************************************************/
	/* Test Product of two Random matrices*/
	int rank_rand = 100;
	Npo = 1000;
	int M = Npo;
	int N = Npo;
	int K = Npo;
	int L = Npo;



  //getting the example configurations from command line
  std::vector<std::unique_ptr<char[]>> argv_data(argc);
  std::vector<char*> argv1(argc);
  for (int i=0; i<argc; i++) {
    argv_data[i].reset(new char[strlen(argv[i])+1]);
    argv1[i] = argv_data[i].get();
    strcpy(argv1[i], argv[i]);
  }
  option long_options[] =
    {{"tst",                     required_argument, 0, 1},
      {"M",                   required_argument, 0, 2},
      {"N",                   required_argument, 0, 3},
      {"K",             required_argument, 0, 4},
      {"L",             required_argument, 0, 5},
      {"rank_rand",          required_argument, 0, 6},
      {NULL, 0, NULL, 0}
    };
  int c, option_index = 0;
  opterr = optind = 0;
  while ((c = getopt_long_only
          (argc, argv1.data(), "",
            long_options, &option_index)) != -1) {

    switch (c) {
    case 1: {
      std::istringstream iss(optarg);
      iss >> tst;
    } break;
    case 2: {
      std::istringstream iss(optarg);
      iss >> M;
    } break;
    case 3: {
      std::istringstream iss(optarg);
      iss >> N;
    } break;
    case 4: {
      std::istringstream iss(optarg);
      iss >> K;
    } break;
    case 5: {
      std::istringstream iss(optarg);
      iss >> L;
    } break;
    case 6: {
      std::istringstream iss(optarg);
      iss >> rank_rand;
    } break;
    default: break;
    }
  }


    if(tst ==1){
		vector<_Complex double> matU(M*rank_rand);
		for (int i=0; i<M*rank_rand; i++)
		matU[i] = (_Complex double)rand() / RAND_MAX;
		MPI_Bcast(matU.data(), M*rank_rand, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		vector<_Complex double> matV(N*rank_rand);
		for (int i=0; i<N*rank_rand; i++)
		matV[i] = (_Complex double)rand() / RAND_MAX;
		MPI_Bcast(matV.data(), N*rank_rand, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		quant_ptr_a=new C_QuantApp(M, N, rank_rand, 0, matU, matV);

		matU.resize(N*rank_rand);
		for (int i=0; i<N*rank_rand; i++)
		matU[i] = (_Complex double)rand() / RAND_MAX;
		MPI_Bcast(matU.data(), N*rank_rand, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		matV.resize(K*rank_rand);
		for (int i=0; i<K*rank_rand; i++)
		matV[i] = (_Complex double)rand() / RAND_MAX;
		MPI_Bcast(matV.data(), K*rank_rand, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		quant_ptr_b=new C_QuantApp(N, K, rank_rand, 0, matU, matV);


		matU.resize(K*rank_rand);
		for (int i=0; i<K*rank_rand; i++)
		matU[i] = (_Complex double)rand() / RAND_MAX;
		MPI_Bcast(matU.data(), K*rank_rand, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		matV.resize(L*rank_rand);
		for (int i=0; i<L*rank_rand; i++)
		matV[i] = (_Complex double)rand() / RAND_MAX;
		MPI_Bcast(matV.data(), L*rank_rand, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		quant_ptr_c=new C_QuantApp(K, L, rank_rand, 0, matU, matV);

	}

	if(tst ==2){
		quant_ptr_a=new C_QuantApp(M, N, 0, 0, 1);
		quant_ptr_b=new C_QuantApp(N, K, 0, 0, 2);
		quant_ptr_c=new C_QuantApp(K, L, 0, 0, 3);
  	if(myrank==master_rank)std::cout<<"M "<<M<<" N "<<N<<" K "<<K<<" L "<<L<<std::endl;
	}

	if(tst ==3){
		quant_ptr_a=new C_QuantApp(M, N, 0, 0, 4);
    quant_ptr_b=new C_QuantApp(N, K, 0, 0, 5);
  	if(myrank==master_rank)std::cout<<"M "<<M<<" N "<<N<<" K "<<K<<std::endl;
	}


	nogeo=1;
	sort_opt=0;

	/*****************************************************************/
	int myseg_m=0;     // local number of unknowns
	int myseg_n=0;     // local number of unknowns
	int myseg_k=0;     // local number of unknowns
	int myseg_l=0;     // local number of unknowns
	int* perms_m = new int[M]; //permutation vector returned by HODLR
	int* perms_n = new int[N]; //permutation vector returned by HODLR
	int* perms_k = new int[K]; //permutation vector returned by HODLR
	int* perms_l = new int[L]; //permutation vector returned by HODLR
	int* groups = new int[size];
	int i_opt;
	double d_opt;
	int cpp=1; //1: use user-defined cpp/c functions for construction

	//quantities for the first holdr
	F2Cptr bmat_dummy;  //hierarchical matrix returned by Fortran code
	F2Cptr option;     //option structure returned by Fortran code
	F2Cptr stats_dummy,stats_a,stats_b,stats_c;      //statistics structure returned by Fortran code
	F2Cptr msh_a,msh_b,msh_c;		   //z_mesh structure returned by Fortran code
	F2Cptr msh0_m, msh0_n, msh0_k, msh0_l;		   //z_mesh structure returned by Fortran code
	F2Cptr kerquant_a,kerquant_b,kerquant_c;   //kernel quantities structure returned by Fortran code
	F2Cptr kerquant_dummy;   //kernel quantities structure returned by Fortran code
	F2Cptr ptree;      //process tree returned by Fortran code

	F2Cptr bf_a,bf_b,bf_c,bf_mv;  //BF returned by Fortran code

	MPI_Fint Fcomm;  // the fortran MPI communicator
	Fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);

	for (int i = 0; i < size; i++)groups[i]=i;
	// create hodlr data structures
	z_c_bpack_createptree(&size, groups, &Fcomm, &ptree);
	z_c_bpack_createoption(&option);
	z_c_bpack_createstats(&stats_dummy);

	// set hodlr options
	z_c_bpack_set_D_option(&option, "tol_comp", tol*1e-1);      // bf_a, bf_b, bf_c uses this tolerance
	z_c_bpack_set_D_option(&option, "tol_rand", tol);           // bf_mv uses this tolerance
	z_c_bpack_set_D_option(&option, "tol_Rdetect", tol*3e-1);   // bf_mv uses this tolerance
	z_c_bpack_set_D_option(&option, "sample_para", sample_para);
	z_c_bpack_set_I_option(&option, "nogeo", nogeo);
	z_c_bpack_set_I_option(&option, "Nmin_leaf", Nmin);
	z_c_bpack_set_I_option(&option, "RecLR_leaf", com_opt);
	z_c_bpack_set_I_option(&option, "xyzsort", sort_opt);
	z_c_bpack_set_I_option(&option, "ErrFillFull", checkerr);
	z_c_bpack_set_I_option(&option, "BACA_Batch", batch);
	z_c_bpack_set_I_option(&option, "LR_BLK_NUM", bnum);
	z_c_bpack_set_I_option(&option, "cpp", cpp);
	z_c_bpack_set_I_option(&option, "LRlevel", lrlevel);
	z_c_bpack_set_I_option(&option, "knn", knn);
	z_c_bpack_set_I_option(&option, "verbosity", 2);
	z_c_bpack_set_I_option(&option, "less_adapt", 1);
	z_c_bpack_set_I_option(&option, "itermax", 10);
	z_c_bpack_set_I_option(&option, "ErrSol", 1);
	z_c_bpack_set_I_option(&option, "elem_extract", 0);// not use block-wise element extraction

  set_option_from_command_line(argc, argv, option);
	z_c_bpack_printoption(&option,&ptree);


    // construct the mesh data for the bf, this should be improved
	int nlevel_m = 0; // 0: tree level, nonzero if a tree is provided
	int* tree_m = new int[(int)pow(2,nlevel_m)]; //user provided array containing size of each leaf node, not used if nlevel=0
	tree_m[0] = M;
	z_c_bpack_construct_init(&M, &Ndim, dat_ptr_m, nns_ptr_m,&nlevel_m, tree_m, perms_m, &myseg_m, &bmat_dummy, &option, &stats_dummy, &msh0_m, &kerquant_dummy, &ptree, &C_FuncDistmn_dummy, &C_FuncNearFar_dummy, quant_ptr_dummy);


	int nlevel_n = 0; // 0: tree level, nonzero if a tree is provided
	int* tree_n = new int[(int)pow(2,nlevel_n)]; //user provided array containing size of each leaf node, not used if nlevel=0
	tree_n[0] = N;
	z_c_bpack_construct_init(&N, &Ndim, dat_ptr_n, nns_ptr_n,&nlevel_n, tree_n, perms_n, &myseg_n, &bmat_dummy, &option, &stats_dummy, &msh0_n, &kerquant_dummy, &ptree, &C_FuncDistmn_dummy, &C_FuncNearFar_dummy, quant_ptr_dummy);

	int nlevel_k = 0; // 0: tree level, nonzero if a tree is provided
	int* tree_k = new int[(int)pow(2,nlevel_k)]; //user provided array containing size of each leaf node, not used if nlevel=0
	tree_k[0] = K;
	z_c_bpack_construct_init(&K, &Ndim, dat_ptr_k, nns_ptr_k,&nlevel_k, tree_k, perms_k, &myseg_k, &bmat_dummy, &option, &stats_dummy, &msh0_k, &kerquant_dummy, &ptree, &C_FuncDistmn_dummy, &C_FuncNearFar_dummy, quant_ptr_dummy);

	int nlevel_l = 0; // 0: tree level, nonzero if a tree is provided
	int* tree_l = new int[(int)pow(2,nlevel_l)]; //user provided array containing size of each leaf node, not used if nlevel=0
	tree_l[0] = L;
	z_c_bpack_construct_init(&L, &Ndim, dat_ptr_l, nns_ptr_l,&nlevel_l, tree_l, perms_l, &myseg_l, &bmat_dummy, &option, &stats_dummy, &msh0_l, &kerquant_dummy, &ptree, &C_FuncDistmn_dummy, &C_FuncNearFar_dummy, quant_ptr_dummy);

// construct the three bfs from entry evaluation
	z_c_bpack_createstats(&stats_a);
	z_c_bf_construct_init(&M, &N, &myseg_m, &myseg_n, nns_ptr_m, nns_ptr_n, &msh0_m, &msh0_n, &bf_a, &option, &stats_a, &msh_a, &kerquant_a, &ptree,&C_FuncDistmn_dummy, &C_FuncNearFar_dummy, quant_ptr_a);
	z_c_bf_construct_element_compute(&bf_a, &option, &stats_a, &msh_a, &kerquant_a, &ptree, &C_FuncBZmn, &C_FuncBZmnBlock, quant_ptr_a); // C_FuncBZmnBlock is not referenced since elem_extract=0

  // {/* test matrix vector multiplication performance*/
  // int nvec=1;
  // int cnt_i = (nvec)*(myseg_n);
  // int cnt_o = (nvec)*(myseg_m);
  // _Complex double* xout = new _Complex double[cnt_o];
  // _Complex double* xin = new _Complex double[cnt_i];
  // for(int i=0;i<cnt_i;i++)
  //   xin[i]=1;
  // char trans='N';
	// z_c_bf_mult(&trans, xin, xout, &myseg_n, &myseg_m, &nvec, &bf_a, &option, &stats_a, &ptree);
  // }

  if(myrank==master_rank)std::cout<<"\nPrinting stats of the first BF: "<<std::endl;
	z_c_bpack_printstats(&stats_a,&ptree);



	z_c_bpack_createstats(&stats_b);
	z_c_bf_construct_init(&N, &K, &myseg_n, &myseg_k, nns_ptr_n, nns_ptr_k, &msh0_n, &msh0_k, &bf_b, &option, &stats_b, &msh_b, &kerquant_b, &ptree,&C_FuncDistmn_dummy, &C_FuncNearFar_dummy, quant_ptr_b);
	z_c_bf_construct_element_compute(&bf_b, &option, &stats_b, &msh_b, &kerquant_b, &ptree, &C_FuncBZmn, &C_FuncBZmnBlock, quant_ptr_b); // C_FuncBZmnBlock is not referenced since elem_extract=0
	if(myrank==master_rank)std::cout<<"\nPrinting stats of the second BF: "<<std::endl;
	z_c_bpack_printstats(&stats_b,&ptree);


if(tst ==3){
  /* test matrix vector multiplication performance*/
  if(myrank==master_rank)
  std::cout<<"\n Test one matrix vector multiplication: "<<std::endl;
  int nvec=1;
  double a=-1.0;
  double b=1.0;
  double s=4.0;
  double h=(b-a)/((double)K);
  vector<_Complex double> xout1((nvec)*(myseg_n));
  for (int ii=0;ii<xout1.size();ii++)
    xout1[ii]={0.0,0.0};
  vector<_Complex double> xout2((nvec)*(myseg_m));
  for (int ii=0;ii<xout2.size();ii++)
    xout2[ii]={0.0,0.0};
  vector<_Complex double> xin((nvec)*(myseg_k));
  for (int ii=0;ii<xin.size();ii++)
    xin[ii]={0.0,0.0};
  vector<_Complex double> xin_glo((nvec)*(K));
  for (int ii=0;ii<xin_glo.size();ii++)
    xin_glo[ii]={0.0,0.0};
  vector<_Complex double> xout_glo((nvec)*(M));
  for (int ii=0;ii<xout_glo.size();ii++)
    xout_glo[ii]={0.0,0.0};
  fstream myfile0;

  // generate the global input vector on all MPI ranks
  for (int nth=0; nth<nvec; nth++){    // this loops over all vectors, currently nvec=1
    for (int i=0; i<myseg_k; i++){
      int i_new_loc = i+1;
      int i_old;
      z_c_bf_new2old_col(&msh0_k,&i_new_loc,&i_old);
      xin_glo.data()[i_old-1+nth*K] = pow(1.0-(a+h*(i_old-1))*(a+h*(i_old-1)),s)+0.0*Im;
      xin.data()[i+nth*myseg_k] = xin_glo.data()[i_old-1+nth*K]; // generate the local input vector
    }
  }
  MPI_Allreduce(MPI_IN_PLACE,xin_glo.data(), K*nvec, MPI_C_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

  if(myrank==master_rank){   // only MPI rank 0 write data to file
    fstream myfile1;
    myfile1.open("sol_u.txt",fstream::out);
    for (int nth=0; nth<nvec; nth++){    // this loops over all vectors, currently nvec=1
      for (int i=0; i<K; i++){
        myfile1 << creal(xin_glo.data()[i+nth*K]) << "  " << cimag(xin_glo.data()[i+nth*K])<< std::endl;
      }
    }
  }

  char trans='N';
	z_c_bf_mult(&trans, xin.data(), xout1.data(), &myseg_n, &myseg_m, &nvec, &bf_b, &option, &stats_b, &ptree);
	z_c_bf_mult(&trans, xout1.data(), xout2.data(), &myseg_n, &myseg_m, &nvec, &bf_a, &option, &stats_a, &ptree);

  // gather the global output on all MPI ranks
  for (int nth=0; nth<nvec; nth++){    // this loops over all vectors, currently nvec=1
    for (int i=0; i<myseg_m; i++){
      int i_new_loc = i+1;
      int i_old;
      z_c_bf_new2old_row(&msh0_m,&i_new_loc,&i_old);
      xout_glo.data()[i_old-1+nth*M] = xout2.data()[i+nth*myseg_m];
    }
  }
  MPI_Allreduce(MPI_IN_PLACE,xout_glo.data(), M*nvec, MPI_C_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

  /* save the matrix vector production results*/
  if(myrank==master_rank){   // only MPI rank 0 write data to file
    fstream myfile;
    myfile.open("vofL.txt",fstream::out);
    for (int nth=0; nth<nvec; nth++){    // this loops over all vectors, currently nvec=1
      for (int i=0; i<M; i++){
        myfile << creal(xout_glo.data()[i+nth*M]) << "  " << cimag(xout_glo.data()[i+nth*M])<< std::endl;
      }
    }
  }
  }


  if(tst ==3){
  /* test TFQMR*/
  if(myrank==master_rank)
  std::cout<<"\n Test TFQMR: "<<std::endl;

  int nvec=1;
  vector<_Complex double> rhs((nvec)*(myseg_m));
  for (int ii=0;ii<rhs.size();ii++)
    rhs[ii]={0.0,0.0};
  vector<_Complex double> rhs_glo((nvec)*(M));
  for (int ii=0;ii<rhs_glo.size();ii++)
    rhs_glo[ii]={0.0,0.0};
  vector<_Complex double> x_glo((nvec)*(K));
  for (int ii=0;ii<x_glo.size();ii++)
    x_glo[ii]={0.0,0.0};
  vector<_Complex double> x((nvec)*(myseg_k));
  for (int ii=0;ii<x.size();ii++)
    x[ii]={0.0,0.0};
  fstream myfile0;

  // generate the global RHS vector on all MPI ranks
  for (int nth=0; nth<nvec; nth++){    // this loops over all vectors, currently nvec=1
    for (int i=0; i<myseg_m; i++){
      int i_new_loc = i+1;
      int i_old;
      z_c_bf_new2old_col(&msh0_m,&i_new_loc,&i_old);
      rhs_glo.data()[i_old-1+nth*M] = 1.0;
      rhs.data()[i+nth*myseg_m] = rhs_glo.data()[i_old-1+nth*M]; // generate the local input vector
    }
  }
  MPI_Allreduce(MPI_IN_PLACE,rhs_glo.data(), M*nvec, MPI_C_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
  if(myrank==master_rank){   // only MPI rank 0 write data to file
    fstream myfile1;
    myfile1.open("rhs.txt",fstream::out);
    for (int nth=0; nth<nvec; nth++){    // this loops over all vectors, currently nvec=1
      for (int i=0; i<M; i++){
        myfile1 << creal(rhs_glo.data()[i+nth*M]) << "  " << cimag(rhs_glo.data()[i+nth*M])<< std::endl;
      }
    }
  }

  C_QuantApp *quant_ptr_H;

	quant_ptr_H=new C_QuantApp(M,K,myseg_n,0,0); // the last two parameters are not used

  quant_ptr_H->bf_a=&bf_a;
	quant_ptr_H->msh_a=&msh_a;
	quant_ptr_H->ptree_a=&ptree;
	quant_ptr_H->stats_a=&stats_a;
	quant_ptr_H->option_a=&option;

	quant_ptr_H->bf_b=&bf_b;
	quant_ptr_H->msh_b=&msh_b;
	quant_ptr_H->ptree_b=&ptree;
	quant_ptr_H->stats_b=&stats_b;
	quant_ptr_H->option_b=&option;


  F2Cptr kerquant_H;
  z_c_bpack_tfqmr_noprecon(x.data(),rhs.data(),&myseg_m,&nvec,&option, &stats_a, &ptree, &kerquant_H, &C_FuncHMatVec, quant_ptr_H);

  // gather the global solution on all MPI ranks
  for (int nth=0; nth<nvec; nth++){    // this loops over all vectors, currently nvec=1
    for (int i=0; i<myseg_k; i++){
      int i_new_loc = i+1;
      int i_old;
      z_c_bf_new2old_row(&msh0_k,&i_new_loc,&i_old);
      x_glo.data()[i_old-1+nth*K] = x.data()[i+nth*myseg_k];
    }
  }
  MPI_Allreduce(MPI_IN_PLACE,x_glo.data(), K*nvec, MPI_C_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

  /* save the tfqmr results*/
  if(myrank==master_rank){   // only MPI rank 0 write data to file
    fstream myfile;
    myfile.open("x.txt",fstream::out);
    for (int nth=0; nth<nvec; nth++){    // this loops over all vectors, currently nvec=1
      for (int i=0; i<K; i++){
        myfile << creal(x_glo.data()[i+nth*K]) << "  " << cimag(x_glo.data()[i+nth*K])<< std::endl;
      }
    }
  }
  }

if(tst ==1 || tst ==2){
	z_c_bpack_createstats(&stats_c);
	z_c_bf_construct_init(&K, &L, &myseg_k, &myseg_l, nns_ptr_k, nns_ptr_l, &msh0_k, &msh0_l, &bf_c, &option, &stats_c, &msh_c, &kerquant_c, &ptree,&C_FuncDistmn_dummy, &C_FuncNearFar_dummy, quant_ptr_c);
	z_c_bf_construct_element_compute(&bf_c, &option, &stats_c, &msh_c, &kerquant_c, &ptree, &C_FuncBZmn, &C_FuncBZmnBlock, quant_ptr_c); // C_FuncBZmnBlock is not referenced since elem_extract=0
	if(myrank==master_rank)std::cout<<"\nPrinting stats of the third BF: "<<std::endl;
	z_c_bpack_printstats(&stats_c,&ptree);


	//////////////////// use resulting bf as matvec to create a new bf

	// quantities for the matvec-based bf construction

	F2Cptr option_mv;     //option structure returned by Fortran code
	F2Cptr stats_mv;      //statistics structure returned by Fortran code
	F2Cptr msh_mv;		   //z_mesh structure returned by Fortran code
	F2Cptr kerquant_mv;   //kernel quantities structure returned by Fortran code
	F2Cptr ptree_mv;      //process tree returned by Fortran code
	C_QuantApp *quant_ptr_mv; //user-defined object

	quant_ptr_mv=new C_QuantApp(M,L,N,K,0);
	quant_ptr_mv->bf_a=&bf_a;
	quant_ptr_mv->msh_a=&msh_a;
	quant_ptr_mv->ptree_a=&ptree;
	quant_ptr_mv->stats_a=&stats_a;
	quant_ptr_mv->option_a=&option;

	quant_ptr_mv->bf_b=&bf_b;
	quant_ptr_mv->msh_b=&msh_b;
	quant_ptr_mv->ptree_b=&ptree;
	quant_ptr_mv->stats_b=&stats_b;
	quant_ptr_mv->option_b=&option;

	quant_ptr_mv->bf_c=&bf_c;
	quant_ptr_mv->msh_c=&msh_c;
	quant_ptr_mv->ptree_c=&ptree;
	quant_ptr_mv->stats_c=&stats_c;
	quant_ptr_mv->option_c=&option;


	z_c_bpack_createptree(&size, groups, &Fcomm, &ptree_mv);
	z_c_bpack_copyoption(&option,&option_mv);
	z_c_bpack_createstats(&stats_mv);

	z_c_bpack_set_I_option(&option_mv, "nogeo", 1); // no geometrical information
	z_c_bpack_set_I_option(&option_mv, "xyzsort", 0);// natural ordering


	z_c_bf_construct_init(&M, &L, &myseg_m, &myseg_l, nns_ptr_m, nns_ptr_l, &msh0_m, &msh0_l, &bf_mv, &option_mv, &stats_mv, &msh_mv, &kerquant_mv, &ptree_mv,&C_FuncDistmn_dummy, &C_FuncNearFar_dummy, quant_ptr_mv);
	z_c_bf_construct_matvec_compute(&bf_mv, &option_mv, &stats_mv, &msh_mv, &kerquant_mv, &ptree_mv, &C_FuncBMatVec, quant_ptr_mv);

	if(myrank==master_rank)std::cout<<"\nPrinting stats of the fourth BF: "<<std::endl;
	z_c_bpack_printstats(&stats_mv,&ptree_mv);

	z_c_bpack_deletestats(&stats_mv);
	z_c_bpack_deleteproctree(&ptree_mv);
	z_c_bpack_deletekernelquant(&kerquant_mv);
	z_c_bpack_deleteoption(&option_mv);
	z_c_bf_deletebf(&bf_mv);

	delete quant_ptr_mv;
	z_c_bpack_deletestats(&stats_dummy);
	z_c_bpack_deletestats(&stats_a);
	z_c_bpack_deletestats(&stats_b);
	z_c_bpack_deletestats(&stats_c);
	z_c_bpack_deleteproctree(&ptree);
	z_c_bpack_deletemesh(&msh_a);
	z_c_bpack_deletemesh(&msh_b);
	z_c_bpack_deletemesh(&msh_c);
	z_c_bpack_deletekernelquant(&kerquant_a);
	z_c_bpack_deletekernelquant(&kerquant_b);
	z_c_bpack_deletekernelquant(&kerquant_c);
	z_c_bpack_delete(&bmat_dummy);
	z_c_bpack_deletekernelquant(&kerquant_dummy);
	z_c_bpack_deleteoption(&option);
	z_c_bf_deletebf(&bf_a);
	z_c_bf_deletebf(&bf_b);
	z_c_bf_deletebf(&bf_c);

	z_c_bpack_deletemesh(&msh0_m);
	z_c_bpack_deletemesh(&msh0_n);
	z_c_bpack_deletemesh(&msh0_k);
	z_c_bpack_deletemesh(&msh0_l);

	delete quant_ptr_a;
	delete quant_ptr_b;
	delete quant_ptr_c;
	delete[] perms_m;
	delete[] perms_n;
	delete[] perms_k;
	delete[] perms_l;
	delete[] tree_m;
	delete[] tree_n;
	delete[] tree_k;
	delete[] tree_l;
}else{
	z_c_bpack_deletestats(&stats_dummy);
	z_c_bpack_deletestats(&stats_a);
	z_c_bpack_deletestats(&stats_b);
	z_c_bpack_deleteproctree(&ptree);
	z_c_bpack_deletemesh(&msh_a);
	z_c_bpack_deletemesh(&msh_b);
	z_c_bpack_deletekernelquant(&kerquant_a);
	z_c_bpack_deletekernelquant(&kerquant_b);
	z_c_bpack_delete(&bmat_dummy);
	z_c_bpack_deletekernelquant(&kerquant_dummy);
	z_c_bpack_deleteoption(&option);
	z_c_bf_deletebf(&bf_a);
	z_c_bf_deletebf(&bf_b);

	z_c_bpack_deletemesh(&msh0_m);
	z_c_bpack_deletemesh(&msh0_n);
	z_c_bpack_deletemesh(&msh0_k);

	delete quant_ptr_a;
	delete quant_ptr_b;
	delete[] perms_m;
	delete[] perms_n;
	delete[] perms_k;
	delete[] tree_m;
	delete[] tree_n;
	delete[] tree_k;
}



	Cblacs_exit(1);
	MPI_Finalize();                                 // Terminate MPI. Once called, no other MPI routines may be called
    return 0;
}
//------------------------------------------------------------------------------


