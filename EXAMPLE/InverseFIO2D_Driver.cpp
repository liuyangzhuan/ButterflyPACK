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
 * @brief This c++ example compresses a 2D Fourier integral operators (FIO) as a butterfly using entry-evaluation-based APIs, then constructs a preconditioner for the inverse FIO using matvec-based APIs. The example works on the double-complex data type.
*/
//------------------------------------------------------------------------------
#include <iostream>
#include <math.h>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <iomanip>
#include <memory>

#include <pthread.h>

#include <cmath>
#include <cassert>
#include <iostream>
#include <random>
#include <vector>
#include <atomic>
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
char const transN='N';
char const transT='T';

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
  int _ker=1;
  int _m_rand;
  int _n_rand;
  int _d=1;

  vector<double> _data_m;
  vector<double> _data;

  std::vector<int> _Hperm;
  std::vector<int> _iHperm;

  std::vector<int> _Hperm_m;
  std::vector<int> _iHperm_m;

  int _nloc = 0;
  int _mloc = 0;

  F2Cptr* bmat;  //hierarchical matrix returned by Fortran code
  F2Cptr* bf_a;  //BF returned by Fortran code
  F2Cptr* stats_a;      //statistics structure returned by Fortran code
  F2Cptr* msh_a;		   //mesh structure returned by Fortran code
  F2Cptr* ptree_a;      //process tree returned by Fortran code
  F2Cptr* option_a;      //option structure returned by Fortran code


  C_QuantApp() = default;

  C_QuantApp(int m_rand, int n_rand, int d, vector<double> data_m, vector<double> data, int ker)
    : _m_rand(m_rand), _n_rand(n_rand), _d(d), _data(move(data)),_data_m(move(data_m)),_ker(ker){
	}
  C_QuantApp(int m_rand, int n_rand, int mloc, int nloc, int ker)
    : _mloc(mloc),_nloc(nloc),_m_rand(m_rand), _n_rand(n_rand),_ker(ker){
	}


  inline void Sample(int m, int n, _Complex double* val){

    if(_ker==1){

        // m still need to convert to the original order, using new2old of rows
        int m_1;
        m_1 = _Hperm_m[m];
        double k1 = _data_m[(m_1-1) * _d];
        double k2 = _data_m[(m_1-1) * _d+1];

        // n still need to convert to the original order, using new2old of cols
        int n_1;
        n_1 = _Hperm[n];
        double x1 = _data[(n_1-1) * _d];
        double x2 = _data[(n_1-1) * _d+1];

        double xk = x1*k1+x2*k2;
        double tmp = (2*M_PI)* (xk);
        *val = cos(tmp)+Im*sin(tmp);

		}else if(_ker==2){
        // m still need to convert to the original order, using new2old of rows
        int m_1;
        m_1 = _Hperm_m[m];
        double k1 = _data_m[(m_1-1) * _d];
        double k2 = _data_m[(m_1-1) * _d+1];

        // n still need to convert to the original order, using new2old of cols
        int n_1;
        n_1 = _Hperm[n];
        double x1 = _data[(n_1-1) * _d];
        double x2 = _data[(n_1-1) * _d+1];

        double xk = x1*k1+x2*k2;
        double sx = (2+sin(2*M_PI*x1)*sin(2*M_PI*x2))/16.0;
        double cx = (2+cos(2*M_PI*x1)*cos(2*M_PI*x2))/16.0;
        double kr = sqrt(pow(sx,2)*pow(k1,2) + pow(cx,2)*pow(k2,2));
        double tmp = (2*M_PI)* (xk + kr);
        *val = cos(tmp)+Im*sin(tmp);
    }
    else{
      cout<<"unsupported kernel"<<endl;
		}
  }
};


// The sampling function wrapper required by the Fortran HODLR code
inline void C_FuncZmn(int *m, int *n, _Complex double *val, C2Fptr quant) {

  C_QuantApp* Q = (C_QuantApp*) quant;
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
inline void C_FuncZmnBlock(int* Ninter, int* Nallrows, int* Nallcols, int64_t* Nalldat_loc, int* allrows, int* allcols, _Complex double* alldat_loc, int* rowidx,int* colidx, int* pgidx, int* Npmap, int* pmaps, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;
}

// The extraction sampling function wrapper required by the Fortran HODLR code
inline void C_FuncBZmnBlock(int* Ninter, int* Nallrows, int* Nallcols, int64_t* Nalldat_loc, int* allrows, int* allcols, _Complex double* alldat_loc, int* rowidx,int* colidx, int* pgidx, int* Npmap, int* pmaps, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;

//   z_c_bf_extractelement(Q->bf,Q->option,Q->msh,Q->stats,Q->ptree,Ninter,Nallrows, Nallcols, Nalldat_loc, allrows,allcols,alldat_loc,rowidx,colidx,pgidx,Npmap,pmaps);

}


// The matvec sampling function wrapper required by the Fortran HODLR code
inline void C_FuncHMatVec(char const *trans, int *nin, int *nout, int *nvec, _Complex double const *xin, _Complex double *xout, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;

  int cnt = (*nvec)*(*nout);
  _Complex double* xin1 = new _Complex double[cnt];
  _Complex double* xbuf1 = new _Complex double[(*nvec)*(Q->_mloc)];



  if(*trans=='N'){
    z_c_bf_mult(trans, xin, xbuf1, nin, &(Q->_mloc), nvec, Q->bf_a,Q->option_a,Q->stats_a,Q->ptree_a);
    for (int ii=0; ii<(*nvec)*(Q->_mloc); ii++)
      xbuf1[ii]=conj(xbuf1[ii]);
    z_c_bf_mult(&transT, xbuf1, xout, &(Q->_mloc), nout, nvec, Q->bf_a,Q->option_a,Q->stats_a,Q->ptree_a);
    for (int ii=0; ii<cnt; ii++)
      xout[ii]=conj(xout[ii]);
  }else if(*trans=='T'){
    for (int ii=0; ii<cnt; ii++)
      xin1[ii]=conj(xin[ii]);
    z_c_bf_mult(&transN, xin1, xbuf1, nin, &(Q->_mloc), nvec, Q->bf_a,Q->option_a,Q->stats_a,Q->ptree_a);
    for (int ii=0; ii<(*nvec)*(Q->_mloc); ii++)
      xbuf1[ii]=conj(xbuf1[ii]);
    z_c_bf_mult(trans, xbuf1, xout, &(Q->_mloc), nout, nvec, Q->bf_a,Q->option_a,Q->stats_a,Q->ptree_a);
  }

  delete[] xbuf1;
  delete[] xin1;

}


// The matvec sampling function wrapper required by the Fortran HODLR code
inline void C_FuncBMatVec(char const *trans, int *nin, int *nout, int *nvec, _Complex double const *xin, _Complex double *xout, C2Fptr quant, _Complex double *a, _Complex double *b) {
  C_QuantApp* Q = (C_QuantApp*) quant;
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
        z_c_bpack_set_D_option(&option0, "tol_comp", opt_d*0.1);
        z_c_bpack_set_D_option(&option0, "tol_rand", opt_d);
        z_c_bpack_set_D_option(&option0, "tol_Rdetect", opt_d*3e-1);
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
        iss >> opt_d;
        z_c_bpack_set_D_option(&option0, "rankrate", opt_d);
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
	int Ndim=2; //data dimension
	double starttime, endtime;
	double* dat_ptr_m, *dat_ptr_n, *dat_ptr_k, *dat_ptr_l;
	int* nns_ptr_m, *nns_ptr_n, *nns_ptr_k, *nns_ptr_l;
	int nogeo;  // 1: no geometrical information passed to hodlr, dat_ptr and Ndim are dummy
	int ker=2 ; // kernel choice

	int Nmin=64; //finest leafsize
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

  int tst = 2;
	int lrlevel=100;

if(myrank==master_rank){
	z_c_bpack_getversionnumber(&v_major,&v_minor,&v_bugfix);
	std::cout<<"ButterflyPACK Version: "<<v_major<<"."<<v_minor<<"."<<v_bugfix<<std::endl;
}



	/*****************************************************************/
	Npo = 1000;
	int M = Npo;
	int N = Npo;



  //getting the example configurations from command line
  std::vector<std::unique_ptr<char[]>> argv_data(argc);
  std::vector<char*> argv1(argc);
  for (int i=0; i<argc; i++) {
    argv_data[i].reset(new char[strlen(argv[i])+1]);
    argv1[i] = argv_data[i].get();
    strcpy(argv1[i], argv[i]);
  }
  option long_options[] =
    {{"M",                   required_argument, 0, 1},
      {"N",                   required_argument, 0, 2},
      {"ker",             required_argument, 0, 3},
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
      iss >> M;
    } break;
    case 2: {
      std::istringstream iss(optarg);
      iss >> N;
    } break;
    case 3: {
      std::istringstream iss(optarg);
      iss >> ker;
    } break;
    default: break;
    }
  }



  if(myrank==master_rank)std::cout<<"M "<<M<<" N "<<N<<std::endl;



	nogeo=0;
	sort_opt=1;

	/*****************************************************************/
	int myseg_m=0;     // local number of unknowns
	int myseg_n=0;     // local number of unknowns
	int* perms_m = new int[M]; //permutation vector returned by HODLR
	int* perms_n = new int[N]; //permutation vector returned by HODLR
	int* groups = new int[size];
	int i_opt;
	double d_opt;
	int cpp=1; //1: use user-defined cpp/c functions for construction
	//quantities for the first holdr
	F2Cptr bmat_dummy;  //hierarchical matrix returned by Fortran code
	F2Cptr option;     //option structure returned by Fortran code
	F2Cptr stats_dummy,stats_a;      //statistics structure returned by Fortran code
	F2Cptr msh_a;		   //z_mesh structure returned by Fortran code
	F2Cptr msh0_m, msh0_n;		   //z_mesh structure returned by Fortran code
	F2Cptr kerquant_a;   //kernel quantities structure returned by Fortran code
	F2Cptr kerquant_dummy;   //kernel quantities structure returned by Fortran code
	F2Cptr ptree;      //process tree returned by Fortran code

	F2Cptr bf_a;  //BF returned by Fortran code

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
	z_c_bpack_set_I_option(&option, "ErrSol", 0);

  set_option_from_command_line(argc, argv, option);
	z_c_bpack_printoption(&option,&ptree);

  double tmp;
  z_c_bpack_getoption(&option, "format", &tmp);
  int format_temp=(int)tmp; // always set format =1 for the BF
  z_c_bpack_set_I_option(&option, "format", 1);


  vector<double> data_geo_m;
  int Ms = int(sqrt(M));
  data_geo_m.resize(Ndim*M);
  for(int ii=0;ii<M;ii++){
    int idx_x = ii%Ms;
    int idx_y = ii/Ms;
    data_geo_m[(ii) * Ndim] = idx_x-Ms/2.0;
    data_geo_m[(ii) * Ndim+1] = idx_y-Ms/2.0;
  }
  vector<double> data_geo_n;
  int Ns = int(sqrt(N));
  data_geo_n.resize(Ndim*N);
  for(int ii=0;ii<N;ii++){
    int idx_x = ii%Ns;
    int idx_y = ii/Ns;
    data_geo_n[(ii) * Ndim] = idx_x/((double)Ns);
    data_geo_n[(ii) * Ndim+1] = idx_y/((double)Ns);
  }
  quant_ptr_a=new C_QuantApp(M, N, Ndim, data_geo_m, data_geo_n, ker);


  // construct the mesh data for the bf, this should be improved
	int nlevel_m = 0; // 0: tree level, nonzero if a tree is provided
	int* tree_m = new int[(int)pow(2,nlevel_m)]; //user provided array containing size of each leaf node, not used if nlevel=0
	tree_m[0] = M;
	z_c_bpack_construct_init(&M, &Ndim, data_geo_m.data(), nns_ptr_m,&nlevel_m, tree_m, perms_m, &myseg_m, &bmat_dummy, &option, &stats_dummy, &msh0_m, &kerquant_dummy, &ptree, &C_FuncDistmn_dummy, &C_FuncNearFar_dummy, quant_ptr_dummy);
  quant_ptr_a->_Hperm_m.resize(M);
  std::copy(perms_m, perms_m + M, quant_ptr_a->_Hperm_m.begin());

	int nlevel_n = 0; // 0: tree level, nonzero if a tree is provided
	int* tree_n = new int[(int)pow(2,nlevel_n)]; //user provided array containing size of each leaf node, not used if nlevel=0
	tree_n[0] = N;
	z_c_bpack_construct_init(&N, &Ndim, data_geo_n.data(), nns_ptr_n,&nlevel_n, tree_n, perms_n, &myseg_n, &bmat_dummy, &option, &stats_dummy, &msh0_n, &kerquant_dummy, &ptree, &C_FuncDistmn_dummy, &C_FuncNearFar_dummy, quant_ptr_dummy);
  quant_ptr_a->_Hperm.resize(N);
  std::copy(perms_n, perms_n + N, quant_ptr_a->_Hperm.begin());

//////////////////////  construct the bf of FIO (A) from entry evaluation
	z_c_bpack_createstats(&stats_a);
	z_c_bf_construct_init(&M, &N, &myseg_m, &myseg_n, nns_ptr_m, nns_ptr_n, &msh0_m, &msh0_n, &bf_a, &option, &stats_a, &msh_a, &kerquant_a, &ptree,&C_FuncDistmn_dummy, &C_FuncNearFar_dummy, quant_ptr_a);
	z_c_bf_construct_element_compute(&bf_a, &option, &stats_a, &msh_a, &kerquant_a, &ptree, &C_FuncBZmn, &C_FuncBZmnBlock, quant_ptr_a); // C_FuncBZmnBlock is not referenced since elem_extract=0
	if(myrank==master_rank)std::cout<<"\nPrinting stats of the first BF: "<<std::endl;
	z_c_bpack_printstats(&stats_a,&ptree);


	//////////////////// use resulting A^*A as matvec to create a new hodlr

	// quantities for the second holdr
	F2Cptr bmat1;  //hierarchical matrix returned by Fortran code
	F2Cptr option1;     //option structure returned by Fortran code
	F2Cptr stats1;      //statistics structure returned by Fortran code
	F2Cptr msh1;		   //d_mesh structure returned by Fortran code
	F2Cptr kerquant1;   //kernel quantities structure returned by Fortran code
	F2Cptr ptree1;      //process tree returned by Fortran code
	C_QuantApp *quant_ptr1; //user-defined object

	quant_ptr1=new C_QuantApp(myseg_m,myseg_n,M, N, ker);
	quant_ptr1->bf_a=&bf_a;
	quant_ptr1->msh_a=&msh_a;
	quant_ptr1->ptree_a=&ptree;
	quant_ptr1->stats_a=&stats_a;
	quant_ptr1->option_a=&option;

	z_c_bpack_createptree(&size, groups, &Fcomm, &ptree1);
	z_c_bpack_copyoption(&option,&option1);
	z_c_bpack_createstats(&stats1);

	z_c_bpack_set_I_option(&option1, "format", format_temp);// HODLR or H format
	z_c_bpack_set_I_option(&option1, "LRlevel", 0);// LR format
	z_c_bpack_set_I_option(&option1, "per_geo", 1);// periodic geometry points
	z_c_bpack_set_D_option(&option1, "period1", 1e0);// period in the first dimension
	z_c_bpack_set_D_option(&option1, "period2", 1e0);// period in the second dimension

  // tol=1e-4;
	// z_c_bpack_set_D_option(&option1, "tol_comp", tol);
	// z_c_bpack_set_D_option(&option1, "tol_rand", tol);           // bf_mv uses this tolerance
	// z_c_bpack_set_D_option(&option1, "tol_Rdetect", tol*3e-1);   // bf_mv uses this tolerance

	// z_c_bpack_set_I_option(&option1, "nogeo", 1); // no geometrical information
	// z_c_bpack_set_I_option(&option1, "xyzsort", 0);// natural ordering


  z_c_bpack_printoption(&option1,&ptree);

	int Npo1 = N;
	int myseg1=0;     // local number of unknowns
	int* perms1 = new int[Npo1]; //permutation vector returned by HODLR
	//tree1 and nlevel1 should be provided by the caller, otherwise natural ordering is used
	int nlevel1 = 0; // 0: tree level, nonzero if a tree is provided
	int* tree1 = new int[(int)pow(2,nlevel1)]; //user provided array containing size of each leaf node, not used if nlevel=0
	tree1[0] = Npo1;
	// int Ndim1=0; //data dimension
	// double* dat_ptr1;

	z_c_bpack_construct_init(&Npo1, &Ndim, data_geo_n.data(), nns_ptr_n,&nlevel1, tree1, perms1, &myseg1, &bmat1, &option1, &stats1, &msh1, &kerquant1, &ptree1, &C_FuncDistmn_dummy, &C_FuncNearFar_dummy, quant_ptr1);
	z_c_bpack_construct_matvec_compute(&bmat1, &option1, &stats1, &msh1, &kerquant1, &ptree1, &C_FuncHMatVec, quant_ptr1);


	// factor hodlr
	z_c_bpack_factor(&bmat1,&option1,&stats1,&ptree,&msh1);

	// solve the system
	int nrhs=1;
	_Complex double* b = new _Complex double[nrhs*myseg_m];
	_Complex double* xtrue = new _Complex double[nrhs*myseg_n];
	_Complex double* xbuf = new _Complex double[nrhs*myseg_n];
	_Complex double* x = new _Complex double[nrhs*myseg_n];

	//////////////////// Generate a true solution xtrue, and its rhs b using A
	for (int i = 0; i < nrhs*myseg_n; i++){
		xtrue[i]=1;
	}
  z_c_bf_mult(&transN, xtrue, b, &myseg_n, &myseg_m, &nrhs, &bf_a, &option, &stats_a, &ptree);


  //////////////////// Generate an approximate solution using (A^*A)^-1A^*b
	for (int i = 0; i < nrhs*myseg_m; i++)
		b[i]=conj(b[i]);
  z_c_bf_mult(&transT, b, xbuf, &myseg_m, &myseg_n, &nrhs, &bf_a, &option, &stats_a, &ptree);
	for (int i = 0; i < nrhs*myseg_n; i++)
		xbuf[i]=conj(xbuf[i]);
	z_c_bpack_solve(x,xbuf,&myseg_n,&nrhs,&bmat1,&option1,&stats1,&ptree);

  double norm1=0,norm2=0,norm1t=0,norm2t=0;
  for (int i = 0; i < nrhs*myseg_n; i++){
		norm1+=pow(cabs(xtrue[i]),2);
		norm2+=pow(cabs(xtrue[i]-x[i]),2);
	}

  MPI_Reduce(&norm1,&norm1t,1,MPI_DOUBLE,MPI_MAX,master_rank,MPI_COMM_WORLD);
  MPI_Reduce(&norm2,&norm2t,1,MPI_DOUBLE,MPI_MAX,master_rank,MPI_COMM_WORLD);


  if(myrank==master_rank)std::cout<<"|x-xtrue|/|xtrue| for inverse FIO: "<< sqrt(norm2t)/sqrt(norm1t) <<std::endl;

	if(myrank==master_rank)std::cout<<"Printing stats of the second Bmat: "<<std::endl;
	z_c_bpack_printstats(&stats1,&ptree1);

	z_c_bpack_deletestats(&stats1);
	z_c_bpack_deleteproctree(&ptree1);
	z_c_bpack_deletemesh(&msh1);
	z_c_bpack_deletekernelquant(&kerquant1);
	z_c_bpack_delete(&bmat1);
	z_c_bpack_deleteoption(&option1);
	delete quant_ptr1;
	delete[] perms1;
	delete[] tree1;


	z_c_bpack_deletestats(&stats_dummy);
	z_c_bpack_deletestats(&stats_a);
	z_c_bpack_deleteproctree(&ptree);
	z_c_bpack_deletemesh(&msh_a);
	z_c_bpack_deletekernelquant(&kerquant_a);
	z_c_bpack_delete(&bmat_dummy);
	z_c_bpack_deletekernelquant(&kerquant_dummy);
	z_c_bpack_deleteoption(&option);
	z_c_bf_deletebf(&bf_a);

	z_c_bpack_deletemesh(&msh0_m);
	z_c_bpack_deletemesh(&msh0_n);

	delete quant_ptr_a;
	delete[] perms_m;
	delete[] perms_n;
	delete[] tree_m;
	delete[] tree_n;


	Cblacs_exit(1);
	MPI_Finalize();                                 // Terminate MPI. Once called, no other MPI routines may be called
    return 0;
}
//------------------------------------------------------------------------------


