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
 * @brief This c++ example uses LR to compute the rank of surface and volumetric Green's function interactions between a pair of source and obvervation geometries, no butterfly algorithm is used.
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
  vector<double> _data_m;
  vector<double> _data_n;
  int _d = 0;
  int _n = 0;
  int _ker=1; // 1. coplannar unit plates 2. parallel unitplates 3. two unit cubes
  int _ppw=5; // points per wavelength
  double _wavelen=1.0; // wavelength
  int _Nperdim=100; //number of points per dimension

  int _rank_rand;
  int _m_rand;
  int _n_rand;
  int _ni1_rand;
  int _ni2_rand;
  std::vector<_Complex double> _MatU;
  std::vector<_Complex double> _MatV;

  std::vector<int> _Hperm_m;
  std::vector<int> _iHperm_m;
  std::vector<int> _Hperm_n;
  std::vector<int> _iHperm_n;
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

  C_QuantApp(double wavelen, int ppw, int ker, int Nperdim, int Ndim)
    : _wavelen(wavelen), _ppw(ppw), _ker(ker),_Nperdim(Nperdim),_d(Ndim){
	}


  inline void Sample(double pos_o[3], double pos_s[3], _Complex double* val){


    double dist = sqrt(pow(pos_o[0]-pos_s[0],2)+pow(pos_o[1]-pos_s[1],2)+pow(pos_o[2]-pos_s[2],2));
    double waven=2*M_PI/_wavelen;
    *val = (cos(-waven*dist)+Im*sin(-waven*dist))/dist;
  }
};


// The sampling function wrapper required by the Fortran HODLR code
inline void C_FuncZmn(int *m, int *n, _Complex double *val, C2Fptr quant) {

  C_QuantApp* Q = (C_QuantApp*) quant;

  double pos_o[3];
  pos_o[0] = Q->_data_m[(*m-1) * Q->_d];
  pos_o[1] = Q->_data_m[(*m-1) * Q->_d+1];
  pos_o[2] = Q->_data_m[(*m-1) * Q->_d+2];

  double pos_s[3];
  pos_s[0] = Q->_data_n[(*n-1) * Q->_d];
  pos_s[1] = Q->_data_n[(*n-1) * Q->_d+1];
  pos_s[2] = Q->_data_n[(*n-1) * Q->_d+2];


  Q->Sample(pos_o,pos_s,val);

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

  // m1 and n1 still need to convert to the original order, using new2old of msh0_m and msh_bf
  int m1_1, n1_1;
  m1_1 = Q->_Hperm_m[m1-1];
  n1_1 = Q->_Hperm_n[n1-1];
  double pos_o[3];
  pos_o[0] = Q->_data_m[(m1_1-1) * Q->_d];
  pos_o[1] = Q->_data_m[(m1_1-1) * Q->_d+1];
  pos_o[2] = Q->_data_m[(m1_1-1) * Q->_d+2];

  double pos_s[3];
  pos_s[0] = Q->_data_n[(n1_1-1) * Q->_d];
  pos_s[1] = Q->_data_n[(n1_1-1) * Q->_d+1];
  pos_s[2] = Q->_data_n[(n1_1-1) * Q->_d+2];


  Q->Sample(pos_o,pos_s,val);
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
	int Ndim=3; //data dimension
	double starttime, endtime;
	double* dat_ptr_m, *dat_ptr_n, *dat_ptr_k, *dat_ptr_l;
	int* nns_ptr_m, *nns_ptr_n, *nns_ptr_k, *nns_ptr_l;
	int nogeo;  // if 1: no geometrical information passed to hodlr, dat_ptr and Ndim are dummy; if 0: geometrical information passed
	int ker=0 ; // kernel choice

	int Nmin=16; //finest leafsize
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

  int tst = 1; // 1. coplannar unit plates 2. parallel unit plates 3. two unit cubes
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
  int ppw=5;
  int Nperdim;



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
      {"ppw",          required_argument, 0, 6},
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
      iss >> ppw;
    } break;
    default: break;
    }
  }


  for(int t=5; t<20;t++){
    double wavelen = 0.25/pow(2,t/2.0);
    double ds = wavelen/ppw;
    vector<double> data_geo_m;
    vector<double> data_geo_n;
    if(tst==1){ // two colinear plate
      Nperdim = ceil(1.0/ds);
      M = Nperdim*Nperdim;
      N = Nperdim*Nperdim;
      data_geo_m.resize(M*Ndim);
      data_geo_n.resize(N*Ndim);

      for(int m=0;m<M;m++){
        int ii = (m+1)%Nperdim+1;
	      int jj = (m+1)/Nperdim;
        data_geo_m[(m) * Ndim]=ii*ds+2;
        data_geo_m[(m) * Ndim+1]=jj*ds;
        data_geo_m[(m) * Ndim+2]=0;
      }
      for(int n=0;n<N;n++){
        int ii = (n+1)%Nperdim+1;
	      int jj = (n+1)/Nperdim;
        data_geo_n[(n) * Ndim]=ii*ds;
        data_geo_n[(n) * Ndim+1]=jj*ds;
        data_geo_n[(n) * Ndim+2]=0;
      }

    }else if(tst==2){ // two parallel plate
      Nperdim = ceil(1.0/ds);
      M = Nperdim*Nperdim;
      N = Nperdim*Nperdim;
      data_geo_m.resize(M*Ndim);
      data_geo_n.resize(N*Ndim);

      for(int m=0;m<M;m++){
        int ii = (m+1)%Nperdim+1;
	      int jj = (m+1)/Nperdim;
        data_geo_m[(m) * Ndim]=0;
        data_geo_m[(m) * Ndim+1]=ii*ds;
        data_geo_m[(m) * Ndim+2]=jj*ds;
      }
      for(int n=0;n<N;n++){
        int ii = (n+1)%Nperdim+1;
	      int jj = (n+1)/Nperdim;
        data_geo_n[(n) * Ndim]=2;
        data_geo_n[(n) * Ndim+1]=ii*ds;
        data_geo_n[(n) * Ndim+2]=jj*ds;
      }
    }else if(tst==3){ // two 3D cubes
      Nperdim = ceil(1.0/ds);
      M = Nperdim*Nperdim*Nperdim;
      N = Nperdim*Nperdim*Nperdim;
      data_geo_m.resize(M*Ndim);
      data_geo_n.resize(N*Ndim);
      for(int m=0;m<M;m++){
        int ii = ((m+1)%(int)pow(Nperdim,2)+1)%Nperdim+1;
        int jj = ((m+1)%(int)pow(Nperdim,2)+1)/Nperdim;
        int kk = (m+1)/pow(Nperdim,2);
        data_geo_m[(m) * Ndim]=ii*ds;
        data_geo_m[(m) * Ndim+1]=jj*ds;
        data_geo_m[(m) * Ndim+2]=kk*ds;
      }
      for(int n=0;n<N;n++){
        int ii = ((n+1)%(int)pow(Nperdim,2)+1)%Nperdim+1;
        int jj = ((n+1)%(int)pow(Nperdim,2)+1)/Nperdim;
        int kk = (n+1)/pow(Nperdim,2);
        data_geo_n[(n) * Ndim]=ii*ds+2;
        data_geo_n[(n) * Ndim+1]=jj*ds;
        data_geo_n[(n) * Ndim+2]=kk*ds;
      }
    }
    quant_ptr_a=new C_QuantApp(wavelen, ppw, tst, Nperdim, Ndim);
    quant_ptr_a->_data_m.resize(Ndim*M);
    quant_ptr_a->_data_m = data_geo_m;
    quant_ptr_a->_data_n.resize(Ndim*N);
    quant_ptr_a->_data_n = data_geo_n;

    if(myrank==master_rank)std::cout<<"tst "<<tst<<" M "<<M<<" N "<<N<<std::endl;

    nogeo=0;
    sort_opt=2;

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
    z_c_bpack_set_I_option(&option, "verbosity", 1);
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
    z_c_bpack_construct_init(&M, &Ndim, data_geo_m.data(), nns_ptr_m,&nlevel_m, tree_m, perms_m, &myseg_m, &bmat_dummy, &option, &stats_dummy, &msh0_m, &kerquant_dummy, &ptree, &C_FuncDistmn_dummy, &C_FuncNearFar_dummy, quant_ptr_dummy);
    quant_ptr_a->_Hperm_m.resize(M);
    std::copy(perms_m, perms_m + M, quant_ptr_a->_Hperm_m.begin());

    int nlevel_n = 0; // 0: tree level, nonzero if a tree is provided
    int* tree_n = new int[(int)pow(2,nlevel_n)]; //user provided array containing size of each leaf node, not used if nlevel=0
    tree_n[0] = N;
    z_c_bpack_construct_init(&N, &Ndim, data_geo_n.data(), nns_ptr_n,&nlevel_n, tree_n, perms_n, &myseg_n, &bmat_dummy, &option, &stats_dummy, &msh0_n, &kerquant_dummy, &ptree, &C_FuncDistmn_dummy, &C_FuncNearFar_dummy, quant_ptr_dummy);
    quant_ptr_a->_Hperm_n.resize(N);
    std::copy(perms_n, perms_n + N, quant_ptr_a->_Hperm_n.begin());

  // construct the three bfs from entry evaluation
    z_c_bpack_createstats(&stats_a);
    z_c_bf_construct_init(&M, &N, &myseg_m, &myseg_n, nns_ptr_m, nns_ptr_n, &msh0_m, &msh0_n, &bf_a, &option, &stats_a, &msh_a, &kerquant_a, &ptree,&C_FuncDistmn_dummy, &C_FuncNearFar_dummy, quant_ptr_a);
    z_c_bf_construct_element_compute(&bf_a, &option, &stats_a, &msh_a, &kerquant_a, &ptree, &C_FuncBZmn, &C_FuncBZmnBlock, quant_ptr_a); // C_FuncBZmnBlock is not referenced since elem_extract=0


    if(myrank==master_rank)std::cout<<"\nPrinting stats of the block: "<<std::endl;
    z_c_bpack_printstats(&stats_a,&ptree);

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

  }



	Cblacs_exit(1);
	MPI_Finalize();                                 // Terminate MPI. Once called, no other MPI routines may be called
    return 0;
}
//------------------------------------------------------------------------------


