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
 * @brief This c++ driver illustrates the c++ interface to ButterflyPACK's Fortran subroutines for compression a kernel ridge regression (KRR) matrix, with the entry-evaluation-based APIs. This file works on the double data type. This file is a simplified version of InterfaceTest.cpp which provides more examples.
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
#include <sstream>
#include <cstring>
#include <getopt.h>
#include <unistd.h>

#include "dC_BPACK_wrapper.h"



//------------------------------------------------------------------------------
using namespace std;


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

// 2-norm distance
inline double dist2(double *x, double *y, int d) {
  double k = 0.;
  for (int i = 0; i < d; i++)
    k += pow(x[i] - y[i], 2.);
  return k;
}

// dot product of two real vectors
inline double dot_product(double* v, double* u, int d)
{
    double result = 0.0;
    for (int i = 0; i < d; i++)
        result += v[i]*u[i];
    return result;
}


// Gauss Kernel
inline double Gauss_kernel(double *x, double *y, int d, double h) {
  double dists;
  dists = dist2(x, y, d);
  if(dists> -log(1e-30)*2.0*pow(h,2.0)){
	return 0;
  }else{
	return exp(-dists / (2. * h * h));
  }
}

//R^4 kernel
inline double K07_kernel(double *x, double *y, int d) {
  double dists;
  dists = dist2(x, y, d);
  return pow(dists,4);
}

// sqrt(R^2+h) kernel
inline double K08_kernel(double *x, double *y, int d, double h) {
  double dists;
  dists = dist2(x, y, d);
  return sqrt(pow(dists,2)+h);
}

// 1/sqrt(R^2+h) kernel
inline double K09_kernel(double *x, double *y, int d, double h) {
  double dists;
  dists = dist2(x, y, d);
  return 1.0/sqrt(pow(dists,2)+h);
}

// Polynomial kernel (X^tY+h)^2
inline double K10_kernel(double *x, double *y, int d, double h) {
  double dotp;
  dotp = dot_product(x, y, d);
  return pow(dotp+h,2);
}




// The object handling kernel parameters and sampling function
class C_QuantApp {
public:
  vector<double> _data;
  int _d = 0;
  int _n = 0;
  double _h = 0.;
  double _l = 0.;
  int _ker=1; //

//   int _rank_rand;
//   int _n_rand;
//   std::vector<double> _MatU;
//   std::vector<double> _MatV;
//   std::vector<double> _MatFull;

//   std::vector<int> _Hperm;
//   std::vector<int> _iHperm;
//   int _nloc = 0;

//   F2Cptr* bmat;  //hierarchical matrix returned by Fortran code
//   F2Cptr* bf;  //BF returned by Fortran code
//   F2Cptr* stats;      //statistics structure returned by Fortran code
//   F2Cptr* msh;		   //mesh structure returned by Fortran code
//   F2Cptr* ptree;      //process tree returned by Fortran code
//   F2Cptr* option;      //option structure returned by Fortran code


  C_QuantApp() = default;

  C_QuantApp(vector<double> data, int d, double h, double l, int ker)
    : _data(move(data)), _d(d), _n(_data.size() / _d),
      _h(h), _l(l),_ker(ker){
    assert(size_t(_n * _d) == _data.size());
	}


//   C_QuantApp(int n, int ker, vector<double> MatFull)
//     : _n(n), _ker(ker), _MatFull(move(MatFull)){
// 	// cout<<_n_rand<<_rank_rand<<_MatU.size()<<endl;
//     assert(size_t(_n * _n) == _MatFull.size());
// 	}



  inline void Sample(int m, int n, double* val){
	switch(_ker){
	case 1: //Gaussian kernel
		*val = Gauss_kernel(&_data[m * _d], &_data[n * _d], _d, _h);
		if (m==n)
		*val += _l;
		break;
	case 2: //R^4 kernel
		*val = K07_kernel(&_data[m * _d], &_data[n * _d], _d);
		break;
	case 3: //sqrt(R^2+h) kernel
		*val = K08_kernel(&_data[m * _d], &_data[n * _d], _d, _h);
		break;
	case 4: //1/sqrt(R^2+h) kernel
		*val = K09_kernel(&_data[m * _d], &_data[n * _d], _d, _h);
		break;
	case 5: //Polynomial kernel (X^tY+h)^2
		*val = K10_kernel(&_data[m * _d], &_data[n * _d], _d, _h);
		break;
	case 6: //Low-rank product of two random matrices
		// *val =0;
		// for (int k = 0; k < _rank_rand; k++){
		// 	*val += _MatU[k*_n_rand+m]*_MatV[k*_n_rand+n];
		// }
		break;
	case 7: //Full matrix
		// *val =_MatFull[n*_n+m];
		// // *val =_MatFull[_Hperm[n]*_n+_Hperm[m]];
		break;
	}
  }
};


// The sampling function wrapper required by the Fortran ButterflyPACK interface
// The function computes one entry A_{m,n}, where m,n starts at 1
inline void C_FuncZmn(int *m, int *n, double *val, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;
  Q->Sample(*m-1,*n-1,val);
}

// The distance function wrapper required by the Fortran ButterflyPACK interface
// not needed in this driver, but need to provide a dummy function to the interface
inline void C_FuncDistmn(int *m, int *n, double *val, C2Fptr quant) {
}

// The compressibility function wrapper required by the Fortran ButterflyPACK interface
// not needed in this driver, but need to provide a dummy function to the interface
inline void C_FuncNearFar(int *m, int *n, int *val, C2Fptr quant) {
}

// The extraction sampling function wrapper required by the Fortran ButterflyPACK interface
inline void C_FuncZmnBlock(int* Ninter, int* Nallrows, int* Nallcols, int64_t* Nalldat_loc, int* allrows, int* allcols, double* alldat_loc, int* rowidx,int* colidx, int* pgidx, int* Npmap, int* pmaps, C2Fptr quant) {
// not needed in this driver, but need to provide a dummy function to the interface
}


// Read a data file into a vector
template<typename T>
vector<T> write_from_file(string filename) {
  vector<T> data;
  ifstream f(filename);
  string l;
  while (getline(f, l)) {
    istringstream sl(l);
    string s;
    while (getline(sl, s, ','))
      data.push_back(stod(s));
  }
  return data;
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
      {{"nmin_leaf",                     required_argument, 0, 1},  // leafsize in the hierarchical partitioning
       {"tol_comp",                   required_argument, 0, 2},  // relative tolerance for matrix construction
       {"tol_rand",                   required_argument, 0, 3},  // relative tolerance for matrix inversion
       {"tol_Rdetect",             required_argument, 0, 4},     // relative tolerance for rank detection during matrix inversion
       {"tol_itersol",             required_argument, 0, 5},     // convergence tolerance for TFQMR iterative solver if precon=2 or 3
       {"n_iter",          required_argument, 0, 6},			 // maximum iteration count for TFQMR
       {"level_check",         required_argument, 0, 7},		 // the level in the hierarchical partitioning where the randomized construction algorithm is tested, set to 10000 by default (no checking)
       {"precon",                  required_argument, 0, 8}, 	// the use mode of butterflypack: 1: as a direct solver 2: as an iterative solver (compress the matrix and pass it to TFQMR without preconditioner), 3: as a preconditioned iterative solver (compress the matrix and invert the matrix and pass them to TFQMR, using approximate matrix inverse as a preconditioner)
       {"xyzsort",      required_argument, 0, 9}, 	// the hierarchical partitioning algorithm: 0: no permutation 1: permutation based on KD-tree 2: permutation based on cobble-like partitioning
       {"lrlevel",     required_argument, 0, 10},   // the level the hierarchical partitioning (top-down numbered) above which butterfly is used and below which low-rank is used
       {"errfillfull",       required_argument, 0, 11}, // errfillfull: a slow (n^2), thorough error checking is performed after the compression of each block
       {"baca_batch",      required_argument, 0, 12}, // block size in batched ACA when reclr_leaf=4 or 5
       {"reclr_leaf",      required_argument, 0, 13}, // low-rank compression algorithms 1:SVD 2:RRQR 3:ACA 4:BACA 5:BACA_improved 6:Pseudo-skeleton 7: ACA with naive parallelization
       {"nogeo",     required_argument, 0, 14}, // whether there is geometry information provided 1: is no geometry (xyzsort can not be 1 or 2), 0: there is geometry
       {"less_adapt",            required_argument, 0, 15}, // 1: improved randomized butterfly construction, default to 1
       {"errsol",           required_argument, 0, 16}, // 1: generating an artificial true solution vector, compute the RHS with compressed matrix, solve the system, and compare with true solution vector
       {"lr_blk_num",                  required_argument, 0, 17}, // sqrt of #of subblocks in H-BACA, default to 1
       {"rank0",  required_argument, 0, 18}, // initial rank guess in the randomized butterfly algorithm, default to 32
       {"rankrate", required_argument, 0, 19}, // increasing ratio of the rank guess in each iteration, default to 2
       {"itermax",               required_argument, 0, 20}, // maximum number of iterations in the randomized butterfly algorithm, default to 10
       {"powiter",  required_argument, 0, 21}, // order of power iteration in the randomized low-rank construction
       {"ilu", required_argument, 0, 22},	// whether symmetric gauss-seidel is used when format=2
       {"nbundle",     required_argument, 0, 23}, // multiply nbundle sets of vectors together in randomized butterfly algorithm for better flop performance, default to 1
       {"near_para",  required_argument, 0, 24}, // admissibility parameter when format=2, strong admissibility typically requires near_para>2.0
       {"format",  required_argument, 0, 25}, // the hiearchial matrix format: 1: HODLR/HODBF 2: H matrix 3: HSS-BF
       {"verbosity", required_argument, 0, 26}, // verbosity for the printing (-1, 0, 1, 2), -1 suppresses everything, 2 prints most details
       {"rmax", required_argument, 0, 27},	// preestimate of the maximum rank for allocating buffers, default to 1000
       {"sample_para", required_argument, 0, 28}, // oversampling factor in the nlogn entry-evaluation-based butterfly algorithm, default to 2
       {"pat_comp",    required_argument, 0, 29}, // pattern of entry-evaluation-based butterfly compression: 1 from right to left, 2 from left to right, 3 from outer to inner
       {"knn",         required_argument, 0, 30}, // nearest neighbouring points used in improved BACA and entry-evaluation-based butterfly compression
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
        d_c_bpack_set_I_option(&option0, "Nmin_leaf", opt_i);
      } break;
      case 2: {
        std::istringstream iss(optarg);
        iss >> opt_d;
        d_c_bpack_set_D_option(&option0, "tol_comp", opt_d);
        d_c_bpack_set_D_option(&option0, "tol_rand", opt_d);
        d_c_bpack_set_D_option(&option0, "tol_Rdetect", opt_d*0.1);
      } break;
      case 3: {
        std::istringstream iss(optarg);
        iss >> opt_d;
        d_c_bpack_set_D_option(&option0, "tol_rand", opt_d);
      } break;
      case 4: {
        std::istringstream iss(optarg);
        iss >> opt_d;
        d_c_bpack_set_D_option(&option0, "tol_Rdetect", opt_d);
      } break;
      case 5: {
        std::istringstream iss(optarg);
        iss >> opt_d;
        d_c_bpack_set_D_option(&option0, "tol_itersol", opt_d);
      } break;
      case 6: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        d_c_bpack_set_I_option(&option0, "n_iter", opt_i);
      } break;
      case 7: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        d_c_bpack_set_I_option(&option0, "level_check", opt_i);
      } break;
      case 8: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        d_c_bpack_set_I_option(&option0, "precon", opt_i);
      } break;
      case 9: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        d_c_bpack_set_I_option(&option0, "xyzsort", opt_i);
      } break;
      case 10: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        d_c_bpack_set_I_option(&option0, "LRlevel", opt_i);
      } break;
      case 11: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        d_c_bpack_set_I_option(&option0, "ErrFillFull", opt_i);
      } break;
      case 12: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        d_c_bpack_set_I_option(&option0, "BACA_Batch", opt_i);
      } break;
      case 13: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        d_c_bpack_set_I_option(&option0, "RecLR_leaf", opt_i);
      } break;
      case 14: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        d_c_bpack_set_I_option(&option0, "nogeo", opt_i);
      } break;
      case 15: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        d_c_bpack_set_I_option(&option0, "less_adapt", opt_i);
      } break;
      case 16: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        d_c_bpack_set_I_option(&option0, "ErrSol", opt_i);
      } break;
      case 17: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        d_c_bpack_set_I_option(&option0, "LR_BLK_NUM", opt_i);
      } break;
      case 18: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        d_c_bpack_set_I_option(&option0, "rank0", opt_i);
      } break;
      case 19: {
        std::istringstream iss(optarg);
        iss >> opt_d;
        d_c_bpack_set_D_option(&option0, "rankrate", opt_d);
      } break;
      case 20: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        d_c_bpack_set_I_option(&option0, "itermax", opt_i);
      } break;
      case 21: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        d_c_bpack_set_I_option(&option0, "powiter", opt_i);
      } break;
      case 22: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        d_c_bpack_set_I_option(&option0, "ILU", opt_i);
      } break;
      case 23: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        d_c_bpack_set_I_option(&option0, "Nbundle", opt_i);
      } break;
      case 24: {
        std::istringstream iss(optarg);
        iss >> opt_d;
        d_c_bpack_set_D_option(&option0, "near_para", opt_d);
      } break;
      case 25: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        d_c_bpack_set_I_option(&option0, "format", opt_i);
      } break;
      case 26: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        d_c_bpack_set_I_option(&option0, "verbosity", opt_i);
      } break;
      case 27: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        d_c_bpack_set_I_option(&option0, "rmax", opt_i);
      } break;
      case 28: {
        std::istringstream iss(optarg);
        iss >> opt_d;
        d_c_bpack_set_D_option(&option0, "sample_para", opt_d);
      } break;
      case 29: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        d_c_bpack_set_I_option(&option0, "pat_comp", opt_i);
      } break;
      case 30: {
        std::istringstream iss(optarg);
        iss >> opt_i;
        d_c_bpack_set_I_option(&option0, "knn", opt_i);
      } break;
      default: break;
      }
    }
  }

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
	int ker=1 ; // kernel choice
	int Npo=10000;   // matrix size
	int Ndim=8; //data dimension
	double starttime, endtime;
	double* dat_ptr;
	int* nns_ptr;
	int nogeo;  // 1: no geometrical information passed to butterflypack, dat_ptr and Ndim are dummy


	int Nmin=200; //finest leafsize
	double tol=1e-4; //compression tolerance
	double sample_para=2.0; //oversampling factor in entry evaluation
	int com_opt=5; //1:SVD 2:RRQR 3:ACA 4:BACA 5:BACA_improved 6:Pseudo-skeleton
	int sort_opt=1; //0:natural order 1:kd-tree 2:cobble-like ordering 3:gram distance-based cobble-like ordering
	int checkerr = 0; //1: check compression quality
	int batch = 16; //batch size for BACA
	int bnum = 1; //sqrt of #of subblocks in H-BACA
	int knn=0; //k nearest neighbours stored per point
	int format=1; //1: HOD-LR/BF  2:H-LR/BF 3: HSSBF/SHNBF 4: HSSBF-MD/SHNBF-MD 5: B-LR/B-BF
	double near_para=0.01; // control the admissibility
	C_QuantApp *quant_ptr;
	int v_major,v_minor,v_bugfix; //version numbers

    int tst = 1;
	int lrlevel=0;

	int nlevel = 0; // 0: tree level, nonzero if a tree is provided
	int* tree = new int[(int)pow(2,nlevel)]; //user provided array containing size of each leaf node, not used if nlevel=0
	string trainfile("../EXAMPLE/KRR_DATA/susy_10Kn");
	// string fullmatfile("../EXAMPLE/FULLMAT_DATA/FHODLR_colmajor_real_double_40000x40000.dat");
	// string leaffile("../EXAMPLE/FULLMAT_DATA/leafs_40000_noheader.dat");
	tree[0] = Npo;

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
      {"N",                   required_argument, 0, 2},
      {"Ndim",                   required_argument, 0, 3},
      {"ker",             required_argument, 0, 4},
      {"h",             required_argument, 0, 5},
      {"lambda",          required_argument, 0, 6},
      {"trainfile",          required_argument, 0, 8},
    //   {"nlevel",          required_argument, 0, 9},
    //   {"leaffile",          required_argument, 0, 10},
    //   {"fullmatfile",          required_argument, 0, 11},
      {"help",          no_argument, 0, 'h'},
      {NULL, 0, NULL, 0}
    };
  int c, option_index = 0;
  opterr = optind = 0;
  while ((c = getopt_long_only
          (argc, argv1.data(), "h",
            long_options, &option_index)) != -1) {

    switch (c) {
    case 1: {
      std::istringstream iss(optarg);
      iss >> tst;
    } break;
    case 2: {
      std::istringstream iss(optarg);
      iss >> Npo;
    } break;
    case 3: {
      std::istringstream iss(optarg);
      iss >> Ndim;
    } break;
    case 4: {
      std::istringstream iss(optarg);
      iss >> ker;
    } break;
    case 5: {
      std::istringstream iss(optarg);
      iss >> h;
    } break;
    case 6: {
      std::istringstream iss(optarg);
      iss >> lambda;
    } break;
    case 7: {
    //   std::istringstream iss(optarg);
    //   iss >> rank_rand;
    } break;
	case 8: {
	  std::istringstream iss(optarg);
  	  iss >> trainfile;
    } break;
	case 9: {
	  std::istringstream iss(optarg);
  	  iss >> nlevel;
    } break;
	case 10: {
	//   std::istringstream iss(optarg);
  	//   iss >> leaffile;
    } break;
	case 11: {
	//   std::istringstream iss(optarg);
  	//   iss >> fullmatfile;
    } break;
	case 'h': {
		if(myrank==master_rank)
		std::cout<<" tst=1: testing data sets with csv formats with ker 1:5 \n "<<std::endl;
	} break;
    default: break;
    }
  }


if(myrank==master_rank){
	d_c_bpack_getversionnumber(&v_major,&v_minor,&v_bugfix);
	std::cout<<"ButterflyPACK Version: "<<v_major<<"."<<v_minor<<"."<<v_bugfix<<std::endl;
}


	/*****************************************************************/
	/* @brief Test Kernels using UCI data sets */
if(tst==1){
    vector<double> data_train = write_from_file<double>(trainfile + "_train.csv");
	assert(Npo == data_train.size() / Ndim);
	quant_ptr=new C_QuantApp(data_train, Ndim, h, lambda,ker);
	dat_ptr = new double[data_train.size()];
	for(int ii=0;ii<data_train.size();ii++)
		dat_ptr[ii] = data_train.data()[ii];
	nogeo=0;
}

// 	/*****************************************************************/
// 	/* Test Full matrices*/
// if(tst==4){
// 	if(ker !=7){
// 		if(myrank==master_rank)std::cout<<"Forcing ker to 7 for tst=4."<<std::endl;
// 		ker = 7;
// 	}
// 	// Npo = 40000;

// 	delete(tree);
// 	// nlevel = 7;
// 	vector<int> t1((int)pow(2,nlevel));
// 	vector<double> matFull(Npo*Npo);
// 	vector<int> perm(Npo);
// 	if(myrank==master_rank){
// 		// vector<double> matFull1(Npo*Npo);
// 		ifstream f(fullmatfile, ios::binary);
// 		f.read((char*)matFull.data(), sizeof(double)*Npo*Npo);

// 		// perm = write_from_file<int>("../EXAMPLE/FULLMAT_DATA/sorder_40000.dat");


// 		// vector<int> perm1=perm;
// 		// for(int ii=0;ii<Npo;ii++)
// 		// 	perm.data()[perm1.data()[ii]]=ii;


// 		t1 = write_from_file<int>(leaffile);

// 		// std::cout<<matFull.data()[Npo*Npo-1]<<" "<<perm.data()[Npo-1]<<" "<<ccc<<std::endl;
// 	}
// 	MPI_Bcast(matFull.data(), Npo*Npo, MPI_DOUBLE, 0, MPI_COMM_WORLD);
// 	// MPI_Bcast(perm.data(), Npo, MPI_INT, 0, MPI_COMM_WORLD);
// 	MPI_Bcast(t1.data(), (int)pow(2,nlevel), MPI_INT, 0, MPI_COMM_WORLD);

// 	tree = new int[(int)pow(2,nlevel)];
// 	for(int ii=0;ii<(int)pow(2,nlevel);ii++)
// 	tree[ii] = t1.data()[ii];


// 	quant_ptr=new C_QuantApp(Npo, ker, matFull);
// 	nogeo=1;
// 	sort_opt=0;
// }



	/*****************************************************************/

	if(myrank==master_rank)std::cout<<"Npo "<<Npo<<" Ndim "<<Ndim<<std::endl;

	int myseg=0;     // local number of unknowns
	int* perms = new int[Npo]; //permutation vector returned by HODLR
	int* groups = new int[size];
	int i_opt;
	double d_opt;
	int cpp=1; //1: use user-defined cpp/c functions for construction

	//quantities for the first holdr
	F2Cptr bmat;  //hierarchical matrix returned by Fortran code
	F2Cptr option;     //option structure returned by Fortran code
	F2Cptr stats;      //statistics structure returned by Fortran code
	F2Cptr msh;		   //d_mesh structure returned by Fortran code
	F2Cptr kerquant;   //kernel quantities structure returned by Fortran code
	F2Cptr ptree;      //process tree returned by Fortran code


	MPI_Fint Fcomm;  // the fortran MPI communicator
	Fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);

	for (int i = 0; i < size; i++)groups[i]=i;
	// create butterflypack data structures
	d_c_bpack_createptree(&size, groups, &Fcomm, &ptree);
	d_c_bpack_createoption(&option);
	d_c_bpack_createstats(&stats);


	// set default butterflypack options
	d_c_bpack_set_D_option(&option, "tol_comp", tol);
	d_c_bpack_set_D_option(&option, "tol_rand", tol);
	d_c_bpack_set_D_option(&option, "tol_Rdetect", tol*1e-1);
	d_c_bpack_set_D_option(&option, "sample_para", sample_para);
	d_c_bpack_set_D_option(&option, "near_para", near_para);
	d_c_bpack_set_I_option(&option, "nogeo", nogeo);
	d_c_bpack_set_I_option(&option, "Nmin_leaf", Nmin);
	d_c_bpack_set_I_option(&option, "RecLR_leaf", com_opt);
	d_c_bpack_set_I_option(&option, "xyzsort", sort_opt);
	d_c_bpack_set_I_option(&option, "ErrFillFull", checkerr);
	d_c_bpack_set_I_option(&option, "BACA_Batch", batch);
	d_c_bpack_set_I_option(&option, "LR_BLK_NUM", bnum);
	d_c_bpack_set_I_option(&option, "cpp", cpp);
	d_c_bpack_set_I_option(&option, "LRlevel", lrlevel);
	d_c_bpack_set_I_option(&option, "knn", knn);
	d_c_bpack_set_I_option(&option, "verbosity", 1);
	d_c_bpack_set_I_option(&option, "less_adapt", 1);
	d_c_bpack_set_I_option(&option, "itermax", 10);
	d_c_bpack_set_I_option(&option, "ErrSol", 0);
	d_c_bpack_set_I_option(&option, "format", format);

	// set command-line butterflypack options
	set_option_from_command_line(argc, argv,option);

	// print out butterflypack options
	d_c_bpack_printoption(&option,&ptree);

    // construct matrix with geometrical points
	d_c_bpack_construct_init(&Npo, &Ndim, dat_ptr, nns_ptr,&nlevel, tree, perms, &myseg, &bmat, &option, &stats, &msh, &kerquant, &ptree, &C_FuncDistmn, &C_FuncNearFar, quant_ptr);
	d_c_bpack_construct_element_compute(&bmat, &option, &stats, &msh, &kerquant, &ptree, &C_FuncZmn, &C_FuncZmnBlock, quant_ptr);

	// factor matrix
	d_c_bpack_factor(&bmat,&option,&stats,&ptree,&msh);

	int nrhs=1;
	double* b = new double[nrhs*myseg];
	double* x = new double[nrhs*myseg];

	// generate a global rhs vector b_glo
    vector<double> x_glo(Npo*nrhs,0.0);
    vector<double> b_glo(Npo*nrhs,0.0);
	b_glo.data()[0]=1.0; // all but 1 element is zero in the RHS

	// map b_glo to the local rhs vector b
	for (int i=0; i<myseg; i++){
      int i_new_loc = i+1;
      int i_old;
      d_c_bpack_new2old(&msh,&i_new_loc,&i_old);
      for (int nth=0; nth<nrhs; nth++){
        b[i+nth*myseg] = b_glo.data()[i_old-1+nth*Npo];
      }
    }

	// solve the system
	d_c_bpack_solve(x,b,&myseg,&nrhs,&bmat,&option,&stats,&ptree);

	// map local solution vector x to the global solution vector x_glo
	for (int i=0; i<myseg; i++){
      int i_new_loc = i+1;
      int i_old;
      d_c_bpack_new2old(&msh,&i_new_loc,&i_old);
      for (int nth=0; nth<nrhs; nth++){
        x_glo.data()[i_old-1+nth*Npo] = x[i+nth*myseg];
      }
  }
  MPI_Allreduce(MPI_IN_PLACE,x_glo.data(), Npo*nrhs, MPI_C_DOUBLE, MPI_SUM, MPI_COMM_WORLD);



	if(myrank==master_rank)std::cout<<"Printing stats of the first hierarchical matrix: "<<std::endl;
	d_c_bpack_printstats(&stats,&ptree);

	d_c_bpack_deletestats(&stats);
	d_c_bpack_deleteproctree(&ptree);
	d_c_bpack_deletemesh(&msh);
	d_c_bpack_deletekernelquant(&kerquant);
	d_c_bpack_delete(&bmat);
	d_c_bpack_deleteoption(&option);

	delete quant_ptr;
	delete[] perms;
	delete[] tree;

	delete[] b;
	delete[] x;


	Cblacs_exit(1);
	MPI_Finalize();                                 // Terminate MPI. Once called, no other MPI routines may be called
    return 0;
}
//------------------------------------------------------------------------------
