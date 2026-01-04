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
 * @brief This c++ driver provides a few examples to illustrate the c++ interface to ButterflyPACK's Fortran subroutines, particularly the entry-evaluation and matvec-based APIs. This file works on the double data type.
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

#include "dBPACK_wrapper.h"



//------------------------------------------------------------------------------
using namespace std;


const double pi = 4.0 * atan(1.0);

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
/**  2-norm distance */
inline double dist2(double *x, double *y, int d) {
  double k = 0.;
  for (int i = 0; i < d; i++)
    k += pow(x[i] - y[i], 2.);
  return k;
}

/**  dot product of two real vectors */
inline double dot_product(double* v, double* u, int d)
{
    double result = 0.0;
    for (int i = 0; i < d; i++)
        result += v[i]*u[i];
    return result;
}


/**  Gauss Kernel */
inline double Gauss_kernel(double *x, double *y, int d, double h) {
  double dists;
  dists = dist2(x, y, d);
  if(dists> -log(1e-30)*2.0*pow(h,2.0)){
	return 0;
  }else{
	return exp(-dists / (2. * h * h));
  }
}

/**  Laplacian/Exponential Kernel */
inline double Laplace_kernel(double *x, double *y, int d, double h) {
  double dists;
  dists = dist2(x, y, d);
  // if(dists> -log(1e-30)*2.0*pow(h,2.0)){
	// return 0;
  // }else{
	return exp(-sqrt(dists)/(h));
  // }
}



/** R^4 kernel */
inline double K07_kernel(double *x, double *y, int d) {
  double dists;
  dists = dist2(x, y, d);
  return pow(dists,4);
}

/**  sqrt(R^2+h) kernel */
inline double K08_kernel(double *x, double *y, int d, double h) {
  double dists;
  dists = dist2(x, y, d);
  return sqrt(dists+h);
}

/**  1/sqrt(R^2+h) kernel */
inline double K09_kernel(double *x, double *y, int d, double h) {
  double dists;
  dists = dist2(x, y, d);
  return 1.0/sqrt(dists+h);
}

/**  Polynomial kernel (X^tY+h)^2 */
inline double K10_kernel(double *x, double *y, int d, double h) {
  double dotp;
  dotp = dot_product(x, y, d);
  return pow(dotp+h,2);
}

// // LR Product of two random matrices
// inline double LR_random(int *m, int *n) {
  // double dotp;
  // dotp = dot_product(x, y, d);
  // return pow(dotp+h,2);
// }


inline double GreenFun(int m, double w, double v0, double tau)
{
    double out;
    double q = -(m-2.0)/2.0;
    double v =q;
    double z = w*tau;
    double t0 = -sqrt(pi)/2*pow(2*tau/w,q)*((sqrt(2/(pi*z))*(sin(z)*cos(-v*pi) + cos(z)*sin(-v*pi)))*sin(q*pi) -sqrt(2/(pi*z))*(cos(z)*cos(-v*pi)-sin(z)*sin(-v*pi))*cos(q*pi));
    out =  v0*t0;
    return out;
}

inline double GreenFun_kernel(double *x, double *y, int d, double w)
{
    double dists = dist2(x, y, d);
    int self = sqrt(dists)<1e-20? 1:0;
    if(self==1){
        return 100.0;
    }else{
        double s0 = 2.0;
        double tau = sqrt(pow(s0,2)* dists);
        double D1 = s0/(2.0*pi);
        return GreenFun(d, w, D1, tau);
    }
}

/**  The object handling kernel parameters and sampling function */
class C_QuantApp {
public:
  vector<double> _data;
  int _d = 0;
  int _n = 0;
  double _h = 0.;
  double _w = 1.5;
  double _l = 0.;
  int _ker=1; //

  int _rank_rand=32;
  int _n_rand;
  std::vector<double> _MatU;
  std::vector<double> _MatV;
  std::vector<double> _MatFull;

  std::vector<int> _Hperm;
  std::vector<int> _iHperm;
  int _nloc = 0;

  F2Cptr* bmat;  //hierarchical matrix returned by Fortran code
  F2Cptr* bf;  //BF returned by Fortran code
  F2Cptr* stats;      //statistics structure returned by Fortran code
  F2Cptr* msh;		   //mesh structure returned by Fortran code
  F2Cptr* ptree;      //process tree returned by Fortran code
  F2Cptr* option;      //option structure returned by Fortran code


  C_QuantApp() = default;

  C_QuantApp(vector<double> data, int d, double h, double l, int ker)
    : _data(move(data)), _d(d), _n(_data.size() / _d),
      _h(h), _l(l),_ker(ker){
    assert(size_t(_n * _d) == _data.size());
	}

  C_QuantApp(int n_rand, int rank_rand, int ker, vector<double> MatU, vector<double> MatV)
    : _n_rand(n_rand), _rank_rand(rank_rand), _ker(ker), _MatU(move(MatU)), _MatV(move(MatV)){
	// cout<<_n_rand<<_rank_rand<<_MatU.size()<<endl;
    assert(size_t(_n_rand * _rank_rand) == _MatU.size());
	}

  C_QuantApp(int n, int ker, vector<double> MatFull, vector<int> perm)
    : _n(n), _ker(ker), _MatFull(move(MatFull)), _Hperm(move(perm)){
	// cout<<_n_rand<<_rank_rand<<_MatU.size()<<endl;
    assert((size_t)_n * (size_t)_n == _MatFull.size());
    // cout<<(size_t)_n * (size_t)_n<<" "<<_MatFull.size()<<endl;
	}

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
		*val =0;
		for (int k = 0; k < _rank_rand; k++){
			*val += _MatU[k*_n_rand+m]*_MatV[k*_n_rand+n];
		}
		break;
	case 7: //Full matrix
		*val =_MatFull[(size_t)n*(size_t)_n+(size_t)m];
		// *val =_MatFull[_Hperm[n]*_n+_Hperm[m]];
		break;
	case 8: //Laplacian kernel
		*val = Laplace_kernel(&_data[m * _d], &_data[n * _d], _d, _h);
		if (m==n)
		*val += _l;
		break;
	case 9: //Laplacian kernel with low-rank update
		*val = Laplace_kernel(&_data[m * _d], &_data[n * _d], _d, _h);
		if (m==n)
		*val += _l;
    //  adding a low-rank update U*U^T of rank R where U(i,j) = j/R/10
		for (int k = 0; k < _rank_rand; k++){
			*val += pow(k+1,2)/(double)pow(_rank_rand,2)/100;
		}
		break;
	case 10: //Wave equation kernel
		*val = GreenFun_kernel(&_data[m * _d], &_data[n * _d], _d, _w);
		break;
	}
  }
};

/** The sampling function wrapper required by the Fortran HODLR code */
inline void C_FuncZmn(int *m, int *n, double *val, C2Fptr quant) {

  C_QuantApp* Q = (C_QuantApp*) quant;
  Q->Sample(*m-1,*n-1,val);
}

/** The sampling function wrapper required by the Fortran HODLR code */
inline void C_FuncBZmn(int *m, int *n, double *val, C2Fptr quant) {

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


/**  The distance function wrapper required by the Fortran HODLR code */
inline void C_FuncDistmn(int *m, int *n, double *val, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;

}

/**  The compressibility function wrapper required by the Fortran HODLR code */
inline void C_FuncNearFar(int *m, int *n, int *val, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;

}

/**  The extraction sampling function wrapper required by the Fortran HODLR code */
inline void C_FuncZmnBlock(int* Ninter, int* Nallrows, int* Nallcols, int64_t* Nalldat_loc, int* allrows, int* allcols, double* alldat_loc, int* rowidx,int* colidx, int* pgidx, int* Npmap, int* pmaps, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;
  d_c_bpack_extractelement(Q->bmat,Q->option,Q->msh,Q->stats,Q->ptree,Ninter,Nallrows, Nallcols, Nalldat_loc, allrows,allcols,alldat_loc,rowidx,colidx,pgidx,Npmap,pmaps);
}

/**  The extraction sampling function wrapper required by the Fortran HODLR code */
inline void C_FuncBZmnBlock(int* Ninter, int* Nallrows, int* Nallcols, int64_t* Nalldat_loc, int* allrows, int* allcols, double* alldat_loc, int* rowidx,int* colidx, int* pgidx, int* Npmap, int* pmaps, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;

  d_c_bf_extractelement(Q->bf,Q->option,Q->msh,Q->stats,Q->ptree,Ninter,Nallrows, Nallcols, Nalldat_loc, allrows,allcols,alldat_loc,rowidx,colidx,pgidx,Npmap,pmaps);

  // // for(int ii=0;ii<*Nallrows;ii++)cout<<allrows[ii]<<endl;
  // d_c_bpack_extractelement(Q->bmat,Q->option,Q->msh,Q->stats,Q->ptree,Ninter,Nallrows, Nallcols, Nalldat_loc, allrows,allcols,alldat_loc,rowidx,colidx,pgidx,Npmap,pmaps);
}


/**  The matvec sampling function wrapper required by the Fortran HODLR code */
inline void C_FuncHMatVec(char const *trans, int *nin, int *nout, int *nvec, double const *xin, double *xout, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;

    #if 0  // user provides a function to do the matvec in the original order
      int Npo = Q->_n;
      vector<double> xin_glo(Npo*(*nvec),0.0);
      vector<double> xout_glo(Npo*(*nvec),0.0);

      // gather xin into xin_glo
      for (int i=0; i<*nin; i++){
        int i_new_loc = i+1;
        int i_old;
        d_c_bpack_new2old(Q->msh,&i_new_loc,&i_old);
        for (int nth=0; nth<(*nvec); nth++){
          xin_glo.data()[i_old-1+nth*Npo] = xin[i+nth*(*nin)];
        }
      }
      MPI_Allreduce(xin_glo.data(),xin_glo.data(), Npo*(*nvec), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      // ****************************
      // call your function here to multiply A with xin_glo (in the original order) to get xout_glo
      // ****************************


      // scatter xout_glo into xout
      MPI_Allreduce(xout_glo.data(),xout_glo.data(), Npo*(*nvec), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      for (int i=0; i<*nout; i++){
        int i_new_loc = i+1;
        int i_old;
        d_c_bpack_new2old(Q->msh,&i_new_loc,&i_old);
        for (int nth=0; nth<(*nvec); nth++){
          xout[i+nth*(*nout)] = xout_glo.data()[i_old-1+nth*Npo];
        }
      }

    #else
      d_c_bpack_mult(trans, xin, xout, nin, nout, nvec, Q->bmat,Q->option,Q->stats,Q->ptree);
    #endif

}


/**  The matvec sampling function wrapper required by the Fortran HODLR code */
inline void C_FuncBMatVec(char const *trans, int *nin, int *nout, int *nvec, double const *xin, double *xout, C2Fptr quant, double *a, double *b) {
  C_QuantApp* Q = (C_QuantApp*) quant;
  int cnt = (*nvec)*(*nout);
  double* xout1 = new double[cnt];

  d_c_bpack_mult(trans, xin, xout1, nin, nout, nvec, Q->bmat,Q->option,Q->stats,Q->ptree);

  for (int ii=0; ii<cnt; ii++){
	xout[ii] = *b*xout[ii] + *a*xout1[ii];
  }
  delete[] xout1;
}


/**  Read a data file into a vector */
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


////////////////////////////////////////////////////////////////////////////////
/**  --------------------------- Main Code Starts Here ------------------------ */

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
	int ker=8 ; // kernel choice
	int Npo=5000;   // matrix size
	int Ndim=1; //data dimension
	int rank_rand=100;  //rank of the random LR product
	double starttime, endtime;
	double* dat_ptr;
	int* nns_ptr;
	int nogeo;  // 1: no geometrical information passed to hodlr, dat_ptr and Ndim are dummy


	int Nmin=200; //finest leafsize
	double tol=1e-4; //compression tolerance
	double sample_para=2.0; //oversampling factor in entry evaluation
	int com_opt=5; //1:SVD 2:RRQR 3:ACA 4:BACA 5:BACA_improved 6:Pseudo-skeleton
	int sort_opt=1; //0:natural order 1:kd-tree 2:cobble-like ordering 3:gram distance-based cobble-like ordering
	int checkerr = 0; //1: check compression quality
	int batch = 100; //batch size for BACA
	int bnum = 1; //sqrt of #of subblocks in H-BACA
	int knn=0; //k nearest neighbours stored per point
	C_QuantApp *quant_ptr;
	int v_major,v_minor,v_bugfix; //version numbers

    int tst = 1;
	int lrlevel=0;

	int nlevel = 0; // 0: tree level, nonzero if a tree is provided
	int* tree = new int[(int)pow(2,nlevel)]; //user provided array containing size of each leaf node, not used if nlevel=0
	string trainfile("../EXAMPLE/KRR_DATA/susy_10Kn");
	string fullmatfile("../EXAMPLE/FULLMAT_DATA/FHODLR_colmajor_real_double_40000x40000.dat");
	string leaffile("../EXAMPLE/FULLMAT_DATA/leafs_40000_noheader.dat");


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
      {"rank_rand",          required_argument, 0, 7},
      {"trainfile",          required_argument, 0, 8},
      {"nlevel",          required_argument, 0, 9},
      {"leaffile",          required_argument, 0, 10},
      {"fullmatfile",          required_argument, 0, 11},
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
      std::istringstream iss(optarg);
      iss >> rank_rand;
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
	  std::istringstream iss(optarg);
  	  iss >> leaffile;
    } break;
	case 11: {
	  std::istringstream iss(optarg);
  	  iss >> fullmatfile;
    } break;
	case 'h': { std::cout<<" tst=1: testing data sets with csv formats with ker 1:5 \n tst=2: testing randomly generated data sets with ker 1:5 \n tst=3: testing a LR product of two random matrices with ker=6 \n tst=4: testing full matrix and leaves stored in file with ker=7 "<<std::endl; } break;
    default: break;
    }
  }




if(myrank==master_rank){
	d_c_bpack_getversionnumber(&v_major,&v_minor,&v_bugfix);
	std::cout<<"ButterflyPACK Version: "<<v_major<<"."<<v_minor<<"."<<v_bugfix<<std::endl;
}
	tree[0] = Npo;

	/*****************************************************************/

  /** tst=1: Test Kernels for Liza's data sets */
if(tst==1){
    vector<double> data_train = write_from_file<double>(trainfile + "_train.csv");
	assert(Npo == data_train.size() / Ndim);
	quant_ptr=new C_QuantApp(data_train, Ndim, h, lambda,ker);
	dat_ptr = new double[data_train.size()];
	for(int ii=0;ii<data_train.size();ii++)
		dat_ptr[ii] = data_train.data()[ii];
	nogeo=0;
}

	/*****************************************************************/
	/** tst=2: Test Kernels for Random point clouds */
if(tst==2){
	vector<double> data_train(Npo*Ndim);
      for (int i=0; i<Npo*Ndim; i++)
	data_train[i] = (double)rand() / RAND_MAX;
	MPI_Bcast(data_train.data(), Npo*Ndim, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	quant_ptr=new C_QuantApp(data_train, Ndim, h, lambda,ker);
	dat_ptr = new double[data_train.size()];
	for(int ii=0;ii<data_train.size();ii++)
		dat_ptr[ii] = data_train.data()[ii];
	nogeo=0;
}

	/*****************************************************************/
	/** tst=3: Test Product of two Random matrices*/
if(tst==3){
	if(ker !=6){
		if(myrank==master_rank)std::cout<<"Forcing ker to 6 for tst=3."<<std::endl;
		ker = 6;
	}
	vector<double> matU(Npo*rank_rand);
      for (int i=0; i<Npo*rank_rand; i++)
	matU[i] = (double)rand() / RAND_MAX;
	MPI_Bcast(matU.data(), Npo*rank_rand, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	vector<double> matV(Npo*rank_rand);
      for (int i=0; i<Npo*rank_rand; i++)
	matV[i] = (double)rand() / RAND_MAX;
	MPI_Bcast(matV.data(), Npo*rank_rand, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	quant_ptr=new C_QuantApp(Npo, rank_rand, ker, matU, matV);
	nogeo=1;
	sort_opt=0;
}

	/*****************************************************************/
	/** tst=4: Test Full matrices*/
if(tst==4){
	if(ker !=7){
		if(myrank==master_rank)std::cout<<"Forcing ker to 7 for tst=4."<<std::endl;
		ker = 7;
	}
	// Npo = 40000;

	delete(tree);
	// nlevel = 7;
	vector<int> t1((int)pow(2,nlevel));
	vector<double> matFull(Npo*(size_t)Npo);
	vector<int> perm(Npo);
	if(myrank==master_rank){
		// vector<double> matFull1(Npo*Npo);
		ifstream f(fullmatfile, ios::binary);
		f.read((char*)matFull.data(), sizeof(double)*Npo*(size_t)Npo);

		// perm = write_from_file<int>("../EXAMPLE/FULLMAT_DATA/sorder_40000.dat");


		// vector<int> perm1=perm;
		// for(int ii=0;ii<Npo;ii++)
		// 	perm.data()[perm1.data()[ii]]=ii;


		t1 = write_from_file<int>(leaffile);

		// std::cout<<matFull.data()[Npo*Npo-1]<<" "<<perm.data()[Npo-1]<<" "<<ccc<<std::endl;
	}
	MPI_Bcast(matFull.data(), Npo*Npo, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// MPI_Bcast(perm.data(), Npo, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(t1.data(), (int)pow(2,nlevel), MPI_INT, 0, MPI_COMM_WORLD);

	tree = new int[(int)pow(2,nlevel)];
	for(int ii=0;ii<(int)pow(2,nlevel);ii++)
	tree[ii] = t1.data()[ii];


	quant_ptr=new C_QuantApp(Npo, ker, matFull,perm);
	nogeo=1;
	sort_opt=0;
}



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
	// create hodlr data structures
	d_c_bpack_createptree(&size, groups, &Fcomm, &ptree);
	d_c_bpack_createoption(&option);
	d_c_bpack_createstats(&stats);



	// set hodlr options
	d_c_bpack_set_D_option(&option, "tol_comp", tol);
	d_c_bpack_set_D_option(&option, "tol_rand", tol);
	d_c_bpack_set_D_option(&option, "tol_Rdetect", tol*1e-1);
	d_c_bpack_set_D_option(&option, "sample_para", sample_para);
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
	d_c_bpack_set_I_option(&option, "verbosity", 2);
	d_c_bpack_set_I_option(&option, "less_adapt", 1);
	d_c_bpack_set_I_option(&option, "itermax", 10);
	d_c_bpack_set_I_option(&option, "ErrSol", 1);


	d_c_bpack_set_option_from_command_line(argc, argv,option);

  F2Cptr option_save;
  d_c_bpack_copyoption(&option,&option_save);
  if(ker ==8 || ker ==9 || ker ==10){
  d_c_bpack_set_I_option(&option, "format", 1);
  d_c_bpack_set_I_option(&option, "LRlevel", 100);
  d_c_bpack_set_I_option(&option, "knn", 100);
  }
	d_c_bpack_printoption(&option,&ptree);

    // construct hodlr with geometrical points
	d_c_bpack_construct_init(&Npo, &Ndim, dat_ptr, nns_ptr,&nlevel, tree, perms, &myseg, &bmat, &option, &stats, &msh, &kerquant, &ptree, &C_FuncDistmn, &C_FuncNearFar, quant_ptr);
	d_c_bpack_construct_element_compute(&bmat, &option, &stats, &msh, &kerquant, &ptree, &C_FuncZmn, &C_FuncZmnBlock, quant_ptr);


  if(ker !=8 && ker !=9 && ker !=10){
	// factor hodlr
	d_c_bpack_factor(&bmat,&option,&stats,&ptree,&msh);

	// solve the system
	int nrhs=1;
	double* b = new double[nrhs*myseg];
	double* x = new double[nrhs*myseg];

	for (int i = 0; i < nrhs*myseg; i++){
		b[i]=1;
	}
	d_c_bpack_solve(x,b,&myseg,&nrhs,&bmat,&option,&stats,&ptree);
  }

	if(myrank==master_rank)std::cout<<"Printing stats of the first HODLR: "<<std::endl;
	d_c_bpack_printstats(&stats,&ptree);


	//////////////////// use resulting hodlr as matvec to create a new holdr

	// quantities for the second holdr
	F2Cptr bmat1;  //hierarchical matrix returned by Fortran code
	F2Cptr option1;     //option structure returned by Fortran code
	F2Cptr stats1;      //statistics structure returned by Fortran code
	F2Cptr msh1;		   //d_mesh structure returned by Fortran code
	F2Cptr kerquant1;   //kernel quantities structure returned by Fortran code
	F2Cptr ptree1;      //process tree returned by Fortran code
	C_QuantApp *quant_ptr1; //user-defined object

	quant_ptr1=new C_QuantApp();
	quant_ptr1->bmat=&bmat;
	quant_ptr1->msh=&msh;
	quant_ptr1->ptree=&ptree;
	quant_ptr1->stats=&stats;
	quant_ptr1->option=&option;
  quant_ptr1->_n=Npo;

	d_c_bpack_createptree(&size, groups, &Fcomm, &ptree1);
	d_c_bpack_copyoption(&option_save,&option1);
  d_c_bpack_printoption(&option1,&ptree);
	d_c_bpack_createstats(&stats1);

	// d_c_bpack_set_I_option(&option1, "nogeo", 1); // no geometrical information
	// d_c_bpack_set_I_option(&option1, "xyzsort", 0);// natural ordering
	// d_c_bpack_set_I_option(&option1, "format", 1);// HODLR format

	int Npo1 = Npo;
	int myseg1=0;     // local number of unknowns
	int* perms1 = new int[Npo1]; //permutation vector returned by HODLR
	//tree1 and nlevel1 should be provided by the caller, otherwise natural ordering is used
	int nlevel1 = 0; // 0: tree level, nonzero if a tree is provided
	int* tree1 = new int[(int)pow(2,nlevel1)]; //user provided array containing size of each leaf node, not used if nlevel=0
	tree1[0] = Npo1;
	int Ndim1=Ndim; //data dimension
	double* dat_ptr1=dat_ptr;

	d_c_bpack_construct_init(&Npo1, &Ndim1, dat_ptr1, nns_ptr,&nlevel1, tree1, perms1, &myseg1, &bmat1, &option1, &stats1, &msh1, &kerquant1, &ptree1, &C_FuncDistmn, &C_FuncNearFar, quant_ptr);
	d_c_bpack_construct_matvec_compute(&bmat1, &option1, &stats1, &msh1, &kerquant1, &ptree1, &C_FuncHMatVec, quant_ptr1);


	if(myrank==master_rank)std::cout<<"Printing stats of the second HODLR: "<<std::endl;
	d_c_bpack_printstats(&stats1,&ptree1);

	d_c_bpack_deletestats(&stats1);
	d_c_bpack_deleteproctree(&ptree1);
	d_c_bpack_deletemesh(&msh1);
	d_c_bpack_deletekernelquant(&kerquant1);
	d_c_bpack_delete(&bmat1);
	d_c_bpack_deleteoption(&option1);
	delete quant_ptr1;
	delete[] perms1;
	delete[] tree1;


  if(ker !=8 && ker !=9 && ker !=10){
	//////////////////// use resulting hodlr as entry extraction to create a new holdr
	quant_ptr1=new C_QuantApp();
	quant_ptr1->bmat=&bmat;
	quant_ptr1->msh=&msh;
	quant_ptr1->ptree=&ptree;
	quant_ptr1->stats=&stats;
	quant_ptr1->option=&option;

	d_c_bpack_createptree(&size, groups, &Fcomm, &ptree1);
	d_c_bpack_copyoption(&option_save,&option1);
	d_c_bpack_createstats(&stats1);

	d_c_bpack_set_I_option(&option1, "nogeo", 1); // no geometrical information
	d_c_bpack_set_I_option(&option1, "xyzsort", 0);// natural ordering
	d_c_bpack_set_I_option(&option1, "format", 1);// HODLR format
	d_c_bpack_set_I_option(&option1, "elem_extract", 1);// use element extraction


	perms1 = new int[Npo1]; //permutation vector returned by HODLR
	//tree1 and nlevel1 should be provided by the caller, otherwise natural ordering is used
	nlevel1 = 0; // 0: tree level, nonzero if a tree is provided
	tree1 = new int[(int)pow(2,nlevel1)]; //user provided array containing size of each leaf node, not used if nlevel=0
	tree1[0] = Npo1;

	// d_c_bpack_construct_matvec_init(&Npo1, &nlevel1, tree1, perms1, &myseg1, &bmat1, &option1, &stats1, &msh1, &kerquant1, &ptree1);
	// d_c_bpack_construct_matvec_compute(&bmat1, &option1, &stats1, &msh1, &kerquant1, &ptree1, &C_FuncHMatVec, quant_ptr1);


	d_c_bpack_construct_init(&Npo1, &Ndim, dat_ptr, nns_ptr,&nlevel1, tree1, perms1, &myseg1, &bmat1, &option1, &stats1, &msh1, &kerquant1, &ptree1, &C_FuncDistmn, &C_FuncNearFar, quant_ptr);
	d_c_bpack_construct_element_compute(&bmat1, &option1, &stats1, &msh1, &kerquant1, &ptree1, &C_FuncZmn, &C_FuncZmnBlock, quant_ptr1);


	if(myrank==master_rank)std::cout<<"Printing stats of the third HODLR: "<<std::endl;
	d_c_bpack_printstats(&stats1,&ptree1);
	d_c_bpack_printstats(&stats,&ptree);

	d_c_bpack_deletestats(&stats1);
	d_c_bpack_deleteproctree(&ptree1);
	d_c_bpack_deletemesh(&msh1);
	d_c_bpack_deletekernelquant(&kerquant1);
	d_c_bpack_delete(&bmat1);
	d_c_bpack_deleteoption(&option1);
	delete quant_ptr1;
	delete[] perms1;
	delete[] tree1;
  }




if(tst==3){
	//////////////////// use resulting hodlr as matvec to create a new bf

	// quantities for the second holdr
	F2Cptr bf,bf2;  //BF returned by Fortran code
	F2Cptr option2;     //option structure returned by Fortran code
	F2Cptr stats2;      //statistics structure returned by Fortran code
	F2Cptr msh2;		   //d_mesh structure returned by Fortran code
	F2Cptr kerquant2;   //kernel quantities structure returned by Fortran code
	F2Cptr ptree2;      //process tree returned by Fortran code
	C_QuantApp *quant_ptr2; //user-defined object

	quant_ptr2=new C_QuantApp();
	quant_ptr2->bmat=&bmat;
	quant_ptr2->msh=&msh;
	quant_ptr2->ptree=&ptree;
	quant_ptr2->stats=&stats;
	quant_ptr2->option=&option;

	d_c_bpack_createptree(&size, groups, &Fcomm, &ptree2);
	d_c_bpack_copyoption(&option_save,&option2);
	d_c_bpack_createstats(&stats2);

	d_c_bpack_set_I_option(&option2, "nogeo", 1); // no geometrical information
	d_c_bpack_set_I_option(&option2, "xyzsort", 0);// natural ordering

	int M = Npo;
	int N = Npo;
	int myrow=0;     // local number of rows
	int mycol=0;     // local number of columns

	d_c_bf_construct_init(&M, &N, &myrow, &mycol, nns_ptr, nns_ptr, &msh, &msh, &bf, &option2, &stats2, &msh1, &kerquant2, &ptree2,&C_FuncDistmn, &C_FuncNearFar, quant_ptr2);
	d_c_bf_construct_matvec_compute(&bf, &option2, &stats2, &msh1, &kerquant2, &ptree2, &C_FuncBMatVec, quant_ptr2);

	if(myrank==master_rank)std::cout<<"Printing stats of the fourth BF: "<<std::endl;
	d_c_bpack_printstats(&stats2,&ptree2);

	d_c_bpack_deletestats(&stats2);
	d_c_bpack_deleteproctree(&ptree2);
	d_c_bpack_deletekernelquant(&kerquant2);
	d_c_bpack_deleteoption(&option2);

	delete quant_ptr2;



	//////////////////// use resulting hodlr as entry extraction to create a new bf
	quant_ptr2=new C_QuantApp();
	quant_ptr2->bf=&bf;
	quant_ptr2->msh=&msh1;
	quant_ptr2->ptree=&ptree;
	quant_ptr2->stats=&stats;
	quant_ptr2->option=&option;
	quant_ptr2->_n=Npo;

	d_c_bpack_createptree(&size, groups, &Fcomm, &ptree2);
	d_c_bpack_copyoption(&option_save,&option2);
	d_c_bpack_createstats(&stats2);

	d_c_bpack_set_I_option(&option2, "nogeo", 1); // no geometrical information
	d_c_bpack_set_I_option(&option2, "xyzsort", 0);// natural ordering
	d_c_bpack_set_I_option(&option2, "elem_extract", 1);// use block-wise element extraction

	d_c_bf_construct_init(&M, &N, &myrow, &mycol, nns_ptr, nns_ptr, &msh, &msh, &bf2, &option2, &stats2, &msh2, &kerquant2, &ptree2,&C_FuncDistmn, &C_FuncNearFar, quant_ptr2);
	d_c_bf_construct_element_compute(&bf2, &option2, &stats2, &msh2, &kerquant2, &ptree2, &C_FuncBZmn, &C_FuncBZmnBlock, quant_ptr2); // C_FuncBZmn is not referenced since elem_extract=1

	if(myrank==master_rank)std::cout<<"Printing stats of the fifth BF: "<<std::endl;
	d_c_bpack_printstats(&stats2,&ptree2);
	d_c_bpack_printstats(&stats,&ptree);

	d_c_bpack_deletestats(&stats2);
	d_c_bpack_deleteproctree(&ptree2);
	d_c_bpack_deletemesh(&msh2);
	d_c_bpack_deletekernelquant(&kerquant2);
	d_c_bf_deletebf(&bf2);
	d_c_bpack_deleteoption(&option2);

	d_c_bf_deletebf(&bf);
	d_c_bpack_deletemesh(&msh1);
	delete quant_ptr2;

}



	d_c_bpack_deletestats(&stats);
	d_c_bpack_deleteproctree(&ptree);
	d_c_bpack_deletemesh(&msh);
	d_c_bpack_deletekernelquant(&kerquant);
	d_c_bpack_delete(&bmat);
	d_c_bpack_deleteoption(&option);

	delete quant_ptr;
	delete[] perms;
	delete[] tree;



	Cblacs_exit(1);
	MPI_Finalize();                                 // Terminate MPI. Once called, no other MPI routines may be called
    return 0;
}
//------------------------------------------------------------------------------
