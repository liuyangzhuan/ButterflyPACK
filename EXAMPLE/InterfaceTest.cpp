//------------------------------------------------------------------------------
#include <iostream>
#include <math.h>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <iomanip>
 
#include <omp.h>
#include <pthread.h>

#include <cmath>
#include <cassert>
#include <iostream>
#include <random>
#include <vector>
#include <atomic>
#include <mpi.h>

#include "C_HODLR_C_Interface.h"



//------------------------------------------------------------------------------
using namespace std;

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

// // LR Product of two random matrices
// inline double LR_random(int *m, int *n) {
  // double dotp;
  // dotp = dot_product(x, y, d);
  // return pow(dotp+h,2);
// }



// The object handling kernel parameters and sampling function
class C_QuantZmn {
public:
  vector<double> _data;
  int _d = 0;
  int _n = 0;
  double _h = 0.;
  double _l = 0.;
  int _ker=1; // 
  
  int _rank_rand;
  int _n_rand;
  std::vector<double> _MatU;
  std::vector<double> _MatV;
  
  std::vector<int> _Hperm;
  std::vector<int> _iHperm;
  int _nloc = 0;

  C_QuantZmn() = default;
  
  C_QuantZmn(vector<double> data, int d, double h, double l, int ker)
    : _data(move(data)), _d(d), _n(_data.size() / _d),
      _h(h), _l(l),_ker(ker){
    assert(size_t(_n * _d) == _data.size());
	}
  
  C_QuantZmn(int n_rand, int rank_rand, int ker, vector<double> MatU, vector<double> MatV)
    : _n_rand(n_rand), _rank_rand(rank_rand), _ker(ker), _MatU(move(MatU)), _MatV(move(MatV)){
	cout<<_n_rand<<_rank_rand<<_MatU.size()<<endl;
    assert(size_t(_n_rand * _rank_rand) == _MatU.size());
	}  
  
  
  inline void Sample(int m, int n, doublecomplex* val){
	switch(_ker){
	case 1: //Gaussian kernel
		val->i=0;
		val->r = Gauss_kernel(&_data[m * _d], &_data[n * _d], _d, _h);
		if (m==n)
		val->r += _l;
		break;
	case 2: //R^4 kernel
		val->i=0;
		val->r = K07_kernel(&_data[m * _d], &_data[n * _d], _d);	
		break;
	case 3: //sqrt(R^2+h) kernel
		val->i=0;
		val->r = K08_kernel(&_data[m * _d], &_data[n * _d], _d, _h);	
		break;
	case 4: //1/sqrt(R^2+h) kernel
		val->i=0;
		val->r = K09_kernel(&_data[m * _d], &_data[n * _d], _d, _h);
		break;
	case 5: //Polynomial kernel (X^tY+h)^2
		val->i=0;
		val->r = K10_kernel(&_data[m * _d], &_data[n * _d], _d, _h);
		break;
	case 6: //Low-rank product of two random matrices
		val->i=0;
		val->r =0;
		for (int k = 0; k < _rank_rand; k++){
			val->r += _MatU[k*_n_rand+m]*_MatV[k*_n_rand+n];
		}
		break;
	}
  } 	
};


// The sampling function wrapper required by the Fortran HODLR code
inline void C_FuncZmn(int *m, int *n, doublecomplex *val, C2Fptr quant) {
	
  C_QuantZmn* Q = (C_QuantZmn*) quant;	
  Q->Sample(*m,*n,val);
}

// Read a data file into a vector
vector<double> write_from_file(string filename) {
  vector<double> data;
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
// --------------------------- Main Code Starts Here ------------------------ //

int main(int argc, char* argv[])
{

    int myrank, size;                     // Store values of processor rank and total no of procs requestedss
    int master_rank = 0;
	MPI_Init(&argc, &argv); 	                            // Initialize MPI, called only once
    MPI_Comm_size(MPI_COMM_WORLD, &size); 	                // Get no of procs
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 	                // Get no of procs
    MPI_Op op;
	double h; //kernel parameter
	double lambda ; //kernel parameter 
	int ker ; // kernel choice 
	int Npo;   // matrix size
	int Ndim; //data dimension
	double starttime, endtime;
	double* dat_ptr; 
	int preorder;  // 0: holdr code will reorder the elements and return a permutation 1: holdr code will use user-supplied tree to create the holdr structure	
	
	
	/*****************************************************************/
	/* Test Kernels for Liza's data sets */
#if 1	
	string filename("./EXAMPLE/SUSY/susy_10Kn");
	h = 0.1;
	lambda = 10.;
	ker = 1;
	Ndim = 8;

    vector<double> data_train = write_from_file(filename + "_train.csv");
	Npo = data_train.size() / Ndim;

	C_QuantZmn quant(data_train, Ndim, h, lambda,ker);	
	dat_ptr = data_train.data();
	preorder=0;
#endif
	
	/*****************************************************************/
	/* Test Kernels for Random point clouds */
#if 0	
	h = 3.267; //0.47;
	lambda = 10.;
	ker = 1;
	Ndim = 6;
	Npo = 4096;
	
	vector<double> data_train(Npo*Ndim);
      for (int i=0; i<Npo*Ndim; i++)
	data_train[i] = (double)rand() / RAND_MAX;
	MPI_Bcast(data_train.data(), Npo*Ndim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	C_QuantZmn quant(data_train, Ndim, h, lambda,ker);	
	dat_ptr = data_train.data();
	preorder=0;
#endif	
	
	/*****************************************************************/
	/* Test Product of two Random matrices*/
#if 0	
	ker = 6;
	int rank_rand = 10;
	Npo = 10000;
	
	vector<double> matU(Npo*rank_rand);
      for (int i=0; i<Npo*rank_rand; i++)
	matU[i] = (double)rand() / RAND_MAX;
	MPI_Bcast(matU.data(), Npo*rank_rand, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	vector<double> matV(Npo*rank_rand);
      for (int i=0; i<Npo*rank_rand; i++)
	matV[i] = (double)rand() / RAND_MAX;
	MPI_Bcast(matV.data(), Npo*rank_rand, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	C_QuantZmn quant(Npo, rank_rand, ker, matU, matV);	
	preorder=1;
#endif	
	/*****************************************************************/

	if(myrank==master_rank)std::cout<<"Npo "<<Npo<<" Ndim "<<Ndim<<std::endl;
	
	
	int Nmin=200; //finest leafsize 
	double tol=1e-4; //compression tolerance
	int nmpi=size;	 //number of active MPI ranks 
	int ninc=1;      // use one rank in HODLR per ninc ranks in HSS
	int myseg=0;     // local number of unknowns
	int* perms = new int[Npo]; //permutation vector returned by HODLR
	int nlevel = 0; // 0: tree level, nonzero if a tree is provided 
	int* tree = new int[(int)pow(2,nlevel)]; //user provided array containing size of each leaf node, not used if nlevel=0
	int* groups = new int[size];
	int i_opt;
	double d_opt;
	
	
	F2Cptr ho_bf_for;  //forward HODLR returned by Fortran code 
	F2Cptr ho_bf_inv;  //factored HODLR returned by Fortran code 
	F2Cptr option;     //option structure returned by Fortran code 
	F2Cptr stats;      //statistics structure returned by Fortran code
	F2Cptr msh;		   //mesh structure returned by Fortran code
	F2Cptr kerquant;   //kernel quantities structure returned by Fortran code
	F2Cptr ptree;      //process tree returned by Fortran code
	
	
	MPI_Fint Fcomm;  // the fortran MPI communicator
	Fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);  
	
	for (int i = 0; i < size; i++)groups[i]=i;
	// create hodlr data structures
	FC_GLOBAL_(c_createptree,C_CREATEPTREE)(&nmpi, groups, &Fcomm, &ptree);
	FC_GLOBAL_(c_createoption,C_CREATEOPTION)(&option);	
	FC_GLOBAL_(c_createstats,C_CREATESTATS)(&stats);		
	
	// set hodlr options
	
	set_D_option(&option, "tol_comp", 1e-4);
	set_I_option(&option, "preorder", 0);
	set_I_option(&option, "Nmin_leaf", 200);
	set_I_option(&option, "RecLR_leaf", 4); //1:SVD 2:RRQR 3:ACA 4:BACA
	
	
	//create hodlr quantities
	FC_GLOBAL_(c_hodlr_construct,C_HODLR_CONSTRUCT)(&Npo, &Ndim, dat_ptr, &nlevel, tree, perms, &myseg, &ho_bf_for, &option, &stats, &msh, &kerquant, &ptree, &C_FuncZmn, &quant, &Fcomm);	
	
    // construct hodlr with geometrical points	
	FC_GLOBAL_(c_hodlr_construct,C_HODLR_CONSTRUCT)(&Npo, &Ndim, dat_ptr, &nlevel, tree, perms, &myseg, &ho_bf_for, &option, &stats, &msh, &kerquant, &ptree, &C_FuncZmn, &quant, &Fcomm);	
	
	// factor hodlr
	FC_GLOBAL_(c_hodlr_factor,C_HODLR_FACTOR)(&ho_bf_for,&ho_bf_inv,&option,&stats,&ptree);

	// solve the system 
	int nrhs=1;
	doublecomplex* b = new doublecomplex[nrhs*myseg];
	doublecomplex* x = new doublecomplex[nrhs*myseg];
	
	for (int i = 0; i < nrhs*myseg; i++){
		b[i].r=1;
		b[i].i=0;
	}	
	FC_GLOBAL_(c_hodlr_solve,C_HODLR_SOLVE)(x,b,&myseg,&nrhs,&ho_bf_for,&ho_bf_inv,&option,&stats,&ptree);

	MPI_Finalize();                                 // Terminate MPI. Once called, no other MPI routines may be called
    return 0;
}
//------------------------------------------------------------------------------
