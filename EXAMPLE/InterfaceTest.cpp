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

#include "dC_HODLR_wrapper.h"



//------------------------------------------------------------------------------
using namespace std;


extern "C" {
      ///////////////////////////////////////////////
      ////// BLACS //////////////////////////////////
      ///////////////////////////////////////////////
      // void Cblacs_get(int, int, int *);
      // void Cblacs_gridinit(int *, const char *, int, int);
      // void Cblacs_gridmap(int *, int *, int, int, int);
      // void Cblacs_gridinfo(int, int *, int *, int *, int *);
      // void Cblacs_gridexit(int);
      void Cblacs_exit(int);
}

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
class C_QuantApp {
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
  
  F2Cptr* ho_bf;  //HODLR returned by Fortran code 
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
	}
  } 	
};


// The sampling function wrapper required by the Fortran HODLR code
inline void C_FuncZmn(int *m, int *n, double *val, C2Fptr quant) {
	
  C_QuantApp* Q = (C_QuantApp*) quant;	
  Q->Sample(*m,*n,val);
}

// The matvec sampling function wrapper required by the Fortran HODLR code
inline void C_FuncMatVec(char const *trans, int *nin, int *nout, int *nvec, double *xin, double *xout, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;	
  d_c_hodlr_mult(trans, xin, xout, nin, nout, nvec, Q->ho_bf,Q->option,Q->stats,Q->ptree);  
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
	int Ndim=1; //data dimension
	double starttime, endtime;
	double* dat_ptr; 
	int nogeo;  // 1: no geometrical information passed to hodlr, dat_ptr and Ndim are dummy	
	

	int Nmin=200; //finest leafsize 
	double tol=1e-4; //compression tolerance
	int com_opt=4; //1:SVD 2:RRQR 3:ACA 4:BACA
	int sort_opt=1; //0:natural order 1:kd-tree 2:cobble-like ordering 3:gram distance-based cobble-like ordering
	int checkerr = 1; //1: check compression quality 
	int batch = 100; //batch size for BACA
	string filename("smalltest.dat");
	C_QuantApp *quant_ptr, *quant_ptr1;
	
	
    int tst = stoi(argv[1]);
	/*****************************************************************/
	/* Test Kernels for Liza's data sets */
if(tst==1){
	filename = string(argv[2]);
	Ndim = stoi(argv[3]);
	ker = stoi(argv[4]);
	h = stof(argv[5]);
	lambda = stof(argv[6]);
	Nmin = stoi(argv[7]);
	tol = stof(argv[8]);
	com_opt = stoi(argv[9]);
	checkerr = stoi(argv[10]);
	batch = stoi(argv[11]);
	
	// Ndim = 8;
	// h = 0.2;
	// lambda = 10.;
	// ker = 1;
	
    vector<double> data_train = write_from_file(filename + "_train.csv");
	Npo = data_train.size() / Ndim;

	quant_ptr=new C_QuantApp(data_train, Ndim, h, lambda,ker);	
	// dat_ptr = data_train.data();
	dat_ptr = new double[data_train.size()];
	for(int ii=0;ii<data_train.size();ii++)
		dat_ptr[ii] = data_train.data()[ii];
	nogeo=0;
}
	
	/*****************************************************************/
	/* Test Kernels for Random point clouds */
if(tst==2){	
	// h = 3.267; //0.47;
	// lambda = 10.;
	// ker = 1;
	// Ndim = 6;
	// Npo = 4096;
	
	
	Npo = stoi(argv[2]);
	Ndim = stoi(argv[3]);
	ker = stoi(argv[4]);
	h = stof(argv[5]);
	lambda = stof(argv[6]);
	Nmin = stoi(argv[7]);
	tol = stof(argv[8]);
	com_opt = stoi(argv[9]);
	checkerr = stoi(argv[10]);
	batch = stoi(argv[11]);
	
	vector<double> data_train(Npo*Ndim);
      for (int i=0; i<Npo*Ndim; i++)
	data_train[i] = (double)rand() / RAND_MAX;
	MPI_Bcast(data_train.data(), Npo*Ndim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	quant_ptr=new C_QuantApp(data_train, Ndim, h, lambda,ker);	
	// dat_ptr = data_train.data();
	dat_ptr = new double[data_train.size()];
	for(int ii=0;ii<data_train.size();ii++)
		dat_ptr[ii] = data_train.data()[ii];	
	nogeo=0;
}	
	
	/*****************************************************************/
	/* Test Product of two Random matrices*/
if(tst==3){
	ker = 6;
	int rank_rand = 10;
	Npo = 10000;
	
	Npo = stoi(argv[2]);
	rank_rand = stoi(argv[3]);
	Nmin = stoi(argv[4]);
	tol = stof(argv[5]);
	com_opt = stoi(argv[6]);
	checkerr = stoi(argv[7]);
	batch = stoi(argv[8]);		
	
	
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

	if(myrank==master_rank)std::cout<<"Npo "<<Npo<<" Ndim "<<Ndim<<std::endl;
	
	int myseg=0;     // local number of unknowns
	int* perms = new int[Npo]; //permutation vector returned by HODLR
	int nlevel = 0; // 0: tree level, nonzero if a tree is provided 
	int* tree = new int[(int)pow(2,nlevel)]; //user provided array containing size of each leaf node, not used if nlevel=0
	tree[0] = Npo;
	int* groups = new int[size];
	int i_opt;
	double d_opt;
	
	//quantities for the first holdr
	F2Cptr ho_bf;  //HODLR returned by Fortran code 
	F2Cptr option;     //option structure returned by Fortran code 
	F2Cptr stats;      //statistics structure returned by Fortran code
	F2Cptr msh;		   //d_mesh structure returned by Fortran code
	F2Cptr kerquant;   //kernel quantities structure returned by Fortran code
	F2Cptr ptree;      //process tree returned by Fortran code

	// quantities for the second holdr
	F2Cptr ho_bf1;  //HODLR returned by Fortran code 
	F2Cptr option1;     //option structure returned by Fortran code 
	F2Cptr stats1;      //statistics structure returned by Fortran code
	F2Cptr msh1;		   //d_mesh structure returned by Fortran code
	F2Cptr kerquant1;   //kernel quantities structure returned by Fortran code
	F2Cptr ptree1;      //process tree returned by Fortran code

	
	
	MPI_Fint Fcomm;  // the fortran MPI communicator
	Fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);  
	
	for (int i = 0; i < size; i++)groups[i]=i;
	// create hodlr data structures
	d_c_hodlr_createptree(&size, groups, &Fcomm, &ptree);
	d_c_hodlr_createoption(&option);	
	d_c_hodlr_createstats(&stats);		
	
	// set hodlr options
	d_c_hodlr_set_D_option(&option, "tol_comp", tol);
	d_c_hodlr_set_I_option(&option, "nogeo", nogeo);
	d_c_hodlr_set_I_option(&option, "Nmin_leaf", Nmin); 
	d_c_hodlr_set_I_option(&option, "RecLR_leaf", com_opt); 
	d_c_hodlr_set_I_option(&option, "xyzsort", sort_opt); 
	d_c_hodlr_set_I_option(&option, "ErrFillFull", checkerr); 
	d_c_hodlr_set_I_option(&option, "BACA_Batch", batch); 
	

    // construct hodlr with geometrical points	
	d_c_hodlr_construct(&Npo, &Ndim, dat_ptr, &nlevel, tree, perms, &myseg, &ho_bf, &option, &stats, &msh, &kerquant, &ptree, &C_FuncZmn, quant_ptr, &Fcomm);	
	
	// factor hodlr
	d_c_hodlr_factor(&ho_bf,&option,&stats,&ptree,&msh);

	// solve the system 
	int nrhs=1;
	double* b = new double[nrhs*myseg];
	double* x = new double[nrhs*myseg];
	
	for (int i = 0; i < nrhs*myseg; i++){
		b[i]=1;
	}	
	d_c_hodlr_solve(x,b,&myseg,&nrhs,&ho_bf,&option,&stats,&ptree);

	
	// use resulting hodlr as matvec to create a new holdr
	
	quant_ptr1=new C_QuantApp();	
	quant_ptr1->ho_bf=&ho_bf;
	quant_ptr1->msh=&msh;
	quant_ptr1->ptree=&ptree;
	quant_ptr1->stats=&stats;
	quant_ptr1->option=&option;
	
	d_c_hodlr_createptree(&size, groups, &Fcomm, &ptree1);
	d_c_hodlr_copyoption(&option,&option1);	
	d_c_hodlr_createstats(&stats1);	
	
	d_c_hodlr_set_I_option(&option1, "nogeo", 1); // no geometrical information
	d_c_hodlr_set_I_option(&option1, "xyzsort", 0);// natural ordering	
	
	
	int Npo1 = Npo; 
	int myseg1=0;     // local number of unknowns
	int* perms1 = new int[Npo1]; //permutation vector returned by HODLR
	//tree1 and nlevel1 should be provided by the caller, otherwise natural ordering is used
	int nlevel1 = 0; // 0: tree level, nonzero if a tree is provided 
	int* tree1 = new int[(int)pow(2,nlevel1)]; //user provided array containing size of each leaf node, not used if nlevel=0
	tree1[0] = Npo1;
	// construct hodlr from blackbox matvec
	d_c_hodlr_construct_matvec(&Npo1, &nlevel1, tree1, perms1, &myseg1, &ho_bf1, &option1, &stats1, &msh1, &kerquant1, &ptree1, &C_FuncMatVec, quant_ptr1, &Fcomm);

	
	d_c_hodlr_deletestats(&stats1);
	d_c_hodlr_deleteproctree(&ptree1);
	d_c_hodlr_deletemesh(&msh1);
	d_c_hodlr_deletekernelquant(&kerquant1);
	d_c_hodlr_deletehobf(&ho_bf1);
	d_c_hodlr_deleteoption(&option1);
	
	
	
	d_c_hodlr_deletestats(&stats);
	d_c_hodlr_deleteproctree(&ptree);
	d_c_hodlr_deletemesh(&msh);
	d_c_hodlr_deletekernelquant(&kerquant);
	d_c_hodlr_deletehobf(&ho_bf);
	d_c_hodlr_deleteoption(&option);
	
	delete quant_ptr;
	delete quant_ptr1;
	delete perms;
	delete perms1;
	delete tree;
	delete tree1;
	
	
	Cblacs_exit(1);
	MPI_Finalize();                                 // Terminate MPI. Once called, no other MPI routines may be called
    return 0;
}
//------------------------------------------------------------------------------
