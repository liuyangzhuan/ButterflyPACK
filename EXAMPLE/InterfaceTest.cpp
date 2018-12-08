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

#include "dC_BPACK_wrapper.h"



//------------------------------------------------------------------------------
using namespace std;


extern "C" {
      ///////////////////////////////////////////////
      ////// BLACS //////////////////////////////////
      ///////////////////////////////////////////////
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
inline void C_FuncHMatVec(char const *trans, int *nin, int *nout, int *nvec, double const *xin, double *xout, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;	
  d_c_bpack_mult(trans, xin, xout, nin, nout, nvec, Q->ho_bf,Q->option,Q->stats,Q->ptree);  
}

// The matvec sampling function wrapper required by the Fortran HODLR code
inline void C_FuncBMatVec(char const *trans, int *nin, int *nout, int *nvec, double const *xin, double *xout, C2Fptr quant, double *a, double *b) {
  C_QuantApp* Q = (C_QuantApp*) quant;	
  int cnt = (*nvec)*(*nout);
  double* xout1 = new double[cnt];
     
  d_c_bpack_mult(trans, xin, xout1, nin, nout, nvec, Q->ho_bf,Q->option,Q->stats,Q->ptree);  
  
  for (int ii=0; ii<cnt; ii++){
	xout[ii] = *b*xout[ii] + *a*xout1[ii];
  }
  delete xout1;
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
	double h=0.1; //kernel parameter
	double lambda=10.0 ; //kernel parameter 
	int ker=1 ; // kernel choice 
	int Npo=5000;   // matrix size
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
	
	C_QuantApp *quant_ptr;
	
	
    int tst = 1;
	if(argc>1)tst = stoi(argv[1]);
	
	/*****************************************************************/
	/* Test Kernels for Liza's data sets */
if(tst==1){
	string filename("../EXAMPLE/KRR_DATA/susy_10Kn");
	Ndim=8;
	if(argc>2)filename = string(argv[2]);
	if(argc>3)Ndim = stoi(argv[3]);
	if(argc>4)ker = stoi(argv[4]);
	if(argc>5)h = stof(argv[5]);
	if(argc>6)lambda = stof(argv[6]);
	if(argc>7)Nmin = stoi(argv[7]);
	if(argc>8)tol = stof(argv[8]);
	if(argc>9)com_opt = stoi(argv[9]);
	if(argc>10)checkerr = stoi(argv[10]);
	if(argc>11)batch = stoi(argv[11]);
	
    vector<double> data_train = write_from_file(filename + "_train.csv");
	Npo = data_train.size() / Ndim;

	quant_ptr=new C_QuantApp(data_train, Ndim, h, lambda,ker);	
	dat_ptr = new double[data_train.size()];
	for(int ii=0;ii<data_train.size();ii++)
		dat_ptr[ii] = data_train.data()[ii];
	nogeo=0;
}
	
	/*****************************************************************/
	/* Test Kernels for Random point clouds */
if(tst==2){	
	Ndim=6;
	if(argc>2)Npo = stoi(argv[2]);
	if(argc>3)Ndim = stoi(argv[3]);
	if(argc>4)ker = stoi(argv[4]);
	if(argc>5)h = stof(argv[5]);
	if(argc>6)lambda = stof(argv[6]);
	if(argc>7)Nmin = stoi(argv[7]);
	if(argc>8)tol = stof(argv[8]);
	if(argc>9)com_opt = stoi(argv[9]);
	if(argc>10)checkerr = stoi(argv[10]);
	if(argc>11)batch = stoi(argv[11]);
	
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
	/* Test Product of two Random matrices*/
if(tst==3){
	ker = 6;
	int rank_rand = 10;
	Npo = 10000;
	
	if(argc>2)Npo = stoi(argv[2]);
	if(argc>3)rank_rand = stoi(argv[3]);
	if(argc>4)Nmin = stoi(argv[4]);
	if(argc>5)tol = stof(argv[5]);
	if(argc>6)com_opt = stoi(argv[6]);
	if(argc>7)checkerr = stoi(argv[7]);
	if(argc>8)batch = stoi(argv[8]);		
	
	
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
	
	
	MPI_Fint Fcomm;  // the fortran MPI communicator
	Fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);  
	
	for (int i = 0; i < size; i++)groups[i]=i;
	// create hodlr data structures
	d_c_bpack_createptree(&size, groups, &Fcomm, &ptree);
	d_c_bpack_createoption(&option);	
	d_c_bpack_createstats(&stats);		
	
	// set hodlr options
	d_c_bpack_set_D_option(&option, "tol_comp", tol);
	d_c_bpack_set_I_option(&option, "nogeo", nogeo);
	d_c_bpack_set_I_option(&option, "Nmin_leaf", Nmin); 
	d_c_bpack_set_I_option(&option, "RecLR_leaf", com_opt); 
	d_c_bpack_set_I_option(&option, "xyzsort", sort_opt); 
	d_c_bpack_set_I_option(&option, "ErrFillFull", checkerr); 
	d_c_bpack_set_I_option(&option, "BACA_Batch", batch); 
	

    // construct hodlr with geometrical points	
	d_c_bpack_construct(&Npo, &Ndim, dat_ptr, &nlevel, tree, perms, &myseg, &ho_bf, &option, &stats, &msh, &kerquant, &ptree, &C_FuncZmn, quant_ptr, &Fcomm);	
	
	// factor hodlr
	d_c_bpack_factor(&ho_bf,&option,&stats,&ptree,&msh);

	// solve the system 
	int nrhs=1;
	double* b = new double[nrhs*myseg];
	double* x = new double[nrhs*myseg];
	
	for (int i = 0; i < nrhs*myseg; i++){
		b[i]=1;
	}	
	d_c_bpack_solve(x,b,&myseg,&nrhs,&ho_bf,&option,&stats,&ptree);

	
	
	
	
	//////////////////// use resulting hodlr as matvec to create a new holdr
	
	// quantities for the second holdr
	F2Cptr ho_bf1;  //HODLR returned by Fortran code 
	F2Cptr option1;     //option structure returned by Fortran code 
	F2Cptr stats1;      //statistics structure returned by Fortran code
	F2Cptr msh1;		   //d_mesh structure returned by Fortran code
	F2Cptr kerquant1;   //kernel quantities structure returned by Fortran code
	F2Cptr ptree1;      //process tree returned by Fortran code
	C_QuantApp *quant_ptr1; //user-defined object
	
	quant_ptr1=new C_QuantApp();	
	quant_ptr1->ho_bf=&ho_bf;
	quant_ptr1->msh=&msh;
	quant_ptr1->ptree=&ptree;
	quant_ptr1->stats=&stats;
	quant_ptr1->option=&option;
	
	d_c_bpack_createptree(&size, groups, &Fcomm, &ptree1);
	d_c_bpack_copyoption(&option,&option1);	
	d_c_bpack_createstats(&stats1);	
	
	d_c_bpack_set_I_option(&option1, "nogeo", 1); // no geometrical information
	d_c_bpack_set_I_option(&option1, "xyzsort", 0);// natural ordering	
	
	int Npo1 = Npo; 
	int myseg1=0;     // local number of unknowns
	int* perms1 = new int[Npo1]; //permutation vector returned by HODLR
	//tree1 and nlevel1 should be provided by the caller, otherwise natural ordering is used
	int nlevel1 = 0; // 0: tree level, nonzero if a tree is provided 
	int* tree1 = new int[(int)pow(2,nlevel1)]; //user provided array containing size of each leaf node, not used if nlevel=0
	tree1[0] = Npo1;
	
	d_c_bpack_construct_matvec_init(&Npo1, &nlevel1, tree1, perms1, &myseg1, &ho_bf1, &option1, &stats1, &msh1, &kerquant1, &ptree1);
	d_c_bpack_construct_matvec_compute(&ho_bf1, &option1, &stats1, &msh1, &kerquant1, &ptree1, &C_FuncHMatVec, quant_ptr1);

	d_c_bpack_deletestats(&stats1);
	d_c_bpack_deleteproctree(&ptree1);
	d_c_bpack_deletemesh(&msh1);
	d_c_bpack_deletekernelquant(&kerquant1);
	d_c_bpack_deletehobf(&ho_bf1);
	d_c_bpack_deleteoption(&option1);
	delete quant_ptr1;
	delete perms1;
	delete tree1;	
	
	
	
	
if(tst==3){	
	//////////////////// use resulting hodlr as matvec to create a new bf
	
	// quantities for the second holdr
	F2Cptr bf;  //BF returned by Fortran code 
	F2Cptr option2;     //option structure returned by Fortran code 
	F2Cptr stats2;      //statistics structure returned by Fortran code
	F2Cptr msh2;		   //d_mesh structure returned by Fortran code
	F2Cptr kerquant2;   //kernel quantities structure returned by Fortran code
	F2Cptr ptree2;      //process tree returned by Fortran code
	C_QuantApp *quant_ptr2; //user-defined object
	
	quant_ptr2=new C_QuantApp();	
	quant_ptr2->ho_bf=&ho_bf;
	quant_ptr2->msh=&msh;
	quant_ptr2->ptree=&ptree;
	quant_ptr2->stats=&stats;
	quant_ptr2->option=&option;
	
	d_c_bpack_createptree(&size, groups, &Fcomm, &ptree2);
	d_c_bpack_copyoption(&option,&option2);	
	d_c_bpack_createstats(&stats2);	
	
	d_c_bpack_set_I_option(&option2, "nogeo", 1); // no geometrical information
	d_c_bpack_set_I_option(&option2, "xyzsort", 0);// natural ordering	
	
	int M = Npo; 
	int N = Npo; 
	int myrow=0;     // local number of rows
	int mycol=0;     // local number of columns
	
	d_c_bf_construct_matvec_init(&M, &N, &myrow, &mycol, &msh, &msh, &bf, &option2, &stats2, &msh2, &kerquant2, &ptree2);
	d_c_bf_construct_matvec_compute(&bf, &option2, &stats2, &msh2, &kerquant2, &ptree2, &C_FuncBMatVec, quant_ptr2);

	d_c_bpack_deletestats(&stats2);
	d_c_bpack_deleteproctree(&ptree2);
	d_c_bpack_deletemesh(&msh2);
	d_c_bpack_deletekernelquant(&kerquant2);
	d_c_bf_deletebf(&bf);
	d_c_bpack_deleteoption(&option2);

	delete quant_ptr2;
	
}	

	
	
	d_c_bpack_deletestats(&stats);
	d_c_bpack_deleteproctree(&ptree);
	d_c_bpack_deletemesh(&msh);
	d_c_bpack_deletekernelquant(&kerquant);
	d_c_bpack_deletehobf(&ho_bf);
	d_c_bpack_deleteoption(&option);
	
	delete quant_ptr;
	delete perms;
	delete tree;

	
	
	Cblacs_exit(1);
	MPI_Finalize();                                 // Terminate MPI. Once called, no other MPI routines may be called
    return 0;
}
//------------------------------------------------------------------------------
