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

#include <pthread.h>

#include <cmath>
#include <cassert>
#include <iostream>
#include <random>
#include <vector>
#include <atomic>
#include <mpi.h>
#include <memory>

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

  C_QuantApp(int n, int ker, vector<double> MatFull)
    : _n(n), _ker(ker), _MatFull(move(MatFull)){
	// cout<<_n_rand<<_rank_rand<<_MatU.size()<<endl;
    assert(size_t(_n * _n) == _MatFull.size());
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
		*val =_MatFull[n*_n+m];
		// *val =_MatFull[_Hperm[n]*_n+_Hperm[m]];
		break;
	}
  }
};


// The sampling function wrapper required by the Fortran HODLR code
inline void C_FuncZmn(int *m, int *n, double *val, C2Fptr quant) {

  C_QuantApp* Q = (C_QuantApp*) quant;
  Q->Sample(*m-1,*n-1,val);
}

// The distance function wrapper required by the Fortran HODLR code
inline void C_FuncDistmn(int *m, int *n, double *val, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;

}

// The compressibility function wrapper required by the Fortran HODLR code
inline void C_FuncNearFar(int *m, int *n, int *val, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;

}

// The extraction sampling function wrapper required by the Fortran HODLR code
inline void C_FuncZmnBlock(int* Ninter, int* Nallrows, int* Nallcols, int* Nalldat_loc, int* allrows, int* allcols, double* alldat_loc, int* rowidx,int* colidx, int* pgidx, int* Npmap, int* pmaps, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;
  d_c_bpack_extractelement(Q->bmat,Q->option,Q->msh,Q->stats,Q->ptree,Ninter,Nallrows, Nallcols, Nalldat_loc, allrows,allcols,alldat_loc,rowidx,colidx,pgidx,Npmap,pmaps);
}

// The extraction sampling function wrapper required by the Fortran HODLR code
inline void C_FuncBZmnBlock(int* Ninter, int* Nallrows, int* Nallcols, int* Nalldat_loc, int* allrows, int* allcols, double* alldat_loc, int* rowidx,int* colidx, int* pgidx, int* Npmap, int* pmaps, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;

  d_c_bf_extractelement(Q->bf,Q->option,Q->msh,Q->stats,Q->ptree,Ninter,Nallrows, Nallcols, Nalldat_loc, allrows,allcols,alldat_loc,rowidx,colidx,pgidx,Npmap,pmaps);

  // // for(int ii=0;ii<*Nallrows;ii++)cout<<allrows[ii]<<endl;
  // d_c_bpack_extractelement(Q->bmat,Q->option,Q->msh,Q->stats,Q->ptree,Ninter,Nallrows, Nallcols, Nalldat_loc, allrows,allcols,alldat_loc,rowidx,colidx,pgidx,Npmap,pmaps);
}


// The matvec sampling function wrapper required by the Fortran HODLR code
inline void C_FuncHMatVec(char const *trans, int *nin, int *nout, int *nvec, double const *xin, double *xout, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;
  d_c_bpack_mult(trans, xin, xout, nin, nout, nvec, Q->bmat,Q->option,Q->stats,Q->ptree);
}


// The matvec sampling function wrapper required by the Fortran HODLR code
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
	int format=1; //1: HOD-LR/BF  2:H-LR/BF
	double near_para=0.01; // control the admissibility
	C_QuantApp *quant_ptr;
	int v_major,v_minor,v_bugfix; //version numbers

    int tst = 4;
	int lrlevel=0;

	int nlevel = 0; // 0: tree level, nonzero if a tree is provided
	int* tree = new int[(int)pow(2,nlevel)]; //user provided array containing size of each leaf node, not used if nlevel=0
	tree[0] = Npo;

	if(argc>1)tst = stoi(argv[1]);


if(myrank==master_rank){
	d_c_bpack_getversionnumber(&v_major,&v_minor,&v_bugfix);
	std::cout<<"ButterflyPACK Version: "<<v_major<<"."<<v_minor<<"."<<v_bugfix<<std::endl;
}


// 	/*****************************************************************/
// 	/* Test Kernels for Liza's data sets */
// if(tst==1){
// 	string filename("../EXAMPLE/KRR_DATA/susy_10Kn");
// 	Ndim=8;
// 	if(argc>2)filename = string(argv[2]);
// 	if(argc>3)Ndim = stoi(argv[3]);
// 	if(argc>4)ker = stoi(argv[4]);
// 	if(argc>5)h = stof(argv[5]);
// 	if(argc>6)lambda = stof(argv[6]);

// 	if(argc>7)Nmin = stoi(argv[7]);
// 	if(argc>8)tol = stof(argv[8]);
// 	if(argc>9)com_opt = stoi(argv[9]);
// 	if(argc>10)checkerr = stoi(argv[10]);
// 	if(argc>11)batch = stoi(argv[11]);
// 	if(argc>12)bnum = stoi(argv[12]);
// 	if(argc>13)knn = stoi(argv[13]);
//     vector<double> data_train = write_from_file<double>(filename + "_train.csv");
// 	Npo = data_train.size() / Ndim;

// 	quant_ptr=new C_QuantApp(data_train, Ndim, h, lambda,ker);
// 	dat_ptr = new double[data_train.size()];
// 	for(int ii=0;ii<data_train.size();ii++)
// 		dat_ptr[ii] = data_train.data()[ii];
// 	nogeo=0;
// }




	/*****************************************************************/
	/* Test Full matrices*/
if(tst==4){
	ker = 7;

	string mat_file("../EXAMPLE/FULLMAT_DATA/FULLMAT.csv");
	string geo_file("../EXAMPLE/FULLMAT_DATA/Geometry.csv");


	vector<double> geo_data = write_from_file<double>(geo_file);
	Ndim=3;
	Npo = geo_data.size()/Ndim;
	dat_ptr = new double[geo_data.size()];
	for(int ii=0;ii<geo_data.size();ii++)
		dat_ptr[ii] = geo_data.data()[ii];
	nogeo=0;

	vector<double> mat_data = write_from_file<double>(mat_file);
	quant_ptr=new C_QuantApp(Npo, ker, mat_data);


	nogeo=0;
	sort_opt=1;
	Nmin = 50;

	// Command line options

	if(argc>2)Nmin = stoi(argv[2]); 		 // leafsize

	if(argc>3)tol = stof(argv[3]);   		 // compression tolerance
	if(argc>4)com_opt = stoi(argv[4]);       // compression algorithm if LR is used
	if(argc>5)checkerr = stoi(argv[5]);      // whether to fully check the compression accuracy
	if(argc>6)batch = stoi(argv[6]);         // the block size in H-BACA if com_opt=4 or 5
	if(argc>7)bnum = stoi(argv[7]);          // number of bottom level blocks in H-BACA if com_opt=4 or 5
	if(argc>8)knn = stoi(argv[8]);           // number of nearest neighbours per point to be found
	if(argc>9)lrlevel = stoi(argv[9]);       // the top "lrlevel" levels of the matrix uses butterfly, the rest uses low-rank
	if(argc>10)sample_para = stof(argv[10]); // oversampling parameter in butterfly compression
	if(argc>11)format = stoi(argv[11]);      // 1: HOD-LR/BF  2:H-LR/BF
	if(argc>12)near_para = stof(argv[12]);   // control the admissibility: e.g. 0.01: weak admissibility 2.01: strong admissibility
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
	d_c_bpack_set_I_option(&option, "ErrSol", 1);
	d_c_bpack_set_I_option(&option, "format", format);

	d_c_bpack_printoption(&option,&ptree);

    // construct hodlr with geometrical points
	d_c_bpack_construct_init(&Npo, &Ndim, dat_ptr, nns_ptr,&nlevel, tree, perms, &myseg, &bmat, &option, &stats, &msh, &kerquant, &ptree, &C_FuncDistmn, &C_FuncNearFar, quant_ptr);
	d_c_bpack_construct_element_compute(&bmat, &option, &stats, &msh, &kerquant, &ptree, &C_FuncZmn, &C_FuncZmnBlock, quant_ptr);

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


	if(myrank==master_rank)std::cout<<"Printing stats of the first HODLR: "<<std::endl;
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



	Cblacs_exit(1);
	MPI_Finalize();                                 // Terminate MPI. Once called, no other MPI routines may be called
    return 0;
}
//------------------------------------------------------------------------------
