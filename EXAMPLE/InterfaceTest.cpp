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

#include "../SRC/CInterface_config.h"



//------------------------------------------------------------------------------
using namespace std;


inline double dist2(double *x, double *y, int d) {
  double k = 0.;
  for (int i = 0; i < d; i++)
    k += pow(x[i] - y[i], 2.);
  return k;
}

inline double Gauss_kernel(double *x, double *y, int d, double h) {
  double dists;
  dists = dist2(x, y, d);
  if(dists> -log(1e-30)*2.0*pow(h,2.0)){
	return 0;
  }else{
	return exp(-dists / (2. * h * h));
  }
}



class C_QuantZmn {

public:
  vector<double> _data;
  int _d = 0;
  int _n = 0;
  double _h = 0.;
  double _l = 0.;
  std::vector<int> _Hperm;
  std::vector<int> _iHperm;
  int _nloc = 0;

  C_QuantZmn() = default;
  
  C_QuantZmn(vector<double> data, int d, double h, double l)
    : _data(move(data)), _d(d), _n(_data.size() / _d),
      _h(h), _l(l){
	  // std::cout<<"_h "<<_h<<"_n "<<_n<<" _d "<<_d<<"size "<<_data.size()<<std::endl;
    assert(size_t(_n * _d) == _data.size());
	}
  
  inline void Sample(int m, int n, int kerchoice, doublecomplex* val){
	// define your functions here, kerchoice=1: Gaussian kernel
	if(kerchoice ==1){ 
	val->i=0;
	val->r = Gauss_kernel(&_data[m * _d], &_data[n * _d], _d, _h);
	if (m==n)
	val->r += _l;
	}
  } 	
};

inline void C_FuncZmn(int *m, int *n, doublecomplex *val, C2Fptr quant) {
	
  C_QuantZmn* Q = (C_QuantZmn*) quant;	
  // val->r=1;
  // val->i=2;
  // std::cout<<"good"<<std::endl;
  // std::cout<<"h "<<Q->_h<<" n "<<Q->_n<<" d "<<Q->_d<<std::endl;
  
  int kerchoice=1; //kerchoice=1: Gaussian kernel
  Q->Sample(*m,*n,kerchoice,val);
}


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

    int myrank, size;                     // Store values of proceeor rank and total no of procs requestedss
    int master_rank = 0;
	MPI_Init(&argc, &argv); 	                            // Initialize MPI, called only once
    MPI_Comm_size(MPI_COMM_WORLD, &size); 	                // Get no of procs
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 	                // Get no of procs
    MPI_Op op;
	double h = 0.1;
	double lambda = 10.;
  
	string filename("smalltest.dat");
	
	int Npo,Ndim;
    double starttime, endtime;
	if (argc > 1)
	filename = string(argv[1]);
	if (argc > 2)
	Ndim = stoi(argv[2]);
	if (argc > 3)
	h = stof(argv[3]);
	if (argc > 4)
	lambda = stof(argv[4]);
	
    vector<double> data_train = write_from_file(filename + "_train.csv");
    Npo = data_train.size() / Ndim;
  
	if(myrank==master_rank)std::cout<<"Npo "<<Npo<<" Ndim "<<Ndim<<std::endl;
	
	C_QuantZmn quant(data_train, Ndim, h, lambda);	
	
	
   starttime = MPI_Wtime();

	int Nmin=200; //finest leafsize 
	double tol=1e-4; //compression tolerance
	int nth=1;   //OMP_NUM_THREADS
	int nmpi=size;	 //number of active MPI ranks 
	int ninc=1;      // use one rank in H per ninc ranks in HSS
	// int nmpi=2;	 //number of active MPI ranks 
	// int ninc=4;
	int myseg=0;
	int* perms = new int[Npo];
	int* tree = new int[Npo];
	int preorder=0;  // 0: holdr code will reorder the elements and return a permutation 1: holdr code will use user-supplied tree to create the holdr structure
	
	F2Cptr ho_bf_for;
	F2Cptr ho_bf_inv;
	F2Cptr option;
	F2Cptr stats;
	F2Cptr msh;
	F2Cptr ker;
	F2Cptr ptree; 
	
	MPI_Fint Fcomm;  // the fortran MPI communicator
	
	char* ttemp = getenv("OMP_NUM_THREADS");
	if(ttemp)
	{
		nth = atoi(ttemp);
	}	
	
	Fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);  
	
    // construct hodlr with geometrical points	
	FC_GLOBAL_(c_hodlr_fill,C_HODLR_FILL)(&Npo, &Ndim, data_train.data(), &Nmin, &tol, &nth, &nmpi, &ninc, &preorder, tree, perms, &myseg, &ho_bf_for, &option, &stats, &msh, &ker, &ptree, &C_FuncZmn, &quant, &Fcomm);	
	
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
	

	
	
	// printf("Myid: %5d Segment size: %5d\n",myrank, myseg);
	// fflush(stdout);
	
   // endtime   = MPI_Wtime();
   // if(myrank==master_rank)printf("Fill time %f seconds\n",endtime-starttime);	
	
	
	 // int Ncol=100;
	// //int Ncol=1;
	
	// double *xin ;
	// double *xout ;
	
	// if(myseg>0){
		// xin = new double[myseg*Ncol];
		// xout = new double[myseg*Ncol];
		// for (int i=0;i<myseg*Ncol;i++)xin[i]=1.0;
	// }		


   // starttime = MPI_Wtime();

	
	// // h_matrix_apply_(&Npo, &Ncol, xin, xout);	
	// FC_GLOBAL_(h_matrix_apply,H_MATRIX_APPLY)(&Npo, &Ncol, xin, xout);		

   // endtime   = MPI_Wtime();
   // if(myrank==master_rank)printf("Apply time %f seconds, per vector %f seconds\n",(endtime-starttime), (endtime-starttime)/Ncol);	

	// // for (int i=0;i<Npo*Ncol;i++)printf("%16.14f\n",xout[i]);
	
	// // for (int i=0;i<Npo*Ncol;i++)std::cout<<setprecision(14)<<xout[i]<<std::endl;
	
	
	MPI_Finalize();                                 // Terminate MPI. Once called, no other MPI routines may be called
    return 0;
}
//------------------------------------------------------------------------------
