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
#include <complex.h>


#include "G3D/rt_nonfinite.h"
#include "G3D/G3D.h"


#include "zC_BPACK_wrapper.h"



//------------------------------------------------------------------------------
using namespace std;
const int pi = 4.0*atan(1.0);



// The object handling kernel parameters and sampling function
class C_QuantApp {
public:
  vector<double> _data;
  int _d = 0;   // data dimension 2 or 3
  int _n = 0;   // size of the matrix
  double _w = pi;   //angular frequency

  F2Cptr* bmat;  //hierarchical matrix returned by Fortran code
  F2Cptr* bf;  //BF returned by Fortran code
  F2Cptr* stats;      //statistics structure returned by Fortran code
  F2Cptr* msh;		   //mesh structure returned by Fortran code
  F2Cptr* ptree;      //process tree returned by Fortran code
  F2Cptr* option;      //option structure returned by Fortran code


  C_QuantApp() = default;

  C_QuantApp(vector<double> data, int d, double w)
    : _data(move(data)), _d(d), _n(_data.size() / _d),_w(w){
    assert(size_t(_n * _d) == _data.size());
	}


  inline void Sample(int m, int n, _Complex double* val){
	  if(m==n){
		  *val = 0;  //self term is not computed yet
	  }else{
		creal_T out =G3D(&_data[m * _d], &_data[n * _d], _w);
		*val = (double)out.re+ I*(double)out.im;
		// if(abs(m-n)<=10)
		// 	std::cout<<m<<" "<<_data[m * _d]<<" "<<_data[m * _d+1]<<" "<<_data[m * _d+2]<<" "<<_data[n * _d]<<" "<<_data[n * _d+1]<<" "<<_data[n * _d+2]<<" "<<n<<" "<<out.re<<""<<out.im<<std::endl;
	  }
  }
};


// The sampling function wrapper required by the Fortran HODLR code
inline void C_FuncZmn(int *m, int *n, _Complex double *val, C2Fptr quant) {

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
inline void C_FuncZmnBlock(int* Ninter, int* Nallrows, int* Nallcols, int* Nalldat_loc, int* allrows, int* allcols, _Complex double* alldat_loc, int* rowidx,int* colidx, int* pgidx, int* Npmap, int* pmaps, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;
  z_c_bpack_extractelement(Q->bmat,Q->option,Q->msh,Q->stats,Q->ptree,Ninter,Nallrows, Nallcols, Nalldat_loc, allrows,allcols,alldat_loc,rowidx,colidx,pgidx,Npmap,pmaps);
}



// The matvec sampling function wrapper required by the Fortran HODLR code
inline void C_FuncHMatVec(char const *trans, int *nin, int *nout, int *nvec, _Complex double const *xin, _Complex double *xout, C2Fptr quant) {
  C_QuantApp* Q = (C_QuantApp*) quant;
  z_c_bpack_mult(trans, xin, xout, nin, nout, nvec, Q->bmat,Q->option,Q->stats,Q->ptree);
}


// The matvec sampling function wrapper required by the Fortran HODLR code
inline void C_FuncBMatVec(char const *trans, int *nin, int *nout, int *nvec, _Complex double const *xin, _Complex double *xout, C2Fptr quant, _Complex double *a, _Complex double *b) {
  C_QuantApp* Q = (C_QuantApp*) quant;
  int cnt = (*nvec)*(*nout);
  _Complex double* xout1 = new _Complex double[cnt];

  z_c_bpack_mult(trans, xin, xout1, nin, nout, nvec, Q->bmat,Q->option,Q->stats,Q->ptree);

  for (int ii=0; ii<cnt; ii++){
	xout[ii] = *b*xout[ii] + *a*xout1[ii];
  }
  delete[] xout1;
}


inline void readoption(F2Cptr* option, int *ii, int argc, char *argv[]){
	string strings("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");
	int flag=1;
	int i_opt;
	double d_opt;
	flag=1;
	while(flag==1){
		(*ii)++;
		if(*ii==argc)
			return;
	 	strings = string(argv[*ii]);
		string r = strings.substr(0, 2);

		if(r.compare("--") == 0){
			if(strings.compare("--nmin_leaf") == 0){
				(*ii)++;
				i_opt = stoi(argv[*ii]);
				z_c_bpack_set_I_option(option, "Nmin_leaf", i_opt);
			}else if(strings.compare("--tol_comp") == 0){
				(*ii)++;
				d_opt = stof(argv[*ii]);
				z_c_bpack_set_D_option(option, "tol_comp", d_opt);
				z_c_bpack_set_D_option(option, "tol_rand", d_opt);
				z_c_bpack_set_D_option(option, "tol_Rdetect", d_opt*0.1);
			}else if(strings.compare("--tol_itersol") == 0){
				(*ii)++;
				d_opt = stof(argv[*ii]);
				z_c_bpack_set_D_option(option, "tol_itersol", d_opt);
			}else if(strings.compare("--precon") == 0){
				(*ii)++;
				i_opt = stoi(argv[*ii]);
				z_c_bpack_set_I_option(option, "precon", i_opt);
			}else if(strings.compare("--xyzsort") == 0){
				(*ii)++;
				i_opt = stoi(argv[*ii]);
				z_c_bpack_set_I_option(option, "xyzsort", i_opt);
			}else if(strings.compare("--lrlevel") == 0){
				(*ii)++;
				i_opt = stoi(argv[*ii]);
				z_c_bpack_set_I_option(option, "LRlevel", i_opt);
			}else if(strings.compare("--errfillfull") == 0){
				(*ii)++;
				i_opt = stoi(argv[*ii]);
				z_c_bpack_set_I_option(option, "ErrFillFull", i_opt);
			}else if(strings.compare("--baca_batch") == 0){
				(*ii)++;
				i_opt = stoi(argv[*ii]);
				z_c_bpack_set_I_option(option, "BACA_Batch", i_opt);
			}else if(strings.compare("--reclr_leaf") == 0){
				(*ii)++;
				i_opt = stoi(argv[*ii]);
				z_c_bpack_set_I_option(option, "RecLR_leaf", i_opt);
			}else if(strings.compare("--nogeo") == 0){
				(*ii)++;
				i_opt = stoi(argv[*ii]);
				z_c_bpack_set_I_option(option, "nogeo", i_opt);
			}else  if(strings.compare("--lr_blk_num") == 0){
				(*ii)++;
				i_opt = stoi(argv[*ii]);
				z_c_bpack_set_I_option(option, "lr_blk_num", i_opt);
			}else  if(strings.compare("--rank0") == 0){
				(*ii)++;
				i_opt = stoi(argv[*ii]);
				z_c_bpack_set_I_option(option, "rank0", i_opt);
			}else  if(strings.compare("--format") == 0){
				(*ii)++;
				i_opt = stoi(argv[*ii]);
				z_c_bpack_set_I_option(option, "format", i_opt);
			}else  if(strings.compare("--verbosity") == 0){
				(*ii)++;
				i_opt = stoi(argv[*ii]);
				z_c_bpack_set_I_option(option, "verbosity", i_opt);
			}else  if(strings.compare("--pat_comp") == 0){
				(*ii)++;
				i_opt = stoi(argv[*ii]);
				z_c_bpack_set_I_option(option, "pat_comp", i_opt);
			}else  if(strings.compare("--rankrate") == 0){
				(*ii)++;
				d_opt = stof(argv[*ii]);
				z_c_bpack_set_D_option(option, "rankrate", d_opt);
			}else  if(strings.compare("--knn") == 0){
				(*ii)++;
				i_opt = stoi(argv[*ii]);
				z_c_bpack_set_I_option(option, "knn", i_opt);
			}else  if(strings.compare("--cpp") == 0){
				(*ii)++;
				i_opt = stoi(argv[*ii]);
				z_c_bpack_set_I_option(option, "cpp", i_opt);
			}else  if(strings.compare("--elem_extract") == 0){
				(*ii)++;
				i_opt = stoi(argv[*ii]);
				z_c_bpack_set_I_option(option, "elem_extract", i_opt);
			}else{
				std::cout<<"Ignoring unknown option: "<<strings<<std::endl;
				(*ii)++;
			}
		}else{
			flag=0;
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
	int ker=1 ; // kernel choice
	int Npo=2500;   // matrix size
	int Ndim=3; //data dimension
	double starttime, endtime;
	double* dat_ptr;
	double w=16*pi; //angular frequency
	int* nns_ptr;
	int nogeo;  // 1: no geometrical information passed to hodlr, dat_ptr and Ndim are dummy


	int Nmin=100; //finest leafsize
	double tol=1e-4; //compression tolerance
	int com_opt=5; //1:SVD 2:RRQR 3:ACA 4:BACA 5:BACA_improved 6:Pseudo-skeleton
	int sort_opt=1; //0:natural order 1:kd-tree 2:cobble-like ordering 3:gram distance-based cobble-like ordering
	int checkerr = 0; //1: check compression quality
	int batch = 100; //batch size for BACA
	int bnum = 1; //sqrt of #of subblocks in H-BACA
	int knn=0; //k nearest neighbours stored per point
	C_QuantApp *quant_ptr;
	int v_major,v_minor,v_bugfix; //version numbers
	double eps=0.000001;

	//quantities for the first holdr
	F2Cptr bmat;  //hierarchical matrix returned by Fortran code
	F2Cptr option;     //option structure returned by Fortran code
	F2Cptr stats;      //statistics structure returned by Fortran code
	F2Cptr msh;		   //d_mesh structure returned by Fortran code
	F2Cptr kerquant;   //kernel quantities structure returned by Fortran code
	F2Cptr ptree;      //process tree returned by Fortran code


	MPI_Fint Fcomm;  // the fortran MPI communicator
	Fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);


if(myrank==master_rank){
	z_c_bpack_getversionnumber(&v_major,&v_minor,&v_bugfix);
	std::cout<<"ButterflyPACK Version: "<<v_major<<"."<<v_minor<<"."<<v_bugfix<<std::endl;
}

	z_c_bpack_createoption(&option);

	string strings("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");
	int ii=1;
	while(ii<argc){
	 	strings = string(argv[ii]);
		if(strings.compare("-quant") == 0){
			int flag=1;
			int i_opt;
			double d_opt;
			flag=1;
			while(flag==1){
				ii++;
				if(ii==argc)
					break;
				strings = string(argv[ii]);
				string r = strings.substr(0, 2);

				if(r.compare("--") == 0){
					if(strings.compare("--nunk") == 0){
						ii++;
						i_opt = stoi(argv[ii]);
						Npo = i_opt;
						Npo=round(pow(round(sqrt(Npo)),2.0));
					}else if(strings.compare("--omega") == 0){
						ii++;
						d_opt = stof(argv[ii]);
						w=d_opt;
					}else{
						std::cout<<"Ignoring unknown quant parameter: "<<strings<<std::endl;
						ii++;
					}
				}else{
					flag=0;
					ii--;
				}
			}
		}else if(strings.compare("-option") == 0){  // options of ButterflyPACK
			readoption(&option, &ii, argc, argv);
		}else{
			if(myrank==master_rank)
				std::cout<<"Ignoring unknown argument: "<<strings<<std::endl;
		}
		ii++;
	}




	/*****************************************************************/
	/* Test the Taylor-expansion-based 2D Green function */

	vector<double> data_geo;
	data_geo.resize(Ndim*Npo);
	double dtheta=(2*pi-eps)/(sqrt(Npo));
	dtheta=(2*pi-eps-dtheta)/(sqrt(Npo));
	double dphi=(2*pi-eps)/(sqrt(Npo));
	double rmax=0.25;
	double center[3];
	center[0]=0.25;
	center[1]=0.25;
	center[2]=0.25;

	double theta;
	double phi;
	
	int iijj=0;
	for(int ii=1;ii<=sqrt(Npo);ii++){
	for(int jj=1;jj<=sqrt(Npo);jj++){
		theta=(ii)*dtheta;
        phi=(jj)*dphi;
		data_geo[iijj * Ndim] =rmax*sin(theta)*cos(phi)+center[0];
		data_geo[iijj * Ndim+1] =rmax*sin(theta)*sin(phi)+center[1];
		data_geo[iijj * Ndim+2] =rmax*cos(theta)+center[2];
		iijj++;
	}
	}

	quant_ptr=new C_QuantApp(data_geo, Ndim, w);
	dat_ptr = new double[data_geo.size()];
	for(int ii=0;ii<data_geo.size();ii++)
		dat_ptr[ii] = data_geo.data()[ii];
	nogeo=0;



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
	int cpp=1; //1: use user-defined cpp/c functions for construction



	for (int i = 0; i < size; i++)groups[i]=i;
	// create hodlr data structures
	z_c_bpack_createptree(&size, groups, &Fcomm, &ptree);
	z_c_bpack_createstats(&stats);



	z_c_bpack_set_I_option(&option, "cpp", cpp);
	z_c_bpack_printoption(&option,&ptree);


    // construct hodlr with geometrical points
	z_c_bpack_construct_init(&Npo, &Ndim, dat_ptr, nns_ptr,&nlevel, tree, perms, &myseg, &bmat, &option, &stats, &msh, &kerquant, &ptree, &C_FuncDistmn, &C_FuncNearFar, quant_ptr);

	z_c_bpack_construct_element_compute(&bmat, &option, &stats, &msh, &kerquant, &ptree, &C_FuncZmn, &C_FuncZmnBlock, quant_ptr);

	// // factor hodlr
	// z_c_bpack_factor(&bmat,&option,&stats,&ptree,&msh);

	// // solve the system
	// int nrhs=1;
	// _Complex double* b = new _Complex double[nrhs*myseg];
	// _Complex double* x = new _Complex double[nrhs*myseg];

	// for (int i = 0; i < nrhs*myseg; i++){
	// 	b[i]=1;
	// }
	// z_c_bpack_solve(x,b,&myseg,&nrhs,&bmat,&option,&stats,&ptree);


	if(myrank==master_rank)std::cout<<"Printing stats of the first HODLR: "<<std::endl;
	z_c_bpack_printstats(&stats,&ptree);



	z_c_bpack_deletestats(&stats);
	z_c_bpack_deleteproctree(&ptree);
	z_c_bpack_deletemesh(&msh);
	z_c_bpack_deletekernelquant(&kerquant);
	z_c_bpack_delete(&bmat);
	z_c_bpack_deleteoption(&option);

	delete quant_ptr;
	delete[] perms;
	delete[] tree;


	MPI_Finalize();                                 // Terminate MPI. Once called, no other MPI routines may be called
    return 0;
}
//------------------------------------------------------------------------------
