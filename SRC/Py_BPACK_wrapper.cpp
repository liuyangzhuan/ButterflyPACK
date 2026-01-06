#include "BPACK_wrapper.h"
// #include <mpi.h>

#include <cassert>
#include <complex>
#include <cstdint>
#include <iostream>
#include <vector>
#include <memory>
#include <cstring>
#include <getopt.h>
#include <iomanip>
#include <stdlib.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

using namespace std;


struct c_bpack_QuantApp{
    py::object py_func;
    py::object py_meta;
};

struct bpack_handle{
    F2Cptr bmat;
    F2Cptr option;
    F2Cptr stats;
    F2Cptr msh;
    F2Cptr ker;
    F2Cptr ptree;
    c_bpack_QuantApp* quant_ptr;
};



#ifdef HAVE_MPI
extern "C" {
      ///////////////////////////////////////////////
      ////// BLACS //////////////////////////////////
      ///////////////////////////////////////////////
      void Cblacs_exit(int);
}
#endif



// The sampling function wrapper required by the Fortran ButterflyPACK interface
// The function computes one entry A_{m,n}, where m,n starts at 1
inline void c_bpack_FuncZmn(int *m, int *n, C_DT *val, C2Fptr quant) {
}

// The distance function wrapper required by the Fortran ButterflyPACK interface
// not needed in this driver, but need to provide a dummy function to the interface
inline void c_bpack_FuncDistmn(int *m, int *n, double *val, C2Fptr quant) {
}

// The compressibility function wrapper required by the Fortran ButterflyPACK interface
// not needed in this driver, but need to provide a dummy function to the interface
inline void c_bpack_FuncNearFar(int *m, int *n, int *val, C2Fptr quant) {
}

// The blocked sampling function wrapper required by the Fortran ButterflyPACK interface
inline void c_bpack_FuncZmnBlock(int* Ninter, int* Nallrows, int* Nallcols, int64_t* Nalldat_loc, int* allrows, int* allcols, C_DT* alldat_loc, int* rowidx,int* colidx, int* pgidx, int* Npmap, int* pmaps, C2Fptr quant) {
    c_bpack_QuantApp* Q = (c_bpack_QuantApp*) quant;
    // int myrank, size;                     // Store values of processor rank and total no of procs requestedss
    // MPI_Comm_size(MPI_COMM_WORLD, &size); 	                // Get no of procs
    // MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 	                // Get no of procs
    int64_t idx_row=0;
    int64_t idx_col=0;
    int64_t idx_val=0;
    int NinterNew=0;
    int nrmax=0;
    int ncmax=0;
    int nvalmax=0;
    vector<int64_t> idx_row_map(*Ninter,0);
    vector<int64_t> idx_col_map(*Ninter,0);
    vector<int64_t> idx_val_map(*Ninter,0);
    vector<int64_t> inter_map(*Ninter,0);

    // The following for loop filters out the blocks that are not local to the current processor
    for (int nn=0;nn<*Ninter;nn++){
      int pp = pgidx[nn];
      int nprow = pmaps[pp];
      int npcol = pmaps[*Npmap+pp];
      int pid = pmaps[(*Npmap)*2+pp];
      int nr = rowidx[nn];
      int nc = colidx[nn];

      if(nprow*npcol==1){
        // if(myrank==pid){
          idx_row_map[NinterNew]=idx_row;
          idx_col_map[NinterNew]=idx_col;
          idx_val_map[NinterNew]=idx_val;
          idx_val+=nr*nc;
          inter_map[NinterNew]=nn;
          NinterNew++;
        // }else{
        // }
        idx_row+=nr;
        idx_col+=nc;
        nrmax = max(nr,nrmax);
        ncmax = max(nc,ncmax);
        nvalmax = max(nc*nr,nvalmax);
      }else{
        std::cout<<"nprow*npcol>1 in c_bpack_FuncZmnBlock"<<std::endl;
        exit(0);
      }
    }
    idx_row_map.resize(NinterNew);
    idx_col_map.resize(NinterNew);
    idx_val_map.resize(NinterNew);
    inter_map.resize(NinterNew);


    int *rows,*cols;
    // #pragma omp parallel private(rows,cols)
    {
    rows= (int*)malloc(nrmax*sizeof(int));
    cols= (int*)malloc(ncmax*sizeof(int));

    // #pragma omp for
    for (int nn1=0;nn1<NinterNew;nn1++){
      int nn =inter_map[nn1];
      int pp = pgidx[nn];
      int nprow = pmaps[pp];
      int npcol = pmaps[*Npmap+pp];
      int pid = pmaps[(*Npmap)*2+pp];
      int nr = rowidx[nn];
      int nc = colidx[nn];

      for (int idxr=0;idxr<nr;idxr++){
        rows[idxr]=allrows[idx_row_map[nn1]+idxr]-1; // 0-based indices

      }
      for (int idxc=0;idxc<nc;idxc++){
        cols[idxc]=allcols[idx_col_map[nn1]+idxc]-1; // 0-based indices
      }



		py::gil_scoped_acquire gil;  // REQUIRED
		// Convert rows/cols to Python lists (small cost, OK per block)
		py::list py_rows;
		py::list py_cols;

		for (int i = 0; i < nr; ++i) py_rows.append(rows[i]);
		for (int j = 0; j < nc; ++j) py_cols.append(cols[j]);

		// Call Python function
		py::object result = Q->py_func(py_rows, py_cols, Q->py_meta);

		// Expect NumPy array (nr x nc, Fortran order)
		py::array_t<C_DT, py::array::f_style | py::array::forcecast> arr(result);
		auto buf = arr.request();

		if (buf.ndim != 2 || buf.shape[0] != nr || buf.shape[1] != nc) {
			throw std::runtime_error("Python returned array with wrong shape");
		}

		// Copy into C++ buffer (already column-major)
		const C_DT* src = static_cast<const C_DT*>(buf.ptr);
		std::memcpy(&(alldat_loc[idx_val_map[nn1]]), src, sizeof(C_DT) * nr * nc);

    }
    free(rows);
    free(cols);
    }
}



#ifdef __cplusplus
extern "C"{
#endif

void py_bpack_logdet(void ** pyobj, C_DT * phase, C_RDT * logabsdet)
{
    bpack_handle* bpack_obj = (bpack_handle*)(*pyobj);
    c_bpack_logdet(phase, logabsdet,&(bpack_obj->option),&(bpack_obj->bmat));
    *pyobj = (void*)bpack_obj;
}


void py_bpack_init_compute(int Npo, int Ndim, double* dat_ptr, void ** pyobj, int* rankmax, int argc, char* argv[])
{
	pybind11::gil_scoped_acquire gil;
    bpack_handle* bpack_obj = (bpack_handle*)malloc(sizeof(bpack_handle));
	// bpack_obj->quant_ptr = (c_bpack_QuantApp*)malloc(sizeof(c_bpack_QuantApp));
	bpack_obj->quant_ptr = new c_bpack_QuantApp();

	// convert the first argc - 2 elements of argv to string and the last two to func and meta
	int c_argc = argc - 2;
    std::vector<char*> cargv;
    cargv.reserve(c_argc);
    for (int i = 0; i < c_argc; ++i) {
        PyObject* obj =
            reinterpret_cast<PyObject*>(argv[i]);
        if (!PyUnicode_Check(obj)) {
            throw std::runtime_error(
                "Command-line argument is not a Python string"
            );
        }
        const char* s = PyUnicode_AsUTF8(obj);
        if (!s) {
            throw std::runtime_error(
                "Failed to convert Python string to UTF-8"
            );
        }
        cargv.push_back(const_cast<char*>(s));
    }
	char** argv_ptr = cargv.data();
    PyObject* func_obj =
        reinterpret_cast<PyObject*>(argv[argc - 2]);
    PyObject* meta_obj =
        reinterpret_cast<PyObject*>(argv[argc - 1]);

    bpack_obj->quant_ptr->py_func =
        py::reinterpret_borrow<py::object>(func_obj);

    bpack_obj->quant_ptr->py_meta =
        py::reinterpret_borrow<py::object>(meta_obj);

    int myrank=0, size=1;                     // Store values of processor rank and total no of procs requestedss
    int master_rank = 0;

#ifdef HAVE_MPI
	// MPI_Init(&c_argc, &argv_ptr); 	                            // Initialize MPI, called only once
    MPI_Comm_size(MPI_COMM_WORLD, &size); 	                // Get no of procs
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 	                // Get no of procs
    // MPI_Op op;
#endif

	double starttime, endtime;
	int* nns_ptr;
	int nogeo=0;  // 1: no geometrical information passed to butterflypack, dat_ptr and Ndim are dummy


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
	c_bpack_QuantApp *quant_ptr=bpack_obj->quant_ptr;
	int v_major,v_minor,v_bugfix; //version numbers

	int lrlevel=0;

	int nlevel = 0; // 0: tree level, nonzero if a tree is provided
	int* tree = new int[(int)pow(2,nlevel)]; //user provided array containing size of each leaf node, not used if nlevel=0
	tree[0] = Npo;


    // if(myrank==master_rank){
    //     c_bpack_getversionnumber(&v_major,&v_minor,&v_bugfix);
    //     std::cout<<"ButterflyPACK Version: "<<v_major<<"."<<v_minor<<"."<<v_bugfix<<std::endl;
    // }

	/*****************************************************************/

	// if(myrank==master_rank)std::cout<<"Npo "<<Npo<<" Ndim "<<Ndim<<std::endl;

	int myseg=0;     // local number of unknowns
	int* perms = new int[Npo]; //permutation vector returned by HODLR
	int* groups = new int[size];
	int i_opt;
	double d_opt;
	int cpp=1; //1: use user-defined cpp/c functions for construction
	int elem_extract=2;

#ifdef HAVE_MPI
	MPI_Fint Fcomm;  // the fortran MPI communicator
	Fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);
#else
    int Fcomm=321;  // dummy fortran MPI communicator
#endif


	for (int i = 0; i < size; i++)groups[i]=i;
	// create butterflypack data structures
	c_bpack_createptree(&size, groups, &Fcomm, &(bpack_obj->ptree));
	c_bpack_createoption(&(bpack_obj->option));
	c_bpack_createstats(&(bpack_obj->stats));


	// set default butterflypack options
	c_bpack_set_D_option(&(bpack_obj->option), "tol_comp", tol);
	c_bpack_set_D_option(&(bpack_obj->option), "tol_rand", tol);
	c_bpack_set_D_option(&(bpack_obj->option), "tol_Rdetect", tol*1e-1);
	c_bpack_set_D_option(&(bpack_obj->option), "sample_para", sample_para);
	c_bpack_set_D_option(&(bpack_obj->option), "near_para", near_para);
	c_bpack_set_I_option(&(bpack_obj->option), "nogeo", nogeo);
	c_bpack_set_I_option(&(bpack_obj->option), "Nmin_leaf", Nmin);
	c_bpack_set_I_option(&(bpack_obj->option), "RecLR_leaf", com_opt);
	c_bpack_set_I_option(&(bpack_obj->option), "xyzsort", sort_opt);
	c_bpack_set_I_option(&(bpack_obj->option), "ErrFillFull", checkerr);
	c_bpack_set_I_option(&(bpack_obj->option), "BACA_Batch", batch);
	c_bpack_set_I_option(&(bpack_obj->option), "LR_BLK_NUM", bnum);
	c_bpack_set_I_option(&(bpack_obj->option), "cpp", cpp);
	c_bpack_set_I_option(&(bpack_obj->option), "LRlevel", lrlevel);
	c_bpack_set_I_option(&(bpack_obj->option), "knn", knn);
	c_bpack_set_I_option(&(bpack_obj->option), "verbosity", 1);
	c_bpack_set_I_option(&(bpack_obj->option), "less_adapt", 1);
	c_bpack_set_I_option(&(bpack_obj->option), "itermax", 10);
	c_bpack_set_I_option(&(bpack_obj->option), "ErrSol", 0);
	c_bpack_set_I_option(&(bpack_obj->option), "format", format);
	c_bpack_set_I_option(&(bpack_obj->option), "elem_extract", elem_extract);

	// set command-line butterflypack options
	c_bpack_set_option_from_command_line(c_argc, argv_ptr,bpack_obj->option);

	// print out butterflypack options
	c_bpack_printoption(&(bpack_obj->option),&(bpack_obj->ptree));

    // construct matrix with geometrical points
	c_bpack_construct_init(&Npo, &Ndim, dat_ptr, nns_ptr,&nlevel, tree, perms, &myseg, &(bpack_obj->bmat), &(bpack_obj->option), &(bpack_obj->stats), &(bpack_obj->msh), &(bpack_obj->ker), &(bpack_obj->ptree), &c_bpack_FuncDistmn, &c_bpack_FuncNearFar, quant_ptr);
	c_bpack_construct_element_compute(&(bpack_obj->bmat), &(bpack_obj->option), &(bpack_obj->stats), &(bpack_obj->msh), &(bpack_obj->ker), &(bpack_obj->ptree), &c_bpack_FuncZmn, &c_bpack_FuncZmnBlock, quant_ptr);

    double vtmp;
    c_bpack_getstats(&(bpack_obj->stats),"Rank_max",&vtmp);
    *rankmax=int(vtmp);


	delete[] perms;
	delete[] tree;
	delete[] groups;



    *pyobj = (void*)bpack_obj;

}




void py_bpack_factor(void ** pyobj)
{
	bpack_handle* bpack_obj = (bpack_handle*)(*pyobj);
	c_bpack_factor(&(bpack_obj->bmat),&(bpack_obj->option),&(bpack_obj->stats),&(bpack_obj->ptree),&(bpack_obj->msh));
	*pyobj = (void*)bpack_obj;
}



void py_bpack_solve(void ** pyobj, int nrhs, C_DT   *xb_global)
{
	bpack_handle* bpack_obj = (bpack_handle*)(*pyobj);
	int idxs, m_loc, fst_row, m, m_loc1, fst_row1;

	c_bpack_localindices(&(bpack_obj->msh), &idxs, &m_loc, &m);

    /* Get the local B and X*/
    C_DT* b= (C_DT*)malloc(m_loc * nrhs*sizeof(C_DT));
    C_DT* x= (C_DT*)malloc(m_loc * nrhs*sizeof(C_DT));

	c_bpack_vector_global2local(&(bpack_obj->ptree), &(bpack_obj->msh), &nrhs, xb_global, b);
	c_bpack_solve(x,b,&m_loc,&nrhs,&(bpack_obj->bmat),&(bpack_obj->option),&(bpack_obj->stats),&(bpack_obj->ptree));
	c_bpack_vector_local2global(&(bpack_obj->ptree), &(bpack_obj->msh), &nrhs, x, xb_global);

	free(b);
	free(x);
	*pyobj = (void*)bpack_obj;
}



void py_bpack_mult(void ** pyobj, int nrhs, char const *trans, C_DT   *xy_global)
{
	bpack_handle* bpack_obj = (bpack_handle*)(*pyobj);
	int idxs, m_loc, fst_row, m, m_loc1, fst_row1;

	c_bpack_localindices(&(bpack_obj->msh), &idxs, &m_loc, &m);

    /* Get the local X and Y*/
    C_DT* y= (C_DT*)malloc(m_loc * nrhs*sizeof(C_DT));
    C_DT* x= (C_DT*)malloc(m_loc * nrhs*sizeof(C_DT));

	c_bpack_vector_global2local(&(bpack_obj->ptree), &(bpack_obj->msh), &nrhs, xy_global, x);
	c_bpack_mult(trans, x, y, &m_loc, &m_loc, &nrhs, &(bpack_obj->bmat),&(bpack_obj->option),&(bpack_obj->stats),&(bpack_obj->ptree));
	c_bpack_vector_local2global(&(bpack_obj->ptree), &(bpack_obj->msh), &nrhs, y, xy_global);

	free(y);
	free(x);
	*pyobj = (void*)bpack_obj;
}





void py_bpack_free(void ** pyobj)
{
	py::gil_scoped_acquire gil;
	bpack_handle* bpack_obj = (bpack_handle*)(*pyobj);

	c_bpack_deletestats(&(bpack_obj->stats));
	c_bpack_deleteproctree(&(bpack_obj->ptree));
	c_bpack_deletemesh(&(bpack_obj->msh));
	c_bpack_deletekernelquant(&(bpack_obj->ker));
	c_bpack_delete(&(bpack_obj->bmat));
	c_bpack_deleteoption(&(bpack_obj->option));
	// free(bpack_obj->quant_ptr);
	delete bpack_obj->quant_ptr;
	free(bpack_obj);
	bpack_obj = NULL;
	*pyobj = (void*)bpack_obj;
}



void py_bpack_terminate()
{
#ifdef HAVE_MPI
	Cblacs_exit(1);
 	MPI_Finalize();                                 // Terminate MPI. Once called, no other MPI routines may be called
#endif
}



#ifdef __cplusplus
}
#endif