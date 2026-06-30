#include "BPACK_wrapper.h"
#include "butterfly_integration.hpp"


#include <mpi.h>
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

using namespace std;


static inline int product(int arr[], int n) {
  int out = 1;  // Initialize product to 1
  for (int i = 0; i < n; i++) {
      out *= arr[i];
  }
  return out;
}

// The command line parser for the example related parameters
void c_bpack_set_option_from_command_line(int argc, const char* const* cargv,F2Cptr option0) {

	struct OptionHelp {
		const char* name;
		const char* description;
	};

	static const OptionHelp option_help_table[] = {
		{"nmin_leaf",       "leafsize in the hierarchical partitioning"},
		{"tol_comp",        "relative tolerance for matrix construction"},
		{"tol_rand",        "relative tolerance for matrix inversion"},
		{"tol_Rdetect",     "relative tolerance for rank detection during matrix inversion"},
		{"tol_itersol",     "convergence tolerance for TFQMR iterative solver if precon=2 or 3"},
		{"n_iter",          "maximum iteration count for TFQMR"},
		{"level_check",     "the level in the hierarchical partitioning where the randomized construction algorithm is tested, set to 10000 by default (no checking)"},
		{"precon",          "the use mode of butterflypack: 1: as a direct solver 2: as an iterative solver (compress the matrix and pass it to TFQMR without preconditioner), 3: as a preconditioned iterative solver (compress the matrix and invert the matrix and pass them to TFQMR, using approximate matrix inverse as a preconditioner)"},
		{"xyzsort",         "the hierarchical partitioning algorithm: 0: no permutation 1: permutation based on KD-tree 2: permutation based on cobble-like partitioning"},
		{"lrlevel",         "the level in the hierarchical partitioning (top-down numbered) above which butterfly is used and below which low-rank is used"},
		{"errfillfull",     "errfillfull: a slow (n^2), thorough error checking is performed after the compression of each block"},
		{"baca_batch",      "block size in batched ACA when reclr_leaf=4 or 5"},
		{"reclr_leaf",      "low-rank compression algorithms 1:SVD 2:RRQR 3:ACA 4:BACA 5:BACA_improved 6:Pseudo-skeleton 7: ACA with naive parallelization"},
		{"nogeo",           "whether there is geometry information provided 1: is no geometry (xyzsort can not be 1 or 2), 0: there is geometry"},
		{"less_adapt",      "1: improved randomized butterfly construction, default to 1"},
		{"errsol",          "1: generating an artificial true solution vector, compute the RHS with compressed matrix, solve the system, and compare with true solution vector"},
		{"lr_blk_num",      "sqrt of #of subblocks in H-BACA, default to 1"},
		{"rank0",           "initial rank guess in the randomized butterfly algorithm, default to 32"},
		{"rankrate",        "increasing ratio of the rank guess in each iteration, default to 2"},
		{"itermax",         "maximum number of iterations in the randomized butterfly algorithm, default to 10"},
		{"powiter",         "order of power iteration in the randomized low-rank construction"},
		{"ilu",             "whether symmetric gauss-seidel is used when format=2"},
		{"nbundle",         "multiply nbundle sets of vectors together in randomized butterfly algorithm for better flop performance, default to 1"},
		{"near_para",       "admissibility parameter when format=2/3/4/5, strong admissibility typically requires near_para>2.0"},
		{"format",          "the hierarchical matrix format: 1: HODLR/HODBF 2: H matrix 3: HSSBF/SHNBF 4: HSSBF_MD/SHNBF_MD 5: block-LR/BF"},
		{"verbosity",       "verbosity for the printing (-1, 0, 1, 2), -1 suppresses everything, 2 prints most details"},
		{"rmax",            "preestimate of the maximum rank for allocating buffers, default to 1000"},
		{"sample_para",     "oversampling factor in the nlogn entry-evaluation-based butterfly algorithm, default to 2"},
		{"pat_comp",        "pattern of entry-evaluation-based butterfly compression: 1 from right to left, 2 from left to right, 3 from outer to inner"},
		{"knn",             "nearest neighbouring points used in improved BACA and entry-evaluation-based butterfly compression"},
		{"knn_near_para",   "admissibility parameter for guiding the nearest neighbouring points search"},
		{"forwardN15flag",  "whether to use N15 or NlogN algorithm for entry-evaluation-based matrix butterfly compression"},
		{"sample_para_outer","oversampling factor for the outtermost factor matrices in the nlogn entry-evaluation-based butterfly algorithm, default to 2"},
		{"elem_extract",    "0: evaluating entries one by one 1: evaluating entries block by block (may requires communication inside the callback function) 2: evaluating entries block by block (no communication allowed inside the callback function)"},
		{"fastsample_tensor","0: uniformly sample each dimension. 1: uniformly sample the rows of the unfolded matrices on top of 0. 2: use translation invariance"},
		{"trans_invariant", "1: reuse HTENSOR blocks for translational-invariant tensor kernels on uniform Cartesian grids"},
		{"htensor_mvp_level_batch", "number of HTENSOR levels grouped per MVP call; 1 keeps level-by-level memory"},
		{"use_zfp",         "whether to use zfp compression"},
		{"use_qtt",         "whether to use qtt compression"},
		{"hextralevel",         "HMAT: extra levels for top partitioning of the H matrix based on MPI counts. BLR: Maxlevel-hextralevel is the level for defining B-LR/B-BF blocks"},
		{"iter_solver",         "The choice of iterative solvers. 1: TFQMR, 2: GMRES or 3: IterativeRefinement)"},
		{"help",            "print this help message"}
	};

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
		{
		{"help",             no_argument,       0, 1000},
		{"nmin_leaf",                     required_argument, 0, 1},
		{"tol_comp",                   required_argument, 0, 2},
		{"tol_rand",                   required_argument, 0, 3},
		{"tol_Rdetect",             required_argument, 0, 4},
		{"tol_itersol",             required_argument, 0, 5},
		{"n_iter",          required_argument, 0, 6},
		{"level_check",         required_argument, 0, 7},
		{"precon",                  required_argument, 0, 8},
		{"xyzsort",                  required_argument, 0, 9},
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
		{"knn_near_para",         required_argument, 0, 31},
		{"forwardN15flag",         required_argument, 0, 32},
		{"sample_para_outer",         required_argument, 0, 33},
		{"elem_extract",         required_argument, 0, 34},
		{"fastsample_tensor",         required_argument, 0, 35},
		{"use_zfp",         required_argument, 0, 36},
		{"use_qtt",         required_argument, 0, 37},
		{"hextralevel",         required_argument, 0, 38},
		{"iter_solver",         required_argument, 0, 39},
		{"trans_invariant",         required_argument, 0, 40},
			{"htensor_mvp_level_batch", required_argument, 0, 41},
		{NULL, 0, NULL, 0}
		};
	int c, option_index = 0;
	// bool unrecognized_options = false;
	opterr = optind = 0;
	while ((c = getopt_long_only
			(argc, argv.data(), "",
			long_options, &option_index)) != -1) {
		switch (c) {
		case 1000: {
		std::cout << "Available ButterflyPACK Command-Line Options:\n";
		for (const auto& opt : option_help_table) {
			std::cout << "  --" << std::setw(20) << std::left << opt.name
						<< " : " << opt.description << "\n";
		}
		} break;
		case 1: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "Nmin_leaf", opt_i);
		} break;
		case 2: {
		std::istringstream iss(optarg);
		iss >> opt_d;
		c_bpack_set_D_option(&option0, "tol_comp", opt_d);
		c_bpack_set_D_option(&option0, "tol_rand", opt_d);
		c_bpack_set_D_option(&option0, "tol_Rdetect", opt_d*0.1);
		} break;
		case 3: {
		std::istringstream iss(optarg);
		iss >> opt_d;
		c_bpack_set_D_option(&option0, "tol_rand", opt_d);
		} break;
		case 4: {
		std::istringstream iss(optarg);
		iss >> opt_d;
		c_bpack_set_D_option(&option0, "tol_Rdetect", opt_d);
		} break;
		case 5: {
		std::istringstream iss(optarg);
		iss >> opt_d;
		c_bpack_set_D_option(&option0, "tol_itersol", opt_d);
		} break;
		case 6: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "n_iter", opt_i);
		} break;
		case 7: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "level_check", opt_i);
		} break;
		case 8: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "precon", opt_i);
		} break;
		case 9: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "xyzsort", opt_i);
		} break;
		case 10: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "LRlevel", opt_i);
		} break;
		case 11: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "ErrFillFull", opt_i);
		} break;
		case 12: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "BACA_Batch", opt_i);
		} break;
		case 13: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "RecLR_leaf", opt_i);
		} break;
		case 14: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "nogeo", opt_i);
		} break;
		case 15: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "less_adapt", opt_i);
		} break;
		case 16: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "ErrSol", opt_i);
		} break;
		case 17: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "LR_BLK_NUM", opt_i);
		} break;
		case 18: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "rank0", opt_i);
		} break;
		case 19: {
		std::istringstream iss(optarg);
		iss >> opt_d;
		c_bpack_set_D_option(&option0, "rankrate", opt_d);
		} break;
		case 20: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "itermax", opt_i);
		} break;
		case 21: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "powiter", opt_i);
		} break;
		case 22: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "ILU", opt_i);
		} break;
		case 23: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "Nbundle", opt_i);
		} break;
		case 24: {
		std::istringstream iss(optarg);
		iss >> opt_d;
		c_bpack_set_D_option(&option0, "near_para", opt_d);
		} break;
		case 25: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "format", opt_i);
		} break;
		case 26: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "verbosity", opt_i);
		} break;
		case 27: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "rmax", opt_i);
		} break;
		case 28: {
		std::istringstream iss(optarg);
		iss >> opt_d;
		c_bpack_set_D_option(&option0, "sample_para", opt_d);
		} break;
		case 29: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "pat_comp", opt_i);
		} break;
		case 30: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "knn", opt_i);
		} break;
		case 31: {
		std::istringstream iss(optarg);
		iss >> opt_d;
		c_bpack_set_D_option(&option0, "knn_near_para", opt_d);
		} break;
		case 32: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "forwardN15flag", opt_i);
		} break;
		case 33: {
		std::istringstream iss(optarg);
		iss >> opt_d;
		c_bpack_set_D_option(&option0, "sample_para_outer", opt_d);
		} break;
		case 34: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "elem_extract", opt_i);
		} break;
		case 35: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "fastsample_tensor", opt_i);
		} break;
		case 40: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "trans_invariant", opt_i);
		} break;
			case 41: {
			std::istringstream iss(optarg);
			iss >> opt_i;
			c_bpack_set_I_option(&option0, "htensor_mvp_level_batch", opt_i);
			} break;
		case 36: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "use_zfp", opt_i);
		} break;
		case 37: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "use_qtt", opt_i);
		} break;
		case 38: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "hextralevel", opt_i);
		} break;
		case 39: {
		std::istringstream iss(optarg);
		iss >> opt_i;
		c_bpack_set_I_option(&option0, "iter_solver", opt_i);
		} break;
		default: break;
		}
	}
	}





void c_bpack_vector_global2local(F2Cptr* ptree, F2Cptr* msh, int* nvec_p, C_DT   *b_global, C_DT   *b)
{
	int idxs, m_loc, fst_row, m, m_loc1, fst_row1;
	int nvec = *nvec_p;

    int iam, size;
    int master_rank = 0;


#ifdef HAVE_MPI

	MPI_Comm comm;
	int fcomm;
	c_bpack_get_comm(ptree, &fcomm);
	comm = MPI_Comm_f2c((MPI_Fint)fcomm);

	MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &iam);
	MPI_Bcast( &nvec, 1, MPI_INT, 0, comm);
	c_bpack_localindices(msh, &idxs, &m_loc, &m);
	fst_row = idxs-1;

	{
	C_DT* b_global_permed=NULL;
	C_DT* btmp= (C_DT*)malloc(m_loc * nvec*sizeof(C_DT));
	int *counts = NULL, *displs = NULL;

	if (iam == 0) {
		counts = (int *) malloc(size*sizeof(int));
		displs = (int *) malloc(size*sizeof(int));
		b_global_permed= (C_DT*)malloc(m * nvec*sizeof(C_DT));
		for (int j = 0; j < nvec; ++j) {
			for (int i = 0; i < m; ++i) {
				int i_old= i+1;
				int i_new;
				c_bpack_old2new(msh,&i_old,&i_new);
				b_global_permed[j + (i_new-1) * nvec] = b_global[j * m + i];
			}
		}
	}
	m_loc1 = m_loc*nvec;
	fst_row1 = fst_row*nvec;
	MPI_Gather(&m_loc1, 1, MPI_INT, counts, 1, MPI_INT, 0, comm);
	MPI_Gather(&fst_row1, 1, MPI_INT, displs, 1, MPI_INT, 0, comm);
	MPI_Scatterv(b_global_permed, counts, displs, C_MPI_DT,  btmp, m_loc * nvec, C_MPI_DT, 0, comm);

	for (int j = 0; j < nvec; ++j)
	{
		for (int i = 0; i < m_loc; ++i)
		{
			b[j * m_loc + i] = btmp[j + i * nvec];
		}
	}
	if (iam == 0) {
		free(counts);
		free(displs);
		free(b_global_permed);
	}
	free(btmp);
	}
#else
	c_bpack_localindices(msh, &idxs, &m_loc, &m);
	for (int j = 0; j < nvec; ++j)
	{
		for (int i = 0; i < m; ++i)
		{
			b[j * m + i] = b_global[j * m + i];
		}
	}
#endif

}



void c_bpack_vector_local2global(F2Cptr* ptree, F2Cptr* msh, int* nvec_p, C_DT   *b, C_DT   *b_global)
{
	int idxs, m_loc, fst_row, m, m_loc1, fst_row1;
	int nvec = *nvec_p;
    int iam, size;
    int master_rank = 0;


#ifdef HAVE_MPI
	MPI_Comm comm;
	int fcomm;
	c_bpack_get_comm(ptree, &fcomm);
	comm = MPI_Comm_f2c((MPI_Fint)fcomm);

	MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &iam);
	MPI_Bcast( &nvec, 1, MPI_INT, 0, comm);
	c_bpack_localindices(msh, &idxs, &m_loc, &m);
	fst_row = idxs-1;

	{
	C_DT* b_global_permed;
	C_DT* btmp= (C_DT*)malloc(m_loc * nvec*sizeof(C_DT));
	for (int j = 0; j < nvec; ++j)
	{
		for (int i = 0; i < m_loc; ++i)
		{
			btmp[j + i * nvec] = b[j * m_loc + i];
		}
	}
	int *counts = NULL, *displs = NULL;
	if (iam == 0) {
		counts = (int *) malloc(size*sizeof(int));
		displs = (int *) malloc(size*sizeof(int));
		b_global_permed= (C_DT*)malloc(m * nvec*sizeof(C_DT));
	}
	m_loc1 = m_loc*nvec;
	fst_row1 = fst_row*nvec;
	MPI_Gather(&m_loc1, 1, MPI_INT, counts, 1, MPI_INT, 0, comm);
	MPI_Gather(&fst_row1, 1, MPI_INT, displs, 1, MPI_INT, 0, comm);
	MPI_Gatherv(btmp, m_loc * nvec, C_MPI_DT, b_global_permed, counts, displs, C_MPI_DT, 0, comm);
	if (iam == 0) {
		for (int j = 0; j < nvec; ++j) {
			for (int i = 0; i < m; ++i) {
				int i_old= i+1;
				int i_new;
				c_bpack_old2new(msh,&i_old,&i_new);
				b_global[j * m + i] = b_global_permed[j + (i_new-1) * nvec];
			}
		}
		free(counts);
		free(displs);
		free(b_global_permed);
	}
	free(btmp);
	}
#else

	c_bpack_localindices(msh, &idxs, &m_loc, &m);
	for (int j = 0; j < nvec; ++j)
	{
		for (int i = 0; i < m; ++i)
		{
			b_global[j * m + i] = b[j * m + i];
		}
	}

#endif

}


void c_bpack_md_vector_local2global(F2Cptr* ptree, F2Cptr* msh, int* ndim_p, int* nvec_p, C_DT   *b, C_DT   *b_global)
{

	int nvec = *nvec_p;
	int Ndim= *ndim_p;
    int iam, size;
    int master_rank = 0;


	int *idxs = (int *) malloc(Ndim*sizeof(int));
	int *m_loc = (int *) malloc(Ndim*sizeof(int));
	int *m = (int *) malloc(Ndim*sizeof(int));

#ifdef HAVE_MPI
	MPI_Comm comm;
	int fcomm;
	c_bpack_get_comm(ptree, &fcomm);
	comm = MPI_Comm_f2c((MPI_Fint)fcomm);

	MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &iam);
	MPI_Bcast( &nvec, 1, MPI_INT, 0, comm);
	c_bpack_md_localindices(&Ndim, msh, idxs, m_loc, m);
	size_t m_scalar = product(m,Ndim);
	int m_loc1 = product(m_loc,Ndim)*nvec;

	int *counts = NULL, *displs = NULL, *m_locs=NULL, *idxss=NULL;
	C_DT *b_global_tmp = NULL;
	if (iam == 0) {
		counts = (int *) malloc(size*sizeof(int));
		displs = (int *) malloc(size*sizeof(int));
		m_locs = (int *) malloc(Ndim*size*sizeof(int));
		idxss = (int *) malloc(Ndim*size*sizeof(int));
		b_global_tmp= (C_DT*)malloc(m_scalar * nvec*sizeof(C_DT));
	}

	MPI_Gather(m_loc, Ndim, MPI_INT, m_locs, Ndim, MPI_INT, 0, comm);
	MPI_Gather(idxs, Ndim, MPI_INT, idxss, Ndim, MPI_INT, 0, comm);
	MPI_Gather(&m_loc1, 1, MPI_INT, counts, 1, MPI_INT, 0, comm);
	if (iam == 0) {
		displs[0] = 0;
		for (int i = 1; i < size; i++) {
			displs[i] = displs[i-1] + counts[i-1];
		}
	}
	MPI_Gatherv(b, m_loc1, C_MPI_DT, b_global, counts, displs, C_MPI_DT, 0, comm);

	if (iam == 0) {
		//Converting from the "stacked" order (rank 0's data, then rank 1's, etc.) to the permuted global ordering
		for (int s=0;s<size;s++){
			int m_loc_1vec = product(&(m_locs[s*Ndim]),Ndim);
			for (int i=0; i<m_loc_1vec; i++){
				int i_new_scalar = i+1;
				int i_new_md[Ndim];
				c_bpack_singleindex_to_multiindex(&Ndim,&(m_locs[s*Ndim]),&i_new_scalar,i_new_md);
				for (int j=0;j<Ndim;j++){
					i_new_md[j]=i_new_md[j]+idxss[s*Ndim+j]-1;
				}
				c_bpack_multiindex_to_singleindex(&Ndim,m,&i_new_scalar,i_new_md);
				for (int nth=0; nth<nvec; nth++){
					b_global_tmp[(i_new_scalar-1)+m_scalar*nth] = b_global[displs[s]+i+m_loc_1vec*nth];
				}
			}
		}

		//Converting from permuted global ordering to the original ordering
		for (int i = 0; i < m_scalar; ++i) {
			int i_old_scalar=i+1;
			int i_new_scalar;
			int i_old_md[Ndim];
			int i_new_md[Ndim];
			c_bpack_singleindex_to_multiindex(&Ndim,m,&i_old_scalar,i_old_md);
			c_bpack_md_old2new(&Ndim, msh, i_old_md, i_new_md);
			c_bpack_multiindex_to_singleindex(&Ndim,m,&i_new_scalar,i_new_md);
			for (int nth=0; nth<nvec; nth++){
				b_global[nth * m_scalar + i] = b_global_tmp[nth * m_scalar + (i_new_scalar-1)];
			}
		}

		free(counts);
		free(displs);
		free(m_locs);
		free(idxss);
		free(b_global_tmp);
	}
#else

	c_bpack_md_localindices(&Ndim, msh, idxs, m_loc, m);
	size_t m_scalar = product(m,Ndim);
	for (int j = 0; j < nvec; ++j)
	{
		for (int i = 0; i < m_scalar; ++i)
		{
			b_global[j * m_scalar + i] = b[j * m_scalar + i];
		}
	}
#endif

free(idxs);
free(m_loc);
free(m);

}

void c_bpack_construct_init(int* Npo, int* Ndim, double* Locations, int* nns, int* nlevel, int* tree, int* perms, 
	int* Npo_loc, F2Cptr* bmat, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,F2Cptr* ker,F2Cptr* ptree, 
	void (*C_FuncDistmn)(int*, int*, double*,C2Fptr), 
	void (*C_FuncNearFar)(int*, int*, int*,C2Fptr), C2Fptr C_QuantApp){
	// C_FuncDistmn: defines distance
	// C_FuncNearFar:

  // correspond to create_uniform_tree
  // arguments in create_uniform_tree:
  //   point_coords: null_ptr (assign uniform grid) or a pointer to an array of points on the grid
  //   num_points: N, the dimension of the matrix
  //   num_levels: number of levels for factorization
  //   global_bounds: bounding min and max of all dimensions
  //   dimension: 2D, or 3D problem
  //   comm: some form of MPI info, gets mpi_rank, mpi_size
  //   reduction_threshold: only uniform reduction pattern is supported, so really there is only one option which is uniform
  //   pattern: only uniform reduction pattern is supported, so really there is only one option which is uniform
  //   returns: tree structure 

  //   can also create HierarchicalFactorization object
  //   * @param N Total number of points in the problem
  //   * @param prop Matrix property (symmetric, hermitian, or nonsymmetric)
  //   * @param kernel_func Kernel evaluator
  //   * @param dim Spatial dimension (2 or 3)
  //   * @param factorization_type Method for factorizing/inverting matrices (default: CHOLESKY)
  //   * @param num_proxy Number of proxy points (-1 uses default 32 for 2D, 256 for 3D)
  //   * @param tol Compression tolerance (default: 1e-6)
  //   * @param proxy_factor Proxy surface radius factor (default: 2.5)

  // ProgramOptions
  // int num_levels = nlevel;
  // int64_t N = Npo;
  // int64_t grid_size = 0;
  // double tolerance = 0.0;
  // fmm::KernelKind kernel_kind = fmm::KernelKind::LAPLACE;
  // NumberKind number_kind = NumberKind::REAL;
  // int dimension = Ndim;
  // int64_t reduction_threshold = 0;
  // int num_proxy = -1;
  // double wave_divisor = 32.0;
  // double length_scale = 0.1;   // Matérn length scale ℓ
  // double nugget = 1e-6;        // Matérn diagonal nugget σ_n²
  // double kappa = 10.0;         // Yukawa screening parameter κ
  // int cond_samples = 0;        // Power iteration samples for condition number estimate (0 = skip)
  
  // Butterflypack end: need some definition of proxy points

  //   Npo: pass into num_points
  //   Ndim: pass into dimension
  //   Locations: pass into point_coords
  //   nns: 
  //   nlevel: pass into num_levels
  //   tree: type difference with tree returned by create_uniform_tree
  //   perms: permutation vector?
  //   Npo_loc: 
  //   bmat: this stores the h2 solver struct, h2 tree, kernel, etc. 
  //   option: 
  //   stats: 
  //   msh 
  //   ker: kernel types from FMM?
  //   ptree: mpi communicator needed, otherwise not relevant
  //   C_FuncDistmn
  //   C_FuncNearFar
  //   C_QuantApp

	
  double tmp;
  c_bpack_getoption(option, "format", &tmp);
  int format=(int)tmp;
  if(format==7){
	// use datatype C_DT
    // construct H2 solver
    // To Do: need to add CoordType and DataType
	//auto H2_solver = std::make_unique<H2<double, C_DT>>();
    butterfly::H2<double, C_DT>*  H2_solver = new H2<double, C_DT>();

    int fcomm;
    c_bpack_get_comm(ptree, &fcomm);
    MPI_Comm mpi_comm = MPI_Comm_f2c((MPI_Fint)fcomm);
    H2_solver->comm = mpi_comm;
    //*bmat = static_cast<F2Cptr>(H2_solver.release());
	*bmat = static_cast<F2Cptr>(H2_solver);
    
	int rank = 0;
    int size = 1;
	MPI_Comm_rank(H2_solver->comm, &rank);
    MPI_Comm_size(H2_solver->comm, &size);

    // // parse options and throw error if any requirements are not satisfied
    // for (int i = 1; i < argc; ++i) {
    //   const std::string arg = argv[i];
    //   if (arg == "--help" || arg == "-h") {
    //     if (rank == 0) {
    //       print_usage(argv[0]);
    //     }
    //     MPI_Finalize();
    //     return 0;
    //   }
    // }
    
    butterfly::ProgramOptions h2_options;
    try {
      h2_options = butterfly::parse_program_options(Npo, Ndim, Locations, option, stats, ker);
    } catch (const std::exception& e) {
        if (rank == 0) {
            std::cerr << "Argument error: " << e.what() << std::endl;
        }
		throw;
    }

    try {
      if (rank == 0) {

        std::cout << "Run configuration: kernel=" << kernel_kind_to_string(h2_options.kernel_kind)
                  << ", number_type=" << number_kind_to_string(h2_options.number_kind)
                  << ", dimension=" << h2_options.dimension
                  << ", reduction_threshold=" << h2_options.reduction_threshold
                  << ", num_proxy=" << h2_options.num_proxy;
        if (h2_options.kernel_kind == fmm::KernelKind::MATERN52) {
          std::cout << ", length_scale=" << h2_options.length_scale
                    << ", nugget=" << h2_options.nugget;
        }
        if (h2_options.kernel_kind == fmm::KernelKind::YUKAWA) {
          std::cout << ", kappa=" << h2_options.kappa;
        }
        std::cout << std::endl;
      }
      
      butterfly::h2_initiate<double, C_DT>(bmat, h2_options, Locations, rank);
    } catch (const std::exception& e) {
        std::cerr << "Error on rank " << rank << ": " << e.what() << std::endl;
        MPI_Abort(H2_solver->comm, 1);
    }

	//temporary
	delete H2_solver;
	H2_solver = nullptr;
  }else{
	  c_bpack_construct_init_fortran(Npo, Ndim, Locations, nns, nlevel, tree, perms, Npo_loc, bmat, option, stats, msh, ker, ptree, C_FuncDistmn, C_FuncNearFar, C_QuantApp);
  }
}


void c_bpack_construct_element_compute(F2Cptr* bmat, F2Cptr* option,F2Cptr* stats,F2Cptr* msh,
	F2Cptr* ker,F2Cptr* ptree, void (*C_FuncZmn)(int*, int*, C_DT*,C2Fptr),
	void (*C_FuncZmnBlock)(int*, int*, int*, int64_t*, int*, int*, C_DT*, int*, int*, int*, int*, int*, C2Fptr), 
	C2Fptr C_QuantApp){  
  // these functions are important to define
  // C_FuncZmn: returns value at (i,j)th element of matrix -- need to update this to work for kernel
  // C_FuncZmnBlock: returns a block, low priority

	
  double tmp;
  c_bpack_getoption(option, "format", &tmp);
  int format=(int)tmp;
  if(format==7){
	H2<double, C_DT>* H2_solver = static_cast<H2<double, C_DT>*>(*bmat);
	H2_solver->kernel = C_FuncZmn;

	// save C_FuncZmnBlock in bmat, low priority right now
  }else{
	c_bpack_construct_element_compute_fortran(bmat, option, stats, msh, ker, ptree, C_FuncZmn, C_FuncZmnBlock, C_QuantApp);
  }
}


void c_bpack_factor(F2Cptr*bmat, F2Cptr*option, F2Cptr*stats, F2Cptr*ptree, F2Cptr*msh){
  
  // Correspond to hierarchical_factorization_parallel, arguments
  // tree: 
  // kernel: kernel from factorizer
  // tolerance: 
  // is_symmetric: bool -- works for general helmholtz and V3D
  // is_hermitian: bool, not supported right now 
  // factorization_method: provided by factorizer, factorization_type
  // unit_proxy_points: 
  // num_proxy: 
  // proxy_radius: 

  // bmat: can contain tree, and kernel function
  // proxy points may have some difficulties, we can discuss more next week
  double tmp;
  c_bpack_getoption(option, "format", &tmp);
  int format=(int)tmp;
  if(format==7){
    try {
      // To Do: need to convert bmat format to the format for hierarchical_factorization_parallel
      H2<double, C_DT>* H2_solver = static_cast<H2<double, C_DT>*>(*bmat);
      // preset because solver can only solve these right now
      is_symmetric = true;
      is_Hermitian = false;
      auto total_start = std::chrono::high_resolution_clock::now();

      fmm::hierarchical_factorization_parallel(
        tree.get(),
        &kernel,
        options.tolerance,
        is_symmetric,
        is_hermitian,
        factorization_method,
        unit_proxy,
        num_proxy,
        2.5,
        true);

      auto total_end = std::chrono::high_resolution_clock::now();
      auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        total_end - total_start);

      if (rank == 0) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Total factorization time: " << total_duration.count() << " ms" << std::endl;
        std::cout << "========================================\n" << std::endl;
      }

    } catch (const std::exception& e) {
        std::cerr << "Error on rank " << rank << ": " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }else{
	c_bpack_factor_fortran(bmat, option, stats, ptree, msh);
  }
}

void c_bpack_solve(C_DT*x, C_DT*b, int*Nloc, int*Nrhs, F2Cptr*bmat, F2Cptr*option, F2Cptr*stats, F2Cptr*ptree){
  // correspond to hierarchical_solve_parallel, arguments:
  //   tree:
  //   rhs: pass in b
  //   solve_data: SolveDataRequest type, accumulates distributes solution / communication requests during the tree level sweep
  //   verbose: printing

  // and gather_solution_to_root, arguments:
  //   tree
  //   solve_data: pass in from hierarchical_solve_parallel
  //   solution: pass into x
  //   aggregated_rhs: 

  // Ax = b
  // x: final solution
  // b: provided rhs
  // Hloc: 
  // Nrhs: number of rhs columns
  // bmat: factored matrix, do we need this? because it's in tree
  // option 
  // stat
  // ptree
	
  // note to self: figure out aggregated_rhs, solution, solve_data; and mpi stuff (mpi stuff probably ask Tianyu)
  // need to redistribute x into H2, and then call mul_parallel, then extract mul_data, the nredistribute to Butterfly
  // discuss about proxy points with Yang next week
  double tmp;
  c_bpack_getoption(option, "format", &tmp);
  int format=(int)tmp;
  if(format==7){
    try {
      //To Do: pass in b to and redistribute to rhs, either do this here or in construct_init
      // pass in tree and rhs to solve_parallel
      std::vector<std::vector<fmm::SolveDataRequest<CoordType, DataType>>> solve_data(
            options.num_levels);
        fmm::hierarchical_solve_parallel(tree.get(), rhs, solve_data, true);

        std::vector<DataType> solution;
        std::vector<DataType> aggregated_rhs;
        const auto gather_verify_start = std::chrono::high_resolution_clock::now();
        fmm::gather_solution_to_root(tree.get(), solve_data, solution, aggregated_rhs);
        const auto gather_verify_end = std::chrono::high_resolution_clock::now();
        const double gather_verify_ms =
            std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
                gather_verify_end - gather_verify_start).count();
        double gather_verify_max_ms = 0.0;
        MPI_Reduce(&gather_verify_ms, &gather_verify_max_ms, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            std::cout << "Gather communication time for solve verification: "
                      << gather_verify_max_ms << " ms" << std::endl;
        }

      // can conduct h2_verification, only if uniform points

      // must call at the end of the program. Idk if solve ends it all or if another function.
      // maybe create another construct_cleanup() function
      // MPI_Finalize();
    } catch (const std::exception& e) {
        std::cerr << "Error on rank " << rank << ": " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

  }else{
	c_bpack_solve_fortran(x, b, Nloc, Nrhs, bmat, option, stats, ptree);
  }
}

void c_bpack_mult(char const * trans, C_DT const * xin, 
	C_DT* xout, int* Ninloc, int* Noutloc, int* Ncol, 
	F2Cptr* bmat,F2Cptr* option,F2Cptr* stats,F2Cptr* ptree){
  
  // F * xin = xout, where F is the approximated matrix from hierarchical decomposition
  // 

  // can call fft_matvec for uniform grid
	
  double tmp;
  c_bpack_getoption(option, "format", &tmp);
  int format=(int)tmp;
  if(format==7){
    H2* h2 = reinterpret_cast<H2*>(*bmat); 
    bool verbose = true;
    // need to pass in rhs from bmat
    // need to pass in h2 tree from bmat

    hierarchical_mul_parallel(tree, rhs, solve_data, verbose); // can only handle matrix vector multiplication right now
	
  }else{
	c_bpack_mult_fortran(trans, xin, xout, Ninloc, Noutloc, Ncol, bmat, option, stats, ptree);
  }
}


void c_bpack_logdet(C_DT* phase, C_RDT* logabsdet, F2Cptr* option, F2Cptr* bmat){
  double tmp;
  c_bpack_getoption(option, "format", &tmp);
  int format=(int)tmp;
  if(format==7){

  }else{
	c_bpack_logdet_fortran(phase, logabsdet, option, bmat);
  }
}