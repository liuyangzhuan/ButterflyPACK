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
		{"use_zfp",         "whether to use zfp compression"},
		{"use_qtt",         "whether to use qtt compression"},
		{"hextralevel",         "HMAT: extra levels for top partitioning of the H matrix based on MPI counts. BLR: Maxlevel-hextralevel is the level for defining B-LR/B-BF blocks"},
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