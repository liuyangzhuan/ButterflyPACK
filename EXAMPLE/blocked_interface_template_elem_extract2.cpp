
// ************* The extraction sampling function that computes a list of matrix blocks *************
// The performance benefits compared to single entry evaluation-based interface are:
// 1. Possible OpenMP loops over the number of blocks
// 2. Possible OpenMP loops over the entries of each block
// 3. The implementation of loop-by-element instead of loop-by-basis for high-order methods, etc.
// To use this interface, one needs to call bpack_set_option(option, "elem_extract", 2);
    // allrows: 1D array containing the global row indices (in original order starting from 1 to N) stacked together
    // allcols: 1D array containing the global column indices (in original order starting from 1 to N) stacked together
    // alldat_loc: 1D array containing the local entry values defined by pmaps (in column major) stacked together
    // colidx: 1D array containing sizes of columns of each block
    // rowidx: 1D array containing sizes of rows of each block
    // pgidx:  1D array containing the process group number of each block, the number starts from 0
    // Npmap:  number of process groups
    // pmaps:  2D array (Npmapx3) containing number of process rows, number of process columns, and starting process ID in each block
inline void C_FuncZmnBlock_BF(int* Ninter, int* Nallrows, int* Nallcols, int64_t* Nalldat_loc, int* allrows, int* allcols, _Complex double* alldat_loc, int* rowidx,int* colidx, int* pgidx, int* Npmap, int* pmaps, C2Fptr quant) {
    C_QuantApp_BF* Q = (C_QuantApp_BF*) quant;
    int myrank, size;                     // Store values of processor rank and total no of procs requestedss
    MPI_Comm_size(MPI_COMM_WORLD, &size); 	                // Get no of procs
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 	                // Get no of procs
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
        if(myrank==pid){
          idx_row_map[NinterNew]=idx_row;
          idx_col_map[NinterNew]=idx_col;
          idx_val_map[NinterNew]=idx_val;
          idx_val+=nr*nc;
          inter_map[NinterNew]=nn;
          NinterNew++;
        }else{
        }
        idx_row+=nr;
        idx_col+=nc;
        nrmax = max(nr,nrmax);
        ncmax = max(nc,ncmax);
        nvalmax = max(nc*nr,nvalmax);
      }else{
        std::cout<<"nprow*npcol>1 in C_FuncZmnBlock_BF"<<std::endl;
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
        rows[idxr]=allrows[idx_row_map[nn1]+idxr];

      }
      for (int idxc=0;idxc<nc;idxc++){
        cols[idxc]=allcols[idx_col_map[nn1]+idxc];

      }

      // The following code computes a single matrix block (intersection) of sizes nr x nc, with row
      // indices in rows (index is 1-based instead of 0-based) and column indices in cols (index is
      // 1-based instead of 0-based), and entry values (nr x nc entries stored in column major) in alldat_loc starting from idx_val_map[nn1]
      your_function_to_compute_one_block(nr, nc, rows, cols, &(alldat_loc[idx_val_map[nn1]]), Q);

    }
    free(rows);
    free(cols);
    }
}



