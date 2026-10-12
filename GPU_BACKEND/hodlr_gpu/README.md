# HODLR GPU backend

GPU backend of the HODLR format (`format=1`, `LRlevel=0`) for symmetric and
unsymmetric matrices, in double and double complex. It covers the
construction, the factorization, the multiply and the solve, with one GPU per
MPI rank or several ranks sharing one GPU.

## Enabling it

- Build with `enable_h2_gpu=ON`. MAGMA 2.10 is recommended; see
  `example_scripts/run_cmake_build_gnu_perlmutter_openblas_sequential_h2gpu.sh`.
- Set the option `HODLR_use_gpu`:
  - `0`: CPU (the default).
  - `1`: GPU.
  - `2`: GPU with FP64 tensor-core GEMMs.
- The option `HODLR_gpu_pieces` (default `4`) sets the pieces per rank of a
  low-rank block shared by several ranks (the top levels), rounded down to a
  power of 2. `1` keeps the CPU's pieces (as in `LR_HBACA_Leaflevel`), each
  compressed by its owner rank. With more, each rank's part is cut further
  (while both sides have at least 2048 rows and columns) and the pieces are
  handed out over the block's ranks, largest estimated cost first, then
  merged back. This balances the top levels when some pieces cost far more
  than others (3D Laplace, 885k points, 8 GPUs, tol 1e-3: construction 22 s
  with `1`, 14 s with `4`). `BPACK_CHECK=hodlr-exact` sets it aside (the
  CPU's split, for comparisons).
- The run stops unless `format=1`, `LRlevel=0` and ZFP is off (`use_zfp` not 1).
- The construction of a level runs on the GPU when all of these hold:
  - `RecLR_leaf` is 4 (BACA) or 5 (BACA without overlap);
  - `LR_BLK_NUM=1` and `forwardN15flag=0`;
  - the level is not `level_check`.

  Any other level is built on the CPU and uploaded. The factorization, the
  multiply and the solve always run on the GPU.
- The entries of the construction come from the matrix's GPU evaluator
  (`doc/gpu_kernels.md`). Without one, the BACA of the low-rank levels still
  runs on the GPU and asks the application for the entries of each step's
  panels through its block callback (`C_FuncZmnBlock` / `FuncZmnBlock`), one
  call per step for all the blocks of the level: the columns it picked on the
  host, then the rows its pivoting picked on the device (their indices
  downloaded first). The entries are uploaded into the panels. This needs
  `elem_extract=1` (the callback collective over the matrix's communicator) or
  `elem_extract=2` (each rank on its own); `elem_extract=0` stops with an
  error. The dense leaves take one call of the callback per rank (one round
  with `elem_extract=1`; the CPU path takes one a leaf).
- A compressed matrix as the entry source of another one: the second matrix
  uses `elem_extract=1`, natural order (`xyzsort=0`, `nogeo=1`; its index i
  is the source's tree index i), and a block callback that calls
  `c_bpack_extractelement` on the source. When the source lives on a subset
  of the second matrix's ranks (a sub-communicator, any size),
  `c_bpack_extractelement_bridge` takes the second matrix's process tree as
  well, with the source's handles null on the ranks outside it: each
  requested block is extracted on the source rank of its owner (or, for an
  owner outside the source, on that rank modulo the source's size) and sent
  to its owner. In Fortran, set `ptree_caller` of `quant_bmat` for
  `Zelem_block_Extraction`. Both work with GPU or CPU matrices on either
  side. `EXAMPLE/Laplace3D_H2Test_Driver.cpp --hodlr_extract_ranks <P>`
  builds the source on the first P ranks.
- Each rank uses GPU number (node-local rank) mod (visible GPUs), the
  node-local rank read from the launcher's environment (`SLURM_LOCALID` and
  the like), so that a matrix on a sub-communicator needs no collective over
  `MPI_COMM_WORLD`. Ranks that share a GPU split its memory equally.
- Entry extraction (`BPACK_ExtractElement`, also behind `c_bpack_extractelement`
  and a compressed matrix used as the entry source of another construction)
  from a matrix built on the GPU runs where the factors are
  (`HODLR_gpu_extract_blocks`). The part of a requested block that falls in a
  low-rank block U V^T, rows R and columns C, is computed by the ranks that
  hold the rows of its larger side on their GPUs (the owners of R if R has
  at least as many rows as C has columns, else the owners of C). Only the
  rows of the smaller side move, the products run as one batched GEMM per
  call, and only the entries go to the rank that asked for them. Dense leaves
  are read on the host by their rank. A HODLR built from another one's
  entries this way (3D Laplace, 32^3 points, 8 ranks with the source on 4):
  construction 73 s → 3.7 s. When a requested block is spread over several
  ranks (`pmaps` grids larger than 1 x 1), the rows are gathered from the GPUs
  and the products run on the host, as in the CPU extraction.
- The truncated SVDs of the construction (the recompression after BACA and
  the TSQR of the merges at the shared levels) use Gram matrices when
  `tol_comp` is at least 1e-4: R the Cholesky factor of X^H X + s I (s about
  1e-15 (k + sqrt(m)) tr(X^H X)), the SVD of the small factors, and the new
  factors X R^-1 (...), without forming Q. This is GEMM work where the
  Householder QR's panels are memory bound (3D Laplace, 262k points, 4 GPUs:
  construction 11.2 s → 8.4 s). X R^-1 is a contraction, so the truncation
  error stays within the threshold, also when X has lower rank than columns;
  singular values below about sqrt(s / tr) of the largest are not resolved,
  hence the tolerance floor (tighter tolerances, a failed Cholesky
  factorization, and `BPACK_CHECK=hodlr-exact` use the Householder QR, which
  matches the CPU at round-off).
- The host loops of the GPU construction (the BACA core QRs) use
  `OMP_NUM_THREADS` threads, at most 8: with an OpenBLAS built with
  `USE_LOCKING`, more threads making small LAPACK calls mostly wait for its
  lock.

- The low-rank factors built on the GPU stay there: the device copies from
  the construction become the forward blocks of the factorization (copied on
  the device, never uploaded from the host), and their host copies are left
  unfilled (3D Laplace, 262k points, 16 GPUs: 0.4 s of 7.1 s of downloads
  saved, and host memory per rank 2.5–3.4 → 1.7 GB). Entry extraction works
  from the GPUs (see above). A CPU multiply or CPU
  factorization of such a matrix stops with an error (TODO; with
  `HODLR_use_gpu > 0` both run on the GPU, so this only happens if a caller
  switches `HODLR_use_gpu` off after the construction). `BPACK_CHECK=hodlr`
  and `ErrFillFull=1` download the host copies during the construction.

## Environment variables

[`doc/environment_variables.md`](../doc/environment_variables.md) describes
all of them; those that concern this backend:

- `BPACK_GPU_HEAP_FRACTION` (default `0.85`): the device memory pool, shared
  with the H2 GPU backend and split equally among ranks sharing a GPU. `H2
  GPU heap exhausted` errors refer to it.
- `BPACK_GPU_EXCHANGE_MB` (default min(2560, GPU memory / 16), divided by the
  ranks sharing the GPU): the device arena for MPI messages (the merges of
  the construction and the exchanges of the factorization, multiply and
  solve), with GPU-aware MPI.
- `BPACK_GPU_AWARE_MPI` (default: `MPICH_GPU_SUPPORT_ENABLED`): `0` sends all
  MPI messages through host buffers.
- `BPACK_CHECK=hodlr` runs the CPU routine next to each GPU step and prints
  the difference; `hodlr-transpose` also checks the transposed products;
  `hodlr-exact` builds the shared blocks as the CPU does (the same random
  numbers) and truncates by Householder QR, so the construction levels
  compare at round-off (with `hodlr`, the Gram truncation agrees with the CPU
  to about 1e-7 at tol 1e-4). The timings of such runs are not meaningful.
- `BPACK_TRACE=hodlr-qr` prints each block and chunk of the recompression
  QR, and when a Gram truncation falls back to the QR; `sync` synchronizes after each step, so that an asynchronous device
  fault is reported by the kernel that caused it.
