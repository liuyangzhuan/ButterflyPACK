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
  - the matrix has a GPU evaluator of its entries (`doc/gpu_kernels.md`);
  - `RecLR_leaf` is 4 (BACA) or 5 (BACA without overlap);
  - `LR_BLK_NUM=1` and `forwardN15flag=0`;
  - the level is not `level_check`.

  Any other level is built on the CPU and uploaded. The factorization, the
  multiply and the solve always run on the GPU.
- Each rank uses GPU number (node-local rank) mod (visible GPUs). Ranks that
  share a GPU split its memory equally.
- Entry extraction (`BPACK_ExtractElement`, also behind `c_bpack_extractelement`
  and a compressed matrix used as the entry source of another construction)
  works on the factors on the GPUs: the rows of U and V that the requested
  entries need are gathered on the GPUs that hold them and sent to the ranks
  that compute the entries. The small products and the layout of the result
  are those of the CPU extraction.
- The host loops of the GPU construction (the BACA core QRs) use
  `OMP_NUM_THREADS` threads, at most 8: with an OpenBLAS built with
  `USE_LOCKING`, more threads making small LAPACK calls mostly wait for its
  lock.

- The low-rank factors built on the GPU stay there: the device copies from
  the construction become the forward blocks of the factorization (copied on
  the device, never uploaded from the host), and their host copies are left
  unfilled (3D Laplace, 262k points, 16 GPUs: 0.4 s of 7.1 s of downloads
  saved, and host memory per rank 2.5–3.4 → 1.7 GB). Entry extraction takes
  the rows it needs from the GPUs (see above). A CPU multiply or CPU
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
  numbers), so the construction levels compare at round-off. The timings of
  such runs are not meaningful.
- `BPACK_TRACE=hodlr-qr` prints each block and chunk of the recompression
  QR; `sync` synchronizes after each step, so that an asynchronous device
  fault is reported by the kernel that caused it.
