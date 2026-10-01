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
  with `1`, 14 s with `4`). It needs `HODLR_GPU_SPLIT=1`, the default.
- The run stops unless `format=1`, `LRlevel=0` and ZFP is off (`use_zfp` not 1).
- The construction of a level runs on the GPU when all of these hold:
  - the matrix has a device kernel (`c_bpack_h2_set_gpu_kernel`);
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

## Environment switches

Each process reads each variable once, at first use, and the value then
applies to every matrix in that process.

### Performance and memory

| Variable | Default | Meaning |
|---|---|---|
| `HODLR_GPU_SPLIT` | `1` | How a block compressed by several ranks is split into pieces. `0` matches the CPU exactly (the block's own group, the same random numbers, and one piece per rank whatever `HODLR_gpu_pieces`), for comparisons. `1`, in the symmetric HODLR, builds each A21 over its node's group (both children's ranks), so the ranks of A12's rows, which the symmetric format never builds, are not idle, and applies `HODLR_gpu_pieces`. |
| `HODLR_GPU_KEEP_FACTORS` | `1` | `1` keeps the device copies of the low-rank factors from the construction until the forward blocks are assembled on the device, which then copies them device-to-device and never uploads the host arrays. `0` frees them after each level: less device memory during the construction, but the factors go down to the host and back up. `0` also turns off `HODLR_GPU_DEFER_HOST`. |
| `HODLR_GPU_DEFER_HOST` | `1` | `1` never fills the host copies of the factors built on the GPU; the device holds the only copy. Entry extraction takes the rows it needs from the GPUs (see above). A CPU multiply or CPU factorization of such a matrix stops with an error (TODO; with `HODLR_use_gpu > 0` both run on the GPU, so this only happens if a caller switches `HODLR_use_gpu` off after the construction). Saves the download time and host memory (3D Laplace, 262k points, 16 GPUs: 0.4 s of 7.1 s, and host memory per rank 2.5–3.4 → 1.7 GB). `0` downloads them during the construction. The download is always done with `HODLR_GPU_CHECK`, `ErrFillFull=1` or `HODLR_GPU_KEEP_FACTORS=0`. |
| `HODLR_GPU_EXCHANGE_MB` | min(2048, GPU memory / 16), divided by the ranks sharing the GPU | Size in MB of the device arena for MPI messages. It is reserved before the main device heap, and only when there is more than one rank and MPI is CUDA-aware. It holds the merges of the construction and the exchanges of the factorization, multiply and solve; a message that doesn't fit goes through pinned host buffers. |

### Validation and debugging

| Variable | Default | Meaning |
|---|---|---|
| `HODLR_GPU_CHECK` | `0` | `1` runs the CPU routine next to each GPU step and prints the difference: each construction level (drawing the same random numbers) and the dense leaves; the upload; the round trip of the factors back to the host (compared, then taken, so the CPU checks that follow use the fetched copies); the entries of the construction check, extracted again with the host factors overwritten by a sentinel so that every row the extraction reads must come from the GPUs; every multiply; the factorization (logdet); and the solves. `2` also checks the transposed products of each multiply. It runs the CPU construction and factorization as well, so the timings are not meaningful. |
| `HODLR_GPU_DEBUG` | `0` | `1` prints each block and chunk of the recompression QR. `2` also synchronizes after each step, so that an asynchronous device fault is reported by the kernel that caused it. |

### Shared GPU runtime (also used by the H2 GPU backend)

| Variable | Default | Meaning |
|---|---|---|
| `H2_GPU_HEAP_FRACTION` | `0.85` | Fraction of the free device memory the backend reserves as its heap at first use, split equally among ranks sharing a GPU. MAGMA's workspaces, the exchange arena and pinned buffers stay outside it. `H2 GPU heap exhausted` errors refer to this heap. |
| `H2_GPU_DEVICE_EXCHANGE` | on when MPI is CUDA-aware | `0` sends all MPI messages through host buffers. Otherwise device buffers go straight to MPI when `MPICH_GPU_SUPPORT_ENABLED=1` (Cray MPICH) or `H2_GPU_AWARE_MPI=1`. |
