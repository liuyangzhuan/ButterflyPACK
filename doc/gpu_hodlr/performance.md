# GPU HODLR: performance summary

Measurements of the GPU backend of the HODLR format (`format=1`, `LRlevel=0`,
option `HODLR_use_gpu`) on Perlmutter, September 2026. The backend runs the
construction, the factorization, the multiply and the solve on the GPUs, for
symmetric and unsymmetric matrices, in double and double complex, on any
number of nodes. Options and environment switches are in
[`GPU_BACKEND/hodlr_gpu/README.md`](../../GPU_BACKEND/hodlr_gpu/README.md); the run script is
[`example_scripts/runit_hodlr_perlmutter_gpu.sh`](../../example_scripts/runit_hodlr_perlmutter_gpu.sh).

## Setup

- **Nodes:** Perlmutter GPU nodes: 4 NVIDIA A100-SXM4 40 GB and one 64-core
  AMD Milan CPU each. CPU baselines run on the same nodes' CPUs.
- **Software:** GNU compilers, sequential OpenBLAS 0.3.30, MAGMA 2.10.0, CUDA
  13.2, Cray MPICH with CUDA-aware MPI.
- **GPU runs:** one MPI rank per GPU (4 per node, 16 OpenMP threads each).
  `HODLR_use_gpu=1` (FP64 GEMMs) in the first table, `HODLR_use_gpu=2` (FP64
  tensor-core GEMMs) in the others.
- **Problems:**
  - 3D Laplace on an n³ grid (`claplace3d_h2`), symmetric HODLR, BACA,
    `Nmin_leaf=64`, tol 1e-4 unless noted.
  - EMSURF: CFIE (α = 0.5) on a sphere of 128k triangles (`cie3d`), unsymmetric,
    BACA without overlap, tol 1e-6.
  - VIE: 64k unknowns (`cvie3d_h2`, `--scaleGreen 0`, h = 0.025), unsymmetric,
    tol 1e-6.
- **Times** are wall times from the solver's statistics. Construction
  includes all entry evaluation. The first run on a freshly allocated node
  costs about 0.9 s more (GPU warm-up); those runs are excluded or repeated.
- **Pieces:** the numbers were measured with `HODLR_gpu_pieces=1` (the
  default during development) unless noted; the default is now 4, which
  mostly helps the top levels on few GPUs (see the last table).

## GPU vs CPU HODLR on the same nodes, same MPI layout

4 nodes, 16 MPI ranks; CPU: 16 threads per rank (all 256 cores); GPU: one
rank per A100.

| Problem | Rank at level 1 | CPU constr | GPU constr | CPU factor | GPU factor | CPU solve | GPU solve |
|---|---|---|---|---|---|---|---|
| Laplace 64³ (262k), sym | 2686 | 1672 s | 6.7 s (250×) | 65 s | 0.65 s (100×) | 0.71 s | 0.074 s |
| EMSURF CFIE 128k, unsym | 2757 | 1239 s | 12.0 s (103×) | 155 s | 1.71 s (91×) | 0.79 s | 0.10 s |
| VIE 64k, unsym | 1865 | 172 s | 3.8 s (46×) | 28.5 s | 0.44 s (65×) | 0.066 s | 0.013 s |

The GPU and CPU runs build the same ranks and reach the same accuracy (entry
and matrix-vector errors of the order of the tolerance).

**Host memory per rank (MaxRSS):** the GPU runs keep the factors on the
devices, so their host memory is 1.5–2.0 GB per rank, against 3.6–5.2 GB
(Laplace), 9.6–10.3 GB (EMSURF) and 3.0 GB (VIE) for the CPU runs.

## GPU vs CPU HODLR, flat MPI on the CPU

3D Laplace, symmetric; the same node count for both; CPU: 64 MPI ranks per
node, one thread each; construction + factorization.

| Problem | Nodes | tol | GPU (4 ranks/node) | CPU (64 ranks/node) | Speedup |
|---|---|---|---|---|---|
| 48³ (110k) | 2 | 1e-3 | 2.0–2.3 s | 31.9 s | ~15× |
| 48³ (110k) | 2 | 1e-4 | 3.5 s | 134 s | 38× |
| 96³ (885k) | 16 | 1e-3 | 6.1 s | 189 s | 31× |
| 96³ (885k) | 16 | 1e-4 | 16.5 s (15.4 s with pieces 4) | 1395 s | 85–90× |

Separately: factorization 7–52× faster, solve 2–6× faster. The two
implementations reach the same ranks and log-determinants (to 1e-6
relative or better).

## GPU HODLR vs GPU H2 (Color)

The same build and the same 3D Laplace problems; symmetric; tol as listed.
H2 compresses and factors in one pass, so its time is its factorization; the
HODLR time is construction + factorization. H2 uses 216-point leaf boxes,
hence grids of 6·2ᵏ.

| N | GPUs | H2, tol 1e-3 | H2, tol 1e-4 | HODLR, tol 1e-3 | HODLR, tol 1e-4 |
|---|---|---|---|---|---|
| 110k (48³) | 1 | 0.81 s | 1.31 s | 7.0 s | 10.1 s |
| 110k | 8 | 0.91 s | 1.40 s | 2.0–2.3 s | 3.5 s |
| 885k (96³) | 8 | 1.88 s | 3.81 s | 22 s (14.0 s with pieces 4) | out of GPU memory |
| 885k | 64 | 1.73 s | 3.81 s | 6.1 s | 16.5 s (15.4 s with pieces 4) |
| 7.1M (192³) | 64 | 3.74 s | 11.8 s | out of GPU memory | out of GPU memory |

At the same tolerance H2 is also more accurate: its log-determinant agrees
with the tol 1e-4 H2 value to about 2e-7 at tol 1e-3, against 3e-4–5e-4
(tol 1e-3) and 4e-6–5e-5 (tol 1e-4) for HODLR. The HODLR solve, once
factored, is faster than H2's (0.02–0.12 s against 0.06–1.6 s).

**Memory limits of the GPU HODLR (40 GB A100s):** the ranks of 3D volume
problems grow quickly with N (level-1 rank 427 → 469 at tol 1e-3 and
2180 → 3458 at tol 1e-4, from 110k to 885k), and so does the storage of the
factors. 885k at tol 1e-4 needs more than 8 GPUs. 7.1M does not fit on 64
GPUs at either tolerance: the construction of a level-1 piece (442k × 442k)
needs its BACA factors grown to 4096 columns, 2 × 14.5 GB next to the old
copies.

## Scaling and balance of the top levels

| Laplace, tol 1e-3 unless noted | pieces 1 | pieces 4 |
|---|---|---|
| 885k on 8 GPUs | 22.5 s (slowest piece 6.8 s, fastest 0.23 s) | 14.0 s (pieces 2.6 s each) |
| 885k on 64 GPUs, tol 1e-4 | 16.5 s | 15.4 s (the level-1 merge over 64 ranks dominates: 5.2 s) |

Strong scaling of the same 885k problem, tol 1e-3: 22.5 s on 8 GPUs, 6.1 s on
64 GPUs (both pieces 1).

## Notes

- **Validation:** `BPACK_CHECK=hodlr` runs the CPU routine next to every GPU
  step (construction level by level, upload, entry extraction, multiply,
  factorization, solve) and prints the differences. In the development
  runs the ranks matched the CPU's, the logdet differences were at rounding
  level, and the multiply and solve differences were about 1e-15.
- **Statistics:** the GPU runs fill the same statistics lines as the CPU runs,
  counted the same way (flops of the same operations with the same
  formulas). Construction flops differ from the CPU's by 0.8–1.4× per level,
  where the GPU performs different operations; the CPU's own factorization
  flop count depends on the number of ranks.
- **Reproducing:** `example_scripts/runit_hodlr_perlmutter_gpu.sh` with
  `CASE=laplace GRID_SIZE=64`, `CASE=cfie MESH=<dir>/sphere_128000`, or
  `CASE=vie SCALE_GREEN=0 VIE_H=0.025` on 4 nodes, with `PIECES=1 USE_GPU=1`
  for the settings of the first table.
