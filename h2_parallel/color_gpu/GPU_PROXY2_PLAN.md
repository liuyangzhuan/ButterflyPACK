# Adaptive ID training rows (H2_ID_proxy 2) on the GPU: plan and milestones

Status (2026-10-01): A1 to A5 done and validated on 8 and 64 ranks (see "Open
points" for what has not been run). Keep this list current.

## What mode 2 does

`H2_ID_proxy 2` adds far rows to a box's ID target, chosen adaptively
(`append_adaptive_id_training_rows`, color_CA/factorization.hpp). The far
field of the box (every point more than `H2_ID_radius` boxes away, in the
Morton order of the tree) is sampled node by node, breadth first from the
root of the tree. For each node: `BACA_Batch` evenly spaced unused far points
are evaluated against the box's points, their residual under the current ID
is formed, the independent residual rows (a row ID with the tolerance) join
the target, and the node counts as converged when the residual norm is below
tolerance x the target's norm. A converged node is confirmed with up to 4
hold-out rows; a node that fails is subdivided (its children are queued), and
the ID is extended by an ID of the residual, or recomputed.

Until 2026-10-01 this ran only on the materialized target (all ring rows,
unsketched) with eager Schur updates: the streamed sketch was switched off for
mode 2 and `H2_lazy_schur != 0` was rejected. The GPU box path needs the
streamed sketch and lazy updates, so mode 2 never reached the device.

## Scope and decisions

- **The adaptive rows are selected against the sketch** (A1). The sketch
  `Y_B` (n x n) stands for the ring rows: the base ID, the residuals and the
  reference norm use it, and the selected rows are appended to it unsketched
  before the final ID. The rows are kernel rows beyond the ID neighborhood,
  which no Schur update reaches, so they need no lazy fill. This is the form
  the GPU runs; the CPU's streamed path is its reference. `H2_use_sketch 2`
  with mode 2 means this form on the CPU too (factorization and
  compression-only), which changes the results of earlier mode-2 runs with
  that option: agreed with the user on 2026-10-01.
- The materialized form stays available with `H2_use_sketch 0` or `1` (CPU
  only, eager updates).
- Symmetric matrices only, as the streamed path.
- Proxy modes need the global point index on every rank: `--distributed64 1`
  rejects them (unchanged).

## Milestones

Note on the times below: the unstructured loop's factor time (the spheres,
"VIE unstructured") still included the logdet and quick verification and the
first-use costs of the device kernels when these runs were made; since the
evening of 2026-10-01 it excludes them, as the full-grid loop does (sphere
9000, proxy 1, 8 ranks: 1.37 s then, 0.65 s now). Compare unstructured times
only within one of the two accountings.

### A1. CPU: mode 2 on the streamed sketch, lazy updates (done)
`gather_id_target_streamed` takes the ID tolerance and appends the adaptive
rows to `Y_B` (Color and CA levels, structured and unstructured loops, the
GPU path's host IDs). The option check that rejected lazy updates with mode 2
is gone. `BPACK_TRACE=id` prints per box `adaptive=frames:..,sampled:..,
holdout:..,appended:..,ids:..,residual_ids:..` (nodes sampled, rows
evaluated, rows added, IDs of the target after the first, IDs of residuals),
also on the materialized path. `claplace3d_h2` takes `--H2_ID_radius`,
`--H2_ID_proxy`, `--H2_ID_proxy_points`.

Results (build/gpu_prep/proxy2/p1, p2; 8 ranks x 16 threads on 2 GPU nodes;
EFIE spheres with cie3d, tol 1e-4, unstructured Color; VIE S2S 48^3 with
cvie3d_h2 `--scaleGreen 1`, tol 1e-3). acc_mvp / acc_forward / factor time:

| case | proxy 0 | proxy 1 | proxy 2 materialized, eager | proxy 2 streamed, eager | proxy 2 streamed, lazy 2 |
|---|---|---|---|---|---|
| sphere 2300 | 3.6e-5 / 9.3e-4 / 3.9 s | 3.7e-5 / 1.0e-3 / 3.8 s | 3.2e-5 / 7.6e-4 / 5.3 s | 1.4e-5 / 3.7e-4 / 6.5 s | 2.0e-5 / 5.5e-4 / 6.1 s |
| sphere 9000 | 2.6e-5 / 1.7e-3 / 11.5 s | 2.0e-5 / 1.6e-3 / 11.7 s | 6.2e-5 / 3.0e-3 / 21.6 s | 4.5e-5 / 2.3e-3 / 15.8 s | 9.5e-5 / 4.2e-3 / 12.5 s |
| sphere 128K | 1.5e-4 / 1.0e-1 / 95.8 s | 1.3e-4 / 1.1e-1 / 94.6 s | 1.1e-4 / 7.2e-2 / 170 s | 1.4e-4 / 9.4e-2 / 132 s | 1.1e-4 / 8.0e-2 / 98.6 s |
| VIE 48^3 | 7.0e-4 / 3.0e-3 / 13.5 s | | 7.1e-4 / 3.0e-3 / 23.2 s | | 6.6e-4 / 3.0e-3 / 13.7 s |

(proxy 0 and 1 with lazy 2; proxy 0 with eager updates: sphere 9000 3.5e-5 /
2.1e-3 / 14.3 s, 128K 1.3e-4 / 8.2e-2 / 130 s. The materialized VIE run used
`H2_XRR_factor 0`, where proxy 0 takes 12.3 s; the other VIE runs factor 1.)
On these cases mode 2 does not change the accuracy beyond the spread between
runs of the sketched IDs (a factor of 2-4 either way). On the streamed sketch
with lazy updates it costs 1-9% over proxy 0 on the three larger cases (55%
on the 2300 sphere, whose small boxes need many nodes); the materialized form
cost 35-90%.

What the selection does per box (the same picture for both forms): one node
(the root), 68 rows evaluated (64 samples + 4 hold-out; 104 with
`BACA_Batch 100`), 15-35 rows added, no further ID, for every box of the
128K sphere and of VIE 48^3. Only small, nearly full-rank boxes need more:
sphere 2300 (leaf n ~ 15) up to 47 nodes and 340 added rows in a few boxes,
sphere 9000 three boxes with 2-9 nodes.

### A2. GPU box path with the rows selected on the host (done; replaced by A3)
An intermediate step: with mode 2 the eliminator declined the device sketch
and took its host-ID path (`run_ids`: host sketch, adaptive rows and ID; T
uploaded; device elimination), on Color levels, structured and unstructured;
CA levels stayed on the host. `run_ids` still does this for a level with a box
too wide for the device sketch's packed lists.

Results (8 ranks on 8 A100, same cases): VIE 48^3 structured and
unstructured: k identical to the CPU in all 576 boxes, logdet agrees to 15
digits (197409.13276733446 vs ...452), same errors; factor time 7.5 s
structured / 8.9 s unstructured against 13.7 / 14.4 s on the CPU, and 1.05 s
for proxy 0 with the device sketch: the host sketch and ID are what is left
(level 3 on rank 0: ID 1.76 s of 2.1 s). EFIE spheres: 2300 identical to the
CPU (k in all 314 boxes, acc_mvp 2.0086e-5). 9000: the leaf level identical
(272 boxes), k differs in 33 of the 56 level-2 boxes. 128K: k differs in 50
of 1160 leaf boxes and in most boxes above. Proxy 1 shows the same pattern
(9000: 11 of 56 at level 2; 128K, September runs: 33 of 1160 leaf boxes): the
EFIE kernel is symmetric only to ~1e-6, so the device and the host differ in
the last digits of the updated blocks. Errors: 9000 7.4e-5 / 3.5e-3 (CPU
9.5e-5 / 4.2e-3), 128K 1.4e-4 / 9.2e-2 (CPU 1.1e-4 / 8.0e-2); 128K factor
time 78.4 s against the CPU's 98.6 s. Proxy 1 on the device sketch is
unchanged (sphere 9000: acc_mvp 3.8236954e-05 as before, 1.24 s).

Laplace 96^3 (claplace3d_h2, `--distributed64 0`, 8 ranks; build/gpu_prep/
proxy2/p3), Color on every level and a replicated CA leaf (`--CA_level 4`):
k and the adaptive statistics identical to the CPU in all 4672 boxes in both,
logdet -9928892.04119997 (Color; the CPU's to 16 digits) and
-9928892.2652035579 (CA leaf; every digit), GPU-vs-host solve check 1.2e-14.
Every box: one node, 68 rows evaluated, 19-26 added. Factor time: CPU 21.0 s
(proxy 0: 21.1 s); GPU 14.8 s with Color, 17.3 s with the CA leaf on the
host, against 1.96 s and 3.15 s for proxy 0 on the device sketch. That gap is
what A3 and A4 are for.

### A3. The selection on the device (done)
`AdaptiveRowSelector` (adaptive_rows.hpp), called by the eliminator's device
sketch (`sketch_wave`) for the wave's boxes together; kernels in
device_adaptive.cu (the residual under a list of IDs with its norm, the row
append) besides the element kernels and the device QRCP (which takes a
per-item tolerance for the residual IDs). Per box on the device: the target
W = [Y; added rows] in a block with room for the first node's rows (it grows
when a box needs more), and a block of the same shape for the in-place ID.
The decisions, the index of unused far points and the node queues stay on the
host, which steers by one small download per round.

A round of a box is one of: *sample* (the next node: its sampled and hold-out
rows through `eval_blocks`, so kind 3 goes through the triangle-pair
evaluator; residuals and norms; the row ID of the residual; the rows
appended; before it, where due, a failed hold-out's rows and a new ID of the
target), *recheck* (a node that converged under an ID with residual IDs: a new
ID, then the residuals again), *extra* (the ID of a residual), *final* (the ID
of the finished target unless it is current). The sampled points take scratch
slots at the end of the level's point table, rewritten every round
(coordinates from `id_source_point_coords`, which h2_initiate now keeps for
mode 2 with H2_use_gpu; kind 3 needs only the ids). The common case is two
rounds per wave (one sample round with the ID of the sketch, the final ID).
With the device sketch back, mode 2 has the device exchange, the kept solve
factors and the device transition like the other modes.

Results (build/gpu_prep/proxy2/a3, a4; 8 ranks on 8 A100, against the CPU
references of A1):

| case | k and adaptive statistics vs CPU | logdet | factor time: proxy 2 device / proxy 2 host rows (A2) / proxy 0 or 1 device / CPU proxy 2 |
|---|---|---|---|
| Laplace 96^3 | identical, 4672 boxes | -9928892.041199971 (CPU ...9729) | 1.68-1.74 / 14.8 / 1.58-1.61 / 21.0 s |
| VIE 48^3 structured | identical, 576 boxes | 197409.13276733446 (CPU ...452) | 1.26 / 7.5 / 1.05 / 13.7 s |
| VIE 48^3 unstructured | identical, 576 boxes | same | 1.70 / 8.9 / - / 14.4 s |
| sphere 2300 | identical, 314 boxes (boxes of up to 47 nodes included) | | 1.27 / 6.6 / - / 6.1 s |
| sphere 9000 | leaf identical; level 2 as A2 (33 of 56 differ in k) | | 1.36 / 11.6 / 1.37 / 12.5 s |
| sphere 128K | as A2 (50 of 1160 leaf boxes differ in k) | | 4.82 / 78.4 / 4.76 / 98.6 s |

GPU-vs-host solve check 1.2e-14 (Laplace). Errors as the CPU's (Laplace
acc_mvp 7.2453079e-04, VIE 6.6024707e-04 in every digit; sphere 9000 7.4e-5 as
A2; 128K 1.45e-4 / 1.0e-1). The selection costs 1-7% of the factor time: two
rounds per wave, "later rounds" 0.03 s of the Laplace leaf level's 0.3 s
sketch. The EFIE differences from the CPU are the device/host ones that proxy 1
has too (A2).

### A4. CA0 levels (done)
The ghost boxes of a replicated CA level are in the waves like the local
ones, so the same selection runs for them; `ca_level_runs` accepts mode 2.
Every copy of a box selects the same rows: the sampled positions depend only
on the tree and the box, the device ID gives the same bits in every launch
variant, and the selection's kernels sum in a fixed order per box
(device_adaptive.cu, compiled with --fmad=false like device_id.cu).

Results (build/gpu_prep/proxy2/a4; 8 ranks, `BPACK_CHECK=solve,replica`):

| case | replica check | vs CPU CA0 | factor time GPU (with the checks) / CPU |
|---|---|---|---|
| Laplace 96^3, CA leaf | 0 of 3536 ghost copies differ | k, statistics identical (4672 boxes); logdet to 15 digits | 2.53 / 20.6 s |
| Laplace 192^3, CA on levels 5 and 4 | 0 of 9488 and 0 of 3536 | identical (37440 boxes); logdet -89253589.210364759 (CPU ...789) | 12.8 / 123 s |
| VIE 96^3 (complex), CA leaf | 0 of 3536 | identical (4672 boxes); logdet to 15 digits | 4.58 / 66.9 s |

Solve checks 4e-15 to 1.2e-14. (Proxy 0, Laplace 192^3 with the same CA
levels, no checks: 9.15 s.)

### A5. Compression-only (--precon 2) (done)
CPU reference: `h2_skeletonize_box` with `H2_use_sketch 2` and mode 2 sketches
the ring rows (the draws of compute_id_sparse_sketch), selects the adaptive
rows against that sketch and appends them unsketched (`gather_id_workspace`
takes `adaptive_rows = false` for it); with `H2_use_sketch 0` or `1` the
materialized form as before. Device: `compress_level_ids`
(compression_gpu.hpp) with the same `AdaptiveRowSelector`, all boxes of a
chunk at once; the scratch points in `LevelPoints`; kind 3 rows through the
triangle-pair evaluator. Mode 2 with `H2_use_sketch 1` keeps the
compression's IDs on the host.

Results (build/gpu_prep/proxy2/a5; 8 ranks; `--precon 2`, spheres with
`--solve 0`): k and the adaptive statistics identical to the CPU in every box
of every case, and H2_CheckError(compression quick) acc_mvp identical in all
printed digits.

| case | boxes | acc_mvp (CPU = GPU) | construction time GPU / CPU |
|---|---|---|---|
| Laplace 96^3 | 4672 | 8.1426819e-04 | 1.35 / 3.76 s (proxy 0 on the device: 1.24 s, 7.9e-4) |
| VIE 48^3, structured and unstructured | 576 | 6.5163054e-04 | 0.39 / 1.02 s |
| sphere 2300 | 314 | 1.6450503e-06 | |
| sphere 9000 | 328 | 1.3720213e-06 | |
| sphere 128K | 1488 | 1.8379402e-07 | 0.81 / 57.8 s (proxy 1 on the device: 0.75 s, 2.0e-7) |

Device matvec check (Laplace) 4e-16 to 7e-16.

### 64 ranks (build/gpu_prep/proxy2/x64; 16 nodes, 64 A100 40 GB)
Laplace 192^3 (`--distributed64 0`), factor time, two runs each:

| | proxy 0 | proxy 2 | CPU proxy 2 |
|---|---|---|---|
| Color | 3.34 / 3.43 s | 3.43 / 3.43 s | 42.7 s |
| CA leaf (`--CA_level 5`) | 3.80 / 3.76 s | 4.05 / 3.94 s | |

Color with proxy 2 against the CPU: k and the adaptive statistics identical in
all 37440 boxes, logdet -89253588.161176607 and acc_mvp 7.6101737e-04 in every
digit; solve checks 4e-15 and 1.3e-14. CA leaf: 0 of 49824 ghost copies
differ. 512K EFIE sphere: factorization 4.98 s with proxy 2 against 5.46 s
with proxy 1 (acc_mvp 7.7e-4 / 6.0e-4); compression-only with proxy 2 0.90 s,
acc_mvp 7.8e-8.

## Open points
- A box that needs many nodes costs one round each; rounds are batched over
  the wave's boxes, so the count is the largest over the boxes (47 on the 2300
  sphere: 95 rounds over the leaf level's waves, still 1.27 s in total). If
  that ever matters: a cap, the rest finished on the host.
- The global index and coordinates are O(N) per rank; with `--distributed64 1`
  there is no global point list (proxy modes are rejected there).
- Not run yet: mode 2 beyond 64 ranks; `H2_ID_radius > 2` together with mode
  2 on the device (static ring rows plus adaptive rows); the Matern and Python
  drivers; a level that falls back to host sketches (a box wider than the
  device sketch's packed lists).
