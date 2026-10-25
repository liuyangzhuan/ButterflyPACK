# Partial factorization of an H2 front (F11 eliminated, F22 retained)

One H2 over a multifrontal front: the separator DOFs (F11) followed by the
update DOFs (F22), one Morton tree, one set of entry callbacks.  `factor`
runs the RS-S factorization but eliminates only F11 DOFs.  The object then
serves what STRUMPACK's front needs from its four blocks today (F11 solves,
F12/F21 applications, the Schur complement CB = F22 - F21 F11^-1 F12 as
entries and matvecs, a release of the CB) without F12/F21 butterflies and
without a sampled S.

Decisions (2026-10-22): nonsymmetric code only (the symmetric path is
untouched; `H2_partial = N` stays bitwise); the production configuration
only (`H2_use_sketch 2`, `H2_lazy_schur 2`, LU, Color, `H2_unstructured 1`,
CPU and GPU); F11 is the first `H2_partial` points in input order; the
cost of serving F22's entries from the children (STRUMPACK's `extract_2d`)
is a separate concern.

## Vocabulary

- **eliminated**: an F11 redundant DOF, factored out as today.
- **sparsified**: an F22 redundant DOF.  Its ID is applied (its coupling to
  the far field and to every live F11 DOF is moved onto its box's skeleton)
  and nothing else: it is not eliminated, it updates nothing, it generates
  no fill.  Its rows and columns are kept as the *sparsified blocks*.
- **skeleton**: as today, per part.

## Why parts are per point, and why F22 is compressed

The separator is a slab through the subdomain and the update set is the
subdomain's boundary surface, so most boxes are pure but every box along
the slab's perimeter holds both kinds of points, at every level.  A box
therefore has an F11 part and an F22 part (skeleton/redundant per part);
"the box is not eliminated" means "its F22 part is not eliminated".

If F22 points only ride up as skeletons (the level-1 state), every box's
F22 share grows with the box: level-2 near blocks carry all of F22 densely
and the root holds an n_upd x n_upd Schur complement, at the dense
multifrontal cost.  So the F22 part of every box is skeletonized as well,
and its redundant DOFs are sparsified rather than eliminated.

## The algorithm

Per level, per box, in the existing wave order:

1. **F11 part: unchanged.**  ID of the box's F11 columns against the far
   field (streamed sketch, lazy fill), X_RR factored, Schur updates on the
   one-hop neighbours' live DOFs (which include F22 DOFs), lazy two-hop
   generators: `compute_and_modify` and `compute_step_two_internal` on the
   F11 columns.

2. **F22 part: ID and sparsify.**  ID of the box's F22 columns against the
   far field **plus every live F11 DOF of the box and of its one-hop
   neighbours**.  Form X_RR, X_SR, X_NR (and X_RS, X_RN) as the F11 step
   does, factor nothing, update nothing, generate nothing; keep them as the
   sparsified blocks.  From then on the sparsified DOFs are treated as
   eliminated for every far-field purpose (sketches, generators, the
   transition) and the F22 skeleton goes to the parent.

   The extended target is what makes this safe.  After the ID the
   sparsified DOFs are decoupled from every F11 DOF they were coupled to
   (far ones by the usual ID, near ones by the extra rows), so no later F11
   elimination, at this or any coarser level, fills them: eliminating R_P
   fills M(x, R_P) M(R_P, R_P)^-1 M(R_P, y) only for x, y coupled to R_P.
   Without the extra rows the fill on the sparsified rows grows with every
   coarser level up to the root's near field.  The extra rows cost at most
   27 x k11 modified near-block rows, about the far ring already sketched.

3. **Root:** partial dense LU of the root matrix on (root S11, root S22):
   LU of the S11 block, L21, U12, and the dense CB block on S22.

**After the factorization** the object holds, for the solves, the F11
factors (T, X_RR LU, X_SR, X_NR, X_RS, X_RN) and the F22 interpolations
T22 (an F22 part applies only `x_R -= T^T x_S` / `x_S -= T x_R`); and for
the CB, T22 at every level, the sparsified blocks at every level and the
root CB block.  Those are an H2 matrix of CB in interpolative form (nested
bases T22, block-sparse couplings over the sparsified DOFs: the extended
sparsification of CB).  Entries and matvecs come from the same up-pass of
e_i / T(:, j) vectors that `extraction.hpp` runs on a compressed H2,
reading a sparsified block at the level where the coarser of the two DOFs
was sparsified, plus the root block.  `release` frees the sparsified
blocks and the root CB block.

**Solves.**  `partial_forward([b1; b2])`: the forward sweep over every
level and the root's L11 part; returns `b2 - F21 F11^-1 b1` on the F22
DOFs and keeps the swept vector.  `partial_backward(x2)`: the root's
`x1 = U11^-1 (y1 - U12 x2)`, the diagonal solves, the backward sweep with
x2 scattered to the F22 slots; returns x1.

## What it gains

F12/F21 butterflies and the sampled S go away (at 400^3 after step 1:
F12 3320 + F21 1353 + Schur 978 of 31498 rank-seconds, and F12/F21's
4.2e9 nnz); CB is as accurate as F11.  The front's F22 kernel entries are
still needed, so the F22 serving cost stays.

## Changes by file

| Where | Change |
|---|---|
| option plumbing (as `H2_s_op`) | `H2_partial`: the number of leading input points forming F11 (0: full factorization); membership by original id through `new2old` |
| `core/tree.hpp` BoxData | skeleton/redundant per part, T22, sparsified blocks, a `sparsified` state |
| `core/factorization.hpp` `compute_and_modify` | a column part, extra target rows, a sparsify mode that stops after X_RR/X_SR/X_NR (and X_RS/X_RN) and skips `compute_step_two_internal` |
| `unstructured/factorization_impl.hpp` | the F22 phase after the level's last wave (as built: see Status); sparsified DOFs leave the live set; the sparsified blocks survive `clear_modified_interaction_matrices`; root partial LU |
| `unstructured/parent_transition.hpp` | children's live lists are S11 and S22 |
| `unstructured/solver_impl.hpp`, `core/solver.hpp` | T-only F22 parts; `partial_forward` / `partial_backward` |
| `unstructured/extraction.hpp` + a CB matvec | CB entries and matvecs from (T22, sparsified blocks, root block); build state `RS_PARTIAL` |
| `SRC/C_BPACK_wrapper.cpp` | `c_bpack_factor` -> partial; the two partial solves; `c_bpack_extractelement` / `c_bpack_mult` on a partial object = CB; `c_bpack_release_schur` |
| `GPU_BACKEND/h2_gpu/level_eliminator.hpp`, `device_solve.hpp` | the device F22 box step (extended sketch rows, stop after the X blocks, sparsified blocks to the host), `build_parent` over S11 and S22, partial `factor_root`, T-only boxes and partial sweeps in the device solve |
| logdet, quick verification | of F11 (the F11 pivots); the partial solves checked against the kernel |

## Milestones (each validated before the next)

- **P1** CPU: parts, F11-only elimination with the F22 parts all-skeleton
  (dense in F22: validates the partial elimination, the root partial LU
  and the partial solves in isolation), CB entries/matvecs from the root
  block.  Checks: `H2_partial = N` bitwise equal to today; at N ~ 10K the
  partial solves and CB against a dense LAPACK CB of the kernel matrix.
- **P2** F22 sparsification with the extended target; CB extraction and
  matvec from the sparsified form.  Same checks; ranks and memory against
  P1 and against F22 as its own H2.
- **P3** GPU.  CPU = GPU to 1e-12..1e-15.
- **P4** release API; multi-rank (the F22 steps ride the owner-deferred
  machinery); `H2_compress_factors` on the sparsified blocks.
- **P5** hand-off to the STRUMPACK session: a `FrontH2` with one object.

Test driver: `Laplace3D_H2Test_Driver --H2_partial n11` on a slab-plus-faces
point set and on a plain cube split, dense checks at small N.

## Status (2026-10-24)

- **P1** done (CPU, 1 and 4 ranks).
- **P2a** done (CPU), in `partial.hpp`:
  - The F22 step runs once per level, after the level's last wave.  The
    owner-deferred replays of step 7 land only after a wave, so a box's
    F22 step inside its wave would read stale neighbour blocks.  Every box
    reads its neighbours' pre-phase layouts (`skeleton_before_f22`), and
    the records X_NR/X_RN are over the neighbours' raw F22 points.
  - The ID target: the far sketch over the pre-phase skeleton
    (`gather_id_target_streamed_nonsymmetric` with `column_slots`), the
    near live F11 rows and the box's own S11 rows, both orientations.
  - CB entries walk the chains t = P e (per level p = r[R22],
    w = r[S'] + T p); the CB product is an up-pass with P, the per-level
    records and the root block, and a down-pass with P^T.
  - Accuracy at tol 1e-6 against dense LAPACK (skew Laplace, one rank):

    | grid, F11 | F22 at root | CB product | CB entries | forward | backward |
    |---|---|---|---|---|---|
    | 16^3, 1024 | 2945 of 3072 | 1.8e-6 | 5.1e-8 | 8.9e-7 | 1.3e-5 |
    | 20^3, 1200 | 5530 of 6800 | 7.2e-7 | 2.6e-7 | 7.0e-7 | 4.3e-6 |
    | 24^3, 2304 | 7217 of 11520 | 5.2e-7 | 1.9e-7 | 4.7e-7 | 6.3e-6 |
    | 32^3, 4096 | 10387 of 28672 | 2.0e-6 | 1.5e-7 | 1.4e-6 | 2.9e-5 |

    (after the P2b fix; the backward error is the conditioning of
    A11, cond1 360-710 here: its residual is at the tolerance)

    `H2_partial = N` stays bitwise equal to the full factorization.
  - CB entries are batched: the requested points' chains go up as blocks
    sharing a support (one box, one child), each term a gathered record
    slice times two GEMMs per pair of blocks in the same or one-hop boxes,
    parallel over row blocks.  At 64^3 (one rank, 16 threads): 5e7
    entries/s for localized 2000 x 2000 requests, 9e5 entries/s for
    scattered ones (these read nearly all of X_RN/X_NR: bandwidth-bound).
  - The CB product runs every pass over a level's boxes in parallel.  It
    reads every record once: 49 GB at 64^3 (levels 3 and 2, F11 = 16384),
    about 1.5 s per product at the 35-40 GB/s that 16 or 32 threads reach
    (6.9 s single-threaded before); 3 right-hand sides cost the same as 1.
- **P2b** done (CPU, any rank count; the F11 waves' Color/CA machinery
  unchanged):
  - The F22 step reads a remote neighbour through its copy (an assisting
    box on a Color level, a ghost on a CA level) and its views of this box,
    else the neighbour's views (a local neighbour), else the kernel (a pair
    no elimination touched).  After the step each rank publishes its
    sparsified boxes' S_old and kept slots to the ranks within three hops;
    the receivers reduce their copies' live sets (`skel_before`,
    `skeleton_before_f22` keep S_old for the records' row lookup).
  - Solves: the F22 steps between the level transfers
    (`scatter_solution_to_children`, `gather_skeleton_to_parent`), the
    outputs summed over ranks.
  - The CB product and CB entries are collective over the level
    communicators: a level's chains or vectors go to the owners of the
    remote one-hop neighbours, each record term is formed by the box that
    owns it, the result is summed and replicated.  Every rank passes the
    same request lists.
  - The bug behind wrong x1 or NaN in about half the 64^3 runs on 4 and 8
    ranks (and small run-to-run differences on one rank): the F22 step's
    far sketch located the box's rows in a fill source into a scratch
    vector that the ring loop then reused, so the box's columns read the
    last ring box's positions and, past its length, stale values from the
    thread's earlier boxes.  Fixed by keeping the box's positions apart.
    The fix also sparsifies more at the coarse levels (64^3: level 3 kept
    86475 -> 85647, level 2 17601 -> 16320) at the same accuracy.  Found by
    dumping each box's F22 target per run (sketch, near and own parts, per
    neighbour and per fill source) and comparing runs: levels that agree
    bitwise isolate the step that does not.  (Take such norms before
    compute_id_complex: its QR overwrites the target, and the overwritten
    rows mix in the sketch.)
  - Also made safe, not a cause: the nonsymmetric owner-deferred mirror
    could create a box's view inside its parallel pass while another
    mirror reads that box's block list.  Every mirror target's views are
    now made before the parallel pass, as the symmetric mirror already
    requires.  Counted in the old code: no view was ever created there
    (64^3 partial, 40^3 full, 8 ranks).
  - Validation (tol 1e-6, skew Laplace; dense LAPACK where N allows, else
    the in-library checks against the kernel):

    | grid, F11 | ranks | F22 at root | CB product | CB entries | forward | backward |
    |---|---|---|---|---|---|---|
    | 20^3, 1200 | 2, 3, 4, 8 | 5528 of 6800 | 6.9e-7 | 2.3e-7 | 6.7e-7 | 4.2e-6 |
    | 32^3, 4096 | 8 | 10405 of 28672 | 2.0e-6 | 1.1e-7 | 1.4e-6 | 3.0e-5 |
    | 48^3, 9216 | 1 / 8 | 14155 / 14150 | 1.1e-6 / 1.2e-6 | | 8.6e-6 / 8.7e-6 | |
    | 64^3, 16384 | 1 / 4, 8 | 16317 / 16320 | 2.4e-6 / 2.3e-6 | | 4.6e-5 | |

    The 20^3 numbers agree to the printed digits on 2 to 8 ranks.  At
    64^3, repeated runs on 1, 4 and 8 ranks each give the same F22 ranks
    and the same errors (before the fix about half of the 4- and 8-rank
    runs were wrong).
    `H2_partial = N` stays bitwise equal to the full factorization.
- **P3** GPU (H2_use_gpu 1), in steps:
  - **P3a** the F11 waves on the device, everything partial on the host.
    The device ID of a box takes its F11 columns only, as
    `compute_id_partial`: the sketch's columns are permuted to the F11
    ones first, the QRCP runs over those, and T's rows of the F22 slots
    are zeroed; the skeleton lists [F11 skeleton in pivot order, F22
    slots] go to the device for the later waves' sketches.  A partial
    level keeps no solve factors on the device; at its end its blocks and
    factors come back to the host, and the F22 step, the transition, the
    root, the partial solves and the CB run on the host as on the CPU
    path.  A box the device never formed a self block for (nothing
    eliminated, no update) takes the kernel's in the F22 step.  The
    level-end download of the blocks goes in chunks (one buffer for all
    of them doubled the level's device footprint).
    Done 2026-10-23: GPU = CPU to the printed digits of every dense check
    (16^3 to 32^3, 1, 2, 4 and 8 ranks), logdets within 1e-15 to 2e-14;
    `H2_partial = N` within 4e-15; the full GPU factorization unchanged.
    64^3, F11 = 16384, A100 40 GB:

    | ranks | CPU | GPU | F22 steps (host) | F11 waves CPU / GPU |
    |---|---|---|---|---|
    | 1 | 211 s | 144 s | 86 s | 84 / 19 s |
    | 4 | 126 s | 74 s | 27 s | |
    | 8 | 147 s | 97 s | 41 s | |

    The host F22 step (BPACK_TRACE=phase, 64^3 one rank, thread-seconds
    at levels 3 / 2): sketch 576 / 625, ID 96 / 477, records 295 / 266,
    the target's near and own rows 6 / 6.
  - **P3b** the F22 step on the device (`LevelEliminator::f22_boxes`),
    done 2026-10-24.  `partial_f22_phase` plans the boxes on the host (S_old,
    the F11 / F22 positions, the neighbours' positions) and hands them to
    the device, which per chunk of boxes:
    - copies the boxes' S_old points to the point table's scratch and
      sketches them with `sketch_wave`'s F22 columns: the columns are those
      copies, a fill source that recorded the box's full point set has its
      P operand gathered at S_old, and the sketch returns after the fill (the
      same draws and sources as the host's `column_slots` gather);
    - assembles each ID target from the sketch, the device near blocks
      K(nb, b) and K(b, nb) and the self block (the kernel where the device
      never formed one), runs the QRCP, and sorts the ID's columns on the
      host as `compute_id_complex`'s callers do;
    - forms T and the records in four batched GEMM passes and downloads
      them.
    The host then publishes the kept slots, the level's blocks come back
    sliced to the new live sets as they are gathered for the download (the
    host slices nothing), and the host transition runs.  A partial level's
    fill sources stay on the device to its end (the F22 step reads them
    all).  The root's W solve and CB product now run a column block a thread
    (they were single-threaded: 25 s of the 62 at 64^3).
    GPU = CPU to the printed digits of every dense check, 1 to 8 ranks; at
    64^3 the same F22 ranks and errors as P3a.  64^3, F11 = 16384, A100 40 GB
    (CPU: 64 threads a node; both with the parallel root):

    | ranks | CPU | GPU, P3b |
    |---|---|---|
    | 1 | 189 s | 40 s |
    | 4 | 105 s | 30 s |
    | 8 | 129 s | 36 s |

    Left at 64^3 on one GPU (40 s) after P3b: the host transitions and the
    next levels' block uploads, the records' download, the device QRCP.
  - **P3c** done 2026-10-24: the partial levels' parents built on the
    device.  After the F22 step the level's self and near blocks are sliced
    on the device to the new live sets (`LevelEliminator::f22_slice`: a
    gather at the kept positions, the old block freed) and the sparsified
    boxes' device skeletons set to them (the remote ones' from the publish);
    the device transition reads a child's rows in a source that recorded
    them at the live set before the step through that set (a gather at the
    kept positions).  Level 1's transition runs on the device too, and the
    root's partial LU runs there (`partial_root_on_device`: A11 = P L U, W,
    CB by MAGMA), its factors to the host.  The F22 records go to the host
    through a background copier (fixed pinned slots, its own stream) while
    the next chunks run: growing one pinned buffer to a chunk's records had
    cost 5.5 s at 64^3 (pinning ~1 s/GB).
    GPU = CPU to the printed digits of every dense check, 1 to 8 ranks.
    64^3, F11 = 16384 (A100 40 GB; CPU 64 threads a node):

    | ranks | CPU | P3a | P3b | P3c |
    |---|---|---|---|---|
    | 1 | 189 s | 144 s | 40 s | 12.4 s |
    | 4 | 105 s | 74 s | 30 s | 8.1 s |
    | 8 | 129 s | 97 s | 36 s | 14.0 s |

    (8 ranks: 2 a GPU; the factor time includes the GPU-aware MPI warm-up,
    about 6.5 s.)  Left on one GPU: the F22 targets' QRCP (4.1 s of 12.4,
    bandwidth-bound: one pivoted Householder pass a column), the F11 waves
    (2.3 s), the F22 sketches (1.5 s).  The records stay host-side (49 GB at
    64^3: more than a 40 GB card), so the partial solves and the CB product
    run on the host.
  - **GPU follow-ups**, from a 64^3 one-GPU profile (2026-10-24; nsys and a
    per-box dump of the F22 targets).  The step's "targets and ID" is all
    QRCP kernel time (level 3 1.61 s, level 2 2.49 s); under a traffic model
    of 16 sum (m-i)(n-i) bytes a box the kernels reach 282 / 371 GB/s, and
    the column update reads every element twice and the reflector once.  The
    traffic is in the tall targets of the boxes in or next to the F11 slab
    (their exact near F11 rows): at level 2, 32 boxes up to 4615 x 1792 with
    ranks 500-715 carry 90% of it, at level 3, 128 boxes 61%.  The records'
    copies ran at 7 GB/s (appended on one thread into fresh vectors; the D2H
    itself at 20-24), so ~2 s of the step waited for them.  X_NR's rows at
    the sparsified neighbours' R positions are never read (the CB product and
    entries take those couplings from the neighbour's X_RN): 17.5 of the
    49.5 GB.  In this order:
    (d) the records' copy (done, below);
    (b) a two-phase F22 step (every box's ID, the publish, then the records),
        the ID launches grouped by size instead of by the records' chunks
        (12 of level 3's 29 QRCP launches ran on 21-28 of the 108 SMs);
    (e) X_NR without the sparsified neighbours' R rows (needs (b));
    (a) a blocked QRCP (dlaqps-style: one read of the trailing matrix a
        step, a GEMM update a panel).
    A QR before the QRCP of the tall targets halves their traffic, but
    cuSOLVER's geqrf runs at 0.3-1.2 TFLOP/s on these shapes (1.1 s for the
    two levels); a smaller far sketch (d = |S_old| + 8 now, mean ranks 167 of
    474 and 255 of 1338) would change the sampling: not pursued.
  - **(d)** done 2026-10-24: the records are `PartialRecord`s (core/tree.hpp),
    not zero-filled when sized; the copier sizes them on its thread and copies
    each 64 MB slot into them in 1 MB pieces over its threads
    (`HostCopier::Segment::into`), and populates the reserved pages of large
    vector destinations before appending to them (MADV_POPULATE_WRITE, every
    copier user).  64^3 on one GPU, runs interleaved on one node: 12.33 /
    12.37 s before, 10.65 / 11.22 s after, and 11.98 s after ten other runs
    on the node: the copies are page-fault bound (9-14 GB/s, vs 7) and the
    faults' cost varies with the node's memory; the step's waits for them
    went from ~2 s to 0.4-1.6 s.  Checks and logdets unchanged (GPU and CPU
    paths, 1 to 8 ranks).  Measured and dropped (means of two interleaved
    runs): zero-filled destinations (a parallel copy into them ran at 9-10
    GB/s), populating plus appending (11.55 s), the owner sizing the records
    on its own thread (it stalled 0.2-1.1 s a level against the copier's page
    faults; 11.23 s), populating the sized ones first (11.48 s).
  - **(b) + (e)** done 2026-10-24: the F22 step in two passes, on both
    paths: every box's ID and own records (T, X_SR, X_RS, X_RR, X_RN), then,
    the kept sets known and published, every box's X_NR over only its
    neighbours' F22 points that stay live (`F22Part::neighbor_nr_*`, which the
    CB product and entries now read X_NR by).  The neighbours' live sets
    before the step are copied first, since the publish reduces the remote
    copies'.  On the device (`f22_ids`, `f22_records`) a level's IDs go in
    four round-robin chunks, each box's QRCP launched with boxes of its size
    (a few large ones on the cooperative path, the others up to one a
    multiprocessor), and each chunk's own records copy to the host behind the
    next chunks.  64^3 on one GPU, runs interleaved with (d) on one node:
    11.31 / 11.01 s -> 10.95 / 10.89 s; records 49.5 -> 32.0 GB (host RSS
    53.1 -> 36.0 GB); targets and ID 1.67 -> 1.35 s at level 3, 2.50 -> 1.96 s
    at level 2.  The copier still bounds level 3 (17.8 GB at ~7.3 GB/s there,
    2.5 s against 1.9 s of IDs).  The CB product took the same time (one
    product ~1.5 s: not its X_NR bytes).  Checks unchanged, GPU and CPU paths,
    1 to 8 ranks.  Dropped on the way: all IDs, then all records (13.1 s: the
    records' copies had no IDs left to overlap with), and huge pages or 32
    threads for the records' copy (no gain; slower).
  - **(a), first part** done 2026-10-24: Nsight Compute on level 2's
    cooperative QRCP launches (64^3) showed them latency-bound, not
    bandwidth-bound: 37% of the DRAM bandwidth, no eligible warp 92% of the
    cycles, 60% of the stalls on a global load - a warp's loop over a column
    kept one load in flight.  The column loops (the reflector's norm and
    scaling, its application) now load four rows a lane ahead of their sums,
    in the same order (the same results bit for bit; the batch test passes),
    and the cooperative kernel stages the step's reflector in shared memory.
    64^3 on one GPU, interleaved on one node: 10.44 / 10.56 s -> 9.20 /
    9.45 s; targets and ID 1.96 -> 1.35 s at level 2, 1.33 -> 0.86 s at
    level 3; level 2's launches at 52% of the DRAM bandwidth.  Eight rows
    ahead gained nothing more, four blocks a multiprocessor 0.06 s.  Tried on
    the copier and dropped: a second thread setting up the records' pages
    ahead of the stores (no net gain: their D2H then waits for it, and the
    tail - X_NR after the publish - stays).
  - **(a), blocked QRCP** done 2026-10-24: `launch_qrcp_blocked`
    (device_id.cu), LAPACK dlaqps / zlaqps on the device: a step reads the
    trailing matrix once (F's column), a panel of 32 steps updates it once
    (A -= V F^H); the step's own rows (the pivot column's update, the
    reflector's norm, V^H x) in chunks of 256 rows over the blocks, each
    block forming the reflector from the chunks' sums: two grid barriers a
    step, as before, one more at a panel's end.  As dlaqps, it ended the
    panel at every cancelled norm downdate, and on these targets one
    cancelled every 4-5 steps (90% of the panel ends; a test of eps^(3/4)
    instead of eps^(1/2) made it every 4): no faster than the unblocked
    kernel.  A cancelled column now keeps an upper bound of its norm
    (n1 sqrt(temp + 256 sqrt(eps))) until the panel's end, and the panel ends
    early only if a bound would be the next pivot: panels of 29-30 steps, 2-3
    early ends a large box.  The F22 IDs use it from 256 columns on (below,
    its extra barriers cost more: level 4 0.04 -> 0.07 s); the primary waves
    keep `launch_qrcp`.  64^3 on one GPU, three runs each interleaved on one
    node: 9.64 / 9.10 / 8.99 s -> 8.70 / 8.67 / 8.62 s; targets and ID 1.35 ->
    0.95 s at level 2, 0.86 -> 0.75 s at level 3.  Checks and logdets as
    before, 1 to 8 ranks; the batch test covers it (tall targets, every
    launch variant bit for bit, ID residuals and ranks as the unblocked
    kernel's).  The batch test also caught a race on the way (a block that
    ended the panel rewrote its candidate while another block of the box
    still read the candidates to decide the same; memcheck then showed the
    other block's panel running past 32 steps): the panel end's candidates
    are a second set.
- **P4** done 2026-10-25 (multi-rank was already P2b and P3):
  - The release (`c_bpack_partial_schur_release`, `partial_schur_release`)
    frees what only the Schur complement reads: the root's CB block, every
    sparsified box's records and its neighbour lists, the entries' leaf
    index; the partial solves keep the F11 factors, T, the root's LU, W and
    A21.  BPACK_TRACE=memory prints each rank's bytes (H2_PARTIAL_RELEASE).
    The driver now checks after the dense checks that the solves give the
    same results bit for bit after a release and that the Schur calls are
    refused.  64^3: 34.0 GB freed of 40.5 GB resident (31.7 records, 2.1
    root block); with glibc's allocator the resident set drops 39.6 ->
    11.7 GB (1 rank), 43.6 -> 12.7 GB (4 ranks); tcmalloc keeps the pages
    for the process's next allocations.
  - X_NR is one record a neighbour (`F22Part::X_NR[one_hop slot]`, m x r
    each): the CB product and entries read it a neighbour at a time.  On the
    device a neighbour's block is formed and copied apart.  Every check the
    same bit for bit.
  - `H2_compress_factors 1` keeps the records compressed
    (core/partial_record_codec.hpp): lossy Bitcomp in the format of the GPU
    kept factors, a step of 0.02 tol_comp of each record's largest entry,
    512 KB chunks.  The device compresses a chunk's records before their copy
    (`f22_submit_records`, in groups of 1 GB of staging; a group without
    staging goes raw, BPACK_TRACE=phase counts them), the host path a box's
    as it forms them; the product decompresses each record where it reads
    it, the entries only the chunks of the entries they read.  T stays raw
    (the solves), and the F11 waves' own device compression stays off (the
    F22 step reads the near blocks raw).  The step: 32^3 against dense
    LAPACK, schur_mult / schur_extract 2.00e-6 / 1.54e-7 raw, 2.67e-6 /
    6.82e-7 at 0.2, 2.04e-6 / 2.32e-7 at 0.05, 2.01e-6 / 1.68e-7 at 0.02.
    64^3 on one GPU, two runs each interleaved on one node: records 31.6 ->
    5.9 GB, resident set 40.5 -> 22 GB, factor time 9.83 / 9.77 -> 8.49 /
    8.56 s (the copier moves 5.9 GB: level 3's wait for it 0.41 -> 0.03 s);
    a CB product 1.16 -> 3.45 s and the driver's entries 0.022 -> 0.080 s
    (host decompression, about 0.5 GB/s a thread); 4 ranks: 4.27 -> 3.85 s,
    48.2 -> 26.4 GB, product 0.37 -> 1.09 s.  GPU and CPU paths agree to the
    quantization (the records' rounding-level differences can land a step
    apart).
- Next: P5 (the STRUMPACK hand-off).
  Later: the records on the device (with compression, or several cards) for
  a device CB product and partial solves; with compression they would fit
  (5.9 GB at 64^3) and the device decompresses at hundreds of GB/s.
