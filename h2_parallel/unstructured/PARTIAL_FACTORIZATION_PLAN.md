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

## Status (2026-10-22)

- **P1** done (CPU, 1 and 4 ranks).
- **P2a** done on one rank (CPU), in `partial.hpp`:
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
    | 16^3, 1024 | 2945 of 3072 | 1.8e-6 | 5.1e-7 | 8.9e-7 | 1.3e-5 |
    | 20^3, 1200 | 5530 of 6800 | 7.2e-7 | 1.4e-6 | 7.0e-7 | 4.3e-6 |
    | 24^3, 2304 | 7217 of 11520 | 5.2e-7 | 8.8e-7 | 4.7e-7 | 6.3e-6 |
    | 32^3, 4096 | 10615 of 28672 | 2.0e-6 | 8.0e-7 | 1.4e-6 | 2.9e-5 |

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
  - Multi-rank partial runs stop with an error until P2b.
- Next: P2b (multi-rank: the skeleton-update exchange, remote view
  slicing, distributed up/down passes), P3 (GPU).
