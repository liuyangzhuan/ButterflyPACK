# Option B — Matrix–Matrix Solve & Multiply for the H2 Solver (nrhs)

**Goal:** Extend the H2 (format‑7) path from single‑vector solve/multiply to **matrix–matrix** (multiple RHS columns, `nrhs`), so one hierarchical sweep handles all columns at once — mirroring how ButterflyPACK's HODLR direct solve already works (`x`,`b` declared `(Ns_loc × num_vectors)`, one `BPACK_Inv_Mult` call with `num_vectors`, BLAS‑3 per block).

**Why "one sweep":** the tree is walked once and each per‑box dense op becomes a `gemm`/`trsm`/`trmm` on an `npoints × nrhs` block instead of `nrhs` separate `gemv`/`trsv`/`trmv` calls. This amortizes the tree traversal *and* batches the MPI messages (one round carrying `nrhs` columns vs. `nrhs` rounds). BLAS‑3 is also compute‑bound and thread‑scalable, unlike memory‑bound BLAS‑2.

---

## Design decisions (settle first)

1. **Layout: column‑major, leading dimension = `num_points`.** Each box's `right_side`/`left_side` becomes length `num_points × nrhs`; element `(i, c)` at index `i + c*num_points`. Matches BLAS (`ld = num_points`), HODLR's `(Ns_loc × num_vectors)`, and the driver's `b_s[i + c*myseg]`.
2. **Carry `nrhs` on `SolveDataRequest`.** Add `int64_t num_points; int64_t nrhs;` fields (set in `initialize`). Today column count is inferred from `right_side.size() == num_points`; storing it lets every helper compute the leading dimension **without adding an argument to every signature**, so helper *call sites* don't change.
3. **Scratch stays gather/scatter.** Skeleton/redundant rows (`x_R`, `x_S`, `neighbor_values`) are *scattered* row subsets, so keep the "gather to contiguous scratch → BLAS → scatter back" pattern — just make the scratch `num_R × nrhs` (column‑major, `ld = num_R`) and the gather/scatter loops 2D over (row, column). Do **not** try to stride the original block.
4. **Input/output buffers: column‑major `(Nloc × nrhs)`** — matches the driver and HODLR. The C ABI already slurps `Nloc*Nrhs` / `Ninloc*Ncol`.

---

## Files involved

- `SRC/C_BPACK_wrapper.cpp` — C entry points (`c_bpack_solve`, `c_bpack_mult`).
- `butterfly_h2/butterfly_integration.hpp` — `hierarchical_solve_parallel`, `hierarchical_mul_parallel`, `gather_local_solution`, `gather_solution_to_root`.
- `COMMUNICATION-AVOIDING-FMM/color/tree.hpp` — `SolveDataRequest`, `PendingSolveUpdates`.
- `COMMUNICATION-AVOIDING-FMM/color/solver.hpp` — per‑box solve helpers (`apply_forward_elimination`, `apply_backward_substitution`, `apply_diagonal_solve`, `gather_skeleton_to_parent`, `scatter_solution_to_children`, `accumulate_replace_then_add`, `merge_pending_solve`, parent↔child transfer serializers).
- `COMMUNICATION-AVOIDING-FMM/color/apply_mul.hpp` — per‑box multiply helpers (`apply_mul_forward_W`, `apply_diagonal_multiply`, `apply_mul_backward_V_with_pending`, `bunch_kaufman_multiply`).
- `COMMUNICATION-AVOIDING-FMM/color/serialization.hpp` — transport (`gather_boxes_solve`, `transport_and_apply_solve_updates_onehop`, `apply_full_update_local`, `apply_skel_update_local`, `bytes_pending`).
- `COMMUNICATION-AVOIDING-FMM/color/blas_declare.hpp` — BLAS prototypes.

---

## Phased implementation

### Phase 1 — Data structures (foundation) — `tree.hpp`
- `SolveDataRequest` (~329–479): add `int64_t num_points; int64_t nrhs;`. Change `initialize(morton, rank, num_points, nrhs)` → `right_side.resize(num_points*nrhs, 0); left_side.resize(num_points*nrhs, 0);`. Update `copy_minimal_from_box`, `copy_factorization_from_box`, `size()`, `is_initialized()` to preserve layout.
- `PendingSolveUpdates` (~94–97): per‑neighbor segment vectors become `seg_points × nrhs`.
- **Keep `nrhs=1` behavior identical first**, verify no regression before touching BLAS.

### Phase 2 — RHS load / solution extract (entry & exit)
In `butterfly_integration.hpp`:
- `hierarchical_solve_parallel` **signature**: add `int nrhs`.
- Box **init** (~1972): `solve_box.initialize(box.morton_index, rank, box.num_points, nrhs);`
- **Leaf RHS load** (~1981–1989) — the only new logic. `Nloc = rhs.size()/nrhs`; then
  `for c: solve_box.right_side[i + c*box.num_points] = rhs[global_idx + c*Nloc];`  (`left_side = right_side` still fine — whole‑vector copy).
- `hierarchical_mul_parallel` **signature** + **init** (~2672) + **input load** (~2676–2682): same column stride on `input_vec`.
- `gather_local_solution` (~1817–1837): emit column‑major output; change sanity check to `g == Nloc*nrhs`; add `nrhs`.

### Phase 3 — Per‑box BLAS (BLAS‑2 → BLAS‑3)

**Inside the two driver functions (the inlined diagonal blocks):**
- Solve inlined diagonal solve (`butterfly_integration.hpp` ~2260–2335): `b_R` → `r × nrhs` (gather 2D at 2262–2265, scatter 2D at 2332–2335); change `int nrhs = 1` → the box `nrhs` at 2269 / 2291 / 2316. LAPACK `dpotrs_`/`zsychol_solve_`/`dgetrs_`/`zgetrs_`/`zsytrs_` already take `nrhs` — **only the `1` changes** (keep `ldb = r`).
- Mul inlined diagonal multiply (`butterfly_integration.hpp` ~2832–2913): `x_R` → `r × nrhs` (gather 2833–2836, scatter 2910–2913); **`trmv → trmm`** at 2847/2851/2856/2860/2872/2876/2885/2889 (`trmm(side='L',uplo,trans,diag, m=r, n=nrhs, alpha=1, A, lda, x_R, ldb=r)`); `laswp` `nrhs=1 → nrhs` at 2880/2892; `bunch_kaufman_multiply` (2901) needs a **column loop** (`for c: bunch_kaufman_multiply(n, A, n, piv, &x_R[c*r])`).

**Helper bodies (call sites unchanged; bodies must change):**
- `solver.hpp` `apply_forward_elimination` (~906–1124): gemv→gemm at 974, 1000, 1043; scratch 2D.
- `solver.hpp` `apply_backward_substitution` (~1142–1478): gemv→gemm at 1224, 1245, 1338, 1427, 1462.
- `solver.hpp` `apply_diagonal_solve` (~2121–2233): `nrhs=1→nrhs` at 2162; `*trs` at 2171/2191/2213.
- `apply_mul.hpp` `apply_mul_forward_W` (~66): gemv→gemm (151/156, 190/195, 212/217).
- `apply_mul.hpp` `apply_diagonal_multiply` (~486): trmv→trmm (537/542/548/552/569/574/586/590); BK 601–607.
- `apply_mul.hpp` `apply_mul_backward_V_with_pending` (~638): gemv→gemm (713/716, 738/741, 830/835).
- `apply_mul.hpp` `bunch_kaufman_multiply` (~388): single‑vector kernel → column loop or batched version.

### Phase 4 — Transport / serialization (distributed correctness)
The serializers are **length‑prefixed**, so they carry wider payloads automatically — but the *appliers assume 1 value/DOF and will silently corrupt* otherwise:
- `serialization.hpp` `apply_full_update_local` (~3496–3511) & `apply_skel_update_local` (~3514–3531): fix the size checks and change `left_side[dof] += upd[i]` → `left_side[dof + c*num_points] += upd[i + c*seg]` per column.
- `solver.hpp` `accumulate_replace_then_add` (~844–868), `merge_pending_solve` (~870–890): per‑column.
- `serialization.hpp` `bytes_pending(PendingSolveUpdates)` (~240–260): × `nrhs`.
- `solver.hpp` parent↔child transfers: `gather_skeleton_to_parent` (~465–468), `scatter_solution_to_children` (~704–706, 831–833) — `parent_offset` stride must move `nrhs` columns of each skeleton row.

### Phase 5 — C wrapper (open the gate) — `SRC/C_BPACK_wrapper.cpp`
- `c_bpack_solve` (~853): remove the `Nrhs != 1` throw (~884–888); pass `*Nrhs` to `hierarchical_solve_parallel` and `gather_local_solution`. RHS already read as `Nloc*Nrhs` (~902).
- `c_bpack_mult` (~933): remove the `Ncol != 1` throw (~957–961); pass `*Ncol`. **Decide `trans`**: the `trans != 'N'` throw (~964–969) — matrix‑matrix transpose needs the apply run transposed; simplest is to keep only `'N'` initially and defer transpose.

### Phase 6 — Verification (test‑only)
- `gather_solution_to_root` (~3340–3519) + `h2_direct_verification`: extend `MPI_Gatherv` counts/displs (keep the index array ×1, scale the solution/rhs data arrays ×`nrhs`) and column‑major final scatter. Only needed to run the ground‑truth residual for `nrhs > 1`.

---

## BLAS declarations — `blas_declare.hpp`

Already present (no action): templated `gemm_` (~411–425) + `dgemm_`/`zgemm_`; templated `trsm_` (~393–406) + `dtrsm_`/`ztrsm_`; multi‑column LAPACK `dgetrs_`/`zgetrs_` (~210–213), `dpotrs_`/`zpotrs_` (~184–187).

**Must add** (for the mul `trmv→trmm`): `extern "C"` `dtrmm_`/`ztrmm_` prototypes, plus a templated `trmm_<DataType>` dispatcher mirroring `gemm_`/`trsm_`.

No build/link/CMake change — BLAS‑3 is the same OpenBLAS already linked.

---

## Key architectural note

The forward/backward sweeps in the two driver functions **delegate** per‑box work to the helpers. Because `nrhs` rides on `SolveDataRequest`, the helper **call sites** (e.g. `apply_forward_elimination(...)` at ~2162, `gather_skeleton_to_parent` ~2198, transport at ~2136/2189, `apply_diagonal_solve` ~2358) **do not change**. Only the **helper bodies** (Phase 3/4) and the **inlined diagonal blocks + input load + signatures** in the two driver functions change.

So the localized edit set *inside* the two driver functions is small:
| Function | Signature | Init | Input load | Inlined diagonal block |
|---|---|---|---|---|
| `hierarchical_solve_parallel` | +`nrhs` | ~1972 | ~1981–1989 | ~2260–2335 (`b_R` r×nrhs, `nrhs` on `*trs`) |
| `hierarchical_mul_parallel` | +`nrhs` | ~2672 | ~2676–2682 | ~2832–2913 (`x_R` r×nrhs, trmv→trmm, laswp nrhs, BK column loop) |

---

## Implementation order & testing

1. **Phase 1 → 2** with `nrhs` plumbed but keep `nrhs=1` correct; verify **no regression** before BLAS.
2. **Phase 3** (single‑rank matrix–matrix working).
3. **Phase 4** (distributed) — test with `-n 8` early; these are the highest‑risk edits.
4. **Phase 5** (expose via C wrapper).
5. **Phase 6** if you want the direct residual for `nrhs>1`.

**Helper:** add `inline int64_t idx(int64_t i, int64_t c, int64_t ld){ return i + c*ld; }` and use everywhere — prevents the most common indexing bug.

**Test ladder:**
1. `nrhs=1` reproduces current results **bit‑for‑bit** (regression gate).
2. `nrhs=3` matches three independent `nrhs=1` solves to tolerance.
3. Cross‑check against HODLR `z_c_bpack_solve` with `nvec=3` on the same operator.

---

## Top risks

1. **Silent serialization corruption** (Phase 4): length‑prefix carries bytes fine, so it compiles/runs, but the DOF‑indexed appliers scatter columns wrong. Test distributed early.
2. **Leading‑dimension mistakes**: skeleton/redundant sub‑blocks are row subsets → use gather‑to‑contiguous‑scratch (`ld = num_R`), never sub‑array striding.
3. **`trans='T'` for mult**: decide up front whether matrix–matrix multiply needs the transpose path; if not, keep the guard.

---

## Reference: how HODLR already does this (validation)

`SRC/BPACK_solve_mul.f90` `BPACK_Solution` (~471): `DT::x(Ns_loc, num_vectors), b(Ns_loc, num_vectors)` — RHS/solution are 2D `(local × num_vectors)` blocks. Direct path (~499) calls `BPACK_Inv_Mult('N', Ns_loc, num_vectors, b, x, ...)` — one batched pass, BLAS‑3 per block (this **is** Option B). Iterative path (~493–497) loops per column because Krylov is inherently per‑RHS. `blackbox_BPACK_MVP(trans, M, N, num_vect, Vin, Vout, ker)` carries `num_vect` for the matvec. Extending the H2 path replicates this proven pattern.
