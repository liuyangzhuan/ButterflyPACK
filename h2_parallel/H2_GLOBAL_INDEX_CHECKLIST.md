# H2 Global-Index Conversion — Change Checklist

Goal: make the FMM factorization (`hierarchical_factorization_parallel` and the
helpers it calls) evaluate kernel entries through ButterflyPACK's index-based
`C_FuncZmn` callback instead of the FMM's coordinate-based kernels, with
`num_proxy = 0`.

Status legend: `[x]` done · `[~]` in progress · `[ ]` not started

---

## Step 0 — Parent propagation of global indices  `[x]`

- [x] `build_parent_level_interactions` (factorization.hpp:~8210): push
      `child.point_indices[skel_idx]` (global) instead of `child.skeleton_indices`
      (local), kept parallel to `point_coords`.

---

## Step 1 — Store `quant` + index kernel adapter  `[~]`

- [x] Add `H2Kernel<CoordType, DataType>` adapter with `evaluate_block_by_index(...)`
      (butterfly_integration.hpp). Holds `kernel` fn ptr + `quant`.
- [x] `H2` stores the adapter by value: `H2Kernel<CoordType, DataType> kernel;`.
- [x] `c_bpack_construct_element_compute` (C_BPACK_wrapper.cpp, format-7): set
      `H2_solver->kernel.kernel = C_FuncZmn;` and `H2_solver->kernel.quant = C_QuantApp;`.
- [ ] (optional) default-init the fn ptr in `H2Kernel`: `... *kernel)(...) = nullptr;`.
- [x] Fix `hierarchical_factorization_parallel` signature (butterfly_integration.hpp:~599):
      kernel param must be `H2Kernel<CoordType, DataType>* kernel`, NOT a raw
      function pointer (helpers call `kernel->evaluate_block_by_index`, and the
      pointer must carry `quant`).
- [x] In `c_bpack_factor` (format-7): construct the factorizer, then call
      `hierarchical_factorization_parallel(H2_solver->tree.get(), &H2_solver->kernel, ...)`
      (pass address of the by-value member). Currently commented out.

---

## Step 2 — Global indices on assisting boxes (`PointDataRequest`)  `[ ]`

- [x] Add `std::vector<int64_t> point_indices;` to `PointDataRequest` (tree.hpp:~302).
      Full, parallel to `coords` (NOT skeleton-only).
- [x] Populate at both pack sites in `exchange_assisting_for_mortons_onehop`
      (serialization.hpp:~1830 and ~1845): `p.point_indices = b->point_indices;`.
- [x] Update the 3 serialization helpers in lockstep, same field order (insert
      after `skel_indices`, before `on_boundary`):
  - [x] `get_serialized_size(PointDataRequest)` (serialization.hpp:~911)
  - [x] `serialize(PointDataRequest, ...)` (serialization.hpp:~928)
  - [x] `deserialize(PointDataRequest, ...)` (serialization.hpp:~964)

---

## Step 3 — Convert kernel calls to index-based  `[ ]`

Drivers (`distributed_routine_*.cpp`) are intentionally NOT preserved.

### 3a. Add index-helper twins (mirror coord helpers)
- [ ] `indices_ptr_maybe_sliced` ↔ `coords_ptr_maybe_sliced` (factorization.hpp:~4126)
- [ ] `indices_for_count_from_box` / `indices_for_count_from_assist`
      ↔ `coords_for_count_from_*` (serialization.hpp:~1523/~1553)
- [ ] `gather_deferred_xnn_indices` ↔ `gather_deferred_xnn_coords` (factorization.hpp:~1661)
- [ ] `extract_skeleton_indices` ↔ `extract_skeleton_coords` (factorization.hpp:~7962)

### 3b. Rewrite `evaluate_block` → `evaluate_block_by_index` in the six helpers
(entry point in parentheses)
- [ ] `gather_id_workspace` (factorization.hpp:~3716)  — (#1 gather_id_workspace)
- [ ] `compute_step_two_internal` (factorization.hpp:~4280)  — (#2 compute_and_modify)
- [ ] `compute_deferred_xnn_base_matrix` (factorization.hpp:~1695)  — (#3 apply_owner_deferred_xnn…)
- [ ] `apply_updates_with_kernel_symmetric` (serialization.hpp:~2050)  — (#4 transport_and_apply…)
- [ ] `extract_child_interaction` (factorization.hpp:~7545)  — (#5 build_parent_level_interactions)
- [ ] `extract_or_evaluate_child_interaction_for_assisting` (factorization.hpp:~8010)  — (#5)

Conversion rule (coord source → index source):
- `box->point_coords` (full)          → `box->point_indices`
- skeleton slice `point_coords[skel]` → `point_indices[skel]`
- assisting full `req.coords`          → `req.point_indices`
- assisting skeleton                    → `req.point_indices[skel]`
- proxy                                 → none (guarded out)

---

## num_proxy guards  `[ ]`

- [ ] Force `num_proxy = 0` for format 7 (in `parse_program_options` / `ProgramOptions`).
- [ ] Throw if `num_proxy > 0` at `gather_id_workspace` entry (and/or
      `hierarchical_factorization_parallel` entry) — proxy blocks have no index equivalent.

---

## Validation  `[ ]`

- [ ] Single-rank `cvie3d`, tiny problem, `num_proxy=0` → leaf + parent, local
      boxes only (no assisting; Step 2 not yet exercised).
- [ ] Multi-rank → exercises ghost + assisting (depends on Step 2).
- [ ] Compare H2 (format 7) solve vs. classic (format ≠ 7) on the same matrix.

---

## Notes / gotchas

- After Step 3, `factorization.hpp` compiles only against `H2Kernel` (the coord
  drivers no longer instantiate).
- Coord temporaries that only fed `evaluate_block` (e.g. `skeleton_coords`) get
  replaced by index temporaries; leave `point_coords` itself populated.
- Serialization triplet (Step 2) must stay byte-identical in field order across
  all three helpers — mismatch → buffer overflow / stream corruption.
- `C_FuncZmn` uses 1-based `int` indices; the adapter converts from 0-based
  `int64_t` (`+1` and `static_cast<int>`).
- Kernel travels as the `H2Kernel` OBJECT (pointer) through the whole
  factorization so `evaluate_block_by_index` is callable and `quant` rides along.
  The bare function pointer only exists at the edge (`H2_solver->kernel.kernel`).
