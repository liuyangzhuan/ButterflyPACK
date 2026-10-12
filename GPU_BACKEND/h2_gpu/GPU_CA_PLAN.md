# GPU CA (communication-avoiding) factorization: plan and milestones

Status: M0 done; M1, M2 and M4b implemented and validated, M4b with
per-level solve residency and the solve's halo in device memory (results
under each). Keep this list current as milestones finish, so deferred
improvements are not lost.

Environment variables: the entries before 2026-10-01 use the names of the
time (`H2_GPU_SOLVE_KEEP`, `H2_CA_REPLICA_CHECK`, `H2_GPU_WAVE_TRACE`, ...);
[`doc/environment_variables.md`](../../doc/environment_variables.md) maps
them to the current ones (`BPACK_CHECK`, `BPACK_TRACE`, ...), and the
switches they name are gone (keeping solve factors is always on).

## Scope and decisions

- Structured (full-grid) points only; symmetric matrices, `H2_XRR_factor=1`.
- **CA0 first**: replicated CA (`H2_CA_owner_component=0`). Component-owner
  CA (CA3) is a later milestone (M5).
- **Host halo exchange first**: reuse the CPU level-entry gather
  (`gather_CA_factorization_data`, or the staged halo with every stage
  waited before elimination), then upload the ghosts. GPU-to-GPU exchange is
  M2.
- **Host solve on CA levels first**: the GPU solve keeps its "CA levels"
  fallback; CA-level factors go to the host as today. GPU solve on CA levels
  is M4.
- CA levels run lazy mode 1 semantics, as on the CPU (explicit near
  updates, lazy far regenerated from generators); a requested lazy 2 is
  treated as 1 there. Color levels keep the GPU Color path (lazy 2).
- Hybrid levels stay as on the CPU: CA on levels `>= CA_level` with at least
  343 boxes per active rank (3D), Color elsewhere.

## Background (why this helps on GPU)

At 64 ranks the GPU Color leaf level is communication bound: 8 parity waves,
each ending in a transport that waits for neighbors (Laplace 192: exchange
1.9 s of 3.8 s total). CA0 exchanges the halo once at level entry, then each
rank eliminates its boundary plus its neighbors' boundary copies (ghosts) in
the global CA order (blue/orange -> purple -> green, 8 parity waves each),
clears the ghosts, and runs the interior, with no communication inside the
level. The price is replicated work: 29% extra boxes for 16^3 bricks (VIE 384
at 64 ranks: 5282 eliminated per rank vs 4096 local), roughly 50-60% for 8^3
bricks. On CPU this made the CA leaf slower than Color (101.5 s vs 69.7 s);
on GPU compute is cheap, so the trade may flip. M0 measures this first.

Why no lazy 2 on CA levels: `generator_near` compresses Color's per-wave
transport of remote near updates into generators. CA0 has no such transport
(each rank applies its replicas' near updates itself), and with
`generator_near` on, the deferred near-update pass skips pairs with a
non-ghost remote endpoint, expecting generators CA never sends. See
`doc/h2_paper_notes/ca_algorithm.md` section 5.

## Milestones

### M0. Baselines and feasibility (runs only)
Laplace 192 (Nmin_leaf 216) at 8 ranks (16^3 leaf bricks) and 64 ranks (8^3),
plus a 16^3-brick case at 64 ranks (Laplace 192 with Nmin_leaf 27):
- CPU Color, CPU CA0 (`--CA_level 2 --H2_CA_staged_halo 2
  --H2_CA_owner_component 0`), GPU Color.
- Record per level: algorithm, replication factor (CA group sizes vs local
  boxes), halo gather time, elimination time, level time, exchange time,
  logdet and acc_mvp.
- Estimate GPU CA0 leaf time = GPU Color box work x replication + halo
  gather, and compare with the GPU Color leaf. Go/no-go for M1 per brick size.

**M0 results (2026-09-28, logs in build/gpu_prep/ca0).** Replication on CA
levels (rank 0, a process-grid corner; interior ranks replicate more):
1.29x for 16^3 bricks (5282 eliminated / 4096 local), 1.86x for 8^3
(954 / 512). CPU halo gathers 28-162 ms.

| case | brick | CPU Color | CPU CA0 | GPU Color level (box / comm / transition) |
|---|---|---|---|---|
| 8 ranks, level 5 | 16^3 | 55.2 s | 53.8 s | 6.73 s (3.5 / 1.4 / 0.86) |
| 8 ranks, level 4 | 8^3 | 30.8 s | 34.8 s | 2.13 s (0.7 / 0.92 / 0.09) |
| 64 ranks Nmin 216, level 5 | 8^3 | 11.4 s | 19.9 s | 1.71 s (0.6 / 0.77 / 0.06) |
| 64 ranks Nmin 27, level 6 | 16^3 | 2.9 s | 3.2 s | 1.60 s (0.6 / 0.33 / 0.43) |
| totals: 8 ranks / 64 ranks Nmin 216 | | 117.3 / 43.0 s | 164.7 / 50.1 s | 10.3 / 4.0 s |

Estimated GPU CA0 per CA level: 5.6 vs 6.7 s, 1.5 vs 2.1 s, 1.3 vs 1.7 s,
1.3 vs 1.6 s: about 15-30% off each CA level (10-15% of the total), an
optimistic ceiling at 64 ranks because rank 0 replicates least (CPU CA0's
64-rank leaf was 1.75x slower than Color).

Resolved: in the 8-rank CPU CA run the Color level 2 after the CA levels took
57.6 s vs 9.2 s. The time was in `H2Kernel::register_points`, called per box
by `register_level_coordinates` in every transport: its
`coordinate_cache.reserve(size() + n)` rehashed the whole cache (2.2M points
after the CA ghosts were added) whenever the bucket count changed, shrinks
included. Fixed to grow geometrically only when needed (butterfly_types.hpp):
level 2 57.96 -> 8.74 s, CPU CA0 total 164.7 -> 115.4 s (Color 115.7 s),
results unchanged.

### M1. GPU CA0 factorization on CA levels
Host halo gather, owner component 0, lazy 1 semantics, host solve on CA
levels.

Options: `H2_CA_owner_component=0`; `H2_CA_owner_serial` does not apply (it
is CA3's serial oracle); `H2_CA_staged_halo` 0 (monolithic gather) or 2
(staged gather, all four stages waited before elimination on the GPU; the
per-group overlap is M3). Box-by-box reference runs on the CPU use staged
halo 0: the CPU's staged overlap merges pending updates later, which sums in
a different order.

Files: no new folder. Device code goes into `h2_gpu/` (a CA mode of the
level eliminator, a new file for the CA glue: ghost upload, CA schedule, CA
sketch semantics; a CA entry in the shared GPU driver); `core/` and
`factorize_CA_level` stay the CPU reference with `#ifdef H2_HAVE_GPU` hooks,
as for the unstructured port.
1. **Determinism (all GPU copies of a box bitwise identical).** Owner and
   ghost replicas sit in different batches on different ranks, so:
   - **Harness first:** run each batched kernel on identical box inputs in
     different batch compositions (alone, other companions, counts across
     the thresholds: <216 vs >=216 boxes, <=108 for the cooperative ID) and
     compare bitwise: sketch, device ID (both variants), MAGMA LU/TRSM/GEMM,
     tensor-core GEMM, gather/symmetrize/transpose kernels.
   - **Device ID:** thread count (1024 below 216 boxes, else 256) and the
     cooperative split (`coop_blocks_per_box(count, max_n)`) depend on the
     batch, but each column's norms and dot products are one warp's
     lane-stride-32 sums and the pivot is an exact max, so skeleton, R and T
     likely do not; the block-wide sum behind the traced sketch norm does
     (fix it so traces compare).
   - **MAGMA `getrf_vbatched`:** its blocking may follow the batch's largest
     matrix; if the harness confirms, factor X_RR with a fixed-blocking
     batched LU of our own (X_RR is at most a few hundred wide) or group by
     exact size.
   - **GEMMs:** MAGMA and the tensor-core kernel sum each output's k products
     in order, and the tensor-core tile is chosen per matrix: expected safe.
   - **Canonical update order (main risk):** a pair block can receive
     updates from several boxes of one wave; apply them in a fixed order
     (waves in CA order, sources by Morton within a wave), not in the host's
     conflict-free batch rounds, which depend on the rank's other targets.
     Lazy-far fill sources in the sketch already follow `one_hop` order.
   - **Runtime check (opt-in):** each rank hashes the skeleton, T and X_RR
     factors of every box it eliminates; ghost copies are compared with the
     owner's; a mismatch names the box and wave.
   - Scope: bitwise equality across GPU ranks only; CPU CA0 stays a
     box-by-box reference (ranks and sketch norms to rounding).
   - **Harness results** (`tests/batch_determinism.cpp`, target
     `h2_gpu_batch_determinism`; 400 boxes, 6 groupings each against the
     box alone, double and complex):
     - device ID, adaptive launch: pivots and rank identical; R and T
       identical for double but **not complex** (up to 2.8e-12); traced
       norm differs (~4e-16). Fixed: `launch_qrcp(..., batch_independent)`
       (one block of 256 threads per box) is identical in every grouping,
       both types, traced norm included. Use it on CA levels.
     - MAGMA `getrf_vbatched`: **not** batch-independent (up to 1.4e-13,
       most boxes). Launching per power-of-two size class with
       `max_n` = the class size is identical in every grouping, both types.
       To do: the eliminator's CA mode launches the X_RR LU that way.
     - TRSM (right, upper), GEMM NN/TN/NT (MAGMA and tensor cores):
       identical in every grouping, both types.
     - Not covered by the harness: the ordered sketch (its list setup);
       covered by the runtime replica check.
2. **Ghosts as writable device boxes.** After the host gather, upload each
   ghost like a local box: point slots, Schur block, entry near blocks
   (reciprocal views restored by the gather's symmetric pruning), and index
   lists. Assisting boxes (outside the band) stay remote records (points
   only).
3. **CA schedule on the device.** Groups blue/orange, purple, green (with
   ghosts), then interior (local only), each in the 8 Morton-parity waves,
   with the CPU's wave numbers (`group * 8 + wave`) in `elimination_wave`.
4. **CA sketch semantics.** Mirror `gather_id_target_streamed` /
   `gather_id_workspace(use_CA_boundary_semantics=true)`: boundary boxes keep
   full rows of eliminated boundary neighbors outside their far map;
   assisting ring boxes always give full rows. Lazy far fills come from the
   device-resident fill sources (as in Color).
5. **Updates without transport.** No generator exchange and no owner pass:
   explicit near updates for every pair of writable boxes (local and ghost),
   plus the local-side updates of pairs with an assisting endpoint (the
   CPU's deferred XNN "mirror" targets, `include_ghosts=true`).
6. **Clear ghosts** on the device before the interior group; interior and
   transition as in Color (local boxes only; device transition to the next
   level when it runs on the GPU).
7. **Factors to the host** for the host solve (already done by the copier).
8. **Validation:** ID traces box by box against CPU CA0 (1 and 8 ranks), an
   opt-in replica check (each ghost copy's skeleton and T hash against its
   owner's), logdet and acc_mvp against CPU CA0, the host solve, and timing
   against GPU Color at 8 and 64 ranks.

**M1 design (per CA level on the GPU).**
0. Host, unchanged: halo gather (staged halo 0, or 2 with every stage
   waited) fills `level.ghost_boxes` with the ghosts' entry state.
1. Begin: device states for local and ghost boxes (ghosts are writable remote
   states: slots, Schur block, entry near blocks, index lists); one stored
   block per pair with a writable endpoint; assisting boxes: points only.
2. Schedule: the CPU's groups (blue/orange, purple, green, interior) x 8
   parity waves, same wave numbers. Per wave:
   - sketch/ID: the streamed-target rules of `gather_id_target_streamed`
     (already ghost-aware): ghosts resolve like local boxes; assisting boxes
     give full rows (never eliminated here); lazy-far fill sources from
     eliminated writable neighbors (resident); batch-independent ID.
   - elimination as in Color, X_RR LU by size class.
   - pair updates: candidates = writable boxes adjacent to the wave; a pair
     of writable boxes computed once, by its lower-Morton endpoint, into its
     single block; a pair with an assisting endpoint as this rank's copy;
     pairs of two eliminated ghosts and far pairs skipped; per-target order
     canonical (sub-batch j = j-th contribution, sources in one_hop order).
   - no transport, generators or device exchange.
3. Before the interior group: free the ghosts' heavy data (Schur, ghost-ghost
   pairs, factor blocks); keep their skeletons and slots (interior sketches
   read ghosts as ring rows) and the local-ghost pairs (transition).
4. Level end: the host's post-elimination assisting gather
   (`gather_CA_assisting_boxes_factorization`), its skeletons registered on
   the device, then the device transition as in Color.
5. Local factors to the host (copier) for the host solve; the next Color
   level adopts the parent blocks on the device.
6. Checks: box traces vs CPU CA0 (staged halo 0), replica hash check,
   logdet/acc_mvp, timing vs GPU Color.

Files (as built): the CA mode of `LevelEliminator` (`level_eliminator.hpp`:
ghost states, writable-box resolution, batch-independent ID and per-class
LU/solves, canonical pair updates and transition fills, no transport); the
CA entry of the shared driver (`factorization_driver.hpp`:
`ca_level_runs`, `start_ca_level`, `refresh_remote_skeletons`,
`check_ca_replicas`); the CA schedule `factorize_CA_level_gpu` next to
`factorize_CA_level` in `butterfly_factorization.hpp`, with the dispatch in
the level loop. `factorize_CA_level` itself is unchanged.

Differences/risks: the CPU computes a pair at its CA-ownership endpoint and
mirrors it, the GPU at its lower-Morton endpoint into one block: results
match the CPU to rounding (identical across GPU ranks). Ghosts add 30-90%
device memory on CA levels (Laplace 192 at 8 ranks already uses 26 GB per
GPU at the leaf): first tests on cases with headroom (64 ranks, Nmin 27).

**M1 results (2026-09-29, logs in build/gpu_prep/ca0/m1 and m1_64).**
Enabled by `--H2_use_gpu` on CA levels with owner component 0 (`GPU CA
level: on`); `H2_CA_REPLICA_CHECK=1` hashes every eliminated box's factors
(skeleton, T, LU and pivots, X_SR, X_RS, X_NR, X_RR_full) and compares each
ghost copy with its owner's, bitwise.

Correctness:
- Laplace 96, 1 and 8 ranks: ID traces identical to CPU CA0 (staged halo
  0) box by box, ghosts included (8208 boxes at 8 ranks: k, ring rows and
  tags, row counts, sketch norms to the 10 printed digits); logdet equal to
  16-17 digits, acc_mvp / acc_forward identical.
- VIE (`--scaleGreen 1`, complex symmetric) n=96, 8 ranks: logdet and acc
  identical to CPU CA0; replica check 0 of 3536 differ.
- Replica check, 0 differ: Laplace 96 at 8 ranks (3536 copies), Laplace 192
  at 8 ranks (9488), 64 ranks Nmin 27 on both CA levels (124704 and 49824).

What the replica check caught (all fixed):
1. MAGMA `trsm_vbatched` right/lower/unit (the second solve for X_RR^{-1})
   is not batch-independent (up to 2e-13; right/upper is). CA levels launch
   both solves per size class with max_m = max_n = the class (harness:
   identical in every grouping).
2. Owner pass: a pair whose lower endpoint is not writable here was
   computed as X_NR_E[hi] temp2_E[lo]^T instead of temp2_E[hi] X_NR_E[lo]^T
   (equal only through the symmetry of X_RR^{-1}); a target's contributions
   came in the candidate's one_hop order, which differs by the endpoint that
   reaches it. Now always the canonical product, and on CA levels by source
   Morton.
3. Device transition two-hop fills: each rank took its local child as i
   (G_i^T X G_j with swapped roles on the two ranks). On CA levels i is now
   the lower child and the sources come by Morton, so a cross-rank parent
   pair's two copies agree bitwise (the next CA level's ghosts read them).
4. The minimum LU size class is 32: classes 8/16 loaded MAGMA kernels of
   their own at first use (a one-time 0.35 s inside the level).

Other fixes on the way: the CA solve applies ghosts again (their factors
must come to the host; see B3); `download_level_blocks` streams through two
pinned ring slots instead of one synchronous copy per block (1-rank CA leaf
transition 3.15 s before); after the last CA level the parent blocks stay on
the device (the next Color level's box path was predicted with the CA
level's capped lazy mode).

Timing (Laplace 192, factorization; GPU Color vs GPU CA0):

| case | CA levels | GPU Color | GPU CA0 |
|---|---|---|---|
| 8 ranks, 16^3 leaf | leaf only (`CA_level 5`) | 9.2 s (leaf 5.97) | 8.7 s (leaf 5.57; about 5.2 without artifacts, see below) |
| 8 ranks, 16^3 leaf | leaf + level 4 | 9.2 s | host OOM in level 4's gather (4 ranks/node) |
| 64 ranks, Nmin 216 (8^3 leaf) | leaf | 4.00 s (leaf 1.79) | 4.25 s (leaf 2.13) |
| 64 ranks, Nmin 27 (16^3 leaf) | leaf + level 5 | 4.46 s (1.61 + 0.73) | 9.5 s staged halo 2 (2.13 + 5.32); 11.9 s staged halo 0 |

Why CA0 does not win yet:
- Replication: 1.29x boxes for 16^3 bricks, 1.86x for 8^3 (owner GEMMs
  715K vs 574K, 110K vs 63K), while the Color exchange it removes is only
  0.4-0.9 s per leaf at 64 ranks after the CUDA-aware exchange work.
- Twice the waves: CA0 runs 4 groups x 8 parities = 32 waves against
  Color's 16, so per-wave latency (sketch, plan, synchronization) doubles
  (64 ranks Nmin 216 leaf: sketch 0.58 vs 0.29 s).
- The host round trip between consecutive CA levels: parent blocks to the
  host, host halo gather (1.6 s staged, 4.4 s monolithic at 64 ranks level
  5), re-upload (1.5 s for 2.8 GB), plus rank skew; the level's device work
  is 0.7 s. At 8 ranks the host copies of both levels exceed the node's
  host memory with 4 ranks per node.
- CA levels also send the whole solve to the host (M4): 2.2 s vs 0.26 s at
  8 ranks.

Leaf level in detail (rank 0 unless noted; rank 0 is a process-grid
corner, so it replicates least):

| | 8 ranks, 16^3 | 64 ranks, 8^3 (Nmin 216) | 64 ranks, 16^3 (Nmin 27) |
|---|---|---|---|
| boxes/rank: local, CA rank 0, CA average | 4096, 5282, 5282 | 512, 954, ~1290 | 4096, 5282, ~6044 |
| boxes per wave: Color / CA | 169-343 / 54-343 | 27-37 / 18-54 | 169-343 / 54-343 (27-point boxes) |
| elimination: Color / CA0 | 4.99 / 4.12 s | 1.69 / 1.84-1.90 s | 1.08 / 1.20 s |
| Color exchange (of which waiting on sizes) | 1.86 (0.93) s | 0.86 (0.47) s | 0.39 (0.23) s |
| level total: Color / CA0 | 5.97 / ~5.2 s | 1.79 / 2.13 s | 1.61 / 2.13 s |

- Two artifacts inflated the 8-rank CA leaf in the table above: 0.33 s of
  first-use MAGMA loading (size classes under 32, since fixed) in the 5.57 s
  run, and the replica check in the 5.80 s run. Without them the
  elimination is 4.12 s (with the check, which runs after the elimination)
  against Color's 4.99 s: CA0 wins the work-bound leaf by about 12%. To be
  confirmed by a clean run. GPU Color's 64-rank Nmin 216 leaf has its own
  first-use cost (LU 0.27 s vs 0.03 s elsewhere).
- Work-bound (8 ranks, 169-343 boxes per wave): the replicated 1.29x work
  costs about +1 s (sketch +0.6, plan +0.4, device +0.55, overlapping),
  less than the 1.86 s exchange it removes.
- 64 ranks: see the wave traces below. (An earlier reading here, that the
  elimination hardly depends on the box count because all ranks' times
  agree within 3%, was wrong: `finish_level` ends in a collective, so every
  rank's elimination time includes waiting for the slowest one.)

**Leaf level, measured per wave (2026-09-29; logs in
build/gpu_prep/ca0/trace and trace2).** Tools: `H2_GPU_WAVE_TRACE=1` prints,
per device level, the spread of the ranks' wave totals, and the per-wave
host phases and device times of rank 0 and of the slowest rank;
`H2_GPU_WARMUP=1` loads the batched kernels before the levels (first-use
loading had cost up to 0.3 s inside a leaf, in Color as well).

| case (elimination, kernels warm) | GPU Color | GPU CA0 |
|---|---|---|
| 8 ranks, 16^3 bricks | 5.00-5.06 s | 4.07-4.12 s |
| 64 ranks, Nmin 216 | 1.25-1.58 s | 1.78-1.86 s |
| 64 ranks, Nmin 27 | 1.03-1.08 s | 1.11-1.15 s |
| the same, copier at 8 threads (now the default at 4 ranks/node): 8 ranks / 64 ranks Nmin 216 | 4.83 / 1.23-1.25 s | 3.79 / 1.57-1.63 s |

1. Load imbalance decides the CA0 leaf at 64 ranks. Ranks eliminate 954
   (rank 0, a process-grid corner) to 1720 boxes (most-ghosted interior
   rank) against 512 local at Nmin 216, and 5282 to 6936 against 4096 at
   Nmin 27; their waves take 836 to 1369 ms, and the level ends with all
   ranks waiting for the slowest. Rank 0 alone would beat Color (waves 0.87
   s against Color's 0.37 s of waves plus 0.65 s of transports). The
   slowest rank is mostly GPU-bound (device busy ~80% of its waves).
2. Host oversubscription: the background copier's 16 threads plus 16
   OpenMP threads (plus the MPI progress thread) on each rank's 32 logical
   CPUs stalled the waves' OpenMP loops (ID store up to 13 ms per wave, 0.1-0.2
   s per level). Fixed: the copier defaults to half the OpenMP threads (at
   most 16; `H2_GPU_COPY_THREADS` overrides): elimination -8% (CA0 8 ranks),
   -12% (CA0 64 ranks), -3% to -9% (Color).
3. Device ID: every launch variant (256 or 1024 threads, cooperative) now
   gives the same bits (compiled with `--fmad=false`, traced norm summed in
   column order; harness: identical in every grouping, both types), so CA
   levels use the adaptive launch: device sketch -11% on the median rank at
   64 ranks Nmin 216, about -25 ms of elimination; replica check still 0
   differ on both CA levels.
4. Per-wave structure. Work-bound (8 ranks): the GPU is busy ~78% of the
   waves; the host path it waits on (ID store, box plan, launch) is ~0.8 of
   3.8 s. The device sketch's list kernels (random draws and counting sorts,
   one warp per box) cost ~3.5 ms per wave on the stream's critical path.
   On a quarter node the serial planning code runs ~2x slower than on a full
   node (not explained by OpenMP spinning, the MPI progress thread or the
   copier alone).

Candidate next steps for the CA0 leaf (not started):
- L1. Sketch lists on a side stream, overlapping the previous wave's
  elimination (~0.1 s per leaf): needs its own buffers (the heap's reuse is
  safe only in one stream) and its own metadata image.
- L2. The box plan off the critical path: its parts that do not depend on
  the ID (neighbor row layouts, block lookups) while the host waits for the
  sketch, and the rest in parallel over boxes.
- L3. The replication imbalance itself: inherent to CA0 (every neighbor's
  boundary is replicated); CA3 (M5) removes most of it.

Where M2-M4 would leave this: M2 removes the round trip (the second CA
level would cost roughly its device work, about Color's time), M4 the host
solve; the leaf itself stays at parity or slower at 64 ranks, from
replication and the doubled waves. A clear win needs cases where the Color
exchange dominates more than here.

### M2. GPU-to-GPU halo exchange
Send the ghosts' entry blocks (Schur, near blocks; one orientation with
symmetric pruning) straight between device buffers with CUDA-aware MPI,
instead of host serialize, send, deserialize and upload. Keep the staged
priorities (shells, blue/orange, purple, green).

Plan (2026-09-29). Only CA levels after the first gain: the leaf's entries
come from the kernel (its halo gather moves no blocks, 20-80 ms). At 64
ranks Nmin 27 the second CA level (level 5) now pays the leaf transition's
download to the host (~0.4 s), the host halo gather (1.6 s staged, 4.4 s
monolithic) and the re-upload in `begin()` (1.5 s for 2.8 GB), against
~0.7 s of device work; at 8 ranks the host copies exhaust host memory.
1. Keep the blocks on the device across a CA -> CA transition:
   `keeps_blocks_of(lvl)` also true when level `lvl` will run as a device CA
   level (the same collective predicate as `ca_level_runs`, plus no process
   reduction between the levels). If `start_ca_level` then declines the level
   (a box too large for the device sketch), the blocks go to the host there
   and the host CA path runs as today.
2. Host halo gather unchanged: the level's host boxes then carry no blocks,
   so `gather_CA_factorization_data` moves only metadata (points, index and
   neighbor lists, assisting boxes); its symmetric pruning and reconstruction
   find nothing to do. (Monolithic gather on device CA levels; the staged
   one's overlap is M3.)
3. Device ghost-block exchange (new, `h2_gpu/ca_halo.hpp`), after the
   host gather and before the eliminator: for every ghost G a peer requests
   from this rank (the gather's `CARequestPlan`), Schur(G) and each pair
   block (G, X), X one hop from G, unless the host rule would prune it
   (G > X and X is local on the receiver or also requested by it; the
   receiver then already has the pair: its own copy, bitwise equal since the
   M1 transition fix, or X's payload). Blocks keep their device layout (rows
   = the higher box). Per peer: one device pack (gather kernel), one
   `MPI_Isend`/`MPI_Irecv` of device memory, one unpack into heap blocks.
   Both sides derive the list from shared data (request sets, the tree's
   one-hop lists, box sizes), so no block headers travel.
4. `begin()` in CA mode adopts the transition's blocks plus the received
   ghost blocks (today CA mode rejects adopted blocks and packs everything
   on the host).
5. Checks: the replica check on both CA levels (bitwise; the same numbers
   as M1's host path move, so logdet and acc stay identical to M1), 64 ranks
   Nmin 27 and 8 ranks with two CA levels (the host-memory failure), timing
   of the second CA level.

**M2 as built (2026-09-29).** `h2_gpu/ca_halo.hpp`, the transition's
keep rule and the CA eliminator's adoption, as planned, with these changes:
- The send rule is tighter than the host's: a pair (G, X) is not sent when
  X is the receiver's own box (either order), nor when X is also a ghost of
  the receiver and X < G, so every block crosses at most once.
- Messages go through the exchange arena (created now also by multi-rank
  CA levels): buffers of the main heap made CUDA-aware MPI ~30x slower (3.5
  s for 1.5 GB). They move in rounds of one arena chunk per message, packed
  from the blocks and unpacked into the ghosts' own blocks, with no
  message-sized buffers (a fragmented heap could not hold them: Laplace 192
  at 8 ranks); the chunk size is agreed across ranks (per-rank sizes
  truncated messages at 64 ranks).
- `H2_GPU_CA_DEVICE_HALO=0` keeps the host path (M1).

Results (8 ranks): Laplace 96 Nmin 27 (CA levels 5 and 4): level 4 3.22 ->
0.75 s (halo 0.12 s for 1.5 GB), factorization 5.52 -> 2.81 s, replica
check 0 differ on both levels, accuracy identical to the host path. Laplace
192 (CA levels 5 and 4, host memory exhausted before): runs, halo 0.38 s for
6.6 GB, level 4 2.2 s, factorization 8.7 s, logdet 5e-16 from CPU CA0 and
the same acc_mvp / acc_forward. VIE (complex) Nmin 27: replica check 0
differ on both levels.

### M3. Overlap the staged halo with elimination on the GPU
The CPU staged mode 2 overlap: wait only for the stage the next group needs,
divert writes into unarrived ghosts to pending state, merge at arrival.

### M4. GPU solve on CA levels
Remove the device solve's "CA levels" fallback: replicate the ghost solves
as `gather_CA_boxes_solve` does (or an owner-based variant), with the
level's factors kept on the device.

Plan (2026-09-29). Today one CA level sends the whole solve and the
factored multiply to the host (8 ranks Laplace 192: 2.2 s vs 0.26 s).
- M4a, owner-based (recommended first): the device solve runs a CA level
  like a Color level but with the CA schedule: waves = the CA groups x 8
  parities (32 on several ranks, 16 on one), the level's local boxes only,
  each wave's updates to other ranks' boxes sent after it and neighbor
  vectors refreshed before each backward wave (the device solve's existing
  tables and transports, built from CA waves instead of Color waves). The
  same elimination order as the factorization, so the same solution (ghost
  factors equal their owners' bitwise); only local factors are needed, kept
  on the device during the factorization as on Color levels. More, smaller
  messages than the host CA solve (32 per sweep vs 1 forward + 3 backward
  gathers), each a vector piece.
- M4b, replicated (the host CA solve's scheme): ghosts' factors on the
  device too; the ghost and assisting vectors gathered once before the
  forward sweep and before each backward group; the CA schedule over local
  and ghost boxes, updates only to local or ghost boxes. Fewer messages,
  redundant ghost work, more code (tables for ghosts, gathers through the
  host or the device).
- Checks: acc_forward / acc_backward against the host CA solve, the
  multiply too, and solve times against GPU Color.

**M4b as built (2026-09-29).** The user chose M4b (CA's advantage over
Color is less communication; an owner-based solve would give it up).
`h2_gpu/device_solve.hpp`: a CA level's `DeviceSolveLevel` holds the
local boxes, then the ghosts (factors uploaded from `level.ghost_boxes`),
and the assisting boxes in the read-only area (their sizes and skeletons
recorded by the factorization after the post-elimination assisting
gather, since the level clears its copy); `build_ca_solve_tables` builds the
CA waves (the factorization's groups x parities) with updates to local or
ghost boxes only; `ca_forward` / `ca_backward` / `ca_mul_forward` /
`ca_mul_backward` follow the host CA solve and multiply, with the host's
`gather_CA_boxes_solve` between groups (forward: once; backward: after each
group but the last, the ghosts the first time). A CA level qualifies when
the device factored it as replicated CA (`device_solve_store().ca_levels`).

Results: `H2_GPU_SOLVE_CHECK=1` (device against host, relative): Laplace 96
one rank solve 3.6e-15, multiply 4.6e-16; 8 ranks Nmin 27 (two CA levels)
solve 4.2e-15, multiply 4.1e-16; VIE complex solve 6.6e-15, multiply
7.0e-16; 64 ranks (Laplace and VIE, Nmin 27) solve 5e-15 to 3.4e-14,
multiply 4e-16 to 2.3e-15; accuracy identical. Solve on the 8-rank Nmin 27
case: 0.13 s.

**After M2 + M4b: GPU CA0 against GPU Color (2026-09-29, Laplace 192,
kernels warm, build/gpu_prep/ca0/m24_64).**

| case | factorization CA0 / Color | solve CA0 / Color |
|---|---|---|
| 8 ranks (CA levels 5, 4) | 8.84 / 9.57 s | 2.44 s (host: device memory) / 0.26 s |
| 64 ranks, Nmin 27 (CA levels 6, 5) | 5.39-5.69 / 4.38-4.44 s | 0.30 / 0.19 s |
| 64 ranks, Nmin 216 (CA level 5) | 4.06-4.10 / 3.52-3.55 s | 0.27 / 0.20 s |

The second CA level's halo now costs 0.35-0.38 s at 64 ranks (was 1.6 s
staged + 1.5 s re-upload + the download). Two open issues followed, both
fixed below: at 8 ranks the device solve did not fit (its decision was for
all levels at once), and the CA solve's gathers went through the host.

**M4 follow-up 1: per-level residency (2026-09-29).** The CA levels' ghost
factors (at 8 ranks Laplace 192: level 5 local 19.8 GB + ghosts 5.7 GB,
level 4 7.2 + 6.2 GB) exceed the ~33 GB heap; keeping them during the
factorization would not reduce the total. `prepare_device_solve` now
decides per level, the same on every rank (each level's largest need and
the smallest room, by Allreduce): every Color level on the device (else
the whole device solve is off, as before), then the CA levels leaf first
while they fit; a CA level that does not fit (or that the device did not
factor as replicated CA) runs its steps on the host through the callback
`run_device_solve(..., host_level)`. `butterfly_solve.hpp`: the host CA
level steps are functions (`CA_level_solve_forward/backward`,
`CA_level_mul_forward/backward`, returning their gather time) and the host
diagonal phases are per-level functions (`diagonal_solve_level`,
`diagonal_mul_level`), used by the host solve and by the host-level
callback (the device sweeps apply a level's diagonal after its forward
sweep, so the callback does too). The setup prints the decision
("GPU solve: CA levels 4 on the host (4: 11.98 GB; device room 31.61 GB
per rank, resident 27.45 GB)").

Results, 8 ranks Laplace 192 (build/gpu_prep/ca0/fix1): device (levels 5,
3, 2) + host (level 4) against the host solve: solve 4.0e-15 / 2.7e-14,
multiply 5.8e-16 / 8.6e-16; acc_forward 4.3446769e-03, bitwise the
host-only run's. Solve 0.94-0.96 s (was 2.44 s on the host; Color 0.26 s),
of which level 4 on the host 0.70-0.72 s. The CPU build (grid 96, CA, 8
ranks) is unchanged.

**M4 follow-up 2: the CA solve's halo in device memory (2026-09-29).**
`build_ca_halo` (device_solve.hpp, at the solve setup, one Alltoall +
Alltoallv of requests per CA level over the tree's ranks): each rank asks
the owners of its ghosts and of the assisting boxes it reads (Morton,
ghost or assisting, points); copy tables for three kinds (0 the ghosts,
1 the ghosts and assisting boxes, 2 the assisting boxes): pack from the
owner's level vector, unpack into the ghost part of the level vector or
the read-only area. `DeviceSolveRun::ca_halo` packs, exchanges through
the MPI exchange arena (CUDA-aware, the Color solve's `exchange`), and
unpacks; no host download, gather or re-upload. The same exchanges as
before: forward the ghosts only (the host gather also moved the assisting
boxes, unused there); backward after each group but the last, kind 1 the
first time, then kind 2; multiply forward kind 1; multiply backward kind
1 at the start and after each group but the last. `H2_GPU_CA_SOLVE_HALO=0`
restores the host gathers. 8 ranks (level 5 on the device): MPI 0.063 s
(was 0.089), vector transfers 0.004 s (was 0.032), pack/unpack 0.024 s.

Results at 64 ranks (Laplace 192, 16 nodes, build/gpu_prep/ca0/fix2, two
repetitions agreeing to 0.01 s): device against host (`H2_GPU_SOLVE_CHECK`,
Nmin 27, both CA levels on the device) solve 5.1e-15 / 3.4e-14, multiply
6.0e-16 / 8.1e-16. Times (the device solve's own timer; one solve, nrhs 1):

| case | solve: CA0 device halo / CA0 host gathers / Color | multiply: same order |
|---|---|---|
| Nmin 27 (CA levels 6, 5) | 0.207 / 0.297 / 0.194 s | 0.251 / 0.368 / 0.235 s |
| Nmin 216 (CA level 5) | 0.221 / 0.270 / 0.199 s | 0.267 / 0.335 / 0.243 s |

Rank 0's MPI time is now the same as Color's (0.103 vs 0.102 s, was 0.19);
CA0 remains 7-11% slower than Color in the solve (its ghost boxes are
eliminated again in each sweep). Factorization in the same job: CA0 5.38-5.39
/ Color 4.38-4.41 s (Nmin 27), 4.08 / 3.54-3.55 s (Nmin 216), as before.

**Level-by-level: GPU CA0 against GPU Color (2026-09-29, Laplace 192,
factorization, rank 0's level times; build/gpu_prep/ca0/{fix1,fix2}, Color
8 ranks from m24_64).** Color comm = the summed per-color exchange times
(waiting included, as the preprint's t_comm); CA0 extra = CA0's elimination
(without its device halo) minus Color's elimination without comm, i.e. the
ghost boxes eliminated again; CA0 comm = device halo + host gathers.

| run | level (boxes/rank) | Color: total (comm) | CA0: total | CA0 extra | CA0 comm |
|---|---|---|---|---|---|
| 8 ranks, Nmin 216 | 5 leaf (4096) | 5.88 (1.43) | 5.02 | +0.44 | 0.05 |
| | 4 (512) | 2.02 (0.81) | 2.21 | +0.51 | 0.40 |
| | 3 (64), 2 (8) | 0.83 (0.35), 0.35 (0.16) | Color | | |
| | factorization | 9.57 | 8.71 | | |
| 64 ranks, Nmin 216 | 5 leaf (512) | 1.37 (0.68) | 1.86 | +1.01 | 0.11 |
| | 4 (64), 3 (8), 2 (8) | 0.65 (0.37), 0.70 (0.37), 0.58 (0.31) | Color | | |
| | factorization | 3.54 | 4.08 | | |
| 64 ranks, Nmin 27 | 6 leaf (4096) | 1.59 (0.35) | 1.73 | +0.31 | 0.09 |
| | 5 (512) | 0.68 (0.32) | 1.51 | +0.67 | 0.45 |
| | 4 (64), 3 (8), 2 (8) | 0.61 (0.35), 0.71 (0.36), 0.59 (0.30) | Color | | |
| | factorization | 4.41 | 5.39 | | |

Reading: CA0 pays for its ghosts in elimination work and saves Color's
per-color exchanges. On 16^3 bricks (4096 boxes/rank) the ghosts are
+13% of Color's work at 8 ranks (2x2x2: each brick shares 3 faces) and
+43% at 64 ranks (4x4x4: interior bricks share 6; the halo sizes of ranks
1 and 35 differ by 2x); with 216-point boxes at 8 ranks CA0 wins the leaf
by 0.86 s. On 8^3 bricks (512 boxes/rank) the ghosts cost 1.5-3.5x Color's
work, plus a device halo on non-leaf levels, more than Color's exchanges
(0.3-0.8 s): CA0 loses there. Levels with 64 or 8 boxes per rank (below
the 343 of `minimum_CA_boxes_per_active_process`) are Color in both, and
at 64 ranks their exchanges (about 1 s) are most of Color's communication.
Restricting CA to the leaf (`--CA_level` = the leaf level) would give 8
ranks' levels 5-4 about 7.0 s (CA0 leaf 5.02 + Color 2.02; not measured),
against 7.9 s for Color and 7.2 s for CA0 on both.

**Laplace 384^3 (tol 1e-3, Nmin 216) on GPUs (2026-09-29,
build/gpu_prep/ca0/g384, g384_512).** 64 ranks (16 nodes, 884K unknowns
per rank, the preprint's P = 64 run): Color and CA0 (leaf only, and leaf
+ level 5) all stopped at the leaf with the device heap exhausted on the
ranks whose bricks share faces with neighbours on every side (in use 30.2
of 33.4 GB, largest free block ~5 MB): at 8 ranks the same per-rank size
peaked at the heap's capacity, and 6 shared faces (vs 3) double the remote
data. The 80 GB A100 nodes (`-C gpu&hbm80g`, heap ~67 GB) would fit;
cancelled after a 4.5 h queue estimate. Rank counts must be 8^k.

512 ranks (128 nodes, 110K unknowns per rank; job 53 s): per level, rank
0, seconds (Color comm as above):

| level (ranks, boxes/rank) | Color total (comm) | CA0 total |
|---|---|---|
| 6 leaf (512, 512) | 1.80 (0.91) | 2.18 CA (of which the initial host gather 0.25) |
| 5 (512, 64) | 0.71 (0.42) | 0.82 (Color) |
| 4 (512, 8; then 512 -> 64 ranks) | 1.02 (0.41) | 1.02 (Color) |
| 3 (64, 8) | 1.56 (1.04) | 1.46 (Color) |
| 2 (8, 8) | 0.89 (0.49) | 0.85 (Color) |
| factorization | 7.03, of which kernel warm-up 0.79 | 6.59, warm-up 0.07 |

Color ran first on fresh nodes: its warm-up (in the factor time) took
0.79 s against 0.07 s. Without it Color is ~6.2 s and CA0 ~6.5 s; the same
Color levels differ by up to 0.1 s between the runs (noise). Color's
exchanges are 3.27 s of its 6.0 s of level time (54%; per level 51, 59,
40, 67, 55%), mostly the size handshake of each round (where a rank waits
for its neighbours; e.g. leaf 0.51 of 0.99 s) at ~1.7 GB/s received per
rank. Weak scaling from 192^3 on 64 ranks (same 110K per rank): 3.54 ->
~6.2 s (+76%; the preprint's CPU runs +38% from 384^3/64 to 768^3/512):
the leaf +0.43 s and one more distributed level (level 3 on 64 ranks,
comm 1.04 s, all inter-node) plus the process reduction. CA0 at the leaf
(8^3 bricks) loses 0.38 s (2.18 vs 1.80 s; the level time includes the
initial gather, which the level's timer covers): its ghost work doubles
Color's (1.68 vs 0.88 s of elimination without comm), and its host-side
initial gather grows with the rank count (0.02 s at 8 ranks, 0.08 s at
64, 0.25 s at 512). Accuracy: acc_forward
6.85e-3 both, acc_mvp 5.0e-4 / 5.3e-4.

**C1a and G1 (2026-09-29, build/gpu_prep/ca0/c1g1).** Two of the CA0
leaf improvements proposed after the 512-rank run (the others: G2 the
initial gather's plan from geometry without global collectives, P1 the box
plan's ID-independent part ahead of the waves, P2 the sketch lists on a side
stream, C1b no host download for kept levels).
- G1: the post-elimination assisting gather moved the assisting boxes'
  points and indices again, over a new request plan (Allgather, Alltoall,
  Alltoallv, then an Alltoall of sizes, all over the level's ranks), though
  only their skeletons change. The level now keeps the plan of its assisting
  gather (`assisting_plan_*` in TreeLevel), and
  `refresh_CA_assisting_skeletons` (core/serialization.hpp) sends only
  the skeletons, point to point with the plan's peers (sizes, then indices);
  the full gather remains the fallback. CPU and GPU CA paths.
- C1a: CA levels keep their factors on the device during the factorization
  like Color levels (`keep_solve_` no longer excludes CA; the ghosts'
  factors too, keyed by Morton index; the solve setup counts kept ghosts as
  resident). The host download remains (as for Color). Within the keep
  budget (60% of the heap, `H2_GPU_SOLVE_KEEP_FRACTION`), decided by all
  ranks of the level.

Results (results unchanged bitwise: replica check 0 differ on every CA
level, solve and multiply checks as before, accuracy identical; the CPU
build's logdet identical):

| case | post-elimination gathers | solve setup |
|---|---|---|
| 8 ranks (levels 5, 4 CA; neither fits the keep budget) | 27, 19 -> 2, 1 ms | 2.85 s, unchanged |
| 64 ranks Nmin 27 (levels 6, 5 CA, both kept) | 15, 31 -> 1, 2 ms | 1.03 -> 0.16 s (upload 0.63 -> 0.07 s) |
| 64 ranks Nmin 216 (leaf CA, not kept: 9.84 GB of factors on the most-ghosted rank) | 39 -> 2 ms | 1.40 s, unchanged |

Factorization times unchanged (4.12 / 5.39 s). The logdet's last digits
vary from run to run, before these changes too (64 ranks Nmin 216: CA
...062607, ...622, ...637; Color ...534012, ...027, ...042): in
`hierarchical_logdet_parallel` the threads' partial sums are added with
`omp atomic`, in the order the threads finish; the factors do not vary.

**C1b, the keep budget, the copier (2026-09-29, build/gpu_prep/ca0/c1b,
c1b2).** G2 (the initial gather's plan from geometry) was set aside: most of
the gather at 512 ranks looked like start-up skew (Color's first exchange
0.18 s there too); the gather's line now breaks it down.
- Keep decision per level before its waves, alike on every rank
  (`decide_solve_keep` in factorization_driver.hpp): an upper bound of the
  factors' bytes (`solve_keep_bound`, device_solve.hpp: per box of n points
  and ntot neighbor points, (1.5 n^2 + ntot n) doubles) within
  `H2_GPU_SOLVE_KEEP_FRACTION` of the heap (default 0.75, was 0.6; 0.8 since
  2026-09-30), so a
  level that cannot keep them no longer starts to and gives up midway (the
  884K-per-rank leaves); a level kept by the bound then needs only the room.
- C1b: a kept box's X_NR (most of the download) gets no host copy when the
  host does not read it during the level (device exchange, a CA level or a
  single rank; not with the replica check). It is copied back from the kept
  factors (`rescue_skipped_xnr`, `materialize_host_factors`) before a fill
  source leaves the device (the transition restores from host copies; the
  level then downloads X_NR again), when the keep fails or another rank's
  did, before a host transition, and before any host solve or multiply
  (the fallback, H2_GPU_SOLVE_CHECK's reference run). The solve setup takes
  a kept box's X_NR rows from the kept factors; the transition's early
  reclaim does not count kept factors.
- HostCopier: a block was downloaded whole and its segments picked on the
  host; now only the ranges its segments read cross PCIe (ranges under 128
  KB apart merged), and "down" counts the bytes moved.

Results (64 ranks: 40 GB nodes; unchanged bitwise: replica check 0 differ,
solve and multiply checks as before with X_NR skipped, i.e. through the
host copies made on demand; accuracy identical):

| 64 ranks, Laplace 192 | before (fix2) | C1b, host copies skipped | + copier ranges |
|---|---|---|---|
| Nmin 216, CA0: factor / leaf / leaf down / solve setup | 4.08 / 1.86 s / 4.99 GB / 1.40 s | 3.70 / 1.43 / 4.99 / 0.14 | 3.68 / 1.42 / 0.64 / 0.14 |
| Nmin 216, Color: factor / leaf | 3.54 / 1.37 s | 3.46 / 1.31 | 3.28 / 1.24 |
| Nmin 27, CA0 / Color: factor | 5.39 / 4.41 s | 4.99 / 4.35 | 4.97 / 4.19 |

The CA0 leaf's 0.43 s came from its hosts: with X_NR no longer copied into
BoxData, the copier (busy 0.65 -> 0.37 s on rank 0, more on the
most-ghosted ranks) stops competing with the waves' host path (rank 0 plan
0.21 -> 0.11 s; slowest rank's elimination 1.61 -> 1.20 s). Color gains
from the copier ranges (its levels were kept already): leaf 1.31 -> 1.24 s,
its exchanges faster too (GPU-direct MPI shares PCIe with the downloads).
The initial gather at 64 ranks Nmin 216: 90 ms = plan 29 (first collective
24, i.e. waiting for ranks) + request sets 4 + payloads 12 + assisting 45.

Node memory matters: `-C gpu` may give 80 GB A100 nodes (one job here); pin
`-C gpu&hbm40g` for comparisons. 8 ranks Laplace 192 (884K per rank) on 80
GB nodes: nothing reclaimed, the CA0 leaf kept (27.8 GB). Per level, CA0 /
Color:

| 8 ranks | leaf (4096 boxes/rank) | level 4 (512) | factorization |
|---|---|---|---|
| 40 GB nodes (c1b) | 5.05 / 5.90 s | 2.29 / 1.82 s | 8.73 / 9.07 s |
| 80 GB nodes (c1b2) | 3.17 / 3.58 s | 2.39 / 1.33 s | 6.89 / 6.16 s |

The 40 GB leaves run under memory pressure (~1000 reclaims, 15-18-chunk
transitions ~0.95 s, 27 GB downloads), Color's more; with room, Color's
exchanges shrink too (level 4 comm 0.71 -> 0.25 s). CA0 still wins the
16^3-brick leaf and loses the 8^3 level: CA on the leaf only (`--CA_level`
= leaf) would give ~5.7 s against Color's 6.2 s on 80 GB nodes.

**P1/P2 baseline: the leaf's waves after C1b (2026-09-29, 40 GB nodes,
CA on the two finest levels, H2_GPU_WAVE_TRACE=1; build/gpu_prep/ca0/p12).**
Per rank (range over ranks), ms:

| leaf | waves | host plan (box plan / owner-pass plan) | launch | sketch wait | device sketch + elimination |
|---|---|---|---|---|---|
| 8 ranks, CA0 (5282 boxes) | 3601-3711 | 1448-1527 (540 / 700) | 253 | 1459-1507 | 1096 + 1734-1857 |
| 8 ranks, Color (4096) | 2626-2731 | 871-1107 (340 / 410) | 417 | 967-1029 | 817 + 1182-1257 |
| 64 ranks Nmin 27, CA0 (5282-6936) | 622-893 | 377-547 (110 / 220) | 84 | 45-55 | 87-115 + 357-522 |
| 64 ranks Nmin 216, CA0 (954-1720) | 614-1055 | 107-198 (20 / 40) | 44 | 389-664 | 309-490 + 268-472 |

The device is idle from the end of a wave's ID to its first elimination
launch: the box plan, the metadata and the launch, and whatever of the
owner-pass plan outlasts the box region's kernels. That is ~0.7-0.8 s of an
8-rank leaf and ~0.3 s of a Nmin-27 leaf; the Nmin-216 leaf at 64 ranks is
device-bound on its slowest rank (~0.1 s idle). The sketch's first mark
(metadata upload and the list kernel) is 207 ms of the 8-rank leaf.

Both plans are mostly structure that does not depend on the wave's ID: the
neighbor row layout (the neighbors' current sizes), the blocks and targets
each pair reaches (hash lookups, allocations), the owner pass's (target,
source, neighbor) tuples in their contribution order; the ID only sets each
box's k and r and hence buffer offsets and GEMM shapes.

**P1 step 1, then stopped (2026-09-29, build/gpu_prep/ca0/p12).** The
host plan's parts now print with the level ("plan boxes [buffers, items, of
which near-block allocations], owner pass [candidates, pairs, join,
batches]"). 8 ranks, CA0 leaf, rank 0: box plan 0.54 s (buffers 0.23, items
0.31 of which near-block allocations 0.17), owner pass 0.71 s (pairs 0.61).
Every owner target belongs to one candidate (its Schur block; a pair with a
higher box or with another rank's box), so the owner pass now finds the
candidates' targets and tasks in parallel, each in the sequential order
(a CA level's order by source included), then creates new targets and joins
them in candidate order: the same batches, results bitwise unchanged
(logdet, replica check 0 differ, solve check 3.988e-15 / 2.7e-14). Owner
pass 0.71 -> 0.26 s (CA0), 0.40 -> 0.13 s (Color); the 8-rank leaf's
elimination only ~0.1 s less (3.98 -> 3.87 s CA0, 4.80 -> 4.71 s Color),
mostly hidden behind the box region before. The rest (box plan, metadata,
launches, ~1.0 s of the 8-rank CA0 leaf while the device waits) was left:
parallel box items (~0.1-0.2 s) or a chunked box region (~0.3-0.4 s) were
judged not worth it now. P2 (sketch lists on a side stream) skipped.

**In-between sizes (2026-09-30, 40 GB nodes, Nmin_leaf 64, CA on the leaf
only; build/gpu_prep/ca0/mid).** 16^3-box leaf bricks with leaf boxes of
4-5.3 points per edge, between the 110K and 884K per rank of before.
Factor time (warm-up included), Color / CA0, and the leaf, s:

| | factor | leaf | level above |
|---|---|---|---|
| 288^3, 64 ranks (373K/rank) | 6.32 / 6.35 | 2.04 / 2.04 | 1.46 / 1.52 |
| 320^3, 64 ranks (512K) | 7.63 / failed | 2.33 / - | 1.69 / - |
| 320^3, 64 ranks, nothing kept | 8.69 / 8.27 | 2.88 / 3.05 | 2.03 / 1.96 |
| 336^3, 64 ranks (593K) | 8.23 / 8.97 | 2.69 / 3.67 | 2.23 / 1.95 |
| 640^3, 512 ranks (512K), nothing kept | 12.0 / 11.7 | 3.09 / 3.35 | 2.09 / 2.10 |
| 672^3, 512 ranks (593K), nothing kept | 12.3 / 12.8 | 3.44 / 3.89 | 2.36 / 2.44 |

CA0's 320^3 run kept its leaf's factors on the device and then ran out of
heap at level 5 (see the fragmentation entry below); at 336^3 its leaf was
not kept (the keep bound did not fit) and it lost. Ties at best: the CA0
leaf only wins where it keeps its factors (X_NR not downloaded).

**80 GB nodes: 384^3 on 64 ranks and 768^3 on 512 (2026-09-30/10-01,
884K per rank, Nmin 216, CA on the leaf only; build/gpu_prep/ca0/g384_80,
g768_80, a binary frozen before the environment changes).** Factor time,
Color keeping its solve factors / Color not keeping / CA0 not keeping (CA0
cannot keep its leaf there: with its ghosts the factors come to ~48 GB per
rank, the size its first solve uploads): 384^3/64 12.0 / 12.8 / 13.5 s;
768^3/512 16.0 / 17.3 / 19.6 s
(Color keeping's first run carried a 0.60 s cold warm-up). Leaf (rank 0
total; elimination shortest/longest over ranks):

| | Color keep | Color no keep | CA0 no keep |
|---|---|---|---|
| 384^3/64 leaf | 5.42 (4.86/4.97) | 6.38 (5.76/5.92) | 6.27 (3.91/5.32) |
| 768^3/512 leaf | 5.75 (4.98/5.12) | 7.05 (5.82/6.54) | 7.20 (3.96/5.93) |
| level above the leaf, 384^3 / 768^3 | 2.36 / 2.55 | 2.69 / 3.20 | 3.42 / 4.82 |

CA0's leaf eliminates 0.6 s faster than Color's without keeping but its
ghost imbalance grows with the rank count (2 s between ranks at 512: more
interior ranks with ghosts on all sides), and the level above it looked
0.7-1.6 s slower: that is mostly the first device exchanges' cost, which
Color pays in its leaf (entry below). Color weak-scales well (leaf 5.42 ->
5.75 s). Same accuracy (acc_forward 9.05e-3 / 9.06e-3).

**Kept factors and heap fragmentation; the remote generators' buffers in
pieces (2026-10-01, build/gpu_prep/ca0/frag, pieces).** A temporary heap
map at each level start (removed) on the 320^3 CA0 run that failed: on its
most-ghosted ranks level 5 started with 15.5 GB free in 619 ranges, the
largest 0.99 GB, between the leaf's kept factors (11.9 GB in 32 blocks, one
per wave) and the parent blocks the device transition left (6.0 GB in 8756
blocks of at most 1 MB); Color's kept blocks sat together (largest free
9.7 GB). The failing request, 2.07 GB, was `apply_remote_generators`' one
buffer for all generators of the first exchange that came through the host
(their uploaded temp2, X_NR, X_RS, ~10 MB each), and the same pass's
product buffer is built the same way. Moving only the kept blocks would
have left 3.4-4.2 GB; reserving the keep bound up front would have cost the
leaf ~12 GB of its working memory (see the next entry). Instead both
buffers now come in pieces of at most 256 MB (`PiecedBuffer`,
device_heap.hpp: the sub-buffers laid out as before, consecutive ones
grouped into separate heap blocks). The 320^3 CA0 run then completes with
its leaf kept: replica check 0 of 124704 differ, solve check 1.0e-14,
logdet as before; factor time 7.35 s against Color's 7.58 s (leaf 2.22 vs
2.46 s). Compaction of kept blocks stays in reserve for a heap that runs
out with most of its free memory in small ranges.

**The keep decision from an estimate (2026-10-01, build/gpu_prep/ca0/bound,
bound2).** The same diagnostic measured `solve_keep_bound` against the
bytes the levels kept: 0.44-0.55 of it at the leaves and level 5, up to 0.69
higher up (r about n / 2, and neighbours eliminated earlier have shrunk to
their skeletons). The bound's T/X_SR/LU term is now n^2 (was 1.5 n^2; still
a bound), and the decision uses bound x ratio: 1 for the first level, then
1.25 x the largest factors/bound measured on the level below (factor sizes
counted per wave whether kept or not); each level's kept factors stay within
0.8 of the heap less what it held at the decision, else the level's keep
fails (reason printed: budget or no free range) and later levels go back to
the bound. A first version with 0.75 for the leaf admitted 336^3's CA0 leaf,
whose ghosts' working set then left no room: its keep failed during the
level and the fallback (each kept box's X_NR copied back to the host) made
the factorization 19.4 s; hence the bound for the first level. 64 ranks, 40
GB nodes, factor time, CA0 / Color:

| | factor | leaf (elimination) | levels kept |
|---|---|---|---|
| 320^3 | 7.33 / 7.58 | 2.24 (1.65) / 2.42 (1.96) | all / all |
| 336^3 | 8.89 / 8.14 | 3.67 (2.93) / 2.83 (2.40) | all but the leaf / all |

320^3's CA0 level 5 is now kept too (solve setup 0.95 -> 0.17 s). Whether
the CA0 leaf keeps its factors decides CA0 against Color at these sizes.

**The first device exchanges, and an exchange warm-up (2026-10-01,
build/gpu_prep/ca0/bound, mpiwarm, warm).** The Color level above a CA0
leaf looked slower than Color's (192^3/64, Nmin 216: 0.80 vs 0.58 s;
384^3/64 3.42 vs 2.69; 768^3/512 4.82 vs 3.20 s) with the same device work
and the same bytes: its exchanges' payloads were slow in the first waves
(wave 0 at 2.8 GB/s, 12 GB/s from wave 4). Color's leaf shows the same ramp
(4.1 -> 14 GB/s): a pair of ranks' first device-memory exchanges are slow,
and whichever level exchanges first pays it, Color's leaf or the level above
a CA leaf (whose waves exchange nothing). A 64 KB message with each one-hop
peer before the first level removes it, as well as sending through the
whole arena does (so it is the first contact, not the arena's registration):
`warm_up_exchange` (factorization_driver.hpp), always on, once per process,
every level's peers, from and into the arena (host buffers without
GPU-aware MPI), timed with the kernel warm-up apart from the factor time
(0.1 s at 64 ranks). 64 ranks, 40 GB nodes, factor time CA0 / Color:

| | before | with the exchange warm-up |
|---|---|---|
| 192^3, Nmin 216, CA on the leaf only | 3.67 / 3.27 | 3.49 / 3.31 (level above the CA leaf 0.83 -> 0.65 s) |
| 320^3, Nmin 64 | 7.33 / 7.58 | 7.04 / 7.39 (level above the CA leaf 1.94 -> 1.73 s; Color's 1.69) |

Unchanged bitwise (8-rank Laplace 192: logdet, replica check 0 differ,
solve check 3.988e-15 / 2.7e-14). With the leaf's factors kept (320^3), CA0
now beats Color by 5%.

**GPU evaluators in place of the kernel kinds (2026-10-01,
build/gpu_prep/eval).** The library no longer has kernels of its own: the
application registers a GPU evaluator of its entries (doc/gpu_kernels.md;
GPU_INTERFACE/, h2_gpu/evaluator.hpp), an entry evaluator (C++, inlined
in the eval and fused sketch kernels, or CUDA source compiled with NVRTC), a
block evaluator, or an entry-list evaluator. The former kinds are examples in
EXAMPLE/gpu/ (Laplace, VIE, EMSURF's triangle pairs as a block evaluator,
CFIE as an entry evaluator) and EXAMPLE/user_block_funcs_george_gpu.py (the GP
kernel, NVRTC and CuPy). Against HEAD built apart, 40 GB nodes: 8-rank
Laplace 192 bitwise as above; 64 ranks, factor time old / new (two runs
each): 192^3 CA0 3.61, 3.51 / 3.75, 3.57, Color 3.18, 3.13 / 3.22, 3.17;
320^3 CA0 7.14, 7.16 / 7.24, 7.18, Color 7.48, 7.41 / 7.47, 7.52 (320^3
logdets bitwise; 192^3 within one ulp, as the old build varies from run to
run). EMSURF EFIE H2, 8 ranks: sphere 12K 1.17, 1.13 / 1.07, 1.06 s; 49K
1.88, 1.85 / 1.80, 1.85 s (accuracy the same to 7 digits). HODLR
construction: EFIE (now a block evaluator, one host-id synchronization per
batch) 1.91, 1.58 / 1.62, 1.63 s; CFIE 2.01, 1.92 / 1.91, 1.90 s.
BPACK_CHECK=kernel: 0 (Laplace), 1e-16 (VIE), 1e-15 self to 1e-13 far
(EFIE triangle pairs); HODLR has the check too now.

**Evaluator and HODLR warm-ups (2026-10-01, build/gpu_prep/eval, warm).**
An evaluator's first call carried its first-use costs into the first level
(Python GP, 100K points, 8 ranks: the CuPy route's leaf 4.3 s against 0.6 s
with NVRTC, all of it CuPy compiling its kernels and starting up; the NVRTC
compile too, less). Evaluator::warm_up now calls each evaluator once (a
block of up to 32 of the rank's points, and its sketch) before its first use,
outside the reported times: H2 factorization (structured and unstructured;
the unstructured path now also warms the batched kernels), H2 compression,
and HODLR, whose new warm-up (GpuState::warm_up) also runs the
factorization's batched LU, solves and GEMM forms per size range and the
construction's QRs, SVDs and cuSOLVER handles. Python GP: evaluator warm-up
NVRTC 0.70 s, CuPy 4.10 s with a cold kernel cache and 0.23 s with a warm
one; the levels then sum to 0.75 s (NVRTC) and 0.80 s (CuPy, warm cache).
HODLR construction old / new (two runs each): Laplace 48^3 4.62, 4.22 /
4.11, 4.14 s; EFIE 1.70, 1.57 / 1.48, 1.49; CFIE 2.03, 1.90 / 1.80, 1.79;
VIE 1.12, 1.15 / 1.01, 1.03 (warm-up 0.15-0.20 s apart). EMSURF H2 factor
1.32, 1.12 / 1.00, 1.04 s. Laplace H2 bitwise as before. Not done: a
first-message warm-up for HODLR's exchanges (not measured).

**Python: split initialization, warm-up at registration (2026-10-01,
build/gpu_prep/eval, p3_*).** The worker now runs py_bpack_init, registers
the evaluator, warms it up with c_bpack_gpu_warm_up (H2: on the GPU the
factorization will use, after begin_operator_build so ranks sharing a GPU
split it first), then py_bpack_compute: the evaluator serves HODLR too, and
its first-use cost counts in the initialization. GP, 100K points, 8 ranks on
one node, three alternating runs of each H2 route: NVRTC compile 0.52-0.62
s, CuPy warm-up 3.25 s (cold kernel cache) / 0.17-0.25 s (warm); levels
NVRTC 497+109+16, 497+122+15, 490+109+15 ms, CuPy 574+109+15, 580+115+30,
576+111+15 ms: the upper levels equal, the leaf about 80 ms (16%) slower with
CuPy. HODLR from Python: GPU construction 0.82-0.96 s (CPU 1.97-2.01 s),
kernel check 1e-16 self / 3e-16 far, accuracy as the CPU's.


**768^3 on 512 ranks again, with the automatic keep, the pieced buffers and
both warm-ups (2026-10-01, job 59166764, 128 A100 80 GB nodes, commit dc52665
frozen; build/gpu_prep/ca0/g768_80 of ButterflyPACK_H2_new_id_proxy2_gpu).**
Color and CA0 (CA on the leaf only), each twice; the first run of each is
slower (cold nodes for Color: kernel warm-up 0.67 s; the first CA run's leaf
8.04 s), so the second runs are compared. CA0 still loses, by 8% (it was 22%
against Color keeping its factors on 2026-09-30):

| second runs | Color | CA0 | 2026-09-30: Color keep / Color no keep / CA0 no keep |
|---|---|---|---|
| factor time | 15.6 s (first run 16.9) | 16.9 s (first run 18.3) | 16.0 / 17.3 / 19.6 s |
| leaf (level 7) total | 5.58 s | 6.66 s | 5.75 / 7.05 / 7.20 s |
| leaf elimination, shortest / longest rank | 4.98 / 5.07 s | 5.46 / 5.62 s | 4.98 / 5.12, 5.82 / 6.54, 3.96 / 5.93 s |
| level 6 (above the leaf) | 2.68 s | 3.71 s | 2.55 / 3.20 / 4.82 s |
| levels 5 / 4 / 3 / 2 | 1.61 / 1.47 / 2.89 / 1.08 s | 1.65 / 1.41 / 2.1 / 1.09 s | |
| solve factors kept | levels 2-7, 31.2 GB per rank | levels 2-6, 10.5 GB; the leaf (37.7 GB) uploaded at the first solve | |
| solve setup; logdet and quick verification | 0.18 s; 2.8 s | 6.96 s; 9.8 s | |
| heap peak at the leaf | 49.6 GB | 40.7 GB | |

Accuracy as before: Color logdet -6968667795.3627348, acc_mvp 5.23e-4,
acc_forward 9.06e-3; CA0 -6968667834.83634, 5.32e-4, 9.05e-3. No errors.
What is left of CA0's deficit (1.3 s): (1) its leaf's factors, with the
ghosts', do not fit the keep budget on 80 GB, so the leaf downloads them
(rank 0: 27.7 GB down against Color's 7.9 GB, 0.65 s of level-end wait) and
is 1.1 s slower than Color's kept leaf; (2) the Color level above the CA
leaf is still 1.0 s slower than Color's own level 6 (exchange 2.65 vs 1.69 s:
sizes 1.51 vs 1.08, payload 0.95 vs 0.40 s for the same bytes), so the 64 KB
exchange warm-up does not remove the first-exchange cost at this size (it
was 1.6-2.3 s before). The leaf elimination's spread over the ranks no
longer shows in "shortest / longest" (5.46 / 5.62 s against 3.96 / 5.93 s):
not looked into (the longest is what the level waits for). Level 3 (64 ranks
after the reduction, 8 boxes each) varies between 2.0 and 2.9 s across runs
without following the method (here Color 2.89 s, CA0 2.1 s; 2026-09-30:
2.68 / 2.16 / 2.62 s): not understood.

**640^3 on 512 ranks, 40 GB nodes (2026-10-01, job 59185735, commit dc52665
frozen; build/gpu_prep/ca0/g640_40): the 512-rank counterpart of 320^3 on 64
ranks (512K per rank, Nmin 64, 125-point leaf boxes, CA on the leaf only).**
Color: 11.8 s (first run, cold nodes) and 11.3 s, every level kept (14.4 GB
per rank), logdet -3937186064.8561087, acc_mvp 9.43e-4, acc_forward 6.61e-3
(12.0 s on 2026-09-30 with nothing kept). **CA0 kept its leaf's factors
(11.94 GB per rank, 51% of the bound) and levels 6, 5 and 4, then every run
(three, the third with the checks) ran out of heap after level 4's
elimination, in the transition to level 3 (the process reduction to 64
ranks): on 5-13 of the receiving ranks a 1.3-1.4 GB request with 22.8 GB of
the 33.4 GB pool in use and the largest free range 1.1-1.3 GB** (10.5 GB free
in pieces: the kept factors, 18.9 GB in blocks of one wave each, sit between
the levels' blocks, as on the 320^3 CA0 run before the pieced buffers; Color
keeps 13.6 GB without trouble because its kept blocks lie together). This is
the symptom compaction of the kept blocks was held in reserve for.

Up to level 4 (second runs, Color / CA0): leaf 2.42 / 2.63 s (elimination
1.96 / 1.86, shortest ranks 1.89 / 1.75; CA0's gather 0.23 and transition
0.47 against Color's 0.46), level 6 1.85 / 1.99 s, level 5 1.35 / 1.37 s,
level 4 elimination 0.89 / 0.86 s. The level above the CA leaf shows no
first-exchange penalty here (exchange 1.12 vs Color's 1.13 s), so the 1.0 s
penalty at 768^3 comes with the larger messages (4.35 GB received per rank
against 2.98 GB), not with the rank count. CA0's elimination advantage at the
leaf is 0.1 s here (0.3 s at 320^3 on 64 ranks: more interior ranks); with
the gather and the larger transition it was 0.3 s behind Color through level
4, so a finished run would have been about a tie.


### M5. CA3 (component-owner) on the GPU
Unique owner per boundary component, the corner -> edge -> face DAG,
FULL/COMPACT/SKELETON routing, deterministic replay (R1/R2/R3). Risk: small
component batches underfill the GPU.  The design below (2026-10-02, after
M9: CA0 and Color tied at the leaf, CA0's leaf device-bound with 1.7x the
boxes) is the plan the user asked for before the code is written; the CPU
feasibility numbers follow it.

**Why CA3 on the device.** Color's leaf cost is 16 transports, each waiting
for the slowest rank's wave (0.46 s of 1.9 s at 320^3/64, 0.6-1.1 s per
level at 640^3/512); CA0 removes them at the price of 1.7x the device work
on interior ranks.  CA3 eliminates every box once and has three dependent
communication stages (blue -> purple -> green) plus the entry fetch and the
return of the results, instead of 16; its overheads (fetch, routing,
install) are bytes, which the device exchange moves at 7-12 GB/s.

**Reuse (what exists and is kept as is).**
- Host, MPI-free: `h2_parallel/core/dataflow.hpp` (process grid, component
  graph by flood fill per colour, owner assignment with its hash, fetch
  plan, `payload_recipients`, `auxiliary_components`, the counter Runtime).
- Host logic of `owner_schedule.hpp`: the per-component wave partition
  (`owner_wave_partition`, a pure function of the graph, hence the same on
  every rank), the synchronous colour schedule, the pair watermarks
  (R1/R2/R3) that make the numbers independent of arrival order.
- Device: the Color eliminator's generator transport (pack on the device per
  destination, explicit receives with transport numbers, device to device;
  the receiver installs the generator boxes and forms its blocks' updates
  on the device straight from the receive buffer, sources in (wave, Morton)
  order, deltas summed in source order) -- this is CA3's COMPACT payload,
  install and passes in one.
- Device: `device_pieces.hpp` / `ca_halo.hpp` (pieces exchange with chunked
  rounds through the exchange arena) for the entry fetch and the FULL
  return; the eliminator's CA mode (writable remote boxes eliminated here,
  LU and solves by size class so every copy factors identically).
- Device solve: home based.  The FULL payload restores the home rank's
  factors on its device (T, X_SR, LU, pivots, X_NR into the kept store), so
  the Color device solve runs unchanged; the CPU's owner solve is not ported.

**Schedule on the device: synchronous colours.** Per boundary colour c
(blue+orange, purple, green): (1) one batch of all owned components of c,
their boxes in the merged component waves (same-colour components never
share a live pair), with the copies of c that already arrived installed in
the same batch; (2) pack and send the batch's payloads; (3) wait for the
remaining copies of c and install them in one batch.  Then the interior as
today (CA0's interior group).  The arrival-driven preemption of the CPU
schedule is a latency device for thin thread teams; on the GPU a colour's
owned work is one wide batch anyway, and the two schedules have the same
numerics (invariant 3 of owner_schedule.hpp).  Later (M5c) the installs of
a colour become arrival driven within the wait.

**Box kinds on the device.** Local home boxes (as today); fetched boxes
(owned components' boxes of other ranks: entry Schur and near blocks and
points, writable, eliminated here, their FULL state goes home afterwards);
installed copies (eliminated elsewhere: generator state for the pairs of
interest, marked eliminated, their temp2^T and X_RR_full kept as fill
sources of the sketches and of the transition, as the Color path keeps
remote generators); skeleton-only boxes (two hops: index sets, the assisting
boxes' skeleton refresh of the Color transport).

**Payloads.** From `payload_recipients`: FULL to the home rank when the
owner is not the home (the generator plus the solve factors), COMPACT to the
ranks holding a one-hop pair of interest (the generator: X_NR rows for the
recipient's boxes, T, X_RR factors, as the Color transport), SKELETON to the
two-hop readers.  Carried by the device exchange with a per-batch peer set.

**Determinism.** Canonical source order (colour, wave, Morton) with pair
watermarks: the owned waves update their pairs in wave order (R1, the
eliminator's owner pass); the installs apply the batch's sources in (wave,
Morton) order (the generator application); the R2/R3 bookkeeping (which
sources a pair has already consumed when its other endpoint became live) is
host logic over box ids, ported from `pair_watermark`.  Check: the FULL
state at the home against the owner's copy, bitwise (as the replica check);
box-by-box skeletons against the CPU CA3 on Laplace; logdet and
H2_CheckError against the CPU CA3; the GPU solve check.

**Milestones.**
- M5a (runs only): CPU CA3 at the GPU's brick size (Laplace 160 Nmin 64 on 8
  ranks; 320^3 on 64 when a slot is free): components, batches, boxes per
  wave, bytes sent, the CPU's split; projection with the GPU's measured
  per-wave costs; go/no-go.
- M5b: the device CA3 level: entry fetch whole at level entry, synchronous
  colours, payloads over the device exchange, home-based solve.  New file
  `h2_gpu/ca_owner.hpp` (schedule glue, fetch and return plans, payload
  framing), the eliminator's CA mode extended by the box kinds and the
  batch entry points, a `ca_level_runs` for owner mode 3.
- M5c: overlap: entry fetch staged by colour, sends overlapped with the
  installs, arrival-driven installs within a colour's wait.
- M5d: 320^3/64 and 640^3/512 against Color and CA0, per-colour breakdown.

**M5a results (2026-10-02, build/gpu_prep/ca3/cpu8.sh; CPU build, Laplace
160 Nmin 64 on 8 ranks: 16^3 bricks of 125-point boxes, CA on the leaf,
level 5).**  Leaf level (rank 0): Color 16.9 s, CA0 15.5 s, CA3 14.5 s
(elimination 13.5 / 12.0 / 10.5-11.0 s); factor time 62.6 / 61.2 / 56.7 s.
CA3's owner schedule at the leaf (rank 0, a corner rank of the 2^3 grid as
all eight are): 15 owned components with 1374 boxes (a brick's boundary is
1352: ownership is balanced to about one boundary per rank), 11 installed
components, 7 batches; boundary waves 24 of 16-144 boxes on rank 0, 42 of
3-180 on the slowest rank; the CPU's split: box region 1.97 s, passes 1.39
s, R2 replay 0.16, final 0.33, pack 0.11, install 0.46 (its thread-sum:
temp2 solves 4.0 s, block rebuilds 1.2 s), idle wait 0.52 s; 2.4 GB sent,
2.0 GB received per rank (the GPU Color leaf at 320^3/64 receives 2.1 GB,
CA0's halo 3 GB).  Accuracy as Color's (acc_mvp 1.81e-3 against 1.58e-3,
acc_forward 4.72e-3 against 4.97e-3); the serial oracle's H2_CheckError
values equal the arrival-driven run's (the logdets differ in the 16th
digit: the reduction's order).

Projection for the GPU leaf at 320^3/64 (slowest rank, from the measured
Color/CA0 device costs of ~0.23 ms per box): boundary 1374 boxes ~0.32 s +
interior 2744 boxes ~0.64 s + installs ~0.1 s + fetch ~0.1 s + 2.4 GB of
payloads ~0.25 s (overlappable) + three dependency waits ~0.3 s: about
1.45-1.7 s against Color's 1.41 s and CA0's 1.38 s of elimination.  Even at
64 ranks, with the upside at 512 ranks, where Color's 16 per-level waits
cost 0.6-1.1 s and CA3 has four.

**M5b as built (2026-10-03; `h2_gpu/owner_plan.hpp`, `h2_gpu/ca_owner.hpp`,
the eliminator's owner mode).**
- The host schedule is reused unchanged: `owner_schedule_setup` builds the
  graph, the assignment, the wave partition and the filters as for the CPU;
  `build_owner_plan` (ca_owner.hpp) turns it into the eliminator's view
  (`OwnerPlanView`: the owned ghosts, the foreign local boxes, each boundary
  box's sequence number colour * num_waves + wave and component, each owned
  box's recipients as tree ranks from `payload_recipients`, per colour the
  owned waves and the send and receive peers, the interior waves by Morton
  parity, the sketch-context hook that sets `owner_current_sketch_comp`).
- The eliminator (`owner=true`): `is_writable` = eliminated here (a local
  box unless foreign, an owned ghost), `is_interest` = local or owned
  ghost (the owner pass's candidates and pair targets; the Schur block
  only of writable boxes), `held_box` for the sources' host boxes.  The
  sketch plan sets the sketch context per box (the level's
  `sketch_eliminated_filter` is the schedule's).  Each owned wave packs
  its boxes' pieces per recipient (`pack_owner_generators`: the generator
  as the Color transport's GenLayout, the X_NR rows of the recipient's
  interest boxes, plus T, X_SR, LU and pivots for the home rank; SKELETON
  recipients get a header entry) into a block of the exchange arena.  An
  owned ghost's factors are neither kept nor copied to the host here (its
  temp2 and X_RR_full are, as restorable fill sources).
- The transport is per wave (M5c, `owner_progress`): a wave's block of
  the exchange arena is sent as it is, a header and the bulk per
  recipient (tags 733/734, pieces of 256 MB), once its device copies are
  done (a stream event, polled at every wave), and freed once sent; the
  receiver probes for headers at every wave and in the colour's wait and
  posts the bulk receives straight into the arena.  At a colour's end
  every recipient gets a marker (tag 735: how many wave messages this
  rank sent it this colour), and the rank waits for its receive peers'
  markers, as many messages and their bulk.  No colour-sized buffer and
  no concatenation copy: the arena holds a colour's arrivals plus the
  few wave blocks in flight (the first version concatenated the colour's
  pieces per recipient and sent at the colour's end, which spilled from
  the arena into the main heap, whose registration by the CUDA-aware MPI
  made the green colour's transport 1.0-1.5 s on some ranks).  The
  arrivals of a colour are installed at its end in sequence order, per
  wave: the skeletons of its full-rank
  and skeleton-only boxes, then `apply_remote_generators` (owner mode: the
  installed box may be a ghost or a foreign local box; the rows and
  updates cover the boxes of interest) for its generators; the FULL
  factors of this rank's boxes land in the kept store (`land_full_factors`:
  the keep budget as `keep_solve_factors`) or in the host boxes, T always
  on the host; the installed generators' temp2 and X_RR_full go to the host
  boxes in the background (the fill-source restore of the transition).
  The colour's receive buffers return to the arena once the installs'
  launches ran (a stream event), as the next receives may land in them.
- The interior runs as CA0's, after the colours; `TreeLevel::ca_home_factors`
  makes the host and device solves and the multiplies run the level as a
  Color level (no replicated ghosts, no owner solve engine).  A device
  transition ahead of an owner level keeps the blocks on the device; the
  owned components' remote boxes' entry blocks then move between the
  devices by the fetch plan (`exchange_ca_owner_blocks`, the halo's pieces
  exchange over caller lists: `exchange_ca_blocks_for`).
- Rules found by the 8-rank tests (every one a block the first version
  left stale or updated when it was never read): a box homed here but
  eliminated elsewhere gets, at its install, its pairs with ALL its
  neighbors (not only the boxes of interest) replaced or created, as its
  owner's elimination made them, and its owner's post-elimination Schur
  block lands with the FULL state (then takes the later sources' updates
  here); an owned ghost's blocks are spent after its elimination here,
  except its pairs with boxes homed here or still to be eliminated here;
  the generator's view of its neighbors is checked only for the boxes of
  interest (a box of no interest it saw eliminated may have no state
  here).  The solves and multiplies run the level in the schedule's
  sequence order (TreeLevel::ca_home_order), not the Morton parities.
  The installed generators stay resident for the level (no host copies to
  restore them from; the receive buffers return to the arena after the
  device updates).

**M5b first results (2026-10-03, build/gpu_prep/ca3/gpu8.sh, frozen_m10h;
8 ranks, Laplace 160 Nmin 64, CA3 on the leaf).**  The device CA3 level
reproduces the CPU CA3 (cpu8_ca3.log): the leaf's compression ratio
0.284024 and skeleton size 35.5031 identical, level 4's 55.5063 identical,
logdet -50156606.410895944 against -50156606.410895959, every
H2_CheckError value identical to all printed digits (acc_mvp 1.8095410e-03,
acc_forward 4.7232536e-03); the device solve against the host solve
2.9e-15 / 1.4e-14.  CA3 on levels 5 and 4 (the owned fetch between the
devices, 1.81 GB received in 2667 blocks in 0.08 s) runs too (logdet
-50156606.734, acc_mvp 1.47e-3, acc_forward 4.65e-3).  CA0 and Color
unchanged.  Leaf time 2.62 s against Color's 1.39 and CA0's 1.38 s: waves
0.76 s (owned 1374 boxes in 24 waves, interior 2744 in 8), the three colour
transports 0.58 s [the slowest rank's colour 0.35, transfers 0.22: 3.04 GB
sent, 2.50 GB received per rank] plus 0.13 s of installs, the transition
0.24 s as Color's, and about 0.45 s around the waves (level start and
finish) still to be found.  At Laplace 192 (884K points per rank) the
owner level exhausts the 40 GB device (peak 33 GB against Color's 16 GB
at 160^3; the colour payloads spill from the exchange arena into the main
heap and fragment it): to be bounded by a wave-by-wave transport (M5c).

**M5b at 320^3 on 64 ranks (2026-10-03, job 59254764, build/gpu_prep/ca3/g320,
frozen_m10k: per-colour balancing, 256 MB payload pieces).**  CA3 on the
leaf (level 6) against Color and CA0 of the same build in the same job:

| | Color | CA0 | CA3 (owner level) |
|---|---|---|---|
| Factor time | 6.09 s | 6.16 s | 8.17-8.39 s |
| leaf level | 1.94 s | 2.02 s | 4.10-4.20 s |
| leaf: waves / transports / end wait | 0.95 + 0.47 (16 transports) | 1.27 (32 waves) / halo 0.12 | 0.73-0.76 / 0.76-0.77 (3 colours: sizes 0.3, transfers 0.3; 2.14 GB sent, 2.80 GB received) / 1.6-1.7 s |
| logdet | -446707781.18238521 | -446707784.16925561 | -446707783.22014749 |
| acc_mvp / acc_forward | 1.17e-3 / 5.37e-3 | 9.67e-4 / 5.49e-3 | 9.22e-4 / 5.23e-3 |
| device solve / check | 0.338 s | | 0.346 s, 5.6e-15 / 1.05e-14 |

Correct (the solve checks, the accuracy as the others'), with 1194 owned
boxes in 14 components and 12 installed components per rank, but 2.2 s
slower per factorization: 1.6 s of it is a wait at the level's end for
some ranks (every rank's finish takes 0.01 s; the level-end collectives
wait for the late ones), the same wait that shows as 0.9 s on 8 ranks at
160^3, where the exchange trace puts it on one node's four ranks at the
green colour; found with the per-rank timeline below.

- Not bitwise across ranks: a pair held on two ranks accumulates its
  contributions with the owner pass's grouping on one and the generator
  application's on the other, as in the Color path (the CPU CA3's R1/R2/R3
  machinery makes pair copies bitwise equal; the device path is instead
  deterministic by its synchronous schedule, and its fill sources are the
  owners' bits on every rank).  Requires CUDA-aware MPI (else the level
  runs on the host).

The level-end wait was found with the per-rank timeline (`BPACK_TRACE=
exchange`: elapsed since a barrier at the level's start at the end of
every colour's waves and transport): one node's four ranks spent 1.0-1.5 s
in the green colour's transport while the others took 0.5 s; 256 MB
pieces did not change it; a 6 GB exchange arena (`BPACK_GPU_EXCHANGE_MB`)
did (leaf 2423 -> 1540 ms, end wait 0.89 -> 0.06 s): the colour's payload
buffers (the waves' blocks, the concatenated send buffers and the receive
buffers, about 3 GB in the green colour) overflowed the 2.5 GB arena into
the main heap, whose pages the CUDA-aware MPI registers at about 1.6 GB/s.

**M5c: the wave-level transport (2026-10-03, frozen_m10l; 8 ranks,
Laplace 160 Nmin 64, CA3 on the leaf, build/gpu_prep/ca3/gpu8.sh).**  The
protocol above (sends from the waves' arena blocks as they are packed,
receives as the headers arrive, colour-end markers) with the default
2.5 GB arena: the leaf 1999 -> 1068 ms (Color 1386, CA0 1382 ms), the
level-end wait 0.9 -> 0.01 s, the colour transports 0.34 -> 0.28 s (blue
0.03, purple 0.08, green 0.17 s), and the 0.45 s "around the waves" of
the first version gone with the spill.  Factor time 4.39 -> 3.30 s.  The
numbers are the previous build's to every printed digit (logdet
-50156606.410895944, acc_mvp 1.8095410e-03, acc_forward 4.7232536e-03).

The same build gave wrong numbers with CA3 on levels 5 and 4 (level 4's
compression ratio 0.18787 against 0.189075, acc_forward 0.69; once an
invalid pivot in a landed LU at the logdet).  The per-colour trace
(`[exchange] ... entries N (full, compact, skeleton); held boxes of the
colour E, without a payload M; arena fallbacks F`) showed the message sets
complete and identical to the colour transport's, but level 4's blue
colour sends 5.4 GB from rank 0 (8^3 bricks: larger boxes, larger
generators) and 21 of its buffers fell back from the arena to the main
heap.  The main heap's blocks are reused in the stream's order (a freed
block may still be read by launches in flight: the pack's X_NR source,
a replaced Schur block, a released fill source), which is safe for
every launch and for a wave block (its copies are launches) but not for
a receive buffer, which MPI writes as the data arrives.  The colour
transport never hit this: its receive buffers were allocated after a
stream synchronisation.  Fix (frozen_m10m): a receive buffer that falls
back to the main heap synchronises the stream before its receives are
posted (the arena's blocks need nothing: they are freed only after
their reads, by the sends' completion or the installs' event).

frozen_m10m (build/gpu_prep/ca3/val7.sh): the leaf 1076 ms, the numbers
of the previous build (logdet -50156606.410895944, the solve checks
2.9e-15 / 1.4e-14, Factor time 3.37 s); CA3 on levels 5 and 4 correct
again (logdet -50156606.734330446, acc_mvp 1.4740937e-03, acc_forward
4.6494769e-03, as frozen_m10k), CA0 and Color at 192 unchanged.  Level 4
of this case (8^3 bricks of 512-point boxes, 290 owned boxes) takes
2.6-3.0 s against the colour transport's 2.5 s: its blue colour alone
sends 5.4 GB from rank 0 and the colour's 2.0-2.4 GB of arrivals per rank
do not fit the 2.5 GB arena beside the wave blocks in flight (up to 21
fallbacks per rank), so the blue transport waits 0.8-1.1 s and the
fallback syncs stretch the next colour's waves; a non-leaf level is off
the target (CA on the leaf only).

**M5c at 320^3 on 64 ranks (2026-10-03, job 59258716, frozen_m10m, the
2.5 GB arena).**  Correct (logdet -446707783.22014749 as frozen_m10k, the
solve checks 5.6e-15 / 1.05e-14), Factor time 7.99-8.44 s against Color
6.12 and CA0 6.22 s in the same job; the leaf 2.9-3.3 s against 1.51 /
1.35-1.42 s.  The wave run's trace: rank 0's colours wait 0.04 / 0.19 /
0.18 s, but over the 64 ranks the purple colour's wait is up to 1.27 s
(median 0.11) and the green's up to 1.98 s (median 0.25): the colours
move up to 1.9 and 2.4 GB per rank each way, and the arena (2.5 GB) does
not hold a colour's arrivals beside the wave blocks in flight (up to 34
fallbacks per rank, 14-20 on the slow ones; the fallbacks' registration
and syncs make the stragglers everybody waits for at the colour's end
and the level's end: "level 6 end: total 1.7-2.0 s").  The device has
the room (CA3's heap peak 16-17 GB of 33 GB at this size; Color's 19 GB
at 640^3/512): frozen_m10n reserves, when the factorization uses the
component-owner level, an exchange arena of a fifth of the GPU (8 GB of
40, `begin_operator_build(..., 0.2)`; `BPACK_GPU_EXCHANGE_MB` still
overrides; doc/environment_variables.md).  On 8 ranks (gpu8.sh m10n) the
leaf is unchanged (1047 ms, the same numbers, the warm-up exchange still
0.09 s) and level 4 of the two-level case drops from 2.6-3.0 s to 1.02 s
(no fallbacks: blue transport 0.21 s, purple 0.03, green 0.05; Factor
time 5.18 -> 3.60 s, the numbers identical); the 192 case (884K per
rank) still exhausts the heap, now 28.4 GB.

**M5c at 320^3 on 64 ranks with the 8 GB arena (job 59259619, frozen_m10n).**
No arena fallback on any rank; the colour waits over the 64 ranks drop
to at most 0.05 / 0.22 / 0.46 s (medians 0.02 / 0.09 / 0.08); the leaf
1.65-1.80 s (from 2.9-3.3 s; Color 1.51-1.57, CA0 1.35-1.42 s), the
level-end wait 0.34-0.50 s; Factor time 6.74-7.06 s against Color 6.12
and CA0 6.22 s of the previous job (the numbers unchanged: logdet
-446707783.22014749).  The slowest ranks' owner-level body (the
timeline's interior mark) ends at 0.84-0.95 s, so the leaf's remaining
excess over Color is outside the colours: the owner level's start (plan
0.02 s, eliminator 0.10 s), the level-end wait, and about 0.6 s of
host work per level outside the elimination and transition timers
("Total level time" 2676 ms against elimination 1796 + transition 291;
Color: 1900 against 1569 + 276): the owner schedule's setup and the
level's cleanup, to be timed next.

**The owner schedule's host setup (2026-10-03, frozen_m10o: per-level
host-phase timers under verbosity 1, `level host phases (ms, this rank /
slowest): schedule setup, start to elimination, elimination to
transition, after the transition`).**  On 8 ranks at 160^3 the CA3 leaf
spends 100 ms from the level's start to the elimination (the schedule
setup 46-52 ms, the halo gather 18 ms) against Color's 13-18 ms, and 11-27
ms after the transition against Color's 20 ms.  The setup is
`owner_schedule_setup`: `build_component_graph` classifies and
flood-fills the boundary boxes of EVERY brick of the process grid, on
every rank (dataflow.hpp; hash maps, a neighbour vector per box, per
component the visible-rank search over its boxes), then validates the
graph and assigns the owners.  Measured on the login node
(h2_parallel/tests/owner_graph_timing.cpp, 16^3 bricks):

| grid | components | boundary boxes | graph | validate | assign |
|---|---|---|---|---|---|
| 2^3 (160^3 on 8) | 117 | 10.8K | 0.038 s | 0.004 s | 0.000 s |
| 4^3 (320^3 on 64) | 665 | 86.5K | 0.358 s | 0.033 s | 0.002 s |
| 8^3 (640^3 on 512) | 4401 | 692K | 3.60 s | 0.32 s | 0.087 s |

Linear in the number of ranks, replicated on each, inside the factor
time: 0.4 s per owner level at 64 ranks (the 0.5-0.6 s gap above), 4 s
at 512, where the leaf itself takes 1.5 s.  The CPU CA3 pays the same.
Confirmed at 320^3/64 (job 59260109, frozen_m10o): the CA3 leaf's host
phases are schedule setup 380 / 408 ms (this rank / slowest), start to
elimination 478 / 563 ms (the setup, the halo gather 59 ms and ~60 ms
more), after the transition 127 ms, against Color's 11 / 13 and 62 / 71
ms; the elimination itself 1658-1734 ms against 1503-1560.

**M5c at 640^3 on 512 ranks (job 59259703, frozen_m10n, build/gpu_prep/
ca3/g640).**  Correct (logdet -3937186080.0625076 in all three CA3 runs,
the solve checks 8.1e-15 / 9.9e-15), but the Factor time is 12.6-12.8 s
against Color 8.57 and CA0 8.35 s: the leaf's "Total level time" is
6.4 s with an elimination of 1.86-1.99 s and a transition of 0.29 s, so
about 4.1 s of host setup per leaf level, as the harness predicted (the
graph 3.6 s, its validation 0.3 s).  The owner level itself: waves
0.71 s, colour transports 0.42 s (blue 0.06, purple 0.14-0.18, green
0.17-0.23), level-end wait 0.56-0.60 s, the halo gather 0.14-0.20 s;
Color's leaf elimination 1.61-1.69 s, CA0's 1.42-1.51 s.  M5e below is
what makes CA3 viable at this scale.

**M5e (plan): the owner graph from the rank's neighbourhood.**  Every
component a rank uses (owned, shared, auxiliary, the sketch filter's and
the halo's boxes) lies in the 3^3 bricks around it: a corner cluster
spans the 8 bricks of its corner, an edge run the 4 of its edge, a face
the 2 of its face, and the two-hop readers' boxes are within two boxes
of the brick.  So: (1) each rank classifies and flood-fills the boundary
boxes of its own brick and its 26 neighbours (36.5K boxes, a constant),
giving the components' PARTS per brick with a key computable by every
rank that sees the brick (brick id, colour, the part's smallest Morton
index) and the links between parts across brick faces; (2) one
Allgatherv over the level's ranks of each rank's own parts (key, colour,
box count, smallest Morton index, the ranks holding its boxes, the ranks
adjacent to them, the lower-colour parts next to it: a few hundred bytes
per rank); (3) every rank unites the parts into the global components
(a union-find over at most 26 P parts), numbers them by (colour, smallest
Morton index) as today, sums their costs, sharers, adjacent and visible
ranks and predecessors from the parts, and assigns the owners (LPT, as
today, identical on every rank); (4) the box lists are attached only for
the components of the rank's neighbourhood (the parts it flood-filled),
complete for the components touching its brick, which is what the wave
partition, the payload recipients and the schedule's filters read.  The
result is a ComponentGraph with the same ids, owners and descriptors as
the global build; dataflow.hpp stays MPI-free by taking the gather as a
callback (owner_schedule.hpp supplies MPI_Allgatherv; the test supplies
the concatenation of all ranks' parts), and owner_graph_timing.cpp
becomes the test that the local graph of every rank equals the global
one on the 2^3 and 4^3 grids.  Expected cost: ~50 ms per level,
independent of the rank count; the CPU CA3 gets it too.

**M5e as built (2026-10-03; dataflow.hpp `LocalParts`,
`local_component_parts`, `assemble_component_graph`,
`build_component_graph_local`; owner_schedule_setup supplies the
Allgatherv over the level's ranks; the user chose it for both codes).**
A part is a component restricted to one brick (the flood fill within the
brick, seeds in Morton order, so its first box is its smallest index);
its key is (brick, smallest index).  A rank's stream describes its own
brick's parts: colour, box count, the ranks adjacent to its boxes and
those whose halo holds one of them, the keys of the same-colour parts
across the brick's faces (the links), of the lower-colour parts within
one hop (preds) and within the two-hop shell (preds2).  The assembly
parses every rank's stream, unites the parts by the links, numbers the
components by (colour, smallest index), sums the costs and unions the
rank sets and predecessors (preds2 minus preds minus self), derives the
successors and the HEFT ranks, attaches the neighbourhood's box lists
(`boxes_complete` when every sharer is in the neighbourhood: so for
every component touching the brick) and the rank's interior node.
`Component::first_morton` carries the global smallest index (the hash
and the ordering use it), `ComponentGraph::local_bricks` marks a local
graph (validate_graph checks the neighbourhood's bricks and the complete
components' boxes); `describe` and `assignment_hash` read the cost
instead of the local list; auxiliary_components and payload_recipients
throw if a boundary box they look at is outside the neighbourhood (the
brick must be at least 4 boxes wide for the two-hop shell of a
component touching the brick to stay inside it).  Test:
h2_parallel/tests/owner_graph_test.cpp (`make h2_owner_graph_test`):
every brick's local graph against the global one on the 2^3 and 4^3
grids of 16^3 bricks, three bricks of the 8^3 grid, the 2D 4^2 grid and
the 8^3-brick case; ids, colours, costs, owners, HEFT ranks, sharers,
adjacent, visible, preds/succs/preds2/succs2, the complete box lists and
the partial ones' membership, box_to_comp, the interior node, the
auxiliary components, the fetch plan, the payload recipients of the
owned boxes, validate_graph and the hash: all equal.  Cost per rank:
8 bricks: parts 0.014 s + assembly 0.002 s (global 0.036 s); 64 bricks:
0.045 + 0.006 s (0.34 s); 512 bricks: 0.061 + 0.107 s (3.48 + 0.32 s of
validation); the stream is 938 int64 per brick (5.3 MB gathered at 512).
Validation on 8 ranks (frozen_m10p, build/gpu_prep/ca3/val_m5e.sh): the
GPU CA3 leaf and levels 5+4 give the previous numbers to every digit
(logdet -50156606.410895944 / -50156606.734330446, the solve checks),
CA0 and Color at 192 unchanged, and the CPU CA3 (arrival-driven, 7
batches) reproduces its reference cpu8_ca3.log exactly (logdet
-50156606.410895959, acc_mvp 1.8095410e-03, acc_forward 4.7232536e-03;
Factor time 56.2 s against 56.7).  The leaf's schedule setup 33-45 ms
on 8 ranks (46-52 before: the neighbourhood is 27 bricks there, against
the grid's 8).  Then the parts step over dense arrays of the
neighbourhood's cells (colour and part index per cell, neighbours by
coordinate offsets, no hash maps or neighbour vectors): 6 / 17 / 34 ms
per rank on the 2^3 / 4^3 / 8^3 grids (14 / 45 / 61 before), the test
unchanged; what remains at 512 bricks is the assembly and the owner
assignment, 0.11 s together, of which assign_owners (the LPT with its
polish over 4401 components and 512 ranks) is about 0.09 s.  The dense
version (frozen_m10q) validated as m10p on 8 ranks: GPU leaf, levels
5+4, CA0/Color at 192 and the CPU CA3 all unchanged; the GPU leaf's
schedule setup 29-39 ms in situ (the setup also validates the local
graph, partitions the waves and derives the auxiliary components: timed
next).

**M5e at scale (2026-10-03; 320^3/64 job 59266621 with frozen_m10q,
640^3/512 job 59266181 with frozen_m10p; Color and CA0 of the same
builds in the same jobs; every check clean, the numbers as before).**

| | Color | CA0 | CA3 before M5e | CA3 with M5e |
|---|---|---|---|---|
| 320^3/64 Factor time | 6.05 s | 6.09 s | 6.70-6.81 s | 6.41-6.47 s |
| 320^3/64 leaf: total / elimination | 1.93 / 1.52-1.61 s | 1.88 / 1.34-1.40 s | 2.64 / 1.66-1.73 s | 2.23-2.30 / 1.58-1.75 s |
| 640^3/512 Factor time | 8.79 s | 8.50 s | 12.6-12.8 s | 8.86-9.05 s |
| 640^3/512 leaf: total / elimination | 2.12 / 1.59-1.70 s | 2.16 / 1.44-1.53 s | 6.4 / 1.86-1.99 s | 2.64-2.68 / 1.70-1.86 s |

The leaf's host phases with M5e: schedule setup 63-139 / slowest
118-164 ms at 64 ranks (the frozen_m10p setup at 512: 232-252 / 287-400
ms), start to elimination 169-244 / 271-313 ms at 64 and 403-420 /
528-615 ms at 512 (the setup, the halo gather 71 / 131 ms, the rest),
after the transition 81-90 / 106-133 ms (Color 62-71 ms: the
transition's imbalance).  The setup's split on 8 ranks (frozen_m10r,
`[owner-schedule] level L setup` under BPACK_TRACE=phase): 34.5 ms =
parts 4.6, gather 9.3 (the Allgather waits for the slowest rank),
assemble 1.5, validate 4.3, assign 0.1, waves 4.6, auxiliary components
and runtime 9.6, filters 0.1; at 512 ranks the assembly (~20 ms), the
owner assignment (~90 ms: the LPT's polish over 4401 components) and the
5 MB gather come on top.  CA3 is now within 0.1-0.4 s of Color at both
sizes: the leaf's remaining excess is the setup, the halo gather and
the level-end wait (0.36-0.51 s against Color's 0.04-0.10 s: the
stragglers of the colour phases).

**M5f: the setup's trims (2026-10-03).**  (1) The owner assignment's
polish (`detail::polish_assignment`) indexes every rank's components in
`ids` order, so the search for a swap or exchange partner on a rank runs
over that rank's components instead of all of them: the same moves in
the same order (the test compares the owners with
`polish_assignment_reference`, kept for it, for both policies), 17 ms
instead of 101 at 512 bricks.  (2) The wave partition covers the
components whose box lists are complete on the rank (those touching its
neighbourhood entirely: every box the rank can see belongs to one; a
first version with the owned and installed components only broke the
CPU owner solve, which reads the sequence numbers of halo boxes), with
num_waves fixed at the parity bound (8 in 3D) so the sequence numbers
agree across ranks whatever subset a rank partitions; the device plan's
`seq_of_box` holds the boxes with a wave.  (3)
`auxiliary_components` scans the brick's two outer layers only (a
deeper box sees, within two hops, this brick's own boundary, whose
components the rank shares).  On 8 ranks the setup is 20.5 ms (parts
4.7, gather 3.1, assemble 1.5, validate 3.9, waves 1.3, auxiliary and
runtime 5.6).  320^3/64 (job 59267460, frozen_m10s): CA3 Factor time
6.17-6.33 s against Color 5.99 and CA0 6.02 s in the same job (6.41-6.47
with frozen_m10q); the leaf 2.11-2.20 s against 1.88 / 1.89 s, its
setup 47-107 ms of which the gather's wait for the slowest rank 13-73
ms (the local work ~20 ms), the level-end wait 0.24-0.40 s (CA0's
0.43 s: the owner level's per-rank duration spreads from 0.40 to 0.91 s
over the 64 ranks, the interior phase alone from 0.12 to 0.37 s).  The
numbers unchanged (logdet -446707783.22014749, the solve checks).
frozen_m10t (the wave partition over every complete component) passes
the 8-rank set: GPU leaf and levels 5+4, CA0 and Color at 192, and the
CPU CA3 exactly its reference (logdet -50156606.410895959); its setup
21.7 ms on 8 ranks (waves 4.8 ms).

The timeline of the 320^3/64 wave run says where the owner level's
imbalance comes from: the corner ranks of the 4^3 grid end the interior
phase in 0.12-0.13 s, the face and interior ranks in 0.31-0.42 s, with
the same 2744 interior boxes each; the installs of the last colour are
still running on the device when the interior starts, and a rank with
26 neighbours installs more components and holds more pairs of
interest than a corner rank with 7.  The owner assignment balances the
owned boxes only; a per-rank pre-load from the install work (the
`fixed_load` of assign_owners, from the halo's size) would move owned
work off the ranks with many neighbours.  Not started (M5g).

640^3/512 with the trims (job 59267797, frozen_m10s): CA3 Factor time
8.65-8.96 s against Color 8.68 and CA0 8.39 s in the same job (8.86-9.05
with M5e alone, 12.6-12.8 before); the leaf 2.52-2.63 s against 2.04 /
2.07 s, its setup 142-223 ms of which the gather's wait 81-163 ms, the
assembly 12 ms, the assignment 13 ms (90 before the indexed polish); the
level-end wait 0.50-0.54 s (CA0's 0.54 s).  Checks clean (logdet
-3937186080.0625076, the solve checks 8.1e-15 / 9.9e-15).

**M5g: the install-aware owner assignment (2026-10-03, user go-ahead).**
The 320^3/64 wave run's timeline, phase by phase against the per-rank
install counts (the `entries` of the exchange trace): the green waves'
host time grows by 0.33 ms per purple generator the rank installed (the
device is still applying the purple colour's installs when the green
waves start), against about 0.17 ms per owned box of the colour's own
elimination; the interior phase is 0.47-0.52 s on every rank.  The
level's per-rank time is thus e * owned + i * installed with i about
2 e, and the installs are tied to the geometry: a corner rank of the
4^3 grid installs 1061 boxes, an interior rank 1958 (their owned boxes
being balanced at ~1194 each), which is the 1.07 against 1.57 s.
`assign_owners(..., install_cost)`: a colour's load of rank p is
install_cost * (boxes of the colour's components adjacent to p) +
(1 - install_cost) * (boxes p owns); with install_cost > 1 owning
relieves, so the greedy gives a component to its most loaded sharer,
the polish works with the signed item cost (both variants take a
`scale`), and the makespan relief runs in the other direction
(`relieve_by_taking`: the rank at the maximum takes a component from a
co-sharer, chains as before).  owner_schedule_setup passes
kDeviceInstallCost = 2.0 for a device level (sync_colors) and 0 for the
CPU schedule, whose numbers stay exact.  The test checks the indexed
polish against the reference and the local graph against the global
one with install_cost 2 as well, and prints the modelled level load
(install = 2 owned): on the 4^3 grid min/max 3566/7032 -> 4538/6128
(the maximum down 13%), on the 8^3 grid 3534/7032 -> 4250/6912 (2%:
the 216 interior ranks' install work alone, 2 * 3152 boxes, is near
the previous maximum, and they can only take components from each
other).  frozen_m10u passes the 8-rank set with the numbers of
frozen_m10t (on the 2^3 grid every brick is a corner, so the model
changes no owner; the CPU CA3 exact).  320^3/64 (job 59271948): WORSE.
CA3 Factor time 6.57-6.66 s against 6.17-6.33 with frozen_m10s (Color
6.04, CA0 6.04 in the job); the leaf 2.47-2.51 s against 2.11-2.20: the
corner rank 0 now owns 538 boxes (16 blue in 7 waves, 90 purple in 6,
432 green in 4) instead of 1194 and waits 0.87-0.93 s in the colour
transports, and the per-rank end of the owner level moved from
0.40 / 0.61 / 0.91 s (min / median / max over the ranks) to 0.57 /
0.95 / 1.20 s: the ranks that took the work got slower by more than the
corners gained, so the model's install cost of 2 owned boxes is wrong
(the regression confounded the installs with the position: an interior
rank's owned boxes have 26 neighbours and cost more than a corner
rank's domain-boundary boxes, which is where its waves' time went) and
tiny waves cost the device their launch latency.  The mechanism stays
(default 0, tested).  The regression over the two assignments (the
per-rank end of the owner level in the wave runs of frozen_m10s and
frozen_m10u against each rank's owned and installed boxes, from the
assignment replayed on the login node): end = 0.085 s + 0.258 ms * owned
+ 0.179 ms * installed (rms 0.16 s; the neighbour-weighted counts fit no
better), i.e. an install costs 0.7 of an owned box, not 2: with the
boxes-only assignment the slowest ranks (0.83-0.91 s) own 1352-1496 and
install 2056-2840 boxes, the corners (0.40 s) own 1194 and install 1028;
with install_cost 2 the slowest (1.12-1.20 s) own 2376-2544.  Modelled
with 0.7 the maximum load drops about 12% at 64 ranks (an interior
rank's installs alone, 0.7 * 2840, exceed a corner's whole load).  Tried
as frozen_m10v (kDeviceInstallCost = 0.7) at 320^3/64 (job 59273058):
also WORSE.  The corner rank 0 now owns every component it shares (2538
boxes in 26 components, installs 0; waves 1.11-1.14 s), the per-rank
end of the owner level 0.56 / 0.64 / 0.99 s (min / median / max)
against 0.40 / 0.61 / 0.91 with the boxes-only balance, the leaf
2.26-2.27 s against 2.11-2.20, Factor time 6.42-6.55 s against
6.17-6.33 (Color 6.10, CA0 6.06 in the job).  Both directions lose: the
per-rank time is not the linear model's (an owned ghost box costs more
than an owned local one: its entry blocks fetched, its FULL state
shipped home, its pairs with the owner's other boxes created; and a
rank's time is not conserved under a transfer), and the 0.16 s rms of
the fit was the warning.  Outcome of M5g: the install cost model stays
in assign_owners (default 0, tested with the reference polish and on
the local graph) and the device level keeps the boxes-only assignment
(kDeviceInstallCost = 0.0); the 640^3 job of install_cost 2 was
cancelled.  The code was then removed at the user's decision
(2026-10-03 evening: assign_owners lost `install_cost`, the polish
functions their `scale`, `relieve_by_taking` went; the indexed polish,
its reference and the per-sharer home counts stay); this record stays.
Validated (frozen_m10y): the 8-rank set as before, and the CPU CA3 at
320^3/64 twice in one job (59306221): logdet -446707783.22014755 and
-446707783.22014749 from the same binary, every H2_CheckError value,
compression ratio and skeleton size identical, the in-situ graph
comparison clean: the 17th-digit spread of the logdet is its final
cross-rank sum, run to run.  The leaf's imbalance (corners 0.4 s, interior ranks 0.9 s)
is the install work tied to the halo's size, and would need the
install itself to get cheaper (apply_remote_generators: 0.18-0.33 ms
per installed generator box) rather than a different owner.

**CPU at 640^3 on 512 ranks (2026-10-03, job 59288977,
build/gpu_prep/ca3/cpu640; 4 ranks x 16 threads per node, the live
build/).**  Color: Factor time 156 s, solve 1.54 s; CA0: 156 s, solve
2.24 s (logdets as the GPU runs').  Both CA3 runs (arrival-driven and
the serial oracle) fail at the leaf with "owner_rebuild_shipped_source:
partner copy missing: box X neighbour G" on many ranks (238 and 12
messages): the receiver, installing X's generator, rebuilds the pair
(X, G) from G's copy, which the owner left out of the payload
(BLOCK_PARTNER) because its `interest_of(receiver, G)` said the
receiver holds no interest in G, while the receiver's interest set says
it does (G in a component it owns).  In every case X and G lie in the
same brick, adjacent to the receiver's, one layer from a brick edge.
The CPU CA3 had only been validated after M5e on 8 ranks (2^3 grid: no
truncated component anywhere); offline, on the global graph, every
one- and two-hop neighbour of every owned box lies inside its owner's
27-brick neighbourhood (scratch owner_reach), and the outer-layer
auxiliary search equals the full scan (scratch aux_check), so the
owner's local graph should know G.  Diagnosis in situ: BPACK_CHECK=
graphcheck (the local graph compared with the whole grid's build at
every owner level, `compare_component_graphs`) and BPACK_CHECK=
globalgraph (the whole grid's build used), plus a guard in
`interest_of` (a boundary box without a component throws), run at
320^3/64 (job 59290622, cpu640/cpu_graph320.sh): the failure
reproduces at 64 ranks, identically with the local graph, with the
in-situ comparison (which finds no difference between the graphs: the
run reaches the same failure) and with the whole grid's build.  So M5e
is not the cause: the CPU CA3 is broken on grids with non-adjacent
bricks by something earlier in this session (the 8-rank reference
passes because every brick of the 2^3 grid touches every other).  The
last known good CPU CA3 on a 4^3 grid is the M5 feasibility of
2026-09-30 (192^3 on 64 ranks, Nmin 27: 16^3 bricks at level 6), before
M7-M9 (the geometric halo plan, the halo send plan with the tree, the
staged halo's pruning "only when blocks exist") and the M5b-M5g host
changes.  Next: the pre-session commit dc52665 built as a CPU reference
in the worktree ../bpack_ref_dc52665 and run at 320^3/64
(cpu640/cpu_ref320.sh); if it passes, bisect over 54d710a (M7), 040fa78
(M8), 4c8892a (M9a-d), 7bdce9f (M9e) and the uncommitted state.

Bisection (job 59291635, CPU drivers of the commits built in the
worktrees ../bpack_ref_<commit>): dc52665, 54d710a, 040fa78 and 7bdce9f
all pass (Factor time 96.7-96.8 s, logdet -446707783.22014749, the GPU
CA3's to every digit); the live tree fails.  Splitting the uncommitted
changes on HEAD (jobs 59293143, 59294401, 59294791): the live tree.hpp,
butterfly_factorization.hpp and butterfly_solve.hpp pass; the live
dataflow.hpp alone (HEAD's schedule header, hence the whole grid's
build) passes; the live owner_schedule.hpp alone, with its local-graph
calls replaced by the whole grid's build, passes, with HEAD's wave
count too; the live schedule and graph headers together (the local
construction in use) fail.  And the earlier "globalgraph" and
"graphcheck" runs were plain local-graph runs: the env parser rejects
unknown BPACK_CHECK values ("BPACK: BPACK_CHECK: unknown value ...
ignored"), which the logs say and I had not read.  So the local graph
is the suspect again, in something the owner_graph_test does not
compare (it finds every brick's local graph equal to the global one on
this very grid, down to the fetch plans and payload recipients; the
test under AddressSanitizer and UBSan is clean).  With the two values
registered in bpack_env.hpp (job 59295572): the local graph fails, the
in-situ comparison finds no difference and then fails, the whole grid's
build PASSES (Factor time 96.6 s, logdet -446707783.22014749).  The
failing pairs decoded on the login node (scratch failing_pairs): X is a
blue box whose owner has G, a purple box of a component owned by the
receiver, in its own local graph with the right component and owner,
yet the owner's `interest_of(receiver, G)` is false on its local graph
and true on the global one.  The cause: the assignment's CUMULATIVE
policy (the CPU's on a fine level: not coarse) breaks ties by the
number of a component's boxes homed at a sharer, counted from
`k.boxes`, which on a local graph is partial for every component away
from the rank; the owners of those components then differ from rank to
rank (each rank's `interest_of` reads its own graph's owners), the
block for the receiver's pair is not shipped, and the receiver finds no
partner copy.  The device level uses the per-colour policy, which never
reads the box lists, and the in-situ comparison and the test compared
the owners under the per-colour policy only; the 2^3 grid has no far
component at all.  Fix: `Component::home_counts` (boxes per sharer,
aligned with `sharers`), filled by both constructions (the parts'
counts per brick) and read by the tie-break (`Component::home_boxes`);
the comparison checks it and the owners under both policies (the test
too).  Test clean.  With the fix the CPU CA3 at 640^3/512 (job
59296059) runs: Factor time 146 s against Color 156 and CA0 156 s, the
solve 1.79 s against 1.54 and 2.24 s, logdet -3937186080.0625076, the
GPU CA3's to every digit (the CPU and the GPU CA3 agree bitwise on the
logdet at 320^3/64 and 640^3/512); its leaf level 18.4 s (elimination
10.3-13.0 s) against Color's 20.8 (16.1-16.2) and CA0's 23.0 s
(11.9-18.0).  The four CPU variants at 640^3 on 512 ranks (4 x 16
threads per node, jobs 59288977 and 59296059):

| CPU 640^3/512 | Factor time | solve | leaf level | logdet |
|---|---|---|---|---|
| Color | 156 s | 1.54 s | 20.8 s | -3937186064.8561087 |
| CA0 | 156 s | 2.24 s | 23.0 s | -3937186089.4963474 |
| CA3 (arrival-driven) | 146 s | 1.79 s | 18.4 s | -3937186080.0625076 |
| CA3 serial oracle | 207 s | 1.78 s | 79.2 s | -3937186080.0625076 |

The serial oracle gives the arrival-driven run's bits (one component per
batch, 79 s for the leaf).  frozen_m10x (the GPU build with the fix)
passes the 8-rank set with the numbers of frozen_m10w, the CPU CA3 at
160^3 exact.

**Where the time goes at 640^3/512, and what CA can win (2026-10-03).**
The factor-time differences between the variants outside the leaf are
run-to-run noise on identical Color levels (GPU, same job: up to 0.26 s
on a level, e.g. level 3 1.34-1.60 s; CPU, two jobs: level 6 22.6-28.0 s,
its rank-0 communication waits 6.7-11.0 s), so only the leaf is
attributable: GPU leaf Color 2.04, CA0 2.06, CA3 2.52-2.63 s (the "CA0
3% faster" of the totals is noise); CPU leaf Color 20.8 (elimination
16.2), CA0 23.0 (11.9-18.0: a 6.2 s wait for the interior ranks' ghost
re-elimination), CA3 18.4 s (10.3-13.0).  Color's neighbour waits on
rank 0 ("sizes" on the GPU, "Comm time" on the CPU) by level:

| level | 7 (leaf) | 6 | 5 | 4 | 3 | 2 | total | factor time |
|---|---|---|---|---|---|---|---|---|
| GPU | 0.59 | 0.88 | 0.75 | 0.64 | 1.01 | 0.22 | 4.1 s | 8.7 s |
| CPU | 5.0 | 9.9 | 9.8 | 11.8 | 17.1 | 8.8 | 62 s | 156 s |

The leaf holds a seventh of them; a leaf-only CA can win at most that.
And rank 0 is a corner rank: its waits are mostly the interior ranks'
larger per-wave work (more neighbours, more remote generators to apply),
which CA does not remove but moves (CA0: ghost re-elimination; CA3:
installs, 1101 generators on rank 0 against Color's 817, 2.81 against
2.11 GB).  How much of the 4.1 s / 62 s is latency (removable by CA)
and how much work skew (removable only by cheaper neighbour work) is
the measurement to make next (a Color run with per-rank exchange
traces).

**The measurement: BPACK_TRACE=skew (2026-10-04;
butterfly_factorization.hpp `LevelSkewProbe`, `level_skew_report`;
doc/environment_variables.md).**  Per rank and level: the elimination's
wall time and the part waiting for or receiving neighbour data (the
transports' sizes + payload timers, serialization excluded), per wave of
the Color loop too; gathered at the level's end on the level's first
rank, which replays the Color waves with zero latency (a wave that
begins with a transport starts when the rank's own and its neighbours'
previous waves are done, any other when its own is; each takes the
rank's measured busy time) and prints compute (the average busy time),
work skew (the replay's excess over it: ranks with more work in a wave,
and the chains of waits they cause) and latency (the measured time's
excess over the replay), plus averages by number of neighbour ranks;
per-rank CSV per level.  Smoke test on 8 ranks at 160^3 (all corner
bricks: no geometric skew; build/gpu_prep/skew/smoke8.sh, frozen_m10z
and frozen_cpu_skew):

| 160^3/8 | leaf: compute / skew / latency | levels 4, 3, 2: latency |
|---|---|---|
| GPU | 882 / 150 / 27 ms | 117, 101, 71 ms |
| CPU | 11.6 / 1.35 / 0.61 s | 1.0, 0.65, 0.53 s |

The 640^3/512 jobs were cancelled at the user's request (128 nodes);
instead 320^3 on 64 ranks (a 4^3 process grid: ranks with 7, 11, 17
and 26 neighbours), GPU Color twice, CA0, CA3 and CPU Color, CA3, CA0
in one 16-node job (skew/s320.sh, job 59320237).

**320^3/64 results (job 59320237, all in one job).**  CPU Color, by level
(compute / work skew / latency): leaf 13.6 / 1.2 / 0.6 s, level 5
18.4 / 2.6 / 0.9, level 4 15.8 / 3.2 / 0.9, level 3 9.5 / 3.7 / 0.6,
level 2 8.2 / 4.3 / 1.0 s: on the CPU the waits are work skew, four
times the latency, and geometric (busy time by neighbour count, level
3: 7.1 s at 7 neighbours, 8.4 at 11, 10.2 at 17, 12.7 s at 26; data
received 1.4 to 6.3 GB).  GPU Color (two runs, first replay without the
send-ready time): leaf 1012-1076 / 223-326 / 274-322 ms, levels 5-2
compute 640 / 405 / 197 / 182, skew 205-211 / 121-150 / 116-118 / 72,
"latency" 476-481 / 417-433 / 331-332 / 319-334 ms; the busy time
rises with the neighbour count there too (level 3: 129 to 319 ms), but
the replay put the last wave's device tail (asynchronous after the host
returns, synchronised at the next transport's start) after the waits:
the latency column is inflated; a send-ready time per transport was
added (`TransportTimers::presend`) and the GPU Color runs repeated.

| 320^3/64, one job | factor | solve | leaf level |
|---|---|---|---|
| CPU Color | 104 s | 1.31 s | 19.7 s |
| CPU CA3 | 95.0 s | 1.29 s | 16.8 s |
| CPU CA0 | 105 s | 1.60 s | 21.9 s |
| GPU Color (2 runs) | 6.03-6.22 s | 0.354 s | 1.89-1.91 s |
| GPU CA0 | 6.06 s | 0.364 s | 1.97 s |
| GPU CA3 | 6.24 s | 0.365 s | 2.20 s |

With the send-ready time in the replay (job 59320536, frozen_m10za and
frozen_cpu_skew2; the send-ready time is small: 20-120 ms per level on
the GPU) the split, summed over levels 6-2 of Color at 320^3/64:

| Color 320^3/64 | compute | work skew | the rest |
|---|---|---|---|
| CPU | 65.3 s | 15.0 s | 4.0 s |
| GPU | 2.46 s | 0.82 s | 1.78 s |

"The rest" is not message latency but the transfer of the neighbours'
data: per level it equals the most-loaded ranks' received bytes at a
near-constant rate (GPU 11-14 GB/s at levels 6-3, 5 GB/s at level 2 with
8 ranks; CPU 9-12 GB/s).  So the user's reading holds on the CPU (the
waits are work skew, nearly four times the transfers); on the GPU the
transfers are the larger part, the skew second, and both grow with the
number of neighbour ranks (level 3: busy 129 to 319 ms and 0.8 to 3.6 GB
received from 7 to 26 neighbours).  CA removes neither: CA0 replaces
the transfers by re-eliminated ghosts (CPU leaf busy 11.9 s at 7
neighbours to 17.3 s at 26), CA3 by installs (10.8 to 12.2 s), and
CA3 moves more data than Color (rank 0: 2.81 against 2.11 GB).

CPU CA3's level 5, a Color level, takes 22.0 s against Color's 25.9: not
the communication (the split is the same) but 0.18M minor page faults
per rank against 1.17M after the Color leaf (0.81M after CA0), about
3 s of "staged" wall time: how the leaf leaves the host memory, with the
runs' TCMALLOC_RELEASE_RATE=1 (the same pattern in the 640^3 jobs:
level 6 22.6 against 27.5 s).  CPU CA0's leaf: busy 11.9 s at 7
neighbours to 17.3 s at 26 (the ghosts re-eliminated).



**M5 feasibility (2026-09-30, runs only; build/gpu_prep/ca3feas).** CPU
Color / CA0 / CA3 (`--H2_CA_owner_component 3 --H2_CA_staged_halo 2
--H2_CA_owner_serial 0`, arrival driven), Laplace 192, CA on the two finest
levels, CPU build on 40 GB GPU nodes (4 ranks x 16 threads per node). CA3's
serial oracle gives the same bits as the arrival-driven run; accuracy as
Color's (acc_forward 4.35e-3 vs 4.11e-3 at Nmin 216). Factorization: 64
ranks Nmin 216 43.0 / 48.5 / 37.4 s, Nmin 27 40.4 / 41.2 / 38.9 s, 8 ranks
117 / 115 / 103 s. The Color levels (identical in all three) vary by up to
1.8 s between the runs (the Color run first in its job was slowest), so per
CA level (rank 0's level times):

| CPU, CA levels | Color | CA0 | CA3 |
|---|---|---|---|
| 64 ranks Nmin 216: leaf (512 boxes/rank) | 11.63 s (comm 4.23) | 19.18 | 10.60 |
| 64 ranks Nmin 27: levels 6 (4096) + 5 (512) | 2.81 + 6.81 | 3.12 + 11.36 | 3.05 + 5.97 |
| 8 ranks: levels 5 (4096) + 4 (512) | 54.69 + 30.68 | 53.13 + 34.77 | 48.69 + 27.23 |

CA3 is 6-11% faster than Color on its levels and far faster than CA0 on
8^3 bricks (no replication). Its leaf at 64 ranks Nmin 216 (rank 0): box
region 2.07 s, passes 1.56, R2 replay 0.46, final 0.27, pack 0.10, install
0.68, idle wait 0.63 s; 1.6-1.8 GB sent/received per rank (Color's GPU leaf
receives 1.69 GB): about 3.1 s of CA3-only work (install, block rebuild,
replay passes) in a 10.6 s level. On the CPU, Color's 4.2 s of "comm" at
the leaf is mostly waiting for slower neighbors color by color, which CA3's
arrival-driven components avoid.

GPU projection (from the GPU Color levels, 40 GB nodes): the GPU leaf at 64
ranks Nmin 216 is 1.24 s (box work ~0.64, exchanges 0.50, transition 0.10);
CA3 would keep the box work, add its install/rebuild/replay on the device
(the CPU's ~30% extra over box work -> ~0.2-0.3 s) and replace 9 exchange
rounds by ~4 stages of about the same bytes (~0.25-0.3 s): roughly even
with Color. At 884K per rank (8 ranks; 384^3 on 64 ranks with 80 GB), the
16^3 leaf's exchanges are larger (1.3 s at 8 ranks) and CA0 already wins
there; CA3 might take a further ~0.3-0.7 s off the leaf and avoid CA0's
8^3-level loss. The CA3-only work runs in small component batches (15
components, 294 boxes at the 64-rank leaf), the GPU's weak spot.

### M6. Hybrid level selection
Choose CA vs Color per level from brick size and measured costs instead of
the fixed `CA_level` / 343-box rule.

### M7. After the 640^3 run: compaction, geometric halo plan, device process reduction
Three items the 640^3 and 768^3 runs on 512 ranks asked for (2026-10-01),
all shared by CA0 and Color.

**M7a. Compaction of the kept solve factors (done 2026-10-02).** At the end
of every level's transition (`ColorGpuDriver::compact_device_memory`, after
the level's own blocks are freed), the blocks that outlive the level, the
solve factors every level kept (one block per wave) and the parents' blocks
of a device transition, slide down the heap (`DeviceHeap::compact`): in
address order, each to the lowest free range that holds it, else down to the
start of the free range right below it (an overlapping move, done through a
256 MB staging block from the top of the heap, or in pieces of its gap).
The pointers into them (`KeptSolveBox`, the parent `DeviceMatrixT`s) are
host-side and follow (`DeviceHeap::relocated`). Nothing is reserved or
estimated: blocks move at their real size after they exist (the alternative,
reserving the keep bound up front, would have held back ~12 GB at the leaf of
320^3). Rank-local, no decision shared; the moves run on the stream after the
level's device work (all device work is on one stream, so no synchronization
is needed beyond the one that times it). With `--verbosity 1` a level that
moved blocks prints `[gpu] level L compaction: N blocks moved (X GB) in T s;
F GB free, largest range R GB`. Validation (2026-10-02, job 59190541,
build/gpu_prep/compact): the 8-rank Laplace 192 regression, CA0 (CA on
levels 5 and 4) with the solve and replica checks and Color, gives the same
logdet and H2_CheckError values as the clean build (../clean: CA0
-89253589.199881583, acc_mvp 8.6255069e-04; Color -89253588.194065362,
8.1095537e-04), replica check 0 differ; per level the compaction moved 7099
/ 861-879 / 111-128 / 23-40 blocks (6.7-7.4 GB) in 33 / 7-16 / 3-11 / 2-18
ms, after which the free memory is one range (26-32 GB). Factor times
unchanged (10.29 / 8.39-8.42 s against 10.34 / 8.40 s).

**640^3 on 512 ranks, 40 GB nodes, with the compaction (2026-10-02, job
59190670, build/gpu_prep/ca0/g640_40c; CA0 twice, Color once, CA0 with the
checks).** CA0 runs through: the kept factors (13.0-13.7 GB per rank at the
leaf and the levels above) are moved in 30-40 ms per level, leaving one free
range of 19.7-20.4 GB, and the level 4 -> 3 reduction allocates without
trouble. CA0 11.86 s (first run, cold) and 11.32 s; Color 11.63 s in the same
job (11.3 s the day before); replica check 1217216 ghost copies, 0 differ;
CA0 logdet -3937186089.4963474, acc_mvp 9.32e-4, acc_forward 6.80e-3; Color
-3937186064.8561087, 9.43e-4, 6.61e-3. Per level (Color / CA0 second run):
leaf 2.51 / 2.58 s (elimination 1.89-2.00 / 1.64-1.74, CA0's gather 0.26 and
transition 0.38 + 0.04 compaction against 0.37 + 0.04), level 6 1.93 / 2.03,
level 5 1.44 / 1.43, level 4 1.43 / 1.40 (its reduction 0.4 / 0.33 s), level 3
2.82 / 2.46, level 2 1.20 / 1.08 s. So at 512K per rank on 512 ranks CA0 and
Color are within noise of each other; what CA0 has over Color at the leaf
(0.2-0.3 s of elimination) it pays in the gather (0.26 s, M7b) and in the
level above (0.1 s). The levels after the reductions (3 and 2) vary by 0.3 s
between runs regardless of the method.

**M7b. Geometric halo plan for the CA level-entry gather (implemented
2026-10-02; validation below).** The
host gather `gather_CA_factorization_data` costs 225-264 ms at 512 ranks
(14 ms at 8): five collectives over the level communicator (an Allgather of
the rank map, Alltoalls of request counts and sizes, Alltoallvs of request
lists, twice for the ghosts and the assisting boxes) and two request-reply
rounds, growing with the rank count rather than the bytes. The ghost set of
a process is a function of the Morton layout alone (`compute_ghost_boxes_and_colors`:
the one-hop neighbors of its region's boundary boxes, their two-hop
neighbors on the faces of their owners, and the orange completions of the
blue ghosts), and so is its assisting set (one and two hops around its
ghosts, not local and not ghosts). Plan: factor both into pure functions of
(dimension, level, active processes, process id), used by the tree
construction for the process itself and by every sender for its neighbors;
each rank then computes which of its own boxes each neighboring region (within
ceil(5 / region width) regions: ghosts lie within 3 box hops of a region,
assisting boxes within 5) needs as ghosts and as assisting boxes, and pushes
both in one message per neighbor over `tree->comm` (a sizes round, then the
payloads; no collective, no request round). The symmetric-pair pruning keeps
working: the sender has the receiver's full ghost set from the same function.
The receiver checks that every ghost and assisting box arrived and that every
owner it expects lies within the candidate radius (a clean error instead of a
hang if the radius reasoning were wrong). The assisting plan kept for the
post-elimination skeleton refresh (`refresh_CA_assisting_skeletons`) becomes
per-peer lists by global rank. The device halo (M2, ca_halo.hpp) takes the
same lists in place of `make_CA_request_plan` and the Alltoall size check.
Expected: tens of ms at 512 ranks, about 0.2 s off the CA leaf of 640^3.

Implementation: `h2_parallel/core/halo_geometry.hpp` (the box colors moved
there from tree_impl.hpp; `RegionLayout`; `ghosts_of_region`,
`ghosts_of_region_in`, `assisting_of_region_in` as closed forms over the
layout: the one-hop shell of the region's cube, the second shell's boxes on
faces of their regions, the orange neighbors of the blue ghosts (blue boxes
lie next to their region's corners), and the assisting boxes as the boxes
within three hops of the cube or two of the ghosts beyond the shell; and the
reference enumerations of the original code). `tests/halo_geometry_test.cpp`
compares the closed forms with the reference for every region of 3D levels
of width 2-16 and 2D levels of width 2-8 (all pass), and the tree's
`compute_ghost_boxes_and_colors` / `compute_assisting_boxes` compare them
again for the rank's own region at every run (a mismatch throws before the
factorization). `gather_CA_factorization_data` runs
`gather_CA_factorization_data_geometric` (one message per neighboring region
over `tree->comm`: a sizes round, then ghosts and assisting boxes together;
the receiver checks that every expected box arrived, in order) and falls
back to the collective form only for non-uniform box distributions; the
assisting plan is `TreeLevel::assisting_plan` (per peer rank); the device
halo (ca_halo.hpp) derives its lists from the same functions and checks the
message sizes point to point; its chunked rounds moved to
`device_pieces.hpp` (shared with M7c). The `--verbosity 1` gather line now
reads `(geometric plan P, serialize S, messages M, unpack U, symmetric halo
R)`. (The unstructured code has no CA levels: `color_unstructured is
Color-only`, so only the structured code is affected.)

Validation on 8 ranks (2026-10-02, job 59191890, build/gpu_prep/halo):
Laplace 192 with CA on levels 5 and 4 on the GPU (solve and replica checks:
9488 + 3536 ghost copies, 0 differ) and on the CPU build, and Color, give
the clean build's results (GPU CA0 logdet -89253589.199881583, Color
-89253588.194065362; CPU CA0 -89253589.199881539 as before). Gather times at
8 ranks: GPU leaf 47 ms (plan 9, serialize 5, messages 29, unpack 4; the old
form 25-40 ms with its first collective absorbing the ranks' skew, which now
shows in the messages), level 4 27 ms (old 14-15); the CPU level 4 gather,
whose ghosts carry their blocks, 3.9 s against 7.0 s. The point is the rank
count: see the 640^3 run below.

**M7c. Device process reduction for Color and CA0 (implemented 2026-10-02;
validation below).** At a process
reduction (512 -> 64 ranks at level 4 -> 3 of 640^3) the transition downloads
the parents' blocks, the host serializes whole boxes to the new owner, which
deserializes, merges and re-uploads them at the next level's start: 0.6 s at
640^3 for the 4 -> 3 reduction (pack 0.03, receive 0.18, unpack 0.23, send
0.34 s, max over ranks), 0.25 s at 3 -> 2, 0.06 s at 2 -> 1, plus the
upload. Plan: keep the blocks on the device across the reduction. In
`transition()`, when the next level has fewer active ranks and the parent
level would run on the device (`level_eliminator_would_run`, rank-independent;
for a CA parent level also every rank's parents within the device sketch's
box limit, all-reduced over the child level's ranks), each sender sends its
`DeviceLevelBlocks` (the Schur block of each of its parents, the near-field
blocks among them and its copies of those with other ranks' parents) to
`level.parent_level_owner`: a host header (kind, key, rows, cols per block),
then the payloads through the exchange arena in chunked rounds packed
straight from and unpacked straight into the blocks (the mechanism of
ca_halo.hpp, factored out), host-staged when MPI is not GPU-aware. A pair
whose two parents' children live on two senders of the same owner is sent
once (by the sender holding the lower parent); a pair across two owners goes
to both, as the parent level's layout has it. The owner allocates the
incoming blocks (resident), merges them into its own `DeviceLevelBlocks`, and
the next level's eliminator adopts them as after a non-reducing transition;
senders free theirs. The host reduction stays for the boxes' metadata (points,
index and neighbor lists), now without block payloads, so
`share_symmetric_level_edges` has nothing to pair. Level 1 then adopts the
blocks on its single rank and the root is assembled on the device. Expected:
most of the ~0.9 s of reductions at 640^3, and no level-start upload after a
reduction.

Implementation: `device_reduction.hpp` (`reduce_parent_blocks`: the
headers over `tree->comm`, tags 4800-4802, the payloads through
`exchange_device_pieces` of `device_pieces.hpp`, collective over the child
level's communicator for the chunk and the rounds), the decision
`ColorGpuDriver::device_reduction_possible` (the parent level would run on
the device: `level_eliminator_would_run` and the CA tests without this
rank's own activity there, uniform region layouts, every rank's parents
within the device sketch's box limit on a CA parent level; all-reduced over
the child level's ranks) and the call in `transition()`; the unstructured
loop now passes its level communicator too. The `--verbosity 1` transition
line says `blocks merged on this device (process reduction)` on an owner and
`sent to their owner's device` on a sender, followed by a `device reduction`
line with the times, bytes and blocks.

Validation on 8 ranks (2026-10-02, job 59192055, build/gpu_prep/reduce):
Laplace 192 Color and CA0 (CA on levels 5 and 4), solve and replica checks,
and VIE 96 unstructured: the reduction 8 -> 1 into level 1 moves 28 blocks
(the 8 level-1 Schur blocks and their 21 pairs, 0.06 GB; 0.15 GB on VIE) in
under 10 ms, the host reduction then carries only metadata ("process
reduction communication time: 0 ms", 11 ms before) and the root is assembled
and factored on the device; results as the clean build's (Color
-89253588.194065362, CA0 -89253589.199881583, 0 replicas differ).

A first 640^3 run (job 59192186) went through the 512 -> 64 reduction on
the device (0.51 GB in 89 blocks per owner in 0.02 s; level 4 1.40 -> 1.00
s) and then hung in the 64 -> 8 reduction: once ranks are idle, the tree
takes the parent level's owners from the idle ranks
(`compute_active_processes`, `use_inactive_owners`: "Active processes: 8
from rank: 4"), and such an owner, inactive at the child level, had no
eliminator and skipped the reduction while its senders waited for it. The
reduction now runs over the transition communicator (the child level's
ranks and the parent level's): every rank of it decides (idle owners as
receivers, with no blocks of their own) and takes part in the piece
exchange. The 8-rank runs could not show this (8 ranks, no idle rank at
the 2 -> 1 reduction); rerun with the fix (job 59194080) they give the same
results as before.

64 ranks, Laplace 256, Nmin 64 (2026-10-02, job 59194081,
build/gpu_prep/reduce; 16 nodes of 40 GB): the reductions 64 -> 8 at level
3 -> 2 (group owners: each receives 0.45 GB in 89 blocks, 0.01 s) and 8 ->
1 (the idle rank 4 as owner: 8 blocks of 0.02 GB per sender) both run on
the device; against the halo build's host reductions in the same job the
reduction communication time falls from 125 ms to 1 ms and the factor time
from 5.66 s to 5.53 / 5.22 s (Color) and 5.76 to 5.50 / 5.48 s (CA0); per
level (Color, device / host reductions): level 3 (ends with 64 -> 8) 0.68 /
0.86 s, level 2 (ends with 8 -> 1, and no longer re-uploads its blocks)
0.54 / 0.70 s, level 1 (adopts the blocks for the root) 1 / 42 ms, the
levels above alike (1.64 / 1.66, 1.22 / 1.24, 0.86 / 0.89 s); results the
same (Color logdet -221223299.04514426 in both; CA0
-221223299.92170054 against ...051 with the host reduction, the lower
parent's copy of a cross pair in place of the host's choice; H2_CheckError
values identical); the solve check passes. The CA leaf gather at 64 ranks:
107-152 ms (plan 17, serialize 3, messages 86-131), the old form 48-108 ms
at 192^3 on 64 ranks: the point-to-point messages wait for the slowest
rank like the first collective did; see the 640^3 numbers.

**640^3 on 512 ranks, 40 GB nodes, with all three (2026-10-02, job 59194082,
build/gpu_prep/ca0/g640_40d; CA0, Color, CA0, Color, CA0 with the checks).**

| | compaction only (job 59190670) | + geometric halo + device reduction |
|---|---|---|
| CA0 factor time | 11.86 (cold) / 11.32 s | 10.83 (cold) / 9.85 s |
| Color factor time | 11.63 s | 10.02 / 9.97 s |
| CA leaf gather | 259-283 ms (plan 91-102, request sets 31, payloads 41, assisting 97-109) | 94-133 ms (geometric plan 17, serialize 3, messages 70-109, unpack 4) |
| level 4 (ends with 512 -> 64) CA0 / Color | 1.40 / 1.43 s | 1.02 / 0.97 s |
| level 3 (ends with 64 -> 8; no upload at its start) | 2.46 / 2.82 s | 1.69 / 1.96 s |
| level 2 (ends with 8 -> 1) | 1.08 / 1.20 s | 0.92 / 0.88 s |
| level 1 (adopts the blocks for the root) | 59 ms | 2 ms |
| process reduction communication | 541-668 ms | 1 ms |

The leaf and levels 6 and 5 are unchanged (CA0 2.50 / 2.04 / 1.39 s, Color
2.53 / 1.92 / 1.43 s). Results identical to the compaction-only run (CA0
logdet -3937186089.4963474, acc_mvp 9.32e-4, acc_forward 6.80e-3; Color
-3937186064.8561087, 9.43e-4, 6.61e-3); replica check 1217216 ghost copies,
0 differ; the solve check passes. The three reductions on the device move
0.51 GB in 89 blocks per owner (512 -> 64), 8 blocks of 0.03 GB per sender
(64 -> 8) and 8 of 0.02 GB (8 -> 1), each in 0.01-0.02 s. With the host
reductions, their serialization, messages and the re-uploads at the
following levels cost 0.6 + 0.3 + 0.1 s plus the level-start uploads; the
levels after the reductions are now 0.4 + 0.8 + 0.2 s faster. The gather's
remaining 70-109 ms are the messages waiting for the slowest rank (plan 17
ms, serialization 3 ms), as the first collective did before. CA0 and Color
are within 1% of each other at this size (9.85 against 9.97 s); CA0's leaf
elimination advantage (1.64-1.74 against 1.89-2.00 s) is spent on the
gather and the level above the leaf, as before.

### M8. After M7: what stays on the host, and the follow-ups the user chose
At 640^3 on 512 ranks (M7c run) the numerical work is all on the device;
the host keeps the per-wave planning (0.3 + 0.12 s at the leaf, hidden
behind the device sketch), the leaf transition plan (0.23 s), the CA
level-entry gather (0.1 s) and, above all, the Color exchange protocol: per
level 0.5-0.8 s of "sizes" (the per-wave handshake, i.e. waiting for the
slowest peer) plus 0.15-0.66 s of payload, about 5.5 s of the 10 s. The
user chose (1) CA on level 6 as well (`--CA_level 6`, an experiment:
build/gpu_prep/ca0/g640_ca6), (3) a cheaper exchange protocol and (4) no
host copies of the kept boxes' factors plus a faster leaf transition plan;
not the reduction threshold.

**(1) CA on levels 7 and 6 at 640^3 on 512 ranks, 40 GB nodes (2026-10-02,
job 59202906): does not fit, and would not pay.** The leaf runs as before
(2.49-2.67 s, 11.94 GB kept); level 6's device halo (M2) then brings 4.18 GB
of ghost blocks per rank (6320 blocks, sent 4.60 GB) in 1.22-1.26 s (MPI
1.16-1.20 s: the arena-staged rounds run at ~4 GB/s here), and the level
exhausts the 33.4 GB heap on dozens of ranks (the leaf's kept factors, the
level's own blocks and the ghosts' leave 0.2 GB). Even with the memory, the
halo alone (1.26 s) costs more than the level-6 exchange it replaces (1.15
s). CA stays on the leaf at this size; the 512-rank CA7 run in the same job
gave 9.84 s (logdet -3937186089.4963474, as before).

**M8a. Exchange without the sizes handshake (implemented 2026-10-02).** A
rank now sends everything at once (its sizes, the host header with a
transport number, the device bulk) and receives each peer's header and
bulk as soon as that peer's sizes are in (`MPI_Waitany` over the size
receives, per-peer receive buffers), so no rank waits for all its peers
before sending and the transfers from early peers overlap the wait for
late ones. The first version probed the headers for their lengths
(`MPI_Improbe`/`MPI_Mrecv`) instead of a sizes round; with Cray MPICH it
matched a later transport's header while an earlier one was pending (the
transport number in the header caught it: "rank 6 at transport 4 received
the header of transport 5 from rank 0"), so the protocol keeps explicit
receives and the transport number stays as a check.

**M8b. Kept boxes without host copies of X_SR and the LU of X_RR
(implemented 2026-10-02).** As for X_NR before: the eliminator skips the
copies for the boxes it keeps (X_RS is dropped: only the host transport
read it), `materialize_host_factors` brings them back for a host solve or
the solve check, `rescue_skipped_xnr(true)` when the keep fails, and the
log-determinant of a kept level comes from the device
(`device_logdet_kept_level`, `launch_logdet_lu`: one block per box over the
LU's diagonal with the LAPACK pivots, as the host's `accumulate_logdet_lu`;
the host sums the levels and boxes that have a host LU). The device solve
takes a kept box's shapes from the kept factors. The 8-rank CA0 logdet came
out identical to the host's.

**M8c. Leaf transition plan (implemented 2026-10-02).** The lazy far
fill's sources of a child pair came from `collect_lazy_far_sources` (26
Morton decodes, two hash lookups and a sort per two-hop pair, ~275K pairs
at the leaf of 640^3); `build_parent` now lists each local child's
eliminated one-hop neighbors with their waves and coordinates once and
filters that list per pair.

Validation of M8a-c on 8 ranks (2026-10-02, job 59208723,
build/gpu_prep/m8, frozen_m8): Laplace 192 Color and CA0 (CA on levels 5
and 4) with the solve and replica checks and without, VIE 96 (complex, CA
leaf, adaptive ID rows) with the solve check. H2_CheckError values as the
clean build's in every run; replica checks 0 differ; solve checks 3e-15 to
3e-14; Color's logdet -89253588.194065422 against -89253588.194065362 (the
device sums the kept levels' boxes in another order), CA0's
-89253589.199881583 unchanged, VIE's equal to 11 digits. The device solve
runs as before (its setup had first refused the kept boxes' missing host
LU; the check now accepts a kept box, and a kept CA level the solve runs on
the host gets its factors back first). Factor times: Color 7.8-8.6 s
(8.4 before), CA0 8.3 s without checks (the leaf exchange on 8 ranks:
1.68 s [sizes 1.14, rest 0.36] for 6.8 GB received; the transition plan's
0.47 s restore dominates its 0.72 s).

**640^3 on 512 ranks, 40 GB nodes, with M8a-c (2026-10-02, job 59208905,
build/gpu_prep/ca0/g640_m8; CA0, Color, CA0, Color, CA0 with the checks).**

| | M7 build (job 59194082) | + M8a-c |
|---|---|---|
| CA0 factor time | 10.83 (cold) / 9.85 s | 9.91 (cold) / 9.13 s |
| Color factor time | 10.02 / 9.97 s | 9.19 / 9.17 s |
| level 6 exchange (Color) | 1.19 s [sizes 0.76, payload 0.31] | 1.04 s [sizes 0.91, rest 0.01] |
| level 5 / 4 exchange | 1.02 / 0.75 s | 0.83 / 0.70 s |
| level 3 exchange (64 ranks) | 1.55 s [0.82 + 0.66] | 1.23 s [1.12 + 0.04] |
| level 2 exchange (8 ranks) | 0.65 s [0.20 + 0.37] | 0.63 s [0.19 + 0.36] |
| leaf transition plan | 0.23-0.24 s | 0.17-0.19 s |

The payload phase of the exchanges has gone into the wait for the slowest
peer (which stays at 0.6-1.1 s per level, the cost of the wave schedule
itself); the 8-rank level 2 does not change (its transfers were already
waiting on nothing). Results identical to the M7 build's (CA0 logdet
-3937186089.4963474, the device sum reproducing the host's here; Color
-3937186064.8561087; H2_CheckError values unchanged); replica check
1217216 ghost copies, 0 differ; the device solve keeps levels 3-7 and sets
up in 0.2 s as before. CA0 and Color are 0.7-0.8 s faster each; CA0 is
again slightly ahead (9.13 against 9.17 s).

**320^3 on 64 ranks, 40 GB nodes, with M7 + M8 (2026-10-02, job 59220006,
build/gpu_prep/ca0/g320_m8; commit 040fa78 as frozen_m8, and the dc52665
build of 2026-10-01 in the same job).** The 64-rank case where CA0 first
beat Color (512K per rank, Nmin 64, CA on the leaf).

| | dc52665 (same job) | 040fa78 |
|---|---|---|
| Color factor time | 7.87 s | 6.74 / 6.72 s (6.69, 6.71 with the checks / the wave trace) |
| CA0 factor time | 7.43 s | 7.26 (cold, first run of the job) / 6.76 s (6.85 with the wave trace) |
| leaf: Color / CA0 | 2.46 / 2.27 s | 2.23 / 2.25 s |
| level 5 (above the leaf) | 1.74 / 1.86 s | 1.65 / 1.64 s |
| level 4 | 1.29 / 1.22 s | 1.12 / 1.15 s |
| level 3 (ends with 64 -> 8) | 1.20 / 0.97 s | 0.80 / 0.83 s |
| level 2 (ends with 8 -> 1) + level 1 | 0.91 / 0.85 s | 0.67 / 0.66 s |

Results identical to the old build's (CA0 logdet -446707784.16925561; Color
-446707781.18238521 against ...533, the device sum of the kept levels);
H2_CheckError values unchanged; replica check 124704 ghost copies, 0
differ; solve checks 5.5e-15 / 1.05e-14 for both. Color gained 1.15 s, CA0
0.67 s: the exchange and reduction work helps every Color level, and the
Color leaf has an exchange where the CA0 leaf has none, so CA0's lead at
this size is gone (a tie within the run-to-run spread: 6.72 against 6.76 s).
CA0's leaf elimination is still shorter (1.54-1.61 against 1.76-1.82 s) but
its gather (0.13-0.23 s; the messages wait for the slowest rank) and the
level's larger start take that back.

Wave traces of both (co64_wave, ca64_wave; times of the slowest rank of each
level, Color run): a Color level's elimination is its 16 (8) waves plus the
transports between them:

| level (boxes per rank) | waves | of which device sketch wait / plan + launch | between waves (transports) | per transport |
|---|---|---|---|---|
| 6, leaf (4096) | 1.03 s | 0.34 / 0.53 s | 0.47 s | 52 ms |
| 5 (512) | 0.60 s | 0.44 / 0.15 s | 0.76 s | 84 ms |
| 4 (64) | 0.38 s | 0.34 / 0.03 s | 0.61 s | 68 ms |
| 3 (8) | 0.24 s | 0.22 / 0.02 s | 0.41 s | 51 ms |
| 2 (8, 8 ranks) | 0.15 s | 0.11 / 0.03 s | 0.43 s | 53 ms |

A transport on the slowest rank is the end of its wave's device elimination
(10-25 ms), the bulk transfer (230-390 MB per transport and rank at 7-12
GB/s: 25-35 ms) and the remote updates; the other ranks wait for it (rank
0: 0.55-0.87 s between waves). So what is left is device latency on small
batches (levels 4-2: one to eight boxes per wave) and the transfers (10.4 GB
received per rank over the levels, about 1.0-1.4 s of the 6.7 s), not host
work, except at the leaf, whose host planning and launches take 0.41-0.53 s
of its 1.8 s elimination. The CA0 leaf (32 waves, no transports): waves
1.02-1.47 s over the ranks (the interior ranks with 6936 boxes against the
corner ranks' 5282: sketch wait 0.62, plan + launch 0.65 s on the slowest).
After a CA leaf the next level's start costs 80 ms more than in the Color
run (start 86 against 5 ms; also at 640^3): not yet looked into.

### M9. Leaf host work and the CA0 start-up costs (2026-10-02, after the 320^3 traces)
The user chose, from the 320^3 breakdown: the leaf's host planning, and the
CA0-only costs (the 80 ms start of the level after a CA leaf, the halo
gather's wait).  CA stays on the leaf level only.

**M9a. Page-locked buffers in the warm-up.** The staging buffers of
`pinned_pool()` (the copiers' ring slots, 64 MB each, the metadata images of
the waves and of their owner passes, the sketches' result buffer) were
page-locked at their first use, about a second per GB: the copier threads'
ring slots held up the leaf's start (0.1-0.2 s at 64 ranks: the level's CUDA
calls wait behind `cudaMallocHost`), and the exchange copier's again at the
first level with exchanges, which after a CA leaf is the next level (the 80
ms).  `ColorGpuDriver::warm_up_buffers` reserves them before the levels,
with the other warm-ups; the warm-up line reads `kernels X s, buffers B s,
exchanges Y s, evaluator Z s` (0.21-0.22 s for the buffers).  This moves a
one-time cost of the process out of the factor time; a first factorization's
wall time is the same.

**M9b. The halo gather's send plan from the tree.** At 64-512 ranks the
gather's "messages" (70-210 ms) were the corner ranks (7 neighboring
regions) waiting for the interior ranks (26) to plan: each region's
orange completion was computed three times per neighbor, 17 ms on rank 0
and about 60 ms inside.  The plan is now made once with the tree
(`compute_assisting_boxes`: `TreeLevel::ca_halo_sends`, in parallel over
the regions, each region's extra ghosts computed once, `RegionExtra`), the
gather packs its per-peer buffers in parallel, and the whole ghost sets
for the symmetric-pair pruning are computed only when a box carries blocks
to prune (none at the leaf).  The gather line prints the slowest rank's
plan + packing time.  halo_geometry_test passes with the refactored
functions.

**M9c. Transition plan.** The chunks' G and P blocks were found through
hash maps keyed by (source, child) Morton pairs, about four million lookups
at a 4096-box leaf; they are now flat tables indexed by the source's dense
state index and the direction from the source to the child.

**M9d. A wave's near blocks after the box region's launch.** Step 5's store
items and the allocation of their blocks (26 per box, half the box plan's
time) are planned after the box region is launched and go with the owner
pass's image; they may reuse the blocks the launch retired.

8 ranks (jobs 59224901 and 59230263, build/gpu_prep/m9; Laplace 192 Color
and CA0 on the leaf with the solve and replica checks; Laplace 160 Nmin 64,
the per-rank leaf of 320^3/64 and 640^3/512, with wave traces): results
bitwise as before (Color logdet -89253588.194065422; the CA0 leaf-only run
9488 ghost copies, 0 differ, solve checks 2.6e-15 / 1.6e-14).  At 160^3
(rank 0; committed build -> M9a+b -> all four): leaf start 0.1-0.2 s (64
ranks) -> 13 -> 6-13 ms; start of the level after the CA leaf 86 -> 6 ms;
gather 10-17 ms with no plan time; leaf transition plan 0.19-0.20 -> 0.13-0.14
s; box plan before the launch 0.14-0.15 -> 0.06-0.07 s (near blocks 0.07-0.08
s after it).  Leaf level 1.43-1.46 -> 1.37-1.38 s (Color), 1.44 -> 1.38-1.39
s (CA0); factor time 3.62-3.67 -> 3.57 s and 3.68-3.79 -> 3.59-3.62 s.  The
wave totals move little (Color 828 -> 794 ms, CA0 957 -> 951 ms): at the leaf
the waves are host-bound, the host's own work (plan 0.25-0.28 s, launches
0.11-0.13 s, sketch plan and launch 0.08 s) alternating with the device
sketch it waits for (0.30-0.45 s), so moving work behind a launch does not
shorten a wave; only less host work does (the block allocations, 0.05-0.06 s
+ the owner pass's, and the staging of the owner batches are what is left
to trim).

**320^3 on 64 ranks, 40 GB nodes (2026-10-02, job 59232015,
build/gpu_prep/ca0/g320_m9): the committed build (frozen_m8), M9a+b
(frozen_m9a) and all of M9 (frozen_m9b) in one job.**  Two conventions,
kept apart from here on: the loop's "total time" line includes the
warm-up (kernels, buffers, exchanges), the driver's "Factor time" line
excludes it; the tables of M7/M8 quoted "total time" (then with 0.1-0.2 s
of warm-up in it, the same in every build compared).

| | frozen_m8 | M9a+b | M9 (all) |
|---|---|---|---|
| Factor time (excl. warm-up), Color / CA0 | 6.48-6.60 / 6.64-6.71 s | 6.26 / 6.34 s | 6.22-6.26 / 6.21-6.23 s |
| loop total (incl. warm-up), Color / CA0 | 6.61-6.72 / 6.76-6.84 s | 6.60 / 6.69 s | 6.58-6.60 / 6.55-6.58 s |
| warm-up | 0.16-0.22 s | 0.42 s | 0.36-0.42 s (buffers 0.21-0.22) |
| leaf level, Color / CA0 | 2.26 / 2.28 s | 1.99 / 2.01 s | 1.92-1.96 / 1.99-2.00 s |
| leaf start, rank 0 / slowest, Color | 221 / 210 ms | 28 / 46 ms | 38 / 24 ms |
| leaf start, CA0 | 167 / 106 ms | 85 / 30 ms | 84 / 24 ms (rank 0 waits in the level's collectives) |
| start of the level above the CA leaf | 82 ms | 6 ms | 6 ms |
| CA halo gather | 138 ms (plan 17, messages 114) | 55 ms (plan 0, slowest rank's plan + packing 6, messages 50) | 64-103 ms (messages 58-98) |
| leaf transition, Color / CA0 | 0.29 / 0.30 s | 0.29 / 0.30 s | 0.23-0.24 / 0.25-0.26 s |
| leaf waves, rank 0 / slowest, Color | 874 / 1066 ms | 806 / 969 ms | 805 / 1014 ms |
| leaf waves, rank 0 / slowest, CA0 | 1024 / 1475 ms | 953 / 1360 ms | 949 / 1376 ms |

Results identical to the committed build's in every run (CA0 logdet
-446707784.16925561, Color -446707781.18238521, check values unchanged,
replica check 124704 copies 0 differ, solve checks 5.5e-15 / 1.0e-14).

What is real and what is accounting: the page-locking (M9a) took 0.18 s
in the Color leaf's start (two copier threads pinning 128 MB each beside the
level's own start) and 0.1 + 0.08 s in the CA0 leaf's and the next level's
starts, and now takes 0.22 s in the warm-up, serially: a first
factorization's wall time does not change from it (the loop totals show
the level sums falling by 0.27-0.45 s while the warm-up grows by 0.2 s); a
process's later factorizations gain the 0.1-0.2 s outright.  The rest is
real: CA0's gather 0.14 -> 0.06-0.10 s (what remains is the transfer and
the ranks' arrival skew, the plan and packing being 6 ms on the slowest
rank), the leaf transition -0.05 s, the leaf waves -0.07 s (rank 0) to
-0.1 s (slowest rank).  Against the committed build the loop's wall time
is 0.03 s shorter for Color and 0.26 s for CA0; the reported factor times
6.22-6.26 (Color) and 6.21-6.23 s (CA0) are 0.22 s of that below the wall
time.  CA0 and Color are even at this size.

### M9e. The sketch's list kernels, block per box (2026-10-02)
Why CA0 ties Color at the leaf (320^3/64, slowest rank): CA0 eliminates
6936 boxes on an interior rank (4096 local + 2840 ghosts for a 16^3 region)
against Color's 4096, 0.4 s more device time, and saves Color's 0.47 s of
transports; its leaf is device-bound (busy about 90% of the wave time), so
the one lever that favors CA0 is the per-box device time (leveraged 1.7 to
1, while Color's transports do not shrink).  An Nsight Systems profile of
rank 0 (8 ranks, Laplace 160 Nmin 64, Color; build/gpu_prep/m9/prof.sh,
rank 0 under nsys through nsys_own.sh; the report has to be imported from
the .qdstrm on local disk, CFS refuses the .nsys-rep) showed the sketch's
three list kernels (the box's std::mt19937_64 draws, the destination lists
of its kernel rows, those of its fill sources) at 176 ms of the run's 2.34 s
of kernel time, 1-6 ms per launch at the leaf: one warp per box, the
counting sort's per-lane histograms in global memory, a serial chain of
thousands of latency-bound steps per box, visible as the sketch's "upload"
phase (75-190 ms per leaf).

Rewritten in device_sketch.cu with one block of 512 threads per box or
source: the generator's twist runs its three dependency phases with one
state word per thread and tempers the 312 outputs in parallel (the same
arithmetic, bitwise the host's draws); the destination lists come from a
block-wide stable partition (each warp owns a contiguous range of the
input, read 32 at a time in order; counts by ballots, a prefix over warps
and owners, then the scatter), so every list keeps the input order and the
sketch kernel's input is unchanged.  The histogram scratch is gone from
SketchListsItem / SourceListsItem and both callers (level_eliminator.hpp,
compression_gpu.hpp).

Profile, rank 0, whole run: draws 61 -> 12 ms, row lists 66 -> 14 ms, source
lists 49 -> 11 ms (max per launch 2.6-6.4 -> 0.5-1.1 ms).  Leaf (160^3,
rank 0): the sketch's "upload" phase 75 -> 27 ms (Color), 113 -> 37 ms
(CA0); leaf waves 794 -> 767 ms (Color), 951 -> 889 ms (CA0).  Results
bitwise unchanged on 8 ranks (Laplace 192 Color and CA0 with the solve and
replica checks, 160^3 both, VIE 96 unstructured, the EFIE sphere 9000 and
its compression-only run with --precon 2 --solve 0: every logdet, check and
solve-check line identical to frozen_m9b's).

**320^3 on 64 ranks, 40 GB nodes (2026-10-02, job 59251416,
build/gpu_prep/ca0/g320_m9c): frozen_m9b against frozen_m9c in one job**
(Factor time = the driver's line, without the warm-up; the leaf = level 6):

| | M9a-d (frozen_m9b) | + M9e |
|---|---|---|
| Factor time, Color | 6.19 / 6.25 s | 5.95 / 5.97 / 6.02 / 6.03 s |
| Factor time, CA0 | 6.19 / 6.33 s | 6.11 (cold) / 6.00 / 6.00 / 6.00 s |
| leaf level, Color / CA0 | 1.89-1.93 / 1.95-1.99 s | 1.83-1.88 / 1.84-1.93 s |
| leaf waves, slowest rank, Color / CA0 | 1014 / 1426 ms | 950 / 1271 ms |
| sketch "upload" (the list kernels), slowest rank, Color / CA0 | 134 / 193 ms | 35 / 60 ms |
| same, rank 0 | 74 / 114 ms | 26 / 36 ms |

Results bitwise unchanged (logdets, check values, replica check 124704
copies 0 differ, solve checks 5.5e-15 / 1.0e-14).  The leaf's slowest rank
gained 0.06 s (Color) and 0.16 s (CA0), the 1.7-to-1 leverage of CA0's
ghost boxes; the factor times 0.2 s (Color) and 0.2-0.3 s (CA0).  CA0 and
Color remain even at this size (6.00 against 5.97-6.03 s): CA0's leaf
waves still take 0.32 s longer than Color's on the slowest rank, against
Color's 0.46 s of transports, and CA0 pays its gather (0.06-0.1 s) and a
larger transition.  What is left of the list kernels is about 1 ms per
wave.

**640^3 on 512 ranks, 40 GB nodes (2026-10-03, job 59251711,
build/gpu_prep/ca0/g640_m9c): commit 040fa78 (frozen_m8, before M9)
against commit 7bdce9f (frozen_m9c, M9a-e) in one job.**

| | 040fa78 | 7bdce9f |
|---|---|---|
| Factor time, CA0 | 8.97 s | 8.54 / 8.49 s |
| Factor time, Color | 9.06 s | 8.57 / 8.64 s |
| leaf level, CA0 / Color | 2.40 / 2.37 s | 2.09-2.12 / 2.01-2.02 s |
| leaf transition, CA0 / Color | 0.31 / 0.29 s | 0.25-0.26 / 0.23 s |
| CA halo gather | 128 ms (plan 17, messages 104) | 91-112 ms (plan 0, slowest rank's plan + packing 6-7, messages 86-107) |
| warm-up | 0.16-0.83 s | 0.37-0.40 s (buffers 0.21) |

Results bitwise the earlier ones (CA0 logdet -3937186089.4963474, Color
-3937186064.8561087; check values unchanged; replica check 1217216 copies 0
differ; solve checks 8.1e-15 / 9.8e-15).  M9 is worth 0.45 s per
factorization at this size (0.3 s of it at the leaf; 0.2 s of the
reported difference is the pinning moved into the warm-up); CA0 and Color
remain even (8.49-8.54 against 8.57-8.64 s).

## Backlog from the GPU Color work (not CA-specific)

- B1. Level-start wait: the first transport of each leaf level (0.4-0.7 s at
  64 ranks) waits for the slowest rank's level setup; trim the setup it
  waits on. (2026-10-01: part of it was the first device exchanges between
  each pair of ranks, now taken out by the exchange warm-up; see M4.)
- B2. Sphere load balance at 64 ranks: per-rank elimination 0 to 0.46 s at
  512K (surface meshes split unevenly across rank bricks).
- B3. CPU CA levels keep their ghosts' data to the end of the factorization
  (Laplace 192 at 8 ranks: 12.6 GB "CA halo" per rank on top of 24.4 GB
  local). Answered during M1: the CA solve applies the ghosts again with
  their factors (X_RR LU, X_SR, X_RS, T, X_NR from `level.ghost_boxes`), so
  the factors must stay; only their blocks can go, and `clear_ghosts`
  already drops those. (Skipping the ghosts' factor download on the GPU
  broke the solve: acc_forward 0.19.) Shrinking this needs a solve that
  receives the ghost factors it needs, or an owner-based solve (M4).
- Dropped by the user: matvec communication overlap.

## Where the factorization's memory went (2026-10-06)

The goal is more boxes per GPU without much more factorization time.  A
level's **resident** heap is the ceiling (its transient part comes and goes
within a wave), and at 192^3 on 8 ranks, 40 GB nodes, the Color leaf's was
26.18 GB of a 33.40 GB heap - two things and nothing else:

| leaf term | baseline | now | how |
|---|---|---|---|
| fill sources (temp2^T, X_RR_full) | 15.51 GB | **2.11 GB** | compressed, done |
| parent blocks (Schur + edges) | 9.47 GB | 9.47 GB | step 2, not started |
| level Schur/edges | 1.14 GB | 1.14 GB | - |
| near, points, lists | 0.06 GB | 0.06 GB | - |
| **resident** | **26.18 GB** | **12.78 GB** | 2.05x |
| transient | 7.02 GB | 7.02 GB | tiling, not started |

logdet -89253588.194065422 -> -89253588.577188373 (4.3e-9), acc_backward
1.3024e-15 -> 1.5636e-15, Factor time 7.71 -> 8.00 s (+3.8%), one binary.
CA0 the same shape: leaf 23.07 -> 13.74 GB, logdet 3.8e-9 off.

### What the fill-source work established

temp2^T is compressed one **slot** at a time - a slot is one one-hop
neighbor's rows - because both readers take whole slots: the sketch's P side
its own slot, its W side one run per ring block, and a run is a slot.  Slots
are 160 KB at the leaf but 2.2 MB at level 2, so a slot is one chunk at the
leaf and several higher up; the per-slot record carries a chunk range.
Ratios 9.4x at the leaf to 15x at level 2, which is well short of X_NR's
18.6x.

Everything is grouped, and every group bound was found by measurement:

- **kFillStagingBudget 64 MB** - the compressor's staging, and the sub-group
  the wave transposes, packs and frees together.  512 MB could not get its
  staging beside the wave's working memory and 14 groups a level fell back.
- **kFillWaveScratch 512 MB** - the sketch decompresses its sources in groups
  of this, because a wave's distinct slots come to 7.79 GB.  Grouping by
  source keeps every reader with the slots it needs and the batches wide.
- the transition puts a compressed source back **whole** into a block of its
  own chunk: it reads rows from all over a source, and the chunk's lifetime
  already exists for the host-restored ones.

### Step 3: the near-field blocks, 9.47 GB of the leaf's 13.77 GB resident

**Outcome (2026-10-08): done, 9.47 -> 3.25 GB, leaf peak -3.60 GB on Color
and -2.94 GB on CA0.**  Read step 3a, then 3c, then "Step 3c reconciled";
step 3b's conclusion is superseded.

**Measured first, because the design turns on it.** The near blocks are not
the parent blocks, so their ratio was sampled directly - half the budget from
the Schur blocks and half from the edges, since a first sample filled itself
from the Schur blocks alone and never reached the edges, which are 89% of the
bytes:

| level | Schur, loose / tight | edges, loose / tight |
|---|---|---|
| 5 (leaf) | 6.3x / 5.2x | 6.6x / 5.4x (CA0 6.4 / 5.4) |
| 4 | 7.0 / 6.1 | 7.1 / 6.2 |
| 3 | 7.0 / 6.0 | 7.7 / 6.8 |
| 2 | 7.8 / 6.6 | 8.2 / 6.7 |

"tight" is the factors' tolerance fraction divided by four.  The edges
compress as well as the diagonal, so the leaf's 10.61 GB of near blocks
(owner near 9.47 + level Schur/edges 1.14) goes to about 1.6 GB loose or
2.0 GB tight.  **Use the tight one**: 0.35 GB for a four times smaller
quantization step on matrices the fill GEMMs read, against factors only the
solve reads.

**The transient is a separate lever and a smaller one.**  Of the leaf's
largest wave's 6.29 GB, only 12% is scratch (A_NS 0.39 + X_NS 0.39); 82% is
the two nt x r outputs (X_NR 2.57 + temp2 2.57), which tiling can stream but
not avoid.  Level 4 is 14% / 81%, CA0 the same shape.  So tiling nt for the
scratch alone is worth 0.78 GB of 6.29 and is not worth doing on its own.

**And a wave writes nearly all of the near blocks** - 10.18 GB of the 10.61
resident, 62.34 GB over 16 waves, 17.23 GB at their largest over 57203
blocks, 8.9 writes each.  So holding them compressed *between* waves saves
nothing: the bounded set has to come from grouping the **targets** within a
wave, not from iterating sources.

#### Step 3a (done): the owner pass runs in bounded target groups

512 MB groups, launched group-major instead of one global sub-batch sweep.
Each target lies wholly in one group and keeps its contributions' order
within that group's sub-batches, so every block sees the same contributions
in the same order; only which GEMMs share a launch changes, which
batch-independence covers.  Group boundaries follow the targets' canonical
order and their shapes, both rank-independent, so replicas group alike.
Verified **bit-identical** at 96^3 on both paths: -9928907.2811287567 and
-9928906.4595639911, the baselines exactly.

#### Step 3c: the design that replaces 3b's store (2026-10-07)

3b failed on the *store*, not the scheme: a re-compressed block's old bytes
were never reclaimed, and reference counting at group-allocation granularity
left a single live block pinning a whole allocation.  (The t* measurement
below rules out only the **compress-once** policy, where a block stays
uncompressed from its first update to its last.  It says nothing against
compressing and decompressing per subset, where a block is uncompressed only
while its own subset runs.)

**The stores.**

- `B_compact`: the compressed blocks, in a few large bump-allocated extents.
  At a chunk's *first* compression it gets exactly its compressed length -
  **no margin** - and keeps that home for the level.
- `B_newtail`: a workspace in fixed-size slots with a free list, for the
  overflow when a re-compressed chunk is longer than its home.  Only the
  *excess* goes here, not a whole new copy.  A chunk's tail is recorded as
  **slot indices, not pointers**, so doubling the workspace (copy slot k to
  slot k) leaves every record valid.  When a chunk is re-compressed its old
  tail slots return to the free list first, so nothing accumulates.
- `Buf_o` and `Buf_c`: two bounded per-subset workspaces, e.g. 128 MB each.
  `Buf_c` must be sized for the compression *output*, not the compressed
  data: bitcomp needs a worst-case buffer, sum of bitcompMaxBuflen >= the
  subset's raw size.  On the way in it holds only the gathered strings
  (about a sixth of that).  **`Buf_c` is also the compression staging**, so
  staging stops being a cost on top of the store.

**The cycle, per subset G_i.**

1. gather G_i's compressed chunks from `B_compact` (+ their tails from
   `B_newtail`) into `Buf_c`, each chunk contiguous - bitcomp requires that,
   and it is why the scattering of tails is invisible to the codec;
2. decompress `Buf_c` -> `Buf_o`;
3. apply G_i's fill updates in `Buf_o`;
4. re-compress `Buf_o` -> `Buf_c`;
5. copy each chunk back: what fits into its `B_compact` home, the excess into
   `B_newtail` slots.

**The three other readers.**  The near blocks are read outside the update
pass as well, and each needs its own window:

- *the per-wave gather* that assembles A_NS and A_NR for the boxes about to
  be eliminated (these are not stored anywhere - they are gathered out of the
  stored blocks into transient arrays, used, discarded).  Read-once: the
  block is gathered and then retired, because the wave replaces it with a
  smaller one.  No re-compression needed.
- *build_parent*, once at the end of the level, forming the coarser level's
  blocks from this level's.  Also read-once, and it already processes parents
  in chunks with its own budget.
- *apply_remote_generators*, which applies updates from boxes eliminated on
  other ranks (817 at the 192^3 leaf against 36893 new targets).  This one
  does need read-modify-write, like the update pass.

**Determinism.**  Nothing here is rank-dependent: the subsets follow the
targets' canonical order and their shapes, each block's contributions keep
their Morton order within its subset, and replicas being bitwise identical
compress to identical lengths, so they take the identical overflow decisions.

**Memory, at the 192^3 leaf.**

| | Color | CA0 |
|---|---|---|
| B_compact (= sum l_first) | 2.96 GB | 4.13 GB |
| B_newtail | ~0.02 | ~0.03 |
| Buf_c + Buf_o | 0.26 | 0.26 |
| **total, replacing 10.61 / 11.25** | **~3.2** | **~4.4** |
| leaf peak 20.24 / 22.51 becomes | **~12.8** | **~15.7** |

**What it does not recover.**  sum l_last is 0.34 GB against sum l_first
2.96 GB - the blocks shrink 8.7x as their boxes reduce to skeletons - so
`B_compact` is mostly slack late in the level.  Recovering that needs an
occasional compaction of `B_compact`; the peak is mid-level, so it is a
follow-on, not a prerequisite.

**Order of work.**  (1) the store: `B_compact` extents with exact-length
chunk homes plus the slotted `B_newtail`, replacing 3b's group allocations
and reference counts, with blocks still materialized into their own
allocations as now.  This alone should remove the CA0 regression, since the
dead weight is what caused it.  (2) the `Buf_o`/`Buf_c` windows and the three
readers.  Each step gated on Color and CA0.

#### Step 3b: built, measured, and **it does not work** (2026-10-07)

> **Superseded on 2026-10-08.**  The conclusion below - that the policy
> cannot pay - was wrong, and the measurements in it were taken against a
> configuration that was holding the whole adopted set raw for an unrelated
> reason.  Read "Step 3c reconciled" after it.  What survives: the last-write
> distribution table, which is why compress-once is not an option.

The machinery is in and verified; the *policy* cannot be made to pay, and the
reason is a property of the elimination, not of the implementation.

- **3b-i**, the owner pass staging and launching one group at a time, is
  **bit-identical** on both paths (-9928907.2811287567, -9928906.4595639911)
  and time-neutral.  It needed more than a loop change: current_slots_in()
  appends a box's skeleton list into whichever image it is handed and returns
  an offset valid only against that image's uploaded base, so each group
  stages its own image with the stamp bumped to force the re-append; and
  `e.c` was baked into the task entries at the join, so allocation moved into
  the group loop with the pointer fixing.
- **3b-ii**, holding a finished group packed, is correct: logdet moves by
  quantization only (5.1e-7 Color, 9.2e-7 CA0), no guard fires, 5.3-6.3x at
  the tighter tolerance, and the codec runs at the measured ~120 GB/s.
- It first appeared to save nothing because a re-packed block's previous
  bytes were never reclaimed - 6.08 GB of raw packed to 1.15 GB where the
  live set was 0.14 GB, about 8x dead weight, matching the 8.9 writes per
  block.  Reference counting per group allocation fixed that for Color
  (leaf 1.91 -> 1.41 GB at 96^3, near blocks 0.73 -> 0.14) but **CA0 got
  worse** (1.99 -> 3.37), because its 32 waves re-pack more and a single live
  block pins its whole group allocation.

**The measurement that settles it.**  When does a near block see its last
write?

| level | 1st quarter | 2nd | 3rd | last | of |
|---|---|---|---|---|---|
| 5 Color | 0% | 11% | 6% | **83%** | 17.23 GB |
| 5 CA0 | 5% | 10% | 16% | **70%** | 22.15 GB |
| 4 Color | 0% | 20% | 11% | 69% | 5.51 GB |
| 2 (both) | 0 | 0 | 0 | 100% | 0.44 GB |

70-83% of the bytes are written into the last quarter of the waves, so they
are necessarily raw when the peak occurs.  No policy - pack per wave, or
pack once after a block's last writer - can shrink the set at that moment.
At 192^3 the whole piece came to: Color peak 20.72 -> 20.24 GB with Factor
7.81 -> 8.63 s; CA0 peak 22.51 -> **25.33** GB with Factor ~8.3 -> 9.87 s.

**Left in place but switched off** (`kHoldNearBlocksPacked`), on the grounds
that 3a and 3b-i are verified and neutral and the only thing wrong is the
policy.  - Which turned out to be the right call for the wrong reason: the
switch is on now, and the thing wrong was elsewhere.

**What this means for the peak.** ~~The leaf's 9.47 GB of near-field blocks is
not reclaimable by compression.~~  Wrong; see the next section.  3.60 GB of it
is.

#### Step 3c reconciled: the pre-pass was holding the raw set (2026-10-08)

**The store works, and holding the near blocks packed pays on both paths.**
192^3, 8 ranks on 2 nodes, A100 40 GB, `--h2_use_sketch 2 --h2_lazy_schur 2`,
`--H2_compress_factors 1`:

| | Color peak | Factor | CA0 peak | Factor |
|---|---|---|---|---|
| store off (`kHoldNearBlocksPacked` false) | 20.72 (13.77 + 6.95) | 7.49 s | 22.51 (15.56 + 6.94) | 8.41 s |
| store on, pre-pass present | 21.01 (13.90 + 7.12) | 8.69 s | 23.48 (16.37 + 7.11) | 10.00 s |
| store on, pre-pass removed | **17.12** (9.46 + 7.66) | 8.51 s | **19.57** (11.91 + 7.66) | 9.97 s |

**-3.60 GB (-17%) on Color, -2.94 GB (-13%) on CA0**, for +14% and +19% of
Factor time.  logdet unchanged from the store-on rows
(-89253733.904870659, -89253641.194709614), accuracy unchanged from store-off
(Color acc_mvp 3.046e-3 -> 2.946e-3, CA0 3.143e-3 -> 2.897e-3; acc_forward
1.554e-2 -> 1.547e-2 and 1.369e-2 -> 1.312e-2), 96^3 gate clean on both with
the round-trip check and the launch-item validator on.

**What was actually holding the memory.**  Piece 4 put every adopted block
back in raw form in a pre-pass before the owner join, because the join
allocated them.  Step 3b-i moved allocation into the group loop and left the
join recording shapes only, which made the pre-pass redundant - and harmful:
it decompressed the entire adopted set up front, 6.16 GB live at the 192^3
leaf, so the whole near-field set was raw at the peak no matter what the
groups did afterwards.  The groups exist to bound that set; the pre-pass
unbounded it one call earlier.  **The fix was deleting 20 lines.**

So both of 3b's negative measurements were measuring the pre-pass.  The
numbers were real; the attribution was not.

**The peak has moved**, and that is where the next lever is.  Leaf resident
at the peak, Color, store off -> on:

| item | off | on |
|---|---|---|
| owner near blocks | 9.47 | 3.25 |
| adopted blocks | 0.00 | 0.00 |
| kept factors | 1.08 | 2.44 |
| fill temp2^T packed | 1.38 | 2.34 |
| fill X_RR_full | 0.64 | 1.13 |
| level Schur/edges | 1.14 | 0.24 |
| total resident | 13.77 | 9.46 |

The near-field items fell by 7.12 GB; the rest rose by 2.81, because the peak
now occurs at a different wave, one with more of the level's factors already
kept.  On CA0 at the peak: kept factors 3.29, fill temp2^T 2.53, owner near
4.31, X_RR_full 1.23, Schur/edges 0.48.

**`kept factors` + `fill temp2^T` is the duplicated temp2 pair** - 4.78 GB on
Color, 5.82 GB on CA0, against the 1.1-1.4 GB estimated when they sat lower
in the ordering.  That backlog item is now the single largest target at the
leaf, and the cheapest: one is the transpose of the other.  The transient,
7.66 GB, is now the larger half of the peak, and 82% of it is the two
n_tot x r outputs.

**Costs and what is left.**  Decompression traffic at the Color leaf is
32.11 GB in 94 group calls, 2.04 GB in 7 generator calls, 0.61 GB in one
build_parent; CA0 48.99 GB in 152 and 1.03 GB in one (no remote generators on
the CA path).  Of the +1 s, roughly 0.64 s is the 512 MB group budget alone -
730 batched GEMM launches at the leaf instead of 120, a cost that is paid
whether or not anything is compressed.  Tuning that budget is the next cheap
win and is deferred by choice; ~2 GB groups would give ~26 groups but add
about 1.5 GB back to the peak.

Reader 1's bounded gather window (512 MB, `materialize_into_window`) is built
and correct and reclaimed nothing here, because what it was built against was
the pre-pass.  It will matter when the gather is the binding term.

#### The original 3b sketch, for the record

Pack a group once its GEMMs are done and free the raw blocks; materialize a
group before its GEMMs; let every other reader materialize on demand.  The
piece 4 machinery carries over - packed_schur_/packed_edges_,
materialize_adopted(), and assert_not_packed() at every site that acts on a
null ptr, which is what found four readers there.  The codec bill is the
62.34 GB of writes per level, so about 125 GB of traffic each way at the
leaf: a few percent of Factor time at the measured rates, but compression is
slower than decompression and that needs measuring, not assuming.

### temp2 stored resident twice: done (2026-10-08), and it was 1.1 GB not 4.8

**Done.**  temp2 is stored once.  Through the level it is the fill source -
temp2^T, r x ntot, compressed a slot (one one-hop neighbor) at a time, which
is the form and the granularity the sketches read.  The wave that keeps a box
*pledges* that copy as the box's solve factor instead of compressing temp2 a
second time, and at the transition, with the heap at its emptiest,
`transpose_kept_xnr()` decompresses each pledged box's slots end to end (they
are consecutive column ranges, so that is the whole r x ntot), transposes into
the ntot x r the solve's kernels read, and compresses again into a block of
the kept store.

| 192^3, 8 ranks | peak before | after | Factor |
|---|---|---|---|
| Color | 17.12 | **16.03** (resident 9.46 -> 8.37) | 8.51 -> 8.64 s |
| CA0 | 19.57 | **18.92** (leaf 19.57 -> 18.07) | 9.97 -> 10.1 s |

The handover costs 0.24 s a run on Color and 0.36 s on CA0, most of it at the
leaf (17.45 GB decompressed, transposed and recompressed in 66 groups,
0.155 s), and it replaces the keep pass's compression, so Factor moves by
0.13 s.  logdet is unchanged to every digit on both paths and at both sizes -
it comes from the LU of X_RR, not from X_NR - and the kept X_NR comes to
exactly the bytes the old path stored (17.45 -> 1.13 GB, 15.4x).  The one
cost is a second quantization of the solve factor, at the same tolerance but
on blocks cut differently: acc_mvp 2.95e-3 -> 3.11e-3 on Color and
2.90e-3 -> 2.98e-3 on CA0, both still inside the band the uncompressed
baseline sits in (3.05e-3, 3.14e-3).

**I had estimated 4.8-5.8 GB and it is 1.1-1.5.**  "kept factors" at the peak
is not all X_NR: X_NR is 92% of what a *level* keeps, but at the peak wave
only the boxes eliminated so far have been kept, and the rest of their
factors (T, X_SR, the LU of X_RR, the pivots) is 1.35 GB of the Color leaf's
2.44.  The duplicated part was the 1.09 GB of X_NR kept up to that wave.  The
lesson from step 3c applies to the categories themselves and not just to
their producers: a resident figure at the peak is not the level's total.

**Reading the reports.**  The keep budget and the per-level "solve factors:
kept" line count a pledged box's X_NR where it actually is, so the figure
printed at the commit is the fill source's slot-wise size - 1.89 GB at the
Color leaf, against the 1.13 GB the whole-matrix form comes to, because
slot-wise deltas are finer.  The handover's own line gives the final size.

**Where the peak is now.**  On Color it is the leaf's waves, 16.03 GB, of
which the transient is 7.66 - 48% of the peak.

CA0 at `--CA_level 4` put the peak at level 4 instead (18.92 GB), but that was
an artefact of the setting: CA_level is the *first* CA level, so 4 at 192^3
made levels 4 and 5 both CA, and level 4 then held the whole level's block set
replicated - parent blocks 6.74 GB plus a 5.99 GB received halo, 12.73 GB a
rank, the halo being installed blocks and not a message buffer
(ca_halo.hpp:279 emplaces them into blocks.schur/edges for the level).  **CA
runs use CA at the leaf only from now on** (the user's decision, 2026-10-08):
`--CA_level 5` at 192^3, one CA level, no halo exchange.  That gives

| 192^3, 8 ranks, CA at the leaf | peak | Factor |
|---|---|---|
| Color | 16.03 (8.37 + 7.66) | 8.60 s |
| CA0 | 18.07 (10.41 + 7.66) | 9.35 s |

and level 4 falls to 11.84 GB, so the leaf binds on both paths.  CA0's extra
2.04 GB of resident over Color is its larger replicated near field (owner near
4.31 vs 3.25, Schur/edges 0.48 vs 0.24) and its kept factors (1.78 vs 1.35).

### The peak after the two 2026-10-08 commits, and what is left

Color leaf, 16.03 GB = resident 8.37 + transient 7.66:

| | GB | |
|---|---|---|
| **transient: the wave's two outputs** | **5.14** | X_NR 2.57 + temp2 2.57, 82% of the largest wave's 6.29 |
| owner near blocks | 3.25 | compressed, 5.9x |
| fill temp2^T packed | 2.34 | also the kept X_NR now |
| transient: codec and gather staging | 1.37 | 7.66 - 6.29 |
| kept factors | 1.35 | X_NR excluded; T, X_SR, LU, pivots, all raw |
| fill X_RR_full | 1.13 | raw |
| transient: scratch A_NS, X_NS | 0.78 | 12% of the wave |
| level Schur/edges, points, misc | 0.66 | |

Next steps, in the order they are worth doing:

1. **The wave's two outputs, 5.14 GB** - the one big lever left.  Compress
   X_NR and temp2 per one-hop neighbor as they are produced and let the owner
   groups materialize the slices they need, which is the bound the groups
   already impose on the near blocks.  temp2's slot store already exists (the
   fill source), and a slot is exactly the row range a target needs, so no new
   granularity problem.  The traffic is the measured "fill source union" per
   wave, 2.63 GB largest / 36.29 GB summed at level 4, so roughly a doubling
   of the update loop's codec work: expect -3 to -4 GB for +0.5 to 1 s.
2. **Compress X_RR_full, 1.13 GB** - the other fill source, and the only one
   still raw.  One contiguous r x r block a box, so it is one slot in the pass
   that already compresses temp2^T; the readers are the same sketches, which
   already read a lossy temp2^T.  About -0.9 GB for little work.
3. **Compress the kept factors other than X_NR, 1.35 GB at the peak** - T and
   X_SR, leaving the LU of X_RR and the pivots exact.  They are now the
   majority of a stored kept block (1.47 GB of 2.60 at the leaf), so this also
   takes about -1.2 GB a level off the solve's resident factors.
4. **One delta per box for the fill source's slots** - the slot-wise store is
   1.89 GB at the leaf where the same data as whole matrices is 1.13, because
   each slot gets a delta from its own max (factor_compress.hpp keeps the scale
   per block by design).  Sharing the box's delta across its slots is about
   -0.9 GB resident *and* makes the handover's second quantization exact, which
   would undo the acc_mvp 2.95e-3 -> 3.11e-3 of the temp2 merge.  But it
   loosens what the fill GEMMs read, which feeds the factorization itself, so
   it needs logdet and acc measured before it is adopted.

1 to 3 would put the Color leaf near 10.5 GB and CA0 near 12.5, which is
2.0x the boxes a rank held this morning.  Reaching the 8x goal needs FP32 on
the stores (deferred) and the near field itself on top of that.

#### Step 1 measured and designed (2026-10-08)

**The number the design turns on.**  Each candidate C needs, from every source
E next to it, *all* of E's X_NR slots, because A runs over E's whole neighbor
list - so whether a window helps at all depends on how many distinct sources
one group touches.  Measured per owner group under `BPACK_CHECK=fillsrc`
(targets are grouped in canonical order, which is Morton order, so a group is
a spatial patch and its sources are that patch's neighborhood):

| 192^3 leaf | largest group's union | with repeats | summed over groups |
|---|---|---|---|
| Color | **0.82 GB** | 5.26 | 37.59 GB in 102 groups |
| CA0 | **0.82 GB** | 5.27 | 53.22 GB in 156 groups |

So a window of about 1 GB replaces the 5.14 GB the wave holds raw, and the
union is 6.4x smaller than materializing per task.  The traffic is 37.6 GB a
level on Color and 53.2 on CA0, against the 32.1 and 49.0 GB the groups
already decompress for the near blocks: roughly a doubling of the update
loop's codec work, 0.1-0.2 s at the measured 262-365 GB/s, plus compressing
X_NR itself.

**The shape.**  Both operands are stored transposed, which makes a slot
contiguous: temp2^T already is (it is the fill source the sketches read), and
X_NR^T becomes a store of the same kind.  The owner GEMM then runs
(Trans, NoTrans) instead of (NoTrans, Trans) - implemented in both the magma
and the tensor-core path (device_gemm.cu:242) - with lda = ldb = r, so the
window needs no transposes and no second staging.

**Lifetimes differ between the two stores, which is what makes this cheap.**
temp2^T must live for the whole level (later waves' sketches read it) and
already does.  X_NR^T is read only by the owner pass of the wave that produced
it, so it is built, used and freed inside the wave: about 0.45 GB resident
instead of 2.57 GB raw.  Freeing both raw outputs is then worth 2.57 + 2.1 GB
against a 0.82 GB window, about -3.9 GB.

**Stages, each gated.**
- A. The X_NR^T store, built and round-trip checked, nothing reading it.  The
  gate is bit-identical logdet, since no reader has changed; the peak goes
  *up* by the store's 0.45 GB until stage C.
- B. The owner pass reads both stores through the window.  logdet moves here:
  the near-field updates would be computed from quantized temp2 and X_NR,
  where today they use the raw values.  That is the same class of change as
  the near blocks and the fill sources, both already lossy, but it is a change
  to the factorization and not just to storage, so acc_forward/acc_backward
  and logdet's relative change are the gate.
- C. Free the raw temp2 and X_NR once the stores are built.  Only in the
  configuration where nothing else reads them - every box pledged and
  `skip_xnr_` - so the host download and a non-pledged keep keep today's path.

#### Step 1 built: all three stages, and what they were worth (2026-10-08)

| 192^3, 8 ranks, CA at the leaf | Color peak | Factor | CA0 peak | Factor |
|---|---|---|---|---|
| before step 1 | 16.03 | 8.60 s | 18.07 | 9.35 s |
| A: the X_NR^T store | 16.03 | 8.93 | 18.07 | 9.51 |
| B: the owner pass reads the stores | 16.02 | 9.09 | 18.08 | 9.73 |
| C: the raw outputs freed | **15.40** | 9.48 | **17.47** | 10.6 s |

Every stage passed its own gate.  A was bit-identical, nothing having read the
new store.  B moved logdet by 1.4e-8 and 3.1e-8 on Color and 8.5e-10 and
2.2e-8 on CA0, with acc_mvp, acc_forward and acc_backward inside the
uncompressed baseline's band - the near-field updates are computed from
quantized operands now.  C was bit-identical to B, since only the lifetimes
changed.  The window came out exactly as measured: 0.82 GB in the largest
group, 37.6 GB over 102 groups.

**And the whole piece is worth 0.62 GB, where the wave's buffers said 3.9.**
The outputs are 82% of what the largest leaf wave allocates, and freeing all
5.14 GB of them moved the transient peak by 0.61.  So the peak is not reached
while they are live: it is reached later in the wave, inside the owner pass,
where they were already dead.  "The largest wave's transient 6.29 GB: outputs
82%" is a sum of what a wave plans to allocate, and I read it as the
composition of the peak.

That is the fourth time in one day: the Transition label, the read union, the
Adopted kind, and now this.  The pattern is always the same - a total, or a
churn, read as an occupancy - and the fix is always the same instrument.  So
the heap now tags **transient** blocks by kind as well, and prints their
composition at the peak beside the resident one (wave work, wave outputs, wave
sources, wave result, codec staging, source window, sketch).  Nothing further
on the transient until that measurement is in.

**What stages A to C leave behind, which is not nothing:** the near-field
updates no longer need the raw outputs at all, so the wave's working set is
decoupled from the owner pass.  Whatever turns out to set the transient peak,
it can now be moved or bounded without the outputs being pinned behind it.

#### The production runs in chunks, and the peak leaves the factorization (2026-10-08)

`kWaveChunks` chunks of a wave's boxes for the production half; the consumers
still see the whole wave.  8 chunks: Color 12.84 -> **11.33 GB**, CA0
14.57 -> **12.46**, Factor +12.5% on both, logdet bit-identical.  **16 chunks
gave exactly the same peak and cost 0.5-0.8 s more**, because by then the
production is no longer what binds; 8 is the setting.

Three things it needed, each of which broke a run first: the generator pack
registers a wave's transport once (split into plan_generators once a wave and
copy_generators a chunk at a time); the LU's INFO array is the wave's, one
entry a box, and every chunk was writing the same first slots; and the
predicate for freeing a chunk's raw output is not the one for letting the
owner pass read the stores - conflating them meant nothing was freed until the
last chunk, which put the peak *up* to 14.01 GB.

**Where the peak is now, with every allocation tagged.**  The wave's own items
are gone from the Color leaf's peak, and what is left is of a different
character on each path:

| Color leaf, 11.33 GB | GB | | CA0 leaf, 12.46 GB | GB |
|---|---|---|---|---|
| **generator products** (transient) | **4.62** | | owner near blocks | 4.33 |
| level Schur/edges | 3.07 | | fill temp2^T packed | 2.23 |
| fill temp2^T raw | 1.69 | | kept factors | 1.42 |
| owner near blocks | 1.07 | | level Schur/edges | 1.13 |
| codec staging | 0.29 | | wave work (transient) | 1.09 |
| adopted blocks | 0.23 | | codec staging + source window | 1.08 |
| the rest | 0.36 | | the rest | 1.18 |

Color's largest item is the `work` PiecedBuffer of
`apply_remote_generators` - the buffers for the products the *other ranks'*
generators drive, which is why CA0 does not have it at all (a CA level sends
no generators).  CA0's peak is its replicated resident set.  So the two paths
now need different work, and neither needs anything more from the wave's
production.

#### The remote generators add in place (2026-10-08)

A box eliminated on another rank sends us its factors, and applying it means
the same near-field update the owner pass does locally: for every pair (C, A)
of its neighbours that live here, N_{C,A} -= temp2_E[C rows] X_NR_E[A rows]^T.
Each is one small GEMM.  They could not be added straight into their blocks,
because several remote boxes contribute to the same block and a batched GEMM
may not have two entries writing one output - so every partial was given a
slot in a staging buffer and a kernel summed the slots afterwards.  Measured,
that buffer was the Color leaf's largest single item: 4.62 GB.

| 192^3, cumulative over a level's 7 generator calls | leaf |
|---|---|
| partial products staged | 49729, in 11.09 GB |
| targets they sum into | 34339, of 7.63 GB |
| partials a target | 1.4 average, **7 at most** |

Only 1.4 partials a target, so most blocks had a whole copy made of a result
that could have been added in place, and no block needed more than 7.  That is
the local owner pass's situation exactly, and its answer applies: sub-batch j
holds the j-th partial of every target, launched in sequence, so no launch
writes a block twice and each partial goes in with beta = 1.  At most 7
sub-batches, about 69% of the work in the first.

Color peak 11.33 -> **10.69 GB**, the transient 5.02 -> 2.52, Factor unchanged
at 10.6 s, and logdet bit-identical on both paths at both sizes.  CA0 is the
control and is identical throughout, a CA level sending no generators.

**Both paths are now resident-dominated and their transients agree** (2.52 GB
each, which is wave work 1.09 + codec staging 0.56 + source window 0.52 +
small).  The Color leaf's 8.17 GB resident is owner near blocks 3.24, fill
temp2^T 2.04 (which is also the kept X_NR), kept factors 0.99, level
Schur/edges 0.88, adopted blocks 0.48, X_RR_full 0.29, owner source slots
0.21.  CA0's 9.93 is the same list with the replicated sizes.

So the transient is bounded by budgets we choose and the resident is what grows
with boxes a rank.  FP32 on the two compressed stores - 3.24 + 2.04 on Color,
more on CA0 - is now the one lever left that is worth more than 1 GB.

#### One quantization step a box for the fill slots (2026-10-08)

A box's temp2^T is compressed a slot at a time - one slot a one-hop neighbor,
which is the granularity the sketches read - and each slot took its step from
its own maximum.  A finer step costs bits, so the slot-wise store was far
bigger than it needed to be.  The slots of a box now share one step, taken
from the largest of their maxima, which holds the box to one relative
tolerance: the same basis the kept X_NR has always had.  Entirely host-side,
since the deltas were already computed on the host from the per-block maxima.

| 192^3 | peak | resident | fill temp2^T | store ratio | Factor |
|---|---|---|---|---|---|
| Color before | 10.69 | 8.17 | 2.04 | 9.4x | 10.6 s |
| Color after | **9.79** | 7.27 | **1.22** | **15.3x** | 11.1 s |
| CA0 before | 12.46 | 9.93 | 2.23 | | 11.8 s |
| CA0 after | **11.48** | 8.96 | **1.34** | | 12.1 s |

**This is the first step that costs accuracy rather than nothing**, and the
user accepted the trade on the measurement:

| acc_forward | uncompressed | before | after |
|---|---|---|---|
| Color | 1.554e-2 | 1.564e-2 | 1.681e-2 |
| CA0 | 1.369e-2 | 1.593e-2 | 1.584e-2 |

So -0.9 GB on both paths for +7.5% of Color's forward error and none of CA0's,
whose degradation was already there.  acc_mvp does not move (Color 3.09e-3 ->
3.12e-3, CA0 3.06e-3 -> 2.92e-3) and acc_backward stays at 1e-15; logdet moves
1.1e-7 relative on Color and 9.1e-8 on CA0, about ten times the perturbation of
the steps before it.  Anything further here should be gated the same way, on
acc_forward rather than on bytes.

#### The near pass joins the chunks, and the wave's work blocks go (2026-10-08)

The one production item left whole was the work block a box - tmp1, tmp2,
X_SR, A_NS, X_NS - because step 5, the near pass, reads A_NS and X_NS and ran
once over the whole wave.  Moving step 5 and its launch inside the chunk loop
lets a chunk's work blocks go with the rest of its raw output.

Safe because **boxes in one wave are never one-hop neighbours** - that is what
makes them a wave - so no box can see a same-wave neighbour's state change
between chunks.  And unlike dividing the wave, it adds no re-packing: step 5
*creates* each near block and writes it once.

| 192^3 | peak | resident | transient | wave work | Factor |
|---|---|---|---|---|---|
| Color before | 9.79 | 7.27 | 2.52 | 1.09 | 11.1 s |
| Color after | **8.74** | 7.17 | **1.58** | **0.00** | 10.8 s |
| CA0 after | **10.44** | 8.86 | 1.58 | 0.00 | 12.0 s |

logdet bit-identical on both paths at both sizes, and slightly faster.  The
one thing this needed was care with `stores`: left outside the loop it would
accumulate, and every chunk's launch would re-apply the earlier chunks'
updates.

**The transient is now 77% the two budgets** - source window 0.64 and codec
staging 0.58 of 1.58 GB - so what is left of it is a memory-for-launches knob
rather than data.  The Color leaf's resident 7.17 is owner near blocks 3.24,
fill temp2^T 1.22, kept factors 0.99, level Schur/edges 0.88, adopted 0.48,
X_RR_full 0.29, owner source slots 0.12.

#### Step 5's blocks go into the near store as they are written (2026-10-09)

`level Schur/edges` was the near-field blocks the store did not cover: step 5
creates one per (eliminated box, neighbour) as A_NS + update, and it sat raw
until the first owner group to update it came by - 0.88 GB of the Color leaf's
resident set at the peak, 1.13 on CA0.  Step 5 now packs its own blocks right
after add_store has written them, in groups bounded like the owner pass's;
all the store's packer reads of a target is block, key, rows and cols, and
every reader already puts a packed block back before touching it.  Nothing in
the same wave writes them again: the owner pass writes pairs of a source's
neighbours, and these are the source's own.

| 192^3 | max peak | leaf peak | Factor | logdet (rel.) | acc_forward |
|---|---|---|---|---|---|
| Color before | 8.23 | 8.23 | 10.9 s | | 1.681e-2 |
| Color after | **7.97** | **7.75** | 11.0 s | 3.3e-8 | 1.598e-2 |
| CA0 before | 9.92 | 9.92 | 12.3 s | | 1.584e-2 |
| CA0 after | **9.18** | 9.18 | 11.5 s | 8.2e-8 | 1.619e-2 |

At the leaf it packs 5.27 -> 1.04 GB over 128430 blocks in 0.11 s (Color),
`level Schur/edges` is gone from both resident lists, and CA0 is 0.8 s faster.
Accuracy is a wash: Color's acc_forward comes back from 1.681e-2 to 1.598e-2,
CA0's moves 1.584e-2 -> 1.619e-2, both inside the band of the step before.

**On Color the peak has left the leaf.**  The leaf fell to 7.75, and level 4,
which had been sitting at 7.72 all along, now binds at 7.97 - 96% resident:

| Color level 4, 7.97 GB | GB |
|---|---|
| kept factors (the leaf's, all of them now) | 2.64 |
| adopted blocks (the leaf's parents, put back lazily) | 2.29 |
| fill temp2^T raw (level 4's own, before packing) | 1.18 |
| parent blocks | 0.94 |
| owner near blocks | 0.54 |
| transient, all of it | 0.28 |

That is the level-end resident floor predicted when the production was
chunked: once a level's transient is small, the next level up inherits the
whole kept set and the whole parent set at once.  CA0's leaf still binds, at
9.18 = 8.15 resident + 1.03.

#### The remote generators' temp2^T is staged a group at a time (2026-10-09)

`fill temp2^T raw` at Color's level 4 was apply_remote_generators' doing: one
raw block per remote generator for its transposed temp2^T, all allocated up
front, all transposed in one launch, and all waiting on a single compression
pass at the end - 1.18 GB live at the level's peak for the sake of a store
that comes to 0.12.  The same stage-everything-then-collapse shape as the
sketch's [P; W] and the generators' products.

Now one staging sized for the largest group (kFillStagingBudget, 64 MB) is
shared: a group's transposes launch as a contiguous sub-range of the list
already in the image, its slots are packed (the packer waits for the staging
to be read), and the next group reuses it.  A generator whose pack finds no
room gets a block of its own copied out before the staging moves on, so the
fallback is what it was; it never fired here.  Without compression nothing
changes, since then the raw block *is* the fill source.

| 192^3 Color | max peak | level 4 | leaf | Factor |
|---|---|---|---|---|
| before | 7.97 (level 4) | 7.97 = 7.68 + 0.28 | 7.75 | 11.0 s |
| after | **7.75** (leaf) | **7.33** = 6.24 + 1.09 | 7.75 | 10.9 s |

logdet is bit-identical on both paths at both sizes - bitcomp is deterministic
per block, and only the staging address changed - and CA0, which sends no
generators, is identical throughout.  Level 4's resident falls 1.44 GB and
the leaf binds again on Color, 0.42 GB above level 4.

**Where the two paths stand**, every item tagged, no transient above 0.4 GB:

| Color leaf, 7.75 GB | GB | | CA0 leaf, 9.18 GB | GB |
|---|---|---|---|---|
| owner near blocks | 3.24 | | owner near blocks | 4.31 |
| fill temp2^T packed | 1.47 | | kept factors | 1.90 |
| kept factors | 1.47 | | fill temp2^T packed | 1.58 |
| transient, all of it | 1.03 | | transient, all of it | 1.03 |
| X_RR_full, adopted, slots, points | 0.53 | | X_RR_full, adopted, slots, points | 0.35 |

Color's level 4 underneath is the leaf's kept factors (2.65) and its parents
put back lazily, 6.24 resident, which no item here touches.  Across the two
days: Color 20.72 -> 7.75 GB (-63%), CA0 22.51 -> 9.18, Factor 7.49 -> 10.9 s.

#### T and X_SR compressed in the kept block (2026-10-09)

A kept block is T (k x r), X_SR (k x r), X_NR, the LU of X_RR and its pivots.
With X_NR compressed, the leaf's stored block on rank 0 was 2.63 GB of which
X_NR 1.14, T and X_SR 0.48 and the LU with pivots 0.99: the two small factors
were 18% of what the solve keeps, raw, and the LU - which must stay exact -
38%.  Now T and X_SR go through the same compressor batch as X_NR
(keep_solve_factor_group adds two sources a box, pledged boxes included),
each with its own chunk list in the level's chunk_strides pool and its own
`t_bytes`/`xsr_bytes` on KeptSolveBox.  The solve's per-wave batch decompresses
T, X_SR and X_NR of a wave's boxes into the parity buffer one after another and
patches the three SolveBox pointers (launch_patch_xnr; a null destination
leaves a pointer alone, which is how a box with X_NR read in place keeps it).
The two host readers of X_SR (materialize_host_factors_level, rescue_skipped_xnr)
decompress with the host codec as they do X_NR; T has its own host copy.  The
LU and pivots are raw.  land_full_factors (owned ghosts of a multi-rank Color
level) keeps them raw as it keeps X_NR raw; it kept no box in these runs.

One bug on the way: build_solve_xnr_wave is shared with the factorization's
own batches (the owner source window, transpose_kept_xnr), which patch no
pointers and leave the new box arrays empty - `data()` of an empty vector is
null, and a 4096-pointer copy from it is "invalid argument".  Those arrays are
zeroed when a batch has none.

**The step is the whole story.**  At X_NR's step (kXnrToleranceFraction 0.2)
the memory came as predicted and the accuracy did not: acc_forward +55% at
96^3 and +34% at 192^3 on both paths, acc_mvp +40% (the device multiply reads
the same kept T), where X_NR at that step had cost 3%.  T and X_SR are 2.5%
of the raw bytes, so a finer step is nearly free - each 10x is 3.3 bits an
element - and rather than pick one, four were measured at 96^3 (acc_forward
against the X_NR-only run, 1.0098824e-2 Color / 9.7304948e-3 CA0):

| T, X_SR step / tol | stored | acc_forward Color | CA0 |
|---|---|---|---|
| 0.2 (X_NR's) | 5.7x | +55% | +56% |
| 0.02 | 4.4x | +0.62% | +0.71% |
| **0.002** | **3.6x** | **+0.012%** | **+0.015%** |
| 0.0002 | 3.0x | -0.0002% | -0.0008% |

kSmallFactorToleranceFraction = 0.002, passed per source into the batch
(`tolerance_of`, one entry a block, in place of a second compressor call).
Against 0.02 it costs 24 MB a rank at 192^3, 0.3% of the peak.  A trap in the
sweep worth writing down: claplace3d_h2 has an RPATH and all of the H2 code is
in libbutterflypack.so.5, so three "binaries" built with three constants were
one executable loading the last-built library - identical results to every
digit, identical md5s.  The sweep was rerun with one library snapshot per
constant preloaded inside the srun wrapper, with 0.0002 repeated as a control
that reproduced the first run exactly.

| 192^3 | peak | resident | kept factors at the peak | kept store, rank 0 | Factor | time_solve | acc_forward |
|---|---|---|---|---|---|---|---|
| Color before | 7.75 | 6.72 | 1.47 | 2.63 | 10.9 s | 0.467 | 1.5978148e-2 |
| Color after | **7.40** | 6.37 | 1.12 | 2.28 | 11.2 s | 0.476 | 1.5978626e-2 |
| CA0 before | 9.18 | 8.15 | 1.90 | 3.49 | 11.6 s | 0.501 | 1.6185732e-2 |
| CA0 after | **8.74** | 7.72 | 1.47 | 3.05 | 11.6 s | 0.508 | 1.6186747e-2 |

logdet is bit-identical on both paths at both sizes (the factorization never
reads the kept T or X_SR; 96^3 heap-checked), acc_mvp is unchanged to five
digits, the solve is 2% slower (a wave's buffer grows 1%).  Level 4 follows:
Color 7.33 -> 6.98 (kept factors 2.65 -> 2.30), CA0 8.11 -> 7.67.  The
per-level kept T and X_SR line in the report ("kept T and X_SR: 0.48 GB raw ->
0.13 GB stored (3.6x); X_NR 17.45 -> 1.14 GB; the LU and pivots 0.99 GB,
uncompressed") is the measurement that the 2026-10-08 list's "X_NR is N% of
the kept block" line was meant to be; that one never printed (its counters
are read after the level's stats are reset) and now carries the same fields.

acc_backward - the residual of the device multiply applied to the device
solve, 1e-15 everywhere until now - reads 4.3e-13 (96^3) and 5.3e-12 (192^3)
on CA0 at 0.002, deterministic (two runs agree to all digits), while 0.02 and
0.0002 give 1e-15 and Color gives 1e-15 at 0.002.  One quantization step off
in one element would show as ~1e-8, so this is rounding-level over many
elements, and it is specific to the replicated CA leaf, where a ghost box's
factors are kept and compressed on every rank that holds it.
With one rank (the 3D tree needs 8^k, so 1 against 8 is the only
replication-free comparison) CA0 at 0.002 is the Color run to every digit
and acc_backward is 1.2e-15; with 8 ranks Color stays at 1e-15 and CA0 does
not.  A ghost box's factors are eliminated on every rank that holds it, from
inputs accumulated in different orders, so two raw copies differ in the last
bit; quantized, they agree exactly unless an element sits within that bit of
a rounding boundary of the step, and then the two ranks' copies are one step
apart there - a chance event per step (none at 0.02 or 0.0002, one or a few
at 0.002), deterministic, and delta * |x_j| / ||b|| in size, which is the
order seen.  Lossy-compressed replicated factors are consistent across ranks
only to the step: X_NR has carried the same property since it was compressed
(a coarser step, 36x the elements), and at 5e-12 against a 1.6e-2
factorization it is noted, not pursued.

**Where the two paths stand**, every item tagged:

| Color leaf, 7.40 GB | GB | | CA0 leaf, 8.74 GB | GB |
|---|---|---|---|---|
| owner near blocks | 3.24 | | owner near blocks | 4.31 |
| fill temp2^T packed | 1.47 | | fill temp2^T packed | 1.58 |
| kept factors | 1.12 | | kept factors | 1.47 |
| transient, all of it | 1.03 | | transient, all of it | 1.03 |
| X_RR_full, adopted, slots, points | 0.53 | | X_RR_full, adopted, slots, points | 0.35 |

Of the kept factors at the leaf's peak, the LU and pivots are now the
majority and exact by design; what is left of the kept block is FP32 on the
three stores (deferred).  Across the two days: Color 20.72 -> 7.40 GB (-64%),
CA0 22.51 -> 8.74 (-61%), Factor 7.49 -> 11.2 s.

#### With and without bitcomp, same binary (2026-10-09)

`--H2_compress_factors 0` against 1, commit 66a1002, 8 ranks on 2 nodes.
96^3, heap-checked (the peak is level 4 on both):

| 96^3 | Color off | Color on | CA0 off | CA0 on |
|---|---|---|---|---|
| peak (resident + transient) | 6.21 (5.54 + 0.67) | **1.84** (1.08 + 0.76) | 10.02 (9.34 + 0.67) | **2.48** (1.62 + 0.86) |
| the resident's two largest | temp2^T raw 3.20, kept 2.13 | Schur/edges 0.61, owner near 0.27 | kept 4.65, temp2^T raw 4.32 | owner near 1.09, temp2^T 0.26 |
| kept store, rank 0, leaf | 2.18 (1.0x) | 0.28 (7.9x) | 4.69 (1.0x) | 0.55 (8.5x) |
| Factor | 1.39 s | 2.23 s (+60%) | 1.34 s | 2.63 s (+96%) |
| GPU solve | 0.097 s | 0.112 s (+16%) | 0.116 s | 0.136 s (+18%) |
| time_mvp | 0.250 s | 0.274 s | 0.285 s | 0.300 s |
| acc_mvp | 7.47e-4 | 1.55e-3 (2.1x) | 7.53e-4 | 1.48e-3 (2.0x) |
| acc_forward | 3.46e-3 | 1.01e-2 (2.9x) | 3.42e-3 | 9.73e-3 (2.8x) |
| logdet | -9928892.040 | -9928903.097 (1.1e-6) | -9928892.266 | -9928897.922 (5.7e-7) |

Without compression the transient is the sketch's [P; W] (0.63 GB) and the
resident is the raw fill source plus the raw kept factors; with it, the
transient is the budgets (source window 0.54) and the peak is 3.4x / 4.0x
lower.  The accuracy cost is the standing one, and the same at both sizes:
the compressed 96^3 run has read 1.0e-2 since the first log, and at 192^3
the uncompressed Color run reads acc_forward 4.20e-3 / acc_mvp 8.11e-4
against 1.60e-2 / 2.95e-3 compressed - 3.8x and 3.6x.  (A "1.554e-2
uncompressed" quoted earlier in this work was a compressed variant; the
2026-10-06 `*192off` logs are the uncompressed ones.)

192^3 on 8 ranks does not run without compression any more: the leaf dies in
its 12th wave of 16 on a 3.85 GB request with 21.7 GB in use and 2.4 GB in
the largest free range - the raw fill source (17.45 GB a level) and the raw
kept X_NR (another 17.45, which only the pledge lets them share) cannot both
sit on the card, and there is no option to decline the keep (decide_solve_keep
is automatic, per the env policy).  16 ranks is not a valid count (the 3D
tree needs 8^k; 64 would be 16 nodes and 512 boxes a rank, a different
regime).  The uncompressed figures that exist for this case are Color only, from
the 2026-10-06 code (A_v192off / A_f192off), when it still fitted with
nothing to spare; today's compressed runs beside them:

| 192^3, 8 ranks | Color off (2026-10-06) | Color on | CA0 on |
|---|---|---|---|
| peak (resident + transient) | 33.21 (26.18 + 7.02) | **7.41** (6.37 + 1.03) | **8.74** (7.72 + 1.03) |
| the resident's two largest | temp2^T raw 14.82, parent blocks 9.47 | owner near 3.24, temp2^T packed 1.47 | owner near 4.31, temp2^T 1.58 |
| kept store, rank 0, leaf | none: keep failed, 18.92 GB uploaded at the first solve (set-up 2.25 s) | 2.28 (8.3x), set-up 0.11 s | 3.05 (8.6x) |
| Factor | 7.71 s (7.54 in its twin) | 11.2 s (+45%) | 11.6 s |
| GPU solve (fwd / bwd, MPI) | 0.257 s (0.105 / 0.133, 0.068) | 0.464 s (0.183 / 0.229, 0.140), +81% | 0.499 s |
| time_mvp | 2.61 s (the multiply of that date) | 0.890 s | 0.885 s |
| acc_mvp | 8.11e-4 | 2.95e-3 (3.6x) | 2.97e-3 |
| acc_forward | 4.20e-3 | 1.60e-2 (3.8x) | 1.62e-2 |
| logdet | -89253588.194065422 | -89253748.104239851 (1.8e-6) | -89253736.570008814 |

No uncompressed CA0 run at 192^3 exists (the 22.51 GB of 2026-10-08 was
compressed with the near store off, like the 20.72 of Color).
(The 20.72 GB / 7.49 s quoted as "before" in the running tally is the
compressed run with the near store off, not an uncompressed one.)

#### FP32 does not shrink the stores (2026-10-09)

Measured before planning around it (h2_c1_work/prec_bench.cu on the dumped
level-4 wave: 37 real X_NR boxes, 435.75 MB as FP64), bitcomp lossy from
FP64 and from an FP32 copy of the same data, the same delta a box, the H2
path's 64K-element chunks:

| step | type | packed MB | vs FP64 raw | bits/el | comp GB/s (own bytes) | dcmp GB/s | max_err |
|---|---|---|---|---|---|---|---|
| 2e-4 (X_NR, near) | FP64 | 22.91 | 19.0x | 3.36 | 428 | 593 | 6.1e-5 |
| 2e-4 | FP32 | 22.51 | 19.4x | 3.31 | 447 | 636 | 6.1e-5 |
| 2e-6 (T, X_SR) | FP64 | 68.37 | 6.37x | 10.0 | 415 | 503 | 4.8e-7 |
| 2e-6 | FP32 | 67.92 | 6.42x | 10.0 | 421 | 480 | 4.9e-7 |

The packed bytes are the quantized integers' bits - round(x / delta) is the
same integer whichever type x came in - so FP32 input gives 1.7% fewer bytes
and the same error.  What FP32 halves is the codec's time an element (the
same GB/s on half the bytes) and every buffer that exists raw.  The "FP32
halves the compressed stores" written above (the 8x-goal remarks of
2026-10-08 and the kept-block remark of 2026-10-09) is therefore wrong, and
FP32 is a smaller item than it was taken for: the transient and the
decompressed working sets (~0.5-0.7 GB of the leaf's 7.41), the raw parent
blocks under level 4 (0.94 -> 0.47), about half of the codec passes' time
(near packing 0.50 s, handover 0.15 s, the windows: perhaps 0.5 s a leaf),
and the GEMMs if TF32 is acceptable; the LU and its solves stay FP64.  The
route to 8x and 64x boxes a rank is host staging of the stores.

#### CA3 on the current code (2026-10-09): correct, fastest, and three reductions short

CA3 (`--H2_CA_owner_component 3`, CA at the leaf only) had not been run
since 2026-10-06, when two CA levels exhausted the heap.  On commit 66a1002
it runs clean on both sizes, 8 ranks:

| | CA3 off 96^3 | CA3 on 96^3 | CA3 on 192^3 | Color on 192^3 | CA0 on 192^3 |
|---|---|---|---|---|---|
| peak (resident + transient) | 7.02 (6.38 + 0.64) | 4.42 (3.27 + 1.16) | **20.67** (14.01 + 6.66) | 7.41 | 8.74 |
| level-4 peak | 4.33 | 2.23 | 10.00 (8.94 + 1.05) | 6.98 | 7.67 |
| kept store, rank 0, leaf | 2.36 (1.0x) | 1.06 (2.2x) | 5.47 of 19.90 (3.6x) | 2.28 (8.3x) | 3.05 (8.6x) |
| Factor | 1.30 s | 2.30 s | **10.5 s** | 11.2 s | 11.6 s |
| GPU solve | 0.117 s | 0.131 s | 0.471 s | 0.464 s | 0.499 s |
| acc_mvp / acc_forward | 7.94e-4 / 3.46e-3 | 1.43e-3 / 9.32e-3 | 3.18e-3 / 1.61e-2 | 2.95e-3 / 1.60e-2 | 2.97e-3 / 1.62e-2 |
| acc_backward | 6.5e-16 | 6.9e-16 | 1.4e-15 | 1.1e-15 | 5.3e-12 |
| logdet | -9928892.253 | -9928899.530 | -89253739.706 | -89253748.104 | -89253736.570 |

Correct (logdets within 1.3e-7 of the other two paths at 192^3, the same
accuracy, no errors), the fastest factorization of the three, and - unlike
CA0 - no replication: 3532 boxes pledged plus 564 landed is the rank's 4096,
and its near store is Color's 3.24 GB.  Its memory is 2.8x Color's because
three of the reductions never reached the owner-component path, and the
leaf's breakdown prices each:

1. **The production is not chunked, 5.6 GB transient** (wave outputs 4.50,
   wave work 1.09): chunk_production requires `owner_ == nullptr`, because
   pack_owner_generators packs a wave's generators once a wave and does not
   separate into plan and copy (the Color transport's split of 2026-10-08).
2. **The remote generators' temp2^T stays raw, 5.81 GB live (6.09 allocated
   over the level)**: the 10-09 group staging runs, but its raw fallback
   fires for essentially every generator here (1121 of ~5.4 MB), where
   Color's level 4 packs all but 0.15 GB.  Why the packer finds nothing to
   pack for them (fill_item's generator-box lookup, or no room) is the first
   thing to establish for the fix.
3. **C2, the landed FULL factors are raw: 4.5 GB in the kept store** (564
   boxes; 5.47 GB stored at 3.6x where the pledged 3532 reach 15x), and it
   sits under every later level (level 4's kept 5.49, level 2's 6.35).
   land_full_factors is called once per entry; the receive buffers live to
   the end of owner_exchange, so a per-sequence batch through the kept-block
   compressor is possible.

With the three in place the CA3 leaf should sit near Color's (the same
stores, no replication): roughly 20.7 -> ~7 GB and level 4 10.0 -> ~5.8,
keeping the 10.5 s.  That would make CA3 the path with the best memory *and*
time, which is a reason to do these before the host staging rather than
after.

#### CA3's three, done (2026-10-09): 20.67 -> 7.91 GB

1. **The owner-component production is chunked.**  pack_owner_generators is
   now plan_owner_generators (once a wave: headers, pieces, the block of the
   arena, the send's messages), copy_owner_generators (a chunk's pieces into
   it) and finish_owner_generators (the ready event, after the loop - a last
   chunk with no boxes is skipped inside it, which also left the wave's X_NR
   origin block unfreed; both now happen after the loop).  chunk_production
   no longer requires `owner_ == nullptr`.
2. **The remote generators' temp2^T is packed on an owner level.**  fill_item
   looked the generator up in the generator table, which has only remote
   boxes; the rank's own boxes landed FULL were not found and took the raw
   fallback, every one of them.  It uses the box found when the generator was
   set up (`g.box`).  Generator boxes allocated raw at the leaf: 6.09 -> 0.28 GB.
3. **The landed FULL factors are compressed.**  The kept-block builder is
   keep_factor_sources over KeptSource (morton, state, the five pointers,
   shapes, pledge): a wave's boxes come through keep_solve_factor_group from
   the result and source blocks, a sequence's FULL entries through
   land_full_factors straight from their payloads (no pledge, no X_RR_full
   rebuild from this LU: that came with the payload), in groups of
   kXnrStagingBudget.  Kept store at the leaf: 5.47 GB (3.6x) -> 2.33 (8.5x).

| 8 ranks | CA3 96^3 before | after | CA3 192^3 before | **after** | Color | CA0 |
|---|---|---|---|---|---|---|
| peak (resident + transient) | 4.42 (3.27 + 1.16) | 1.97 (1.10 + 0.86) | 20.67 (14.01 + 6.66) | **7.91** (6.87 + 1.04) | 7.41 | 8.74 |
| level-4 peak | 2.23 | 1.45 | 10.00 | 6.86 | 6.98 | 7.67 |
| kept store, rank 0, leaf | 1.06 (2.2x) | 0.29 (8.2x) | 5.47 (3.6x) | 2.33 (8.5x) | 2.28 (8.3x) | 3.05 (8.6x) |
| Factor | 2.30 s | 2.66 s | 10.5 s | 11.2 s | 11.1 s | 11.5 s |
| GPU solve | 0.131 s | 0.138 s | 0.471 s | 0.496 s | 0.475 s | 0.505 s |
| acc_mvp / acc_forward | 1.43e-3 / 9.32e-3 | 1.43e-3 / 9.29e-3 | 3.18e-3 / 1.61e-2 | 2.97e-3 / 1.54e-2 | 2.95e-3 / 1.60e-2 | 2.97e-3 / 1.62e-2 |
| logdet | -9928899.530 | -9928899.446 | -89253739.706 | -89253742.145 | -89253748.104 | -89253736.570 |

Color and CA0 are bit-identical to the runs before (logdet, kept bytes,
peaks, accuracies), which gates the refactor.  CA3's logdet moved 8.4e-9 /
2.7e-8 because the later waves' sketches now read the generators' packed
temp2^T, as Color's did when its fill sources were compressed; its accuracy
is unchanged (96^3) or better (192^3), inside the three paths' spread.  The
time went where the chunking put it on Color: Factor 10.5 -> 11.2 s, now
Color's.  What is still CA3's own: the landed boxes' kept X_NR is compressed
from the payload rather than pledged (their fill pack comes later in the
sequence than the landing), about 0.2 GB of kept factors (1.33 against
Color's 1.12), and their X_RR_full arrives raw with the payload rather than
rebuilt (fill X_RR_full 0.45 against 0.28).

**The three paths at the 192^3 leaf**: Color 7.41 GB / 11.1 s, CA3 7.91 /
11.2, CA0 8.74 / 11.5.  Host staging now applies to all three alike.

#### Host staging (a), done (2026-10-09): the kept factors leave the device between levels

At each level's start (start_level, start_owner_level, start_ca_level) every
deeper level's kept store goes to host memory (stage_kept_level): the kept-LU
logdet is summed on the device first and kept per box, so
device_logdet_kept_level adds it later in the same order and the logdet stays
bit-identical; the blocks are copied end to end with page-aligned starts, then
freed; the boxes keep the device addresses the blocks had.  The first solve's
set-up brings each level back in one block and rebases the pointers
(upload_kept_level); a host-copy reader in between goes through kept_fetch;
the solve's eligibility counts the staged bytes as needed again.  Pageable
memory, synchronous copies: pinning ran at 1.9 GB/s (1.2-1.6 s a leaf store
where the copy took 0.09), a pageable copy stages at ~7 GB/s out and ~4.6 GB/s
back, once.

| 192^3, 8 ranks | Color | CA3 | CA0 |
|---|---|---|---|
| leaf peak | 7.41 -> 7.41 | 7.91 -> 7.91 | 8.74 -> 8.74 |
| level 4 peak | 6.98 -> **4.71** | 6.86 -> **4.54** | 7.67 -> **4.64** |
| level 3 / level 2 peak | 5.30 / 4.72 -> 2.31 / 1.54 | 5.31 / 4.74 -> 2.29 / 1.52 | 6.06 / 5.46 -> 2.32 / 1.52 |
| staged, rank 0 | 3.29 GB | 3.26 GB | 3.99 GB |
| staging copies | 0.45 s (+0.10 alloc) | 0.41 s | 0.52 s (+0.08) |
| Factor | 11.1 -> 11.9 s | 11.2 -> 11.8 s | 11.5 -> 12.6 s |
| solve set-up, first solve | 0.11 -> 0.35 s | 0.11 -> 0.35 s | 0.13 -> 0.41 s |
| time_solve | 0.476 -> 0.473 | 0.496 -> 0.504 | 0.505 -> 0.515 |

logdet and the accuracies are bit-identical on all three paths at both sizes
(96^3 heap-checked).  The leaf is now the only level that binds, by 2.7-4.1
GB; under it nothing changed, which is what (b) and (c) are for.  Factor's
+0.6-1.1 s is the copies (0.4-0.6 s) plus the stream sync the staging puts on
the critical path at each level start; an asynchronous copy with
stream-ordered frees would overlap the next level's host planning, a
follow-up if the time matters.

**The bug, and why ten runs were needed to see it.**  The first gate failed
the solve's plan with bitcomp -2 (invalid data) on a level-3 wave, on a
wave that moved from run to run, while every check said the staging was
right: the copy was byte-identical to the device, the accounting matched the
heap to the byte, every chunk planned alone at its original address, every
chunk's head matched the host copy after the upload, the uploaded block was
intact before and after the failed plan, no pointer lay outside it, and
bitcomp (offline, on real chunks) accepted any pointer order and any
relocation.  All of those compared the staged copy *through the same
mapping* as the solve, so a mapping error passed them all.  What caught it
was an independent record: the boxes' pointers at staging against those at
upload - five of 64 boxes had T, X_SR and the LU rewritten by one offset.
The driver's between-level compaction relocates every kept store's pointers
for the blocks the heap moved; a staged level's pointers are the addresses
its freed blocks had, and those ranges were reused by blocks that then moved.
The compaction now skips staged levels.

#### Host staging (b), done (2026-10-09): the packed fill store lives on the host

The slot-compressed temp2^T of every box (the fill store, 1.47-1.58 GB at the
leaf) now lives in pinned host memory.  No staging of its own was needed:
bitcomp only asks that its compressed input be device-visible, and pinned
host memory is (UVA), so the later readers - the sketches of later waves,
the transition, the handover - decompress straight from the host copy,
reading the small compressed input over PCIe while writing the output on the
device.  The packer copies each packed block to a process-wide arena
(fill_host_arena, bump allocation, reset at the level's end) on a stream of
its own and records an event; the device copy stays until the wave's owner
pass has read it and goes at the wave's end once the copy is done
(evict_fill_copies).  Every reader takes its base from packed_of(st): the
device copy while it exists, else the host one.  Counted at the leaf with
BPACK_TRACE=wave: the owner window reads 128609 slot chunks from the device
and none from the host, the sketches 342873 from the host and 6586 from the
device, and the sketch costs 1.49 -> 1.59 s for it - the PCIe-fed
decompression is nearly free.

What was not free was the pinning.  cudaMallocHost pins at 1.9 GB/s: 0.85 s
for the leaf's 1.6 GB, spent inside the packer, which stalled the host while
the device waited (+0.55 s in "plan boxes" and the same in the owner pass's
device bracket - idle, not work).  Measured on a node (pin_bench.cu):
transparent huge pages are "always" there, malloc'd memory touched and
registered with cudaHostRegister pins 1.6 GB in 0.07 s at the same 26 GB/s.
The arena does that.

| 192^3, 8 ranks | Color | CA3 | CA0 |
|---|---|---|---|
| leaf peak | 7.41 -> **6.10** (4.43 + 1.67) | 7.91 -> **6.57** (4.90 + 1.67) | 8.74 -> **7.33** (5.66 + 1.67) |
| fill temp2^T packed on the device, at the peak | 1.47 -> 0.12 | 1.51 -> 0.12 | 1.58 -> 0.12 |
| fill store on the host, at most | 1.47 GB | 1.51 GB | 1.58 GB |
| Factor | 11.9 -> 12.3 s | 11.8 -> 12.4 s | 12.6 -> 12.7 s |
| time_solve | 0.473 -> 0.472 | 0.504 -> 0.509 | 0.515 -> 0.516 |

logdet and the accuracies are bit-identical on all three paths at both sizes
(96^3 heap-checked).  The leaf's peak moved to a moment with a larger
transient (1.03 -> 1.67: the source window and the sketch vectors), which is
where the budgets come back into play.  What the leaf holds now, Color:
owner near blocks 3.24, transient 1.67, kept factors 0.57, X_RR_full 0.30,
owner source slots 0.15, fill store 0.12.  The near store is three quarters
of the resident, which is (c).

#### Host staging (c), the near store: review before building (2026-10-09)

**What the store is** (step 3c, as built).  The level's Schur complement:
one block per one-hop pair and the Schur block of each box, 57203 blocks of
17.23 GB raw at their largest on the Color leaf, held packed at 5.3x.  A
block's packed bytes are a run of *chunk homes* in B_compact (extents of
32 MB doubling to 512 MB, bump-allocated, a home sized at the chunk's first
compressed length and never moved) plus a *tail slot* (64 KB, in 16 MB
extents) for what a re-compression adds beyond the home.  Resident at the
leaf peak: 3.24 GB Color and CA3, 4.31 CA0 - three quarters of what is
resident now.

**Who writes it.**  The owner pass, a group at a time (kNearGroupBudget
256 MB of raw targets): the group's raw blocks are compressed into the
compressor's device staging and copy_spans scatters the chunks into their
homes and tails, then the raw blocks are freed.  Per leaf level: 451501
block packs in 326 groups, 10.46 GB of packed bytes written (CA0 15.1 GB
in 550 groups), 0.47 s, 8.9 writes a block; 83% of the blocks go final only
in the last quarter of the waves, so nothing can leave early.

**Who reads it.**  One path, materialize_adopted: allocate the raw blocks,
gather each chunk's home and tail into a contiguous "store gathers"
staging (copy_spans), decompress in one batch.  Per leaf level, Color:
update groups 34.9 GB in 181 batches, generators 2.2, build_parent 0.8 -
38 GB decompressed from 7.2 GB of packed input (CA0: 54 GB from 10.3).  The
transition's build_parent then packs the parent blocks as contiguous runs
(packed_owned), which the next level adopts: "adopted blocks" 1.87 GB
under level 4.

**The proposal.**  Put the extents - B_compact's, B_newtail's, and
build_parent's packed_owned runs - in pinned host memory (the arena of (b)
is the model: malloc + THP + cudaHostRegister, device-visible), and change
nothing in the structure: compact_reserve and tail_grow take from the
arena; copy_spans scatters a group's packed chunks into the host homes
(zero-copy writes over PCIe); materialize_adopted's gather reads the host
homes into the device staging (zero-copy reads), and the decompression is
as now.  The level's region resets when its eliminator goes; the parent
runs live in the next level's region, since that level reads them.  The
between-level compaction must leave host pointers alone (the lesson of
(a)).  About 150 lines, all in the allocator and the two copy sites.

**What it buys.**  Leaf, Color: resident 4.43 -> ~1.2, peak 6.10 -> ~2.9 GB;
CA3 6.57 -> ~3.3; CA0 7.33 -> ~3.0.  Level 4 with the adopted runs on the
host: 4.68 -> ~2.3.  Against two days ago that is 20.72 -> ~2.9 on Color,
7x, and the 8x-boxes goal at 192^3 (3 GB x 8 = 24 GB of the card) comes
within reach.

**What it costs, estimated from measured rates.**  Writes: 10.5 GB of
scattered 4-byte-word kernel writes to host memory a leaf level (CA0 15
GB).  Reads: 7.2 GB (CA0 10.3) gathered the same way.  At a PCIe-4 x16
link's ~20 GB/s each way that is 0.35-0.5 s of reads and 0.5-0.75 s of
writes a leaf level, +0.9-1.3 s on an 12.3 s factorization if none of it
overlaps - and today none would: the packer synchronizes after each
group's scatter so the compressor's staging can be reused, and the gather
precedes each group's decompression.  Double-buffering the staging and
running the scatter on the copy stream would hide most of the writes; a
follow-up, not the first version.  The solve is untouched.

**The unknown to measure before building** (the rule from the bitcomp
benchmark: the kernel where it runs, not a guess): the bandwidth of
copy_spans's access pattern - 4-byte words, 64-512 KB spans, coalesced per
warp - against mapped host memory, device-to-host and host-to-device, on
one GPU and with four ranks a node at once.  PCIe is efficient at 128-byte
transactions and poor below; if the pattern lands well under 10 GB/s, the
alternative is a contiguous device staging per group with one bulk copy
each way and a host-side scatter (the background copier's thread has
spare time).  Thirty minutes with pin_bench.cu's harness.

**Measured (bw_bench.cu, 2026-10-09).**  copy_spans's pattern against mapped
host memory runs at the link's rate: 26.3-26.4 GB/s writing to the host and
26.7 GB/s reading from it, for spans of 64, 256 and 512 KB, sequential or
scattered, 4-byte or 16-byte words - the same as cudaMemcpy - and all four
ranks of a node get it at once, each GPU on its own link.  So the estimate
above holds at its better end: 0.40 s of writes and 0.27 s of reads a Color
leaf level (CA0 0.57 and 0.39), +0.7-1.0 s unoverlapped, and no host-side
scatter is needed.

**The honest limit.**  Host memory per rank is 64 GB on these nodes.  At
8x the boxes a rank the three staged stores come to about 8.5 GB x 8 = 68
GB a rank - the edge of what fits, with the near store 26 GB of it.  At
64x they do not fit on the host either; that regime needs a smaller
representation or more ranks, not staging.

#### Host staging (c), done (2026-10-09): the near store in pinned host memory

Built as reviewed, nothing restructured: compact_reserve and tail_grow take
their extents from a pinned host arena - two arenas, by level parity, so the
parent runs build_parent writes for the level above live in that level's
arena and each arena resets when its level's eliminator goes; copy_spans
scatters a group's packed chunks into the host homes and gathers them back
into the device staging in front of each decompression; the between-level
compaction and the frees skip host pointers (heap_.owns).  The chunk-home
structure, the groups and the readers are as they were.

Gate: six runs, logdet bit-identical to (b) on every one (Color
-89253748.104239851, CA3 -89253742.144751042, CA0 -89253736.570008814 at
192^3; -9928903.0965428073 / -9928899.4463859536 / -9928897.9223398492 at
96^3), acc_forward and acc_backward unchanged, heap-checked at 96^3, 0 real
errors.  Peaks as "peak (resident + transient)" where the moment matters:

| 192^3, 8 ranks | leaf peak, (b) -> (c) | level 4 | level 3 | level 2 | Factor | time_solve |
|---|---|---|---|---|---|---|
| Color | 6.10 (4.43+1.67) -> **5.01** (3.46+1.55) | 4.68 -> 3.20 | 2.30 -> 1.90 | 1.54 -> 1.41 | 12.3 -> 14.6 s | 0.472 -> 0.482 |
| CA3 | 6.57 (4.90+1.67) -> **4.14** (2.45+1.69) | 4.51 -> 3.06 | 2.27 -> 1.88 | 1.52 -> 1.40 | 12.4 -> 14.5 | 0.509 -> 0.499 |
| CA0 | 7.33 (5.66+1.67) -> **3.01** (1.34+1.67) | 4.62 -> 3.14 | 2.31 -> 1.91 | 1.52 -> 1.40 | 12.7 -> 15.1 | 0.516 -> 0.507 |
| Color, no bitcomp (2026-10-06, A_f192off / A_v192off) | 33.21 (26.18+7.02) | 16.04 | 10.88 | 8.39 | 7.54 s (7.71 in its twin) | 0.263 (18.92 GB uploaded at the first solve, set-up 2.25 s) |

The uncompressed row is the only one there is at 192^3: on 8 ranks the
uncompressed leaf has died on the heap since the kept store exists (see
"With and without bitcomp"), and CA0 and CA3 never ran uncompressed at this
size.  Its accuracy is the standing comparison, acc_forward 4.20e-3 against
1.60e-2 compressed (acc_mvp 8.11e-4 vs 2.95e-3); logdet -89253588.194 vs
-89253748.104 (1.8e-6).  Against it, (c) is 6.6x less peak memory on Color
for 1.9x the Factor time and 1.8x the solve.

| 96^3, 8 ranks, heap-checked | level 4, (b) -> (c) | level 3 | level 2 | Factor |
|---|---|---|---|---|
| Color | 1.82 (1.06+0.76) -> 1.56 (0.80+0.76) | 1.23 -> 1.05 | 0.84 -> 0.78 | 2.24 -> 2.60 s |
| CA3 | 1.79 (0.94+0.86) -> 1.47 (0.66+0.81) | 1.16 -> 1.00 | 0.84 -> 0.77 | 2.67 -> 3.13 |
| CA0 | 2.21 (1.35+0.86) -> 1.12 (0.26+0.86) | 1.19 -> 1.02 | 0.84 -> 0.78 | 2.83 -> 3.40 |

On the host at 192^3: near store 3.24 GB at most on the Color and CA3 leaf,
4.31 on CA0, 2.26-2.29 at level 4; the pinned arenas come to 7.78 GB a rank
on Color and CA3 and 8.86 on CA0, fill store included.  Against two days
ago the Color leaf is 20.72 -> 5.01 GB, 4.1x; CA0 22.51 -> 3.01, 7.5x.

**What the peaks are now.**  CA0 landed where the review said (3.0 GB):
resident 1.34, of which 0.91 is the leaf's own kept factors, which (a)
stages out only at the next level's start.  CA3 at 4.14 is 0.8 above the
estimate: resident 2.45 = fill X_RR_full 0.47 + temp2^T packed 0.33 + the
level's own Schur/edge blocks 0.62 + kept factors 0.67.  Color at 5.01 is
2.1 above it: with the store gone from the peak moment, the next-highest
moment is a wave chunk early in the level holding 2.88 GB of *level
Schur/edges* - step 5's fresh near blocks of the chunk (A_NS + update, c x k
per neighbor of every box in the chunk, ResidentKind::LevelSchur), allocated
raw when the chunk is planned and packed only by the owner pass's groups
behind it.  That is the next item on Color, and its knob is the wave chunk
(kWaveChunks = 8 a wave, regardless of bytes).  Level 4 is 3.1-3.2 on all
three paths, with the materialized adopted blocks 1.76-1.87 GB of it - the
raw parent blocks a wave's read set restores - and the transient half
unchanged at 1.05-1.09.

**What it costs, measured: +2.2-2.4 s, not the +0.7-1.0 estimated.**  The
transfers themselves run at the link's rate.  On the Color leaf "near
blocks held packed" went 0.47 -> 1.37 s for 10.46 GB written (0.40 s of
that at 26 GB/s) and "adopted blocks materialized" 0.30 -> 0.55 for 7.2 GB
read (0.27 s); level 4: 0.12 -> 0.47 and 0.09 -> 0.19.  Two things were
left out of the estimate:

1. *Pinning.*  The arenas grow 256 MB at a time, synchronously, inside the
   packing: posix_memalign + THP fault-in + cudaHostRegister at 0.15 s a
   GB.  The near arenas pin 6.2 GB of new memory a rank (7.78 less the fill
   store's 1.61): "pinning 0.20 s so far" became 1.15 s (CA3 1.08, CA0
   1.27), all on the critical path - about 0.95 s of the 2.3.  It is paid
   once a process, but the benchmark factorizes once.
2. *The traffic at the other sites*: the plan site's store gathers for the
   box buffers (+0.22 s in "plan boxes [buffers]"), build_parent's parent
   runs written into the next level's arena in the transition (+0.22 s),
   level 4's reads and writes above (+0.45) - about 0.9 s more at the link
   rate, none of it overlapped: the packer still synchronizes after each
   group's scatter so the compressor can reuse its staging, and the gather
   precedes each decompression.

Where it goes next, cheapest first: (c2) pin ahead of demand - a worker
thread pins the next chunk while the current one fills, and pins further
ahead whenever a request had to wait; ~40 lines in PinnedArena, and it hides
most of the 0.95 s.  Then, if wanted, the scatter on a side stream
overlapping the next group's GEMMs (a dedicated span list, an event before
the compressor reuses its staging; the groups of a wave partition its
targets, so the next group's gather has to wait for the scatter only across
a wave boundary) - worth about 0.5 s at 192^3, more lines than (c2) for
less.  The solve is untouched (0.48-0.51 s).

#### (c2), done (2026-10-09): the arenas pin ahead of demand

PinnedArena pins its next chunk on a worker thread (std::async) as soon as
one is taken, and the driver warms all three arenas at each level's start
so the first growth is ready before it is asked for; when the chunk a
request takes is not ready yet, the arena pins twice as far ahead from then
on (kAheadMax 4).  The report line now says how much of the pinning the
levels actually waited for.  Gate: six runs, logdets, accuracies and every
level's peak identical to (c) (the arenas hand out the same addresses in
the same order), 0 errors.

| Factor, (c) -> (c2) | 192^3 | 96^3 | pinning / waited at 192^3 |
|---|---|---|---|
| Color | 14.6 -> 14.3 s | 2.60 -> 2.51 | 1.24 s / 0.36 s |
| CA3 | 14.5 -> 14.0 | 3.13 -> 3.03 | 1.34 / 0.42 |
| CA0 | 15.1 -> 14.9 | 3.40 -> 3.33 | 1.52 / 0.58 |

Less than the 0.95 s because 0.36-0.58 s is still waited for, and the
report shows where: 0.34 of Color's 0.36 is accrued by the leaf.  The near
store's extents double from 32 MB up to kCompactMax = 512 MB, and a request
larger than the arena's 256 MB chunk cannot take a chunk pinned ahead - it
pins synchronously, as before.  A 3.24 GB leaf store has five or six 512 MB
extents at 0.08 s each: the 0.34 s.

**(c3), measured and reverted (2026-10-09): 512 MB chunks.**  With the
chunk at the largest extent every extent takes a chunk pinned ahead, and the
levels' waits did fall - to 0.00 s on Color, 0.13 CA3, 0.07 CA0 at 192^3 -
but the pinning work rose from 1.24 / 1.34 / 1.52 s to 1.91 / 1.95 / 2.06 s
(8.59 GB pinned instead of 7.78, and concurrent pinners slow each other),
and the Factor time did not follow the waits: 14.3 -> 14.2 s Color, 14.0 ->
14.2 CA3, 14.9 -> 15.0 CA0, and at 96^3 it rose 2.51 -> 2.78, 3.03 -> 3.29,
3.33 -> 3.51 while the pinned total went 1.61 -> 2.15 GB.  Logdets and peaks
identical in all six runs.  So pinning on a worker thread is not hidden:
what the main thread stops waiting for, it loses as slower progress while
the worker pins - at 96^3 one for one (+0.30 s of pinning, +0.27 s of
Factor).  The 256 MB chunk stays; the remaining 0.36-0.58 s are the price of
registering 7.8 GB in a process that factorizes once.

**The mechanism, measured (pin_bench2.cu, one A100, 2026-10-09).**  A worker
thread pins 2 GB in 256 MB chunks the way PinnedArena does while the main
thread launches.  The THP fault-in (touch) alone takes 0.156 s and the main
thread keeps 100% of its rate (423k tiny launches/s, or 156 6-ms kernels/s).
cudaHostRegister on touched memory takes 0.093 s, and while it runs the main
thread gets 2% of its launch rate and 12% of its kernel rate: the
registration blocks the process's other CUDA calls almost entirely.  Both
together, 0.250 s at 60-66%.  So of the 0.15 s a GB, the 0.08 of fault-in
hides behind the GPU work and the 0.046 of registration does not, whichever
thread does it: 7.8 GB registered is 0.36 s of the factorization - what (c2)
left, and what (c3) could not take away.  What would: registering less (the
arenas are within 10% of the stores' maxima), or registering while the GPU
holds a long queue of launched work, which the group-synchronous owner pass
never has.  A process that factorizes twice pays it once.

**Where (c) stands after (c2)**, against (b): leaf peak 6.10 -> 5.01 Color,
6.57 -> 4.14 CA3, 7.33 -> 3.01 CA0; level 4 4.5-4.7 -> 3.1-3.2; Factor
+2.0 / +1.6 / +2.2 s at 192^3 (14-17%), of which ~0.36 s is registration,
~1.0 s PCIe traffic at the link rate (the scatter's 0.5 s could overlap the
GEMMs on a side stream, the rest precedes what reads it), and the rest the
plan site's and the transition's reads and writes; +0.27 / +0.36 / +0.50 s
at 96^3.  Solve unchanged.

#### Asynchronous host staging, (c) (2026-10-09): the scatter under the next group, the gather prefetched

Built as proposed.  The near store's scatter (device staging -> host homes)
runs on a near copy stream of its own: the compression synchronizes the
compute stream, so the staging is complete; the scatter is issued with its
own span list and the raw blocks go back at once; the compressor holds its
staging (hold_staging) until the scatter's event, so the next compress
cannot free it under the copy.  A reader waits only for the scatter that
last wrote one of its blocks (PackedBlock::written_seq against the
scatters the host knows complete), which inside a wave is never - a wave's
groups partition its targets - and across a wave boundary is one GPU-side
event wait.  The group loop prefetches the next group's gather (host homes
-> device scratch) on the near copy stream under the current group's
GEMMs, after the compute stream passes its last decompression (an event,
so scratch the heap hands out again is safe), and decompresses from it
when the group's turn comes; materialize_adopted is now plan_gather +
materialize_from_gather, the other sites unchanged.  The level's end waits
for the last scatter before its region resets.

Gate: six runs, logdets, accuracies and every level's peak identical to
(c2), 0 errors; at the 192^3 leaf 166 of 189 group gathers were prefetched
on Color (160 of 204 CA3, 256 of 288 CA0), 50 of 75 at level 4.

| (c2) -> async | Factor 192^3 | Factor 96^3 | leaf "held packed" | leaf "materialized" | leaf group loop |
|---|---|---|---|---|---|
| Color | 14.3 -> 14.1 s | 2.51 -> 2.47 | 1.24 -> 0.87 s | 0.557 -> 0.500 s | 2.59 -> 2.33 s |
| CA3 | 14.0 -> 13.7 | 3.03 -> 3.00 | 1.28 -> 0.91 | 0.586 -> 0.536 | 2.61 -> 2.32 |
| CA0 | 14.9 -> 14.5 | 3.33 -> 3.28 | 1.66 -> 1.11 | 0.736 -> 0.762 | 3.51 -> 3.26 |

The scatter is off the critical path entirely: the packing fell by 0.37 /
0.37 / 0.55 s where (c)'s estimate of its PCIe time was 0.40 / 0.40 / 0.57,
and what it still costs over (b)'s 0.47 / 0.51 / 0.56 is the arenas'
pinning waits (0.37-0.52 s by the leaf's end).  The gather prefetch bought
0.06 s of the 0.25 expected.  The likely reason: the gather is a kernel
(copy_spans), not a copy-engine transfer, and under a batched GEMM that
fills the SMs its blocks are scheduled only as GEMM blocks retire, so it
runs slowly until the GEMMs are done and the next group waits for it
anyway.  Two things to try, in order: the near copy stream at the highest
priority (its blocks go ahead of the GEMMs' as SMs free up; one line), and
cudaMemcpyBatchAsync (CUDA 13.2) for the gather and the scatter, which uses
the copy engines and no SMs at all.

#### Asynchronous host staging, (a) (2026-10-09): measured, reverted

Built: the kept factors' host copy in a pinned arena of its own
(kept_host_arena, pinned ahead like the others), the copies issued on a
kept copy stream at the level's start after the blocks' last writes, the
device blocks returned once the copies had landed (finish_kept_staging,
after the eliminator's construction, before the waves; the heads compared
there), the solve's upload from pinned memory.  Gate: logdets and
accuracies identical, 0 errors.  What it measured, 192^3:

| (c) async -> (a) async | Factor | level-4 start, 2.3-3.0 GB staged | solve set-up | pinned arenas | level 4 peak |
|---|---|---|---|---|---|
| Color | 14.1 -> 14.1 s | copy 0.32 s -> allocation 0.44, copy hidden | 0.35 -> 0.23 s | 7.78 -> 10.05 GB | 3.20 -> 3.20 |
| CA3 | 13.7 -> 13.9 | 0.32 -> 0.33 | 0.34 -> 0.23 | 7.78 -> 10.11 | 3.06 -> 3.06 |
| CA0 | 14.5 -> 14.7 | 0.42 -> 0.43 | 0.39 -> 0.26 | 8.86 -> 11.90 | 3.14 -> **3.51** |

Three things, all foreseeable from what was already measured.  The window
the copy can hide under is the eliminator's construction, 0.01 s ("begin")
- the kept store is complete only at the level's end (the pledged X_NR is
handed over in the transition) and the next level's first wave is its
sketch - so the copy is not hidden, it is only faster (26 GB/s pinned vs
13 pageable), and that saving is paid back by pinning the arena: a 2.3 GB
request is bigger than a chunk pinned ahead and pins synchronously, 0.33-
0.44 s, where the pageable copy cost 0.32-0.42.  Pinned memory a rank goes
up by the kept store, 2.3-3.0 GB.  And on a CA level the start is not
empty: the device halo exchange allocates its received pieces while the
leaving blocks are still on the card, 3.14 -> 3.51 GB at CA0's level 4 -
the one place the overlap existed, it cost memory.  Reverted; the solve
set-up's 0.12 s is not worth the rest.

What would make (a) asynchronous for real: copy each kept block as it is
produced, under the waves, the way the fill store's producer already does
- the T, X_SR, LU and pivots blocks are final per group during the level,
the pledged X_NR blocks per group of the transition's handover - with the
device copy freed when the copy lands.  That hides the whole copy, needs
no level-start window, and takes the level's own kept factors off the
card: 0.67 GB of CA3's and 0.91 GB of CA0's leaf peak (Color's peak moment
holds 0.08).  It changes the kept store's lifecycle (every mid-level
reader of a kept block: rescue_skipped_xnr, the logdet, the handover) and
is a day's work with its own gate; not started.  A cheaper partial, a
shared pool of pinned chunks so the kept store takes the leaf's near
arena's chunks instead of pinning its own, would recover ~0.3 s and the
extra 2.3-3.0 GB of pinned memory, but not hide the copy.

#### The near copy stream by the copy engines (2026-10-09): measured, reverted

Why the prefetched gather bought 0.06 s of 0.25, measured
(h2_c1_work/gather_bench.cu, one A100, driver 580.178): 2000 spans of 4-256
KB, 124 MB, between pinned host memory and the device on a stream of their
own.  The copy kernel runs at 26.7 GB/s alone and at **2.6 GB/s** under a
kernel that fills the SMs on another stream - its blocks run only as the
other kernel's retire - and the highest stream priority changes nothing
(it orders pending blocks; it does not preempt resident ones).
cudaMemcpyBatchAsync (CUDA 12.8+, the copy engines) runs at 25 GB/s alone
and under the load and slows the loaded kernel by nothing - and costs
0.8-1.0 ms of host time per call of 2000 entries, 0.5 us an entry, where
the kernel's issue is a list upload and a launch.

Built (the scatter and the prefetched gather through cudaMemcpyBatchAsync,
guarded by CUDART_VERSION, the kernel otherwise) and gated against the
asynchronous (c) logs: the three 96^3 runs and CA3 at 192^3 bit-identical;
Color and CA0 at 192^3 died in the leaf's fourth wave with "unspecified
launch failure", a device fault surfaced by the fill store's event query.
Not diagnosed, because the surviving run disqualifies the path on its own:
CA3's leaf packing 0.91 -> 1.14 s and materialization 0.54 -> 0.77, +0.46 s
a leaf level, which is the host cost of its 606 calls of ~2000 entries.
The owner pass's groups hold ~2000 blocks of ~20 KB, one chunk each, so the
entry count is what it is.  The copy engines are the right hardware at the
wrong granularity.  Reverted to the kernel version (0dfedec); the stream
priority went with it.  (If this is taken up again: the first suspect for
the fault is the lifetime of the call's host arrays - they were
temporaries - and the way to the entry count is coalescing adjacent homes:
the homes are allocated in the targets' canonical order and the groups
walk that order every wave, so consecutive blocks of a group are often
adjacent in B_compact and in the gather scratch alike.)

**Where asynchronous host staging stands, 2026-10-09.**  (b)'s producer
copies were asynchronous from the start.  (c)'s scatter is hidden
entirely: the packing fell by its whole PCIe time.  (c)'s gather prefetch
overlaps only what leaves SMs free.  (a) cannot be hidden at the level's
start (a 0.01 s window) and a pinned copy pays its pinning; the version
that would hide it copies each kept block as it is produced.  The
registration of the arenas (0.36-0.6 s) cannot be hidden by any thread.
Factor at 192^3 against the synchronous (c2): Color 14.3 -> 14.1 s, CA3
14.0 -> 13.7, CA0 14.9 -> 14.5; peaks unchanged.

#### Adopted blocks, steps A and B (2026-10-09): the generator pass re-packs, the Schur blocks lazy

Level 4 (Color on every path, CA being at the leaf only) peaked at 3.06-3.20
GB with 1.76-1.87 GB of "adopted blocks": the raw copies of this level's
packed near and Schur blocks.  Three things made them, all named in the
piece-4 notes: begin() materialized every adopted Schur block eagerly
(0.46-0.48 GB); the remote-generator pass materialized the union of what
its generators read in one batch (5.1-5.4 GB decompressed a level in 7
passes) and left the updated blocks raw until some later reader happened to
retire them; and a block a wave materialized but did not read was never
retired.

**Step A.**  The generator pass re-packs what it put back, right after its
GEMMs, in kNearGroupBudget groups as the owner pass does ("adopted blocks
re-packed after the generator passes": 9468 blocks, 4.2-4.4 GB raw at level
4; 7717 and 18948 at the CA3 and Color leaf, whose generator passes read
packed near blocks too).  Level 4: 3.20 / 3.06 / 3.14 -> 2.80 / 2.67 / 2.74
GB, with the peak moment moved to a generator pass's own union (adopted
2.16-2.29 GB resident, the decoder's staging 0.27 the only transient); the
leaf too, 5.01 -> 4.88 Color and 4.14 -> 3.97 CA3; levels 3 and 2 by 0.3.
Not bit-identical: every re-pack is one more quantization of those blocks,
and logdet moved by 2.7e-7 at 192^3 with acc_forward 1.60e-2 -> 1.64e-2
Color and 1.54 -> 1.70 CA3 (and 1.01 -> 0.98 at 96^3 Color: noise, not
drift).

**A'.**  Those re-packs at a 4x finer step (pack_near_group's tol_scale
0.25): the owner pass packs the same blocks again at its own step, so the
finer bytes are transient.  Logdet back within 1e-8 to 3e-8 of (c)'s and
acc_forward within 2%:

| (c) async -> A+A'+B | level 4 peak | leaf peak | logdet | acc_forward | Factor |
|---|---|---|---|---|---|
| Color 192^3 | 3.20 -> **2.68** GB | 5.01 -> 4.88 | -89253748.104 -> -89253747.301 | 1.60e-2 -> 1.58e-2 | 14.1 -> 14.1 s |
| CA3 192^3 | 3.06 -> **2.55** | 4.14 -> 3.97 | -89253742.145 -> -89253739.943 | 1.54e-2 -> 1.57e-2 | 13.7 -> 13.7 |
| CA0 192^3 | 3.14 -> **2.61** | 3.01 -> 3.01 | -89253736.570 -> -89253734.923 | 1.62e-2 -> 1.64e-2 | 14.5 -> 14.6 |
| Color 96^3 | 1.56 -> 1.53 | - | -9928903.097 -> -9928902.514 | 1.01e-2 -> 1.03e-2 | 2.47 -> 2.51 |
| CA3 96^3 | 1.47 -> 1.47 | - | -9928899.446 -> -9928899.926 | 9.29e-3 -> 9.22e-3 | 3.00 -> 2.98 |
| CA0 96^3 | 1.12 -> 1.12 | - | -9928897.922 -> -9928898.059 | 9.73e-3 -> 9.67e-3 | 3.28 -> 3.29 |

**Step B.**  The adopted Schur blocks stay packed at the level's start, like
the near blocks.  What made that unsafe before - writers treating a null
pointer as "no block" and filling from the kernel - is covered now: the
wave's gather reads a packed Schur block through its window, the
generator pass and build_parent name it to the materialization that
precedes them, and every other null-pointer site asserts (assert_not_packed),
so a missed site throws instead of silently answering from the kernel.
After the box's elimination the skeleton block takes the slot and the stale
packed entry is never consulted again, as with a replaced near block (the
packer resets it when it next packs the slot).  No guard fired on the six
runs.  Level 4 -0.12 / -0.13 GB: the eager set's share at the peak moment,
a late pass; it was 0.46-0.48 at begin().  The 2.1 GB still resident there
is one generator pass's union, which is step C.

#### Adopted blocks, step C (2026-10-09): the generator pass in bounded chunks

The remote-generator pass now runs in chunks of generators whose adopted
read set fits kNearGroupBudget raw: each chunk's blocks are put back
(materialize_adopted), its generators' updates built into batches of their
own (the per-generator body is a lambda, process_gen), launched, and the
chunk's blocks packed again at the finer step before the next chunk's are
put back.  The generators' own temp2^T - their transposes into the staging
and the slot packing - runs once per pass in an image of its own, before
the chunks; `replaced` persists across chunks (a replaced block must still
not be accumulated into), the retired blocks go back per chunk.  A block
two chunks read is packed twice, at the finer step.

Gate: six runs, 0 errors, no guard fired.  Level 4 at 192^3: 2.68 / 2.55 /
2.61 -> **1.57 / 1.53 / 1.55 GB** (resident 0.30 / 0.22 / 0.29 with no
adopted blocks in it; the peak moment is now a wave's, transient 1.27-1.31);
level 3 1.59 -> 1.10, level 2 1.33 -> 1.07-1.18; the leaf unchanged (4.88 /
3.97 / 3.01).  At 96^3 nothing moves: level 4 holds no adopted set there.
Accuracy: the chunk-shared blocks' extra packs show as logdet 3e-8 to 8e-8
relative against step B (Color -89253747.301 -> -89253740.435, CA3
-89253739.943 -> -89253734.641, CA0 -89253734.923 -> -89253731.853) and
acc_forward 1.58e-2 -> 1.62e-2, 1.57 -> 1.59, 1.64 -> 1.64; against (c)
that is within 3% on every path and size.  Factor 14.1 -> 14.3 s Color,
13.7 -> 13.9 CA3, 14.6 -> 14.6 CA0: the 93-95 materialization batches at
level 4 where there were 75, and an image upload per chunk.

**Where the peaks stand (2026-10-09, against (b) two commits of host staging
ago and against two days ago):**

| 192^3, 8 ranks | leaf | level 4 | level 3 | level 2 | rank's peak |
|---|---|---|---|---|---|
| Color | 4.88 GB (step 5's fresh near blocks of a wave chunk, 2.95) | 1.57 | 1.10 | 1.07 | 4.88 (was 20.72 on 2026-10-07, 6.10 at (b)) |
| CA3 | 3.97 (kept factors 0.67, Schur/edges 0.59, X_RR_full 0.47, temp2^T 0.33) | 1.53 | 1.09 | 1.15 | 3.97 (6.57 at (b)) |
| CA0 | 3.01 (kept factors 0.91) | 1.55 | 1.10 | 1.18 | 3.01 (22.51 on 2026-10-08, 7.33 at (b)) |

The leaf is the rank's peak on every path now, and level 4 is no longer a
floor under it.  Next levers, by size: Color's 2.95 GB of fresh near blocks
in one wave chunk (kWaveChunks = 8 a wave regardless of bytes: chunk by
bytes); the CA paths' kept factors (tier 2 of the producer-side (a), 0.67-
0.91 GB); CA3's fill X_RR_full 0.47.

#### The generator passes' kernel-made targets (2026-10-09): the Color leaf's 2.95 GB

A correction first: the 2.95 GB of "level Schur/edges" at the Color leaf's
peak was not step 5's fresh near blocks - those are packed as they are
written (3.93 -> 0.78 GB over 128 groups a leaf, ~31 MB raw a chunk).  The
only other allocations of that kind are the targets the remote-generator
pass makes from the kernel when a pair has no stored block yet
(ensure_schur, ensure_edge, replace_edge), and step A re-packed only what
the pass had *materialized*, not what it made.  At the Color leaf the
passes sum into 34987 targets of 7.75 GB over the level, most of them
fresh in the early waves; CA3's leaf had 0.59 GB of them at its peak
moment; CA0's leaf, a replicated level with no generator pass, none.

**Built in two commits.**  The three lambdas record the blocks they
allocate and the chunk packs them with the adopted ones it put back (at
the finer step; "packed after the generator passes: 33115 adopted blocks
(5.48 GB raw), 9344 targets made from the kernel (2.95 GB raw)" at the
Color leaf, 14379 of 4.40 GB at CA3's).  That alone took CA3's leaf 3.97
-> 3.39 GB and the 96^3 level 4 peaks 1.47 -> 1.17 (CA3) and 1.53 -> 1.09
(Color), but Color's leaf only 4.88 -> 3.83: its peak moment still held
the 2.95 GB, now with a 0.39 GB transient - the moment inside one chunk
after the targets are allocated and before they are packed.  The chunk
budget counted only the adopted blocks a chunk reads, and the leaf's first
pass reads almost none, so one chunk took every generator.  The second
commit counts what a generator would make from the kernel (a target with
no stored block, and every replace_edge) in the budget as well.

| (c) async -> now | leaf peak | level 4 | logdet | acc_forward | Factor |
|---|---|---|---|---|---|
| Color 192^3 | 5.01 -> 4.88 -> 3.83 -> **2.87** GB | 3.20 -> 1.58 | -89253748.104 -> -89253743.343 (5e-8) | 1.60e-2 -> 1.64e-2 | 14.1 -> 14.3 s |
| CA3 192^3 | 4.14 -> 3.97 -> **3.39** | 3.06 -> 1.51 | -89253742.145 -> -89253737.442 (5e-8) | 1.54e-2 -> 1.57e-2 | 13.7 -> 14.1 |
| CA0 192^3 | 3.01 -> **3.01** | 3.14 -> 1.55 | -89253736.570 -> -89253729.441 (8e-8) | 1.62e-2 -> 1.64e-2 | 14.5 -> 14.5 |
| Color 96^3 | level 4 1.56 -> **1.04** | - | -9928903.097 -> -9928902.363 | 1.01e-2 -> 0.99e-2 | 2.47 -> 2.63 |
| CA3 96^3 | level 4 1.47 -> **1.17** | - | -9928899.446 -> -9928901.097 | 9.29e-3 -> 1.01e-2 | 3.00 -> 3.22 |
| CA0 96^3 | level 4 1.12 -> 1.12 | - | -9928897.922 -> -9928898.027 | 9.73e-3 -> 9.72e-3 | 3.28 -> 3.32 |

Six runs, 0 errors, no guard fired.  The Color leaf's resident at its new
peak moment, a late wave's owner pass: X_RR_full 0.30, temp2^T 0.12, kept
factors 0.57, owner source slots 0.15, under a 1.67 GB transient.

**Where the peaks stand (2026-10-09, end of day).**  Against the uncompressed
Color run of 2026-10-06 (33.21 GB) and two days ago (20.72 Color, 22.51 CA0):

| 192^3, 8 ranks | leaf (the rank's peak) | level 4 | level 3 | level 2 | Factor | solve |
|---|---|---|---|---|---|---|
| Color | **2.87** GB: 1.19 resident (kept 0.57, X_RR_full 0.30) + 1.67 transient | 1.58 | 1.09 | 1.08 | 14.3 s | 0.48 |
| CA3 | **3.39**: 1.69 resident (kept 0.67, Schur 0.59 gone, X_RR_full 0.47, temp2^T 0.33) + 1.70 | 1.51 | 1.09 | 1.06 | 14.1 | 0.49 |
| CA0 | **3.01**: 1.34 resident (kept 0.91) + 1.67 | 1.55 | 1.09 | 1.09 | 14.5 | 0.51 |

Every path's peak is its leaf, and the transient half (the sketch's and
the owner pass's bounded budgets, 1.67-1.70) is now more than half of it.
Next levers: the kept factors at the leaf (0.57-0.91: tier 2 of the
producer-side (a)); CA3's fill X_RR_full 0.47; and the transient's two
budgets, the owner window (kOwnerWindowBytes 1 GB, 0.54 used) and the
sketch's fill slots, which are the first items that are not stores at all.

#### Memory pressure (2026-10-10): pack and stage only when the card is short of room

**The question.**  Against bitcomp off (7.54 s uncompressed, 7.49 s with
the fill sources compressed and everything else on the card), the packed
and staged factorization took 14.3 s at 192^3 on 8 ranks, 1.9x.  The
peak was 2.9-3.4 GB on a 40 GB card: the packing and the staging were
paying their time whether or not the room was needed.

**What the time was** (Color leaf, from the breakdown and the owner-pass
timers added for this): the owner pass 0.68 -> 2.25 s (GEMMs 0.7, pack
0.76, materialize 0.51, window 0.17, stage 0.13, 1396 batches instead of
120); the wave's production chunks 0.26 -> 1.17 s (the "plan boxes
[buffers]" number is a timer artifact: its start is set once before the
chunk loop and reset after each chunk's launches, so it measures the
previous chunk's tail - the compression of its fill sources, kept factors
and X_NR^T with their syncs); the transition 0.53 -> 0.87; the kept
staging 0.45; the arenas' pinning 0.5; level 4 +1.7.  The heap's
transient scan was suspected and acquitted by a counter (1.7M range
steps over 18k takes, milliseconds).  A 60 MB host image of the wave's
interpolation matrices, allocated per chunk and never read on the
device-sketch path, was removed; it was not the cost either.

**Budgets x4 and chunks 8 -> 4 first** (the "window budget" idea): 14.3 ->
13.0 s Color, 14.5 -> 13.1 CA0, peak 5.8-5.9 GB.  Most of it was the
chunking.

**Then the policy.**  A block may be raw or packed at any moment and every
reader handles both, so each pack site asks whether the heap is above a
hold limit now - the capacity less the transient budgets less 10% - and
packs only then: an owner group's targets, step 5's fresh blocks, a
generator chunk's blocks; the fill store's device copies are evicted (and
the host copies made, on demand) only then; the wave's production runs in
chunks only when its outputs would not fit beside what the heap holds; the
transition packs the parents only then; the solve factors are staged at
the next level's start only once a level has packed or chunked (or the
heap is above the limit).  The budgets are x2 of the originals (window 2
GB, fill window 1 GB, group 512 MB, scratch 512 MB, codec stagings 1 GB /
128 MB, handover 512 MB), the chunks 4 a wave, and the transient reserve
behind the hold limit is their sum, 6.5 GB.  Scaling the budgets x4 "with
room" was tried and removed: a wave's groups are formed while the heap is
still under the limit, and an 8 GB window taken then left a 13.8 GB heap
without a gather window for the next chunk.  The gather window halves
until it fits, and a wave that gets none puts its packed blocks back one
by one (9220 at the CA3 leaf in a 11 GB heap) instead of tripping the
guards.  BPACK_GPU_HEAP_FRACTION (an existing variable) shrinks the heap
and forces the regime for a test.

**Gate, default heap (six runs, 0 errors)**, against the packed state
(dd75ba4) and the uncompressed runs:

| default heap | Factor packed -> now | uncompressed | acc_forward packed -> now (uncompressed) | leaf peak | level 4 |
|---|---|---|---|---|---|
| Color 192^3 | 14.3 -> **8.51 s** | 7.54 | 1.64e-2 -> 8.0e-3 (4.2e-3) | 18.3 GB | 10.1 |
| CA3 192^3 | 14.1 -> **8.65** | - | 1.57e-2 -> 9.2e-3 | 18.9 | 7.7 |
| CA0 192^3 | 14.5 -> **8.33** | - | 1.64e-2 -> 8.2e-3 | 19.2 | 10.8 |
| Color 96^3 | 2.63 -> **1.43** | 1.34 | 9.9e-3 -> 5.4e-3 (3.5e-3) | 2.16 | - |
| CA3 96^3 | 3.22 -> **1.73** | - | 1.0e-2 -> 5.5e-3 | 2.11 | - |
| CA0 96^3 | 3.32 -> **1.68** | 1.39 | 9.7e-3 -> 5.2e-3 (3.4e-3) | 2.93 | - |

Nothing packed, no wave chunked, no parent packed, no factor staged (CA3's
leaf transition packed 3 of 9 chunks, briefly over the limit, and its
factors were then staged).  Logdet Color 192^3 -89253593.153 against the
uncompressed -89253588.194 (5.6e-8): only the fill sources and the kept
factors are quantized now.  The solve: 0.38-0.41 s (set-up 0.10).

**Gate, forced pressure** (0 errors unless said): 96^3 at 15% of the heap
(hold limit 2.5-3.0 GB), heap-checked, all three paths: parents packed,
factors staged, 1.6-1.8 s, logdets in the packed band.  192^3 Color at 35%
(13.8 GB, limit 5.4): 20 of 32 owner groups packed, 13 of 16 waves
chunked, parents packed in 15 of 15 chunks, factors staged, near store 3.25
GB on the host, peak 9.6 GB, 11.7 s, acc_forward 1.53e-2.  CA3 at 35%
(11 GB, limit 3.6): 82 of 110 groups packed, 1667 fill sources copied on
demand, one chunk without a window, 11.4 GB peak, 0 errors.  CA0 at 55% and
45% (limits 12.5 and 8.9 GB): 4 and 20 groups packed, parents packed in 8
of 10 and 11 of 11 chunks, peaks 14.5 and 11.0 GB, 10.0 and 10.8 s, 0
errors.

**Known issue, pre-existing path.**  CA0 at 35% fails at level 4's first
wave with "near block shape (hi): block 279 x 360, box 512 n 279, neighbor
79 c 361 (windowed), packed entry live": an adopted block one column short
of a *remote* neighbor's current size, the case begin()'s shape check
skips.  It appears only when the small heap also makes decide_solve_keep
decline the leaf's factors (the raw bound / 4 against 80% of the heap), so
the X_NR goes to the host and the remote skeleton sizes travel another
path; with the factors kept (45%, 55%) the same mixed transitions pass.
Not pursued: the keep-declined path is not where the 8x regime should
run, and the keep rule itself is the next 8x item - at 8x the boxes a
rank the raw bound / 4 exceeds 80% of a 33 GB heap and the leaf would
decline on the full card (its X_NR then going to the host uncompressed).
The rule should compare the compressed size (the measured 8.3-8.6x) and
the kept store should leave the device as it is produced (tier 2).

**Where it stands.**  At 192^3 on 8 ranks the slowdown against bitcomp off
is 1.13x (8.51 vs 7.54 s) with a 18-19 GB peak; the memory regime engages
by itself when the card is short, at the costs measured before (the near
store 5x, staging, 11-14 s at this size), and in between the two mix per
group, per wave and per transition chunk.

#### Toward 8x boxes a rank (2026-10-10): the keep estimate, the hold limit, chunks by bytes

Three device-side edits (e624528) ahead of the 384^3 probe: the keep
decision compares the *compressed* size of a level's factors (the first
level at an assumed 8x, every later one at the ratio the previous level
measured, allreduced as the least compression of any rank) instead of the
raw bound / 4, since at 8x the boxes the raw bound exceeds 80% of the heap
and the leaf would decline on the full card; the hold limit reserves the
kept store it expects (`kept_expected_`), so the packing starts early
enough to leave it room; and the wave's production chunks are counted by
bytes (kWaveChunkBytes 2 GB) rather than a fixed four, so a 8x wave gets
8x the chunks.  192^3 unchanged on the three paths (logdets, acc, peaks,
Factor 8.3-8.8 s).  The 384^3 probe (8 ranks, 32768 leaf boxes a rank)
then died in the leaf's interior waves at ~60 GB of host RSS a rank, not
at level 4's start as expected: the host, not the card, was the limit.

#### Host memory (2026-10-10): the copies the restore never read, and the pinned pool

**What the host held.**  A per-field count of every host box at each
level's end (report_host_boxes, BPACK_TRACE=wave), 192^3 Color, rank 0:

| at the leaf's end | boxes | GB | of which |
|---|---|---|---|
| leaf local boxes | 4096 | 1.27 | X_RR_full 0.99, T 0.24, points and lists 0.04 |
| leaf generator boxes (other ranks') | 817 | 5.15 | X_NR (temp2) 4.94, X_RR_full 0.21 |
| level 4 (at its end) | 512 + 217 | 0.40 + 3.47 | X_RR_full 0.32, temp2 3.32 |

The pinned arenas (near store, fill store) were 4.3-5.9 GB at 35% of the
heap; RSS 8.4 GB at the leaf's end on the full card.  So the host side of
a level was 6.4 GB of box copies and 5 GB of stores, and at 384^3 (the
generator set grows with the surface, x4; the local boxes x8) the probe's
RSS reached 44 GB at the leaf's end plus 34 GB of arenas that never
shrank - the arenas kept every level's peak, and the kept factors of the
leaf (pageable, 18 GB at 384^3) came on top at the next level's start.

**What read the copies.**  The remote generators' temp2 and the local
boxes' X_RR_full are the host side of the transition's restore: a fill
source freed during the waves (early_free_sources_, multi-rank levels) is
put back from them.  But with the fill sources stored compressed the
transition puts a source back from its packed temp2^T
(collect_fill_source, decompression at 365 GB/s) and rebuilds X_RR_full
from the kept LU (xrr_rebuilt); the host copies were a fallback that no
run has read: 0 reclaims and "restore 0.00 GB" in every run, pressure runs
included.  T stays (the matvec reads it on the host).

**What changed.**
1. `packed_restore_` (= device sketch and H2_compress_factors): the
   exchange no longer copies a received generator's temp2 to the host
   (its X_RR_full, r x r, still comes), and the wave download skips a
   kept box's X_RR_full (host_skips_xrr; not on a level whose generators
   go through the host, where emit_generators reads it).  A reclaim under
   it frees a source's X_RR_full only when the LU or a host copy can
   rebuild it, frees the raw temp2^T only when the packed form exists,
   and makes no host copy (it used to materialize the whole level's X_NR
   on the host first - at 8x that alone would be the raw X_NR again).
2. The transition's packed put-back, when the source's X_RR_full is on
   neither side of the card (a remote generator whose block was reclaimed;
   a box whose keep failed), takes it from the host copy or rebuilds it on
   the host from the rescued LU (host_xrr_full_from_lu: P L U with the
   pivots undone last to first, as launch_row_swaps does) through the
   chunk's restore block as a source with no temp2 rows (`xrr_only`).
   drop_kept nulls the states' LU pointers (they pointed into the freed
   blocks: a latent read of freed memory after a mid-level keep failure)
   and first gives every box that dropped X_RR_full for its LU the block
   back (restore_dropped_xrr, host rebuild, the rare path only).
3. One PinnedChunkPool for the process: the arenas (fill store, the two
   near stores) and now the staged kept factors take 256 MB chunks from
   it and give them back at their reset; the driver trims the idle set at
   each level's start, after the kept factors took what they could, to
   half of what is idle (floor 1 GB).  So the leaf's near-store chunks
   serve the leaf's kept factors at the next level's start (no pinning,
   and the kept staging is a pinned copy now: 26 GB/s instead of 13),
   and the host holds the live bytes plus one cache.  Unpinning runs at
   1.6 GB in 0.03 s.

**Gate (192^3, 8 ranks, the corrected build).**

| default heap | Factor before -> now | logdet | acc_forward | host RSS at the leaf's end | leaf download |
|---|---|---|---|---|---|
| Color | 8.55 -> **7.88 s** (bitcomp off 7.54) | -89253593.152662069, identical | 8.0039e-3, identical | 8.40 -> **2.16 GB** | 6.38 -> 0.45 GB |
| CA0 | 8.33 -> 8.96 | -89253606.998565048, as e624528 left it | 9.07e-3, as e624528 | 2.33 GB | |
| CA3 | 8.65 -> 8.59 | -89253625.765752852, as e624528 left it | 8.93e-3, as e624528 | 4.77 GB | |
| Color at 35% | 12.1 -> **10.9** | -89253779.87903136, identical | 1.5442e-2, identical | 8.47 GB | |
| CA3 at 35% | 15.6 -> 11.4 | -89253762.743 (packed band) | | 10.8 GB | |

The leaf's download is what the time was: 6.38 GB a rank of host copies
through the copier, 0.45 GB now, and Color is within 4.5% of the
uncompressed run (7.88 against 7.54).  The host RSS at every level's end
is the tree and T: 2.16 GB.  The CA logdets are the ones e624528 produced
(its lower hold limit - the kept store reserved - packs 2 of CA0's and 4
of CA3's transition chunks where the reference packed 0 and 3), unchanged
by this step; CA0 stages its leaf factors now for the same reason.

**Two stale-object lessons** that cost an afternoon: an object whose
compile had started before a header edit and finished after it carries
the old source with a new timestamp, so make keeps it; and two
`*BPACK_gpu_wrapper` objects do not list the eliminator's header in their
dependencies at all.  Both now show as a run that "rebuilt X_RR_full from
the kept LU" at every level (that path announces itself).

**The 384^3 probe after this.**  The leaf's 16 waves run to the end (57
s of elimination; kept factors 18.6 GB a rank: X_NR 145.9 -> 9.7 GB at
15.1x, T and X_SR 3.9 -> 1.1, the LU and pivots 7.8 uncompressed) with
host RSS at most 44.3 GB a rank and 37-38 GB of pinned chunks at the last
wave (the near store and the fill store, 8x their 192^3 size).  It then
died in the leaf's transition: the heap exhausted at 33.1 of 33.4 GB with
a 236 KB largest free range.  The chunk budget (2 GB, or a quarter of
the free room) counted the G/P buffers, the host-path restores and the
parents' blocks but not the packed put-backs - and at this size every
source is packed and they are the level's raw temp2^T, 146 GB a rank, so
a chunk sized by its parents alone decompressed far more than the room
it had.

#### The leaf's transition at 8x (2026-10-10): the chunk budget, and a leak the chunks exposed

Counting the packed put-backs (packed_restore_bytes) in the chunk
budget was not enough: the 35% runs at 192^3 then died in the transition
too ("no room to put back fill source ... heap used 11922 of 13115 MB,
largest free 5 MB, chunk budget 24 MB"), and the per-chunk trace at
384^3 (BPACK_TRACE=wave) showed why.  The X_RR_full rebuilt from the
kept LU for a chunk's sources (3 r x r a source: the product and its L
and U) went into blocks freed only when the whole transition was done,
so over the leaf's chunks they accumulated - "restored" 0.14 GB after
chunk 1, 14.1 GB after chunk 100 - and the budget, a quarter of the free
room, collapsed to 50 MB and the chunks to one parent each.  At 192^3
the leaf's 10 chunks never showed it.  They are freed with the chunk now
(the chunk's P products are their only reader) and counted in the budget
(xrr_rebuild_bytes).

| 192^3, 8 ranks | Factor | logdet | acc_forward | leaf transition |
|---|---|---|---|---|
| Color | 7.94 s | -89253593.152662069, identical | 8.0039e-3, identical | 0.53 s, 38 chunks (was 10) |
| Color at 35% | 11.0 | -89253771.650 (packed band) | 1.59e-2 | 0.86 s, 38 chunks |
| CA3 at 35% | 11.6 | -89253759.096 (packed band) | 1.54e-2 | 0.96 s, 37 chunks |

**384^3 on 8 ranks (32768 leaf boxes a rank, 8x the 192^3 run) factors
to the end.**  Per rank: leaf elimination 57.1 s and transition 8.2 s
(334 chunks of ~12 parents: budget 2.15 GB = buffers 0.57 + packed
put-backs 1.33 + X_RR_full rebuilt 0.14 + parent blocks 0.24; the heap
at 17 GB with the kept factors 8.9 GB of it); level 5 12.7 + 2.0 s;
level 4 5.1 + 0.7; level 3 2.8 + 0.3; level 2 1.2; the root LU on the
device.  Factorization total 92.3 s, of which the data exchange 17.5 s.
Memory: the heap under its hold limit throughout (leaf limit 11.9 GB:
718 of 776 owner groups packed, 15 of 16 waves chunked, every transition
chunk's parents packed); host RSS high-water 51.5 GB a rank, pinned
chunks at most 38.1 GB (the near store and the fill store at the leaf;
30 GB of staged kept factors at level 3, 18.4 of them the leaf's); the
trims gave back 8.6, 4.3 and 1.6 GB at the starts of levels 5, 4 and 3
(0.41, 0.07, 0.05 s).  The kept factors: leaf 18.6 GB a rank (X_NR
145.9 -> 9.7 GB, 15.1x), level 5 6.1, level 4 1.7, level 3 0.5.  logdet
-792567616.26243639.  The run was then OOM-killed on the host in the
solve's set-up, after the factorization: the solve at 8x is the next
item (the kept factors alone are 27 GB a rank compressed, and a host
fall-back decompresses the X_NR).

The 38 chunks at 192^3 on the full card (10 before, +0.04 s) are the
rebuilt X_RR_full now counted against a 2 GB cap the free room would
allow to be 4x larger; the cap can follow the room if the chunk count
ever matters.

#### The solve at 8x (2026-10-10): streamed kept factors, and waves in segments

After the 384^3 factorization the run was OOM-killed on the host.  The
solve's set-up (prepare_device_solve) wants every Color level's kept
factors resident for the sweeps; at 384^3 that is 27 GB a rank (leaf
18.7, level 5 6.3, level 4 1.8, level 3 0.5) against a budget of 24.2
GB, so it declared the device solve off and materialized the factors on
the host - the X_NR uncompressed, 146 GB a rank.  Two changes, both
decided by room and collectively, no new options:

1. **Streamed levels.**  The staged levels that do not fit beside the
   rest are streamed, largest first: their box pointers are rebased at
   the set-up to one device home (the largest streamed level, 19.6 GB
   here), and at each sweep step the level's bytes come into it from
   their pinned chunks (stream_kept_in, 26.8 GB/s), ahead of the step's
   kernels on the same stream.  A solve streams each level twice (both
   sweeps), the multiply too.  The set-up streams a level once while it
   plans its tables (the planner reads the compressed headers).  The
   pinned chunks of a streamed level stay (they are its factors now);
   release_prepared puts the pointers back and frees the home.
2. **Waves in segments.**  The sweeps decompress a wave's X_NR into a
   scratch before its kernels, and the scratch was sized by the widest
   wave: 9 GB raw at the 384^3 leaf.  A wave is now cut into segments
   whose raw bytes fit one buffer of a scratch all levels share (two
   buffers of a segment; a segment at most 2 GB, an eighth of the room
   when that is less), and a segment is the unit the sweeps decompress
   and launch - the kernels take an order range, so a wave's segments
   run one after another with the wave's accumulation and messages
   unchanged after the last.  The scratch is counted in the solve's
   budget up front.

**384^3 on 8 ranks, end to end** (st5 run): Factor 91.8 s; solve
3.89 s (rank 0: forward 0.82, backward 0.92; the streaming
54 GB in 2.02 s at 26.8 GB/s, 4 levels twice); multiply 4.17 s; H2_CheckError acc_mvp 1.27e-3,
acc_forward ACC0.82, acc_backward ACC0.92.  Host RSS high-water
51.6 GB a rank.  192^3 unchanged on the three paths (Color
7.96 s, logdet identical, solve 0.371 s, multiply 0.442; CA0 8.23 s with
the reference logdet -89253593.818 and acc_forward 8.15e-3; CA3 8.57 s,
logdet -89253593.484 and acc_forward 7.97e-3, its leaf transition
packing 4 of 37 chunks now where it packed 3 of 9).  The 4.3 GB scratch
is reserved on every card now, which the 33 GB card does not notice.

**What the per-rank accounting showed.**  Ranks 1-7 had 8.8 GB of kept
factors resident at the set-up where rank 0 had none: the kept-factor
staging at a level's start is decided per rank (a level that packed or
chunked, or a card above the hold limit), and only rank 0's level 5 had
seen pressure, so levels 3-5 stayed on the cards of the other ranks.
Each rank streams what it has on the host and sizes its home by that,
so the decision is consistent; it also makes the budget tight (planned
4.4 + home 19.6 against 24.2 GB on the tightest rank).  At 64x a level
would not fit the home either: the next step there is streaming by
segment (a wave's share of the kept store), not by level.

#### 768^3 on 64 ranks, the first run (2026-10-11): what each path hit

16 nodes, 64 ranks, 32768 leaf boxes a rank (8x the 192^3 run; CA at the
leaf, level 7).  Color factored in 386 s with the data exchange 176 s of
it; levels 3 and 2 took 71 and 82 s: both ran "host sketches (boxes over
2047 points)" - a sketch list entry packed 20 row bits and 11 destination
bits, so a box of over 2047 points (level-4 skeletons average 179, eight
of them a level-3 box) sent the whole level to host sketches, and with
no device sketch the generators went through the host (level 3: 62 s of
exchange, 26.6 s of size handshakes and 30.9 s assisting; level 2: 53 s;
the ID on the host 6.4 and 25.5 s).  Its solve set-up then wanted the
leaf's 19.7 GB home beside 9.9 GB the two unkept reduced levels and the
scratch need, against 23.3 GB on the tightest rank (which still held 9.9
GB of kept factors: the staging decision is per rank), declared the
device solve off and was OOM-killed materializing the factors on the
host.  CA3 exhausted the 28.4 GB heap in the leaf's owner pass (green
colour, requests of 2.5-7.4 GB against 2.4-6.9 GB largest free).  CA0's
leaf ran (74.3 s against Color's 66.4; 36594 boxes a rank on a corner
brick, +11.7%) and its transition too; then one rank found no room for
the kept X_NR hand-over and another was OOM-killed at 52.3 GB.

#### Three fixes for 768^3 (2026-10-11): the entry width, the ring, and a diagnosis retracted

**Sketch entries: 18 row bits, 13 destination bits** (8191 points a box,
262143 rows; a box's rows are its neighbours' points, at most 27 times
its own, so the bounds agree).  One constant in bpack_gpu_kernels.cuh;
the kernels and the host bounds follow the masks.  192^3 bit-identical.

**The solve streams by wave through a byte ring.**  Streaming a level
through a level-sized home (19.7 GB at the 384^3 leaf) was the budget
blocker.  Now a staged level is planned over the solve's waves (the
Color binning, or the CA groups x parities): the blocks each wave reads,
the blocks numbered by first use and placed whole in a ring grown until
no two blocks of a wave or of two consecutive waves overlap (the prefetch
is one wave deep); a box's pointers are rebased to its blocks' places
once at the set-up (the staged form kept for the host readers); a wave's
blocks are copied in on the decompression stream ahead of its step, the
next wave's under this one's kernels, and the ring remembers what it
holds.  The diagonal solves of a streamed level run with each wave (its
LU is only in the ring then; legal, since a box's redundant part is
updated only by earlier waves - the first version ran every level's
diagonal per wave and cost 10% of the 192^3 solve on the small kernels,
so the other levels keep the one-shot diagonal).  The set-up streams a
wave in while it plans the wave's segments: bitcomp's batch plan reads
the compressed headers.  A first version sized the ring by block count
times the largest block and came to 25.7 GB at the 384^3 leaf (a wave
reads ~26 blocks: its keep groups and the X_NR hand-over groups); by
bytes it is 10.2 GB.

| | before | now |
|---|---|---|
| 384^3 / 8 ranks: streamed levels, home | 3,4,5,6 through 19.6 GB | the leaf only, through a 10.2 GB ring (levels 3-5 resident) |
| solve / multiply | 3.89 / 4.17 s (2.0 s of copies in each) | **2.55 / 2.77 s** (16.5 GB streamed a solve, hidden) |
| acc_forward / acc_backward | 1.58e-2 / 1.0e-15 | identical |
| 192^3 three paths | | bit-identical; solves 0.37-0.41 s as before |

**The CA leaf's "per-wave tail": measured and retracted.**  The
768^3 numbers read as a fixed cost per wave (the chunk tails 21.4 s over
CA0's 32 waves against 12.0 s over Color's 16), and batching the small
boundary waves' compression tails was proposed.  A per-wave host-time
line (BPACK_TRACE=wave) at the 384^3 CA0 leaf on 8 ranks says otherwise:

| rank 0, 384^3 CA0 leaf | boxes a wave | host time a wave | per box |
|---|---|---|---|
| blue/orange, 8 waves | 34-61 | 0.03-0.04 s | ~0.9 ms |
| purple, 8 waves | 195-286 | 0.16-0.22 s | 0.8 ms |
| green, 8 waves | 588-980 | 1.09-1.35 s | 1.4 ms |
| interior, 8 waves | 3375 | 8.1 then 4.1-5.8 s | 1.5 ms |

The small waves are the cheaper ones per box.  The leaf: 61.0 s for
36594 boxes (1.67 ms a box) against Color's 56.7 s for 32768 (1.73 ms,
of which the exchange 11.7 s; 1.37 ms of compute).  So the CA leaf's
extra cost is the ghost layer itself - 11.7% more boxes on a corner
brick, ~20% on an interior one, with 23% more near-block targets
(392762 against 320341) and 26% more device time - not the wave count,
and no tail batching is built.  CA0's leaf beats Color's only where
Color's exchange wait is over ~45% of its compute; at 64 ranks it is 41%
(19.4 of 47 s) and the level time follows the interior bricks.

**CA0 at 8x fails in the normal path.**  The 384^3 CA0 run (8 ranks,
its leaf kept and staged) died at level 5's start with "near block shape
(hi): block 359 x 350, box 4256 n 359, neighbor 749 c 349 (windowed),
packed entry live" on every rank but 0: an adopted block one column
short of a remote neighbour's current size.  The same message was
recorded at 192^3 only with the keep declined at 35% of the heap, and
read as a keep-declined path; at 8x the packing is always on, and this
is where it shows.  The replicated elimination of a ghost box has to pick
the same skeleton on every rank that holds it; the lossy quantization of
its inputs (a whole delta where the raw values straddle a bin) makes the
ID's rank decision flip for one box in thousands, and the parent block
built here disagrees with the owner's by a column.  CA with packing needs
the ghosts' skeletons taken from their owners (as the assisting boxes'
already are), which is a CA-path change, not a memory one.

#### 768^3 on 64 ranks, the rerun (2026-10-11): one level fixed, the next one blocks

With the 13-bit destinations, Color's levels 3 and 2 take the device
path: level 3 eliminates in 4.6 s where it took 71 s, level 6 in 18 s
(30 before; its exchange 5.3 s).  The run then dies in the 64 -> 8
process reduction from level 3 into level 2, which now runs on the
device as well (the first run's went through the host because those
levels had no device sketch): the eight receivers (ranks 7, 14, ...,
56) allocate the merged level-2 blocks raw, one resident block a piece,
on cards already full - rank 7's own level-3 transition found 31.1 of
31.9 GB in use, and ranks 28 and 56 asked for single blocks of 10.4-10.9
GB, which matches no level-2 dimension (a parent of eight 313-point
skeletons gives 50 MB blocks).  Two things to settle before the next
run: what those blocks are, and whether the device sketch is right for
boxes over 2047 points - at level 3 it kept 313 skeleton points on
average where the host sketch kept 196 on the same level-4 input (the
destination-range split, parts of ~557 against the shared memory, is
exercised from 557 points on at every size, so the difference is not the
split as such).  CA3 and CA0 fail as in the first run (the owner pass's
2.5-7.4 GB requests; the kept X_NR hand-over on 13 ranks and the host
OOM on the interior brick of rank 35).  A 192^3 test with 2744-point
leaves was mis-sized (an octree's box sizes step by 8x: 216, 1728,
13824 points; a box between 2048 and 8191 points exists only at the
upper levels of a large problem) and timed out on host sketches.

#### 768^3 on 64 ranks runs end to end (2026-10-11): the staging decision made collective

The diagnostics of the rerun said three things.  Every rank's eight
level-3 boxes have 600-2130 points and the device sketch keeps very
large skeletons on the late boxes (k up to 1475 of 2001 points; the
level's average 311 against the host sketch's 196 in the first run),
monotone in the elimination order - so an accuracy difference, not a
corrupt sketch; the rule is the same on both sides (d = the box's
points, 4 nonzeros a column).  The 10.7 GB allocations on the ranks
that failed were not the reduction's: they came before it, in those
ranks' own level-3 transition - the G and P buffers of one parent's fill
GEMMs (r x k a pair, with r and k both ~1000-1500 at these skeletons),
which the chunk loop takes whole for its first parent.  And those ranks'
cards held ~27 GB of kept factors that rank 0 had staged: the staging
at a level's start was each rank's own decision (its pressure or its
card above the hold limit), and the ranks that saw none kept levels 6-4
resident.

**The decision is the level's now** (an allreduce of the flag over the
level's ranks in stage_solve_factors), the reduction tests the receivers'
room (what its senders hold against the card, with a margin) and takes
the host reduction otherwise, and the per-rank card and box lines are
under BPACK_TRACE=wave.  With that, 768^3 Color on 64 ranks runs to the
end:

| 768^3, 64 ranks, Color | first run (host sketches at 3, 2) | now |
|---|---|---|
| factorization | 386 s (exchange 176) | **119.7 s** (exchange 37.6) |
| levels 7 / 6 / 5 / 4 / 3 / 2 | 66.7 / 30.0 / 8.0 / 5.9 / 70.7 / 82.5 s | 65.6 / 16.2 / 7.5 / 5.5 / 2.9 / 3.8 s |
| 64 -> 8 reduction into level 2 | through the host | on the device: receivers get 3.0-3.6 GB with 0.2-0.7 GB in use |
| cards at level 3's start, heaviest ranks | 33 GB in use | 3.7 GB (parent blocks) |
| solve / multiply | OOM-killed at the set-up | 3.39 / 3.81 s, the leaf through a 10.2 GB ring |
| accuracy | - | acc_mvp 1.25e-3, acc_forward 2.14e-2, acc_backward 1.3e-15 |
| host RSS high-water | | 51.6 GB a rank (54.3 at level 1) |

The logdet of this run is reported with a negative sign (-1,
log|det| 6968677636.29) where the first run's, with host sketches at
levels 3 and 2, had +1 and 6968679026.42: the two approximations differ
by 2e-7 in log|det| and in the inertia's parity, which makes the
level-3 skeleton question one of correctness, not cost.

Open on this path: the level-3 skeleton inflation (a 2.5x larger
level-2 block set and 1.6x more solve work at those levels than the host
sketch would give; the cause - the lossy fill rows, or the device
sketch's sampling of them - is not established), and the transition's
per-parent buffer floor (10.7 GB for one parent here; fits once the
cards are empty, but it is the next limit up).  CA3 and CA0 at 768^3
fail as before, in their leaves.

#### The rank growth and the sign of the determinant (2026-10-11): the fill store's quantization

The three 768^3 Color runs disagreed above level 5 and nowhere else:

| 768^3, 64 ranks | level 7 | 6 | 5 | 4 | 3 | 2 | logdet |
|---|---|---|---|---|---|---|---|
| first run (host sketches at 3 and 2) | 44.27 | 72.02 | 123.86 | 178.72 | 196.01 | 83.75 | +1, -6968679026.42 |
| rerun (level 3 on the device) | 44.27 | 72.02 | 123.86 | 178.34 | 313.13 | - | - |
| end to end (A_l768c) | 44.27 | 72.02 | 123.96 | 166.95 | 219.96 | 219.81 | -1, -6968677636.29 |
| end to end again (A_dg768o, same library) | 44.27 | 72.02 | 123.76 | 167.00 | 220.40 | 221.06 | -1, -6968677635.06 |

The skeleton count was run-to-run noise that grew by a decade a level (0.1%
at level 5, 7% at 4, 40% at 3), and within a level it grew with the colour:
at level 3 the first colour kept 2.4% of its points and the last 39%
(per-rank box lines, `BPACK_TRACE=wave`).  The only lossy input the top
levels have is the fill store - the memory-pressure packs (owner groups,
parents, near blocks) all read 0 at levels 5 and above - and it is packed
at 2 * 0.2 * tol_comp * max|temp2^T of the box|.  A later box's sketch is
its kernel rows less W_E P_E over its eliminated neighbours E, with both
factors built from E's packed temp2^T, so the quantization noise lands in
Y at a fraction of the fill, and the fill is of the order of the kernel
rows for the late colours; the ID at tol keeps whatever noise stands above
tol |R_00|.  Hence the growth with the colour, the sensitivity to which
blocks happened to be packed, and a determinant whose approximate
factorization is no longer positive definite.

**Measured** (lib_dg: the fill store and the kept X_NR 100x finer at
levels <= 5, nothing else changed; two runs):

| 768^3, 64 ranks | 5 | 4 | 3 | 2 | logdet | acc_mvp | acc_forward | factor | solve |
|---|---|---|---|---|---|---|---|---|---|
| committed (A_dg768o) | 123.76 | 167.00 | 220.40 | 221.06 | -1, -6968677635.06 | 1.21e-3 | 1.98e-2 | 140.2 s | 3.53 s |
| finer at <= 5, run a | 101.95 | 109.00 | 111.09 | 76.19 | +1, -6968677260.85 | 5.60e-4 | 1.40e-2 | 134.2 s | 2.78 s |
| finer at <= 5, run b | 101.91 | 109.08 | 110.79 | 76.47 | +1, -6968677256.72 | 5.76e-4 | 1.40e-2 | | 2.67 s |

Levels 4, 3 and 2 took 3.69 / 1.27 / 0.94 s against 5.78 / 3.03 / 3.87.  The
per-level logdet lines (new under `BPACK_TRACE=wave`) show the sign is the
parity of 180k sign flips (negative pivots and row swaps) over the levels;
the -1 was that parity, not one level.

**The rule, against the uncompressed run** (192^3 on 8 ranks; the lossless
Color figures are the 2026-10-06 `A_v192off` log; `BPACK_DG` was a
run-time knob of the sweep library only):

| fill store step | level 5 | 4 | 3 | 2 | logdet | acc_mvp | acc_forward |
|---|---|---|---|---|---|---|---|
| uncompressed (2026-10-06) | 42.533 | 61.540 | 76.025 | 65.313 | -89253588.194 | 8.1e-4 | 4.2e-3 |
| a box, 0.2 (committed) | 42.550 | 62.409 | 87.283 | 114.641 | -89253593.153 | 1.89e-3 | 8.00e-3 |
| a slot, 0.2 | 42.531 | 61.755 | 79.266 | 86.000 | -89253589.701 | 2.30e-3 | 7.97e-3 |
| a box, 0.002 | 42.517 | 61.511 | 75.842 | 64.594 | -89253588.235 | 1.96e-3 | 7.46e-3 |
| a slot, 0.002 | 42.521 | 61.542 | 75.834 | 65.563 | -89253588.208 | 1.90e-3 | 7.48e-3 |
| a box, 0.0002 | 42.531 | 61.502 | 75.873 | 65.281 | -89253588.202 | 1.80e-3 | 7.25e-3 |

The fraction matters, the grouping much less: 0.002 brings the skeletons
back to the uncompressed ones within 1% at level 2 and 0.2% below it, and
the logdet to 5e-10 (0.2: 5e-8).  The remaining acc gap to the uncompressed
run (acc_mvp 1.9e-3 against 8.1e-4) is the kept X_NR's own step, which the
sweep did not touch.

The leaf's fill store does not need it: at 384^3 on 64 ranks, the leaf
packed finer or not gave the leaf 43.449 / 43.475 and level 5 65.603 /
65.600, and the leaf's store is the one that costs memory (packed 100x
finer at every level, the 768^3 run was OOM-killed on the host at the
leaf).  Nor can every level above it have it at 8x: with the rule "fine
above the leaf" (lib_fix) 768^3 was OOM-killed on the host in level 6,
whose fill store is packed under pressure and goes to the host - its
pinned chunks read 47-50 GB a rank at waves 8-14 against 42-43, the RSS
63 GB against 60, four ranks a 256 GB node.  The same rule at 384^3 on 8
ranks (the same level geometry, less of it interior) fitted: level 5's
store 60.4 -> 10.4 GB (5.8x) against 3.9 (15.3x), host +1 GB, the level
+0.85 s; the levels above it 43.83 67.86 98.86 92.21 71.69 against 43.83
69.84 111.57 140.54 164.03, factor 91.1 s against 93.9, solve 1.70 s
against 2.60, acc_forward 1.58e-2 both.  192^3 under it, all three paths:
Color 42.55 61.55 75.86 65.00 (logdet -89253588.46, acc_forward 7.8e-3),
CA0 42.28 59.87 75.93 65.16, CA3 42.21 59.22 75.53 65.20, every sign +1.

So the step is decided a level at a time (`decide_fill_precision`, at
the level's first pack, when the keep rule's estimate is in the hold
limit): `kFillToleranceFraction = 0.002` when the store's bound at the
fine step (27 n^2 an n-point box, over `kFillFineRatio` = 6) fits in
half the room under the hold limit, else the kept X_NR's 0.2; the leaf
never; the kept X_NR unchanged.  At 768^3 on 64 ranks that is coarse at
levels 7 and 6 (bounds 55 and 18 GB at the fine step against rooms of
1.3 and 14 GB) and fine from level 5 up (5.9 against 12); at 1x
everything above the leaf is fine.  Under it:

| | level 7 | 6 | 5 | 4 | 3 | 2 | logdet | acc_mvp | acc_forward | solve |
|---|---|---|---|---|---|---|---|---|---|---|
| 768^3, 64 ranks, committed | 44.27 | 72.02 | 123.76 | 167.00 | 220.40 | 221.06 | -1, -6968677635.06 | 1.21e-3 | 1.98e-2 | 3.53 s |
| 768^3, 64 ranks, the rule | 44.27 | 72.02 | 123.57 | 128.27 | 107.90 | 76.59 | +1, -6968677341.91 | 1.04e-3 | 1.71e-2 | 2.98 s |
| 384^3, 8 ranks, committed | 43.83 | 69.84 | 111.57 | 140.54 | 164.03 | | +1, -792567616.16 | 1.27e-3 | 1.58e-2 | 2.60 s |
| 384^3, 8 ranks, the rule | 43.83 | 69.84 | 109.14 | 91.88 | 71.86 | | +1, -792567587.49 | 1.14e-3 | 1.50e-2 | 2.42 s |

(The 384^3 rows list levels 6 to 2.)  192^3 on 8 ranks under the rule is
the "fine above the leaf" result above, all three paths.  Level 5 of
768^3 stays at 123.6 because level 6's store is coarse, as level 5 of
384^3 at 8x (69.84 against 65.60 at 1x); the level above a coarse store
pays a few percent, the ones above that converge.

**The lossless 384^3 reference** (64 ranks, `H2_compress_factors 0`, on
80 GB cards - the leaf's raw kept factors do not fit a 40 GB card even
at 1x): 43.449 65.607 88.071 97.268 71.906, logdet +1 -792566483.68,
acc_mvp 5.2e-4, acc_forward 7.25e-3.  The fine-packed 64-rank run
(lib_dg: fill store and kept X_NR at 0.002 above the leaf) reads 43.475
65.600 88.04 97.34 72.00, -792566485.70, 5.1e-4, 7.5e-3 - within 0.1%
at every level and 3e-9 in the logdet; the committed packing read
43.475 66.50 101.80 153.95 168.53 and -792566539.70.  So the fine step
is the uncompressed elimination to the digits the quick verification
resolves, and the kept X_NR's step is what is left of the accuracy gap.

**Factor time, and a caution.**  The rule's clean 768^3 run took 137.7 s
against the committed night run's 119.7, and three repeats in one job
took 226 / 193 / 193 s - none of it in the levels' work: the per-level
eliminations, exchanges and transitions equal the night run's within
noise, and the rest is one rank at a time stalled in `cudaHostRegister`
(the pool's pin-ahead): "pinning 62.25 s so far" at level 6 against
6.1 s, the leaf's "after the transition" 60-97 s against 1.9, level 6's
waves 8-10 (where 11 chunks are pinned) 4-14 s of launch on one rank
against 0.9.  The committed library's own morning run shows it too
(leaf transition 24 s on rank 0, max over ranks 28.5 s against 8.1 /
9.9 at night), on different ranks every job, so it is the host memory
state of the nodes, not the code: with 37-50 GB a rank pinned (150-200
GB of a 256 GB node) the kernel's page-pinning slows by 10x when the
node's memory is fragmented.  Compare factor times only within one job,
and read the "pinning N s so far" line first.  The pool's cap and its
trim policy are the lever if this recurs (open).

**CA under the rule (2026-10-11, same day):** unchanged.  384^3 on 8
ranks at 8x: CA0's leaf runs (60.7 s against Color's 56.8) and the level
above dies on the same blocks as before ("near block shape (hi): block
278 x 357, box 16384, neighbor 2343 c 356" on rank 4, identical to the
2026-10-10 run), so the ghost/owner skeleton mismatch is deterministic
and not the packing noise; CA3 exhausts the heap in the leaf's owner
pass (2-3.6 GB requests).  768^3 on 64 ranks: CA0's leaf runs (73.9 s
against Color's 65.6, no leaf speedup at this size) and dies at the
hand-over as before; CA3 exhausts the heap in the owner pass (1-7 GB
requests).  The CA items stand as listed: ghosts' skeletons from their
owners, the owner pass chunked to the room, the interior ranks' host
memory.

#### 768^3 on 80 GB cards (2026-10-11): the leaf comparison is confounded by the heap

The one configuration where the analysis says CA0 should win the leaf -
8x a rank with nothing packed - run on 16 hbm80g nodes (57.5 GB heap,
hold limit 44.6 GB at the leaf against 11.9 on the 40 GB cards):

| leaf, 768^3 / 64 ranks | Color, 40 GB (night) | Color, 80 GB | CA0, 40 GB | CA0, 80 GB |
|---|---|---|---|---|
| elimination | 65.6 s | **144.2 s** | 73.9 s | 89.4 s |
| exchange | 19.4 s | 57.3 s | 0.1 s | 0.1 s |
| host plan boxes [buffers] | 11.6 [11.1] s | 40.2 [39.7] s | 21.4 [20.8] s | 33.9 [33.4] s |
| sketch meta build | 1.0 s | 11.4 s | 3.8 s | 9.3 s |
| owner groups packed / left raw | 718 / 58 | 367 / 409 | 822 / 101 | 408 / 515 |
| waves chunked | 15 of 16 | 7 of 16 | 23 of 32 | 11 of 32 |
| heap free ranges at the level end | 81 | **8224** | 97 | 6046 |
| transient takes, ranges scanned from the top | 145k, 3.7e8 | 142k, **2.5e9** | 161k, 1.1e9 | 159k, 2.1e9 |
| wave time over ranks min / median / max | | 75 / 120 / 132 s | 47 / 62 / 74 s | 54 / 68 / 89 s |

Both leaves got slower on the bigger card, Color's by 2.2x, and the
numbers say why: with the hold limit at 44.6 GB the level packs less
(409 groups left raw against 58), the raw blocks fragment the heap
(8224 free ranges against 81), and a transient take walks the free
ranges from the top - 2.5 billion range visits over the leaf's 142k
takes, 17500 a take.  That is host time in the allocator: "plan boxes
[buffers]" 39.7 s against 11.1, the sketch's meta build 11.4 against
1.0, and the device phases stretched over the launch gaps.  Color pays
it twice, because its waves then wait for the slowest rank's allocator
(exchange 57 s, the plan phase 35 / 69 / 79 s over ranks); CA0's leaf
has no such wait and so "wins" the leaf, 89 against 144 s, for a reason
that has nothing to do with the ghost layer or the exchange.  The
comparison is not meaningful until the allocator's take is not O(free
ranges); the heap already keeps a size index (best fit) for the
resident blocks, and the transient takes should not scan.

**The heap's transient take through size classes (2026-10-11).**  The
free ranges are now also indexed by size class (floor of log2), and a
per-wave take looks at the request's own class for the highest range of
at least the size and at the top of each larger class, where every
range fits: the same placement (the highest range that fits), bounded by
one class's population instead of the free-range count.  192^3 is
bit-identical on all three paths (logdet -89253588.461341903), its leaf
4.41 -> 3.64 s and the scan 3.4e7 -> 8e5 ranges; 384^3 on 8 ranks at 8x
within its run-to-run band (logdet -792567587.529 against .486), its
leaf 57.3 -> 48.3 s (plan boxes 11.8 -> 8.8 s, the scan 3.7e8 -> 7.5e6),
the factorization 91.3 -> 81.2 s; 768^3 on 64 ranks at 40 GB: the leaf
65.6 -> 63.1 s with the exchange wait 19.4 -> 10.6 s (the slowest rank's
plan phase was the allocator's variance: the waves' host time over the
ranks 39.9 / 42.7 / 53.5 s against 41.3 / 50.4 / 58.3 before) and the
wave buffers' allocation 11.1 -> 8.4 s.  So even on the 40 GB cards the
walk had cost about 3 s of compute and 9 s of waiting at the leaf, and
Color's exchange wait at 64 ranks is 20% of its compute (10.6 of 52.5 s),
not 41%.
(The afternoon's 40 GB total stayed at 190.8 s on the same pinning
stalls as the morning's: "after the transition" 69 s at the leaf on a
rank other than 0, level 6 26 s.)

**The 80 GB comparison redone with the size-class take (2026-10-11):**

| leaf, 768^3 / 64 ranks, 80 GB | Color | CA0 |
|---|---|---|
| elimination | 52.1 s | **49.0 s** |
| exchange | 15.0 s | 0.1 s |
| boxes a rank (min / median / max) | 32768 | 36594 / 39512 / 41176 |
| waves' host time over ranks (min / median / max) | 33.4 / 35.9 / 39.0 s | 36.3 / 41.5 / 48.8 s |
| device elimination over ranks | 20.6 / 21.7 / 23.9 s | 22.9 / 26.8 / 32.7 s |
| wave buffers' allocation | 6.2 s | 6.8 s |
| owner groups packed / left raw, waves chunked | 367 / 409, 7 of 16 | 408 / 515, 11 of 32 |
| ranges the allocator visited | 2.1e7 | 8.5e6 |

With the allocator out of the way, CA0's leaf beats Color's by 6%: its
slowest rank (an interior brick, 26% more boxes than Color's) takes
48.8 s of compute where Color's slowest takes 39.0 s and then 15.0 s of
exchange, i.e. the ghost layer's 26% costs less than Color's waiting
(41% of its compute here, 20% on the 40 GB cards where the leaf's
packing evens the ranks out).  That is the result expected on
2026-10-11 05:03 and it holds only where neither path packs: on the 40
GB cards CA0 loses by 7-13%.  Color's whole factorization on 80 GB is
93.0 s against 119.7 at 40 GB (the leaf 52.1 against 65.6, level 6 13.7
against 16.2), solve 1.85 s, acc_forward 1.46e-2, host 30 GB.  CA0
still dies one level up (the ghost mismatch) and CA3 in its owner pass
(5.2 GB requests, 31-36 GB in use), so the leaf is the only level the
comparison covers.

#### Why CA0 dies above its leaf (2026-10-11): the replicas are not the owners' copies

`BPACK_CHECK=replica` at 192^3 on 8 ranks, CA0 at the leaf, the
committed library: at the default heap 6217 of 9488 ghost copies at
level 5 and 3179 of 3536 at level 4 differ bitwise from their owners'
(in T, LU, X_SR, X_RS, X_NR and X_RR_full, a few in X_NR alone; the
skeleton the same), and the run completes; at a 35% heap 6671 of 9488
differ, one of them (ghost 623 on rank 1) in the skeleton itself with
the same count, and the next level's halo dies ("rank 0 expects
177212535 elements from rank 1, which sends 177200035").  The 384^3
check OOM-killed the host (the check keeps the X_NR host copies).

**Correction, same day: a regression from the compression, not a
property of the design.**  The check read 0 differ on every CA level up
to 2026-10-02 (the batch-independence harness covers every batched
kernel), and it had not been run with `H2_compress_factors 1` since.  At
96^3 on 8 ranks, CA0, the same library: compression off, 0 of 3536
ghost copies differ; on, 1948 of 3536 differ, already at the leaf, in T,
LU, X_SR, X_RS, X_NR and X_RR_full.  Some path the compression added
depends on a wave's or a rank's composition (candidates: the sketch's
fill groups formed by a scratch budget over the wave, which can split a
box's sources over different GEMMs on its owner and its replica; the
production chunks and source groups formed by the room; a quantization
step shared by a group).  The fix is to make each of those a function of
the box alone and re-verify with the check (0 differ), which keeps CA0's
leaf free of any exchange; taking the ghosts' skeletons from their
owners instead would add a per-wave dependency on the neighbours, the
thing CA0 exists to avoid.

#### CA0's replicas bitwise again (2026-10-11): two rank-dependent paths the compression added

1. **The sketch's fill.**  The fill sources a wave reads are decompressed
   in groups cut by a scratch budget over the whole wave, and Y -= W P was
   one product a (group, box): a box's sources fell into groups by what
   else its wave held, so its fill was summed in a different order on a
   ghost's replica than on its owner.  Now, with compression, a term a
   source in Morton order (the groups are Morton ranges; raw sources are
   placed by Morton too), a box's terms in one group as passes.  96^3
   CA0: 1948 of 3536 ghost copies differed, now 0; 192^3: 6217 of 9488
   and 3179 of 3536, now 0 and 0.  Uncompressed is unchanged (one group,
   one product a box).
2. **Holding near blocks packed under pressure.**  Decided by each rank
   from its own heap for every block, so a ghost's blocks were quantized on
   one rank and raw on its owner (192^3 at a 35% heap: 5219 of 9488 differ,
   skeletons among them - the crash).  Now a block that exists on more
   than one rank (both its boxes ghosts here or held by another rank as
   ghosts or assisting boxes, from `ca_halo_sends`) is held packed after
   every write when the level decides so at its start (every rank alike: a
   MAX of "my near set would not fit raw" - at the leaf half the full-size
   blocks of the rank's boxes and their one-hop pairs, above it a tenth of
   the card beyond what is used), and a block private to the rank under the
   rank's own pressure, as before.  A private block cannot reach a
   replicated box's elimination: the other copy's rank does not have it.
   Measured on the way: a level-wide decision taken a wave late (a
   non-blocking MAX at every wave) could not follow the green group's first
   wave, which grows the near set by ~20 GB at 8x (the heap went 6 -> 32
   GB, out of memory); packing the shared blocks always cost 1.65 s of a
   4.05 s leaf at 192^3, where nothing needs packing.

| CA0 at the leaf | before | now |
|---|---|---|
| 192^3, 8 ranks: replica check | 6217 / 3179 differ | 0 / 0 |
| 192^3, 8 ranks, 35% heap: the leaf | 5219 differ, halo count mismatch above | 0 differ (level 4 then out of device memory at its start: the 35% heap's limit) |
| 192^3, 8 ranks: levels, total | 4.05 2.28 0.65 0.28 s, 9.04 s | 4.11 2.28 0.67 0.29 s, 9.12 s |
| 384^3, 8 ranks (8x) | dies above the leaf (near block shape) | **end to end**: leaf 46.2 s (Color 48.3), total 83.6 s (Color 81.2), acc_forward 1.50e-2 |

At 384^3 the Color level above CA0's leaf takes 15.3 s against 11.9 in
the Color run, its exchange 4.3 s against 1.1: open.

#### CA0 at 768^3 on 40 GB: three memory steps after the leaf (2026-10-11)

With the replicas fixed, CA0 at 768^3 on 64 ranks got past its leaf
(53.8 s against Color's 63.1) and lost an interior rank to the node's
memory during the leaf's transition: the job gets 241 GB a node (the
cgroup limit; MemTotal shows 263), ~60 GB a rank, and CA0's interior
ranks (41176 boxes against 32768) held 49 GB pinned at the leaf's end and
grew 7 GB more while the transition packed the parents into the next
level's host region (the per-chunk trace now prints each rank's RSS).

1. **The near store compacted at the level end.**  A block's chunk homes
   keep the length of its first compressed write, when it was largest,
   and the blocks shrink to their skeletons as their boxes go, so the
   store held 28.09 GB with 2.15 GB live at the 384^3 CA0 leaf (25.40 /
   1.70 Color; 12.69 / 2.16 at level 5).  finish() copies the live
   blocks out as contiguous runs (build_parent's form, which every reader
   takes) through pageable memory and back into the region's freed
   chunks, in 0.5 s; a first version copied into a freshly pinned region
   and a slow-pinning rank held the level end 21 s.  Results bitwise
   unchanged (384^3 CA0 logdet, 192^3 replica check).
2. **The parents' host copies check the job's room.**  `host_room_bytes()`
   reads the nearest cgroup with a memory limit (Slurm's job, shared by
   the node's ranks) less its usage, plus the pinned pool's idle chunks;
   below a 12 GB reserve a group of parents stays raw on the device.  It
   no longer triggers at 768^3 once (1) frees the store.
3. **The kept X_NR hand-over after the level's blocks are freed.**  It
   needs room for the boxes' compressed X_NR (12 GB on an interior rank)
   and ran first in release_level_data, with 8.4 GB free; the near and
   Schur blocks, the sketch buffers and the adopted blocks now go before
   it (it reads only the fill sources).

768^3 CA0 on 64 ranks, 40 GB, end to end: factor 109.0 s (leaf 56.2,
then 19.8 / 7.6 / 5.6 / 2.7 / 0.9 s), logdet +1 -6968676797.12, acc_mvp
1.05e-3, acc_forward 1.68e-2, acc_backward 6.6e-16, solve 3.44 s,
multiply 3.67 s.  Run to run the pressure decisions follow the heap's
state, which the background copier's frees make timing-dependent, so
the logdet of the same configuration moves in the ninth digit (384^3
Color: -792567588.01 three times, -792567590.66 once).

**Color and CA0 at 768^3 in one job (40 GB, 16 nodes, the same nodes):**

| 768^3, 64 ranks | Color | CA0 |
|---|---|---|
| leaf | 53.1 s (waves 40.1 / 43.0 / 44.5 s over ranks, exchange 10.0) | 56.0 s (waves 45.4 / 50.4 / 53.8, no exchange) |
| levels 6 / 5 / 4 / 3 / 2 | 17.1 / 7.6 / 5.7 / 2.7 / 0.9 s | 19.8 / 7.7 / 5.6 / 2.6 / 0.9 s |
| factor | **102.2 s** | 108.2 s |
| acc_forward / solve | 1.66e-2 / 2.84 s | 1.73e-2 / 3.44 s |

With the heap's walk gone and no pinning stall in either run, Color's
leaf is the faster on the 40 GB cards: its slowest rank waits 10 s in
the exchange, CA0's interior brick computes 9 s longer than Color's
slowest (26% more boxes).  The earlier "CA0 53.8 against Color 63.1"
compared runs on different nodes, Color's with a pinning stall.  On the
80 GB cards (above) CA0 won its leaf by 6%.

#### CA3 at 8x on 40 GB, and the pinned pool (2026-10-11, evening)

The heap-exhausted message now names the request's kind and lists what
holds the heap.  CA3 at 384^3 on 8 ranks, step by step:

1. *Exchange overflow, 2-7 GB single requests*: an owner-wave payload too
   large for the exchange arena fell back to one main-heap block.  Now a
   main-heap block if a range fits, else pinned host memory (the device
   reads it at the same address; the copies that read payloads take
   cudaMemcpyDefault).
2. *Wave outputs, 7.3 GB*: the temp2 blocks were grouped by a bound on the
   source block, which with the device exchange grows by only X_RR_full a
   box, so one temp2 block took a whole unchunked wave.  Both bounded.
3. *Keep failed, no contiguous range* (1.1 GB with 6 GB free in 9000
   ranges): a keep group that finds no single range now splits in halves
   (keep_sources_or_halves) - with the keep failed, the chunked
   production stopped and the next wave asked for 17 GB of outputs.

With these CA3's leaf ran (50.6 s; Color 49.0, CA0 47.0), its transition
and the next level, then the node's host memory ran out (ranks at 57-59
GB of a ~60 GB share).  In the next run the received payloads waiting for
their install took 17.8 GB of the main heap and the leaf's own blocks
found no room; bounding the payloads to the hold limit sent them to the
host instead (and the node ran out) and sent 33 GB of Color's leaf
buffers to the host too (leaf 48.9 -> 64.2 s), so that bound is not in.
**What blocks CA3 at 8x is its install at the colour's end**: every
payload of a colour is held until the colour ends, ~18 GB a rank for the
green colour at 8x, and neither the 40 GB card nor the host share has
it.  Installing the received payloads at each wave boundary (they touch
only other owners' boxes of the current colour, which no wave of the
colour on this rank reads) would bound them to about a wave's worth
(~2 GB); owner_exchange's install is ~270 lines entangled with the
colour's bookkeeping.  Not done.

**The pinned pool.**  Idle chunks now stay pinned between levels and are
unpinned only when the job is short of host memory on its node
(`host_room_bytes() < kHostReserveBytes`); unpinning half at every level
start made the next level pin again, the step that runs 10x slower on a
fragmented node.  The chunks are mmap'd rather than malloc'd: trimming
24 GB of idle chunks had left the RSS unchanged (the allocator kept the
freed memory).  Color and CA0 at 384^3 unchanged bitwise (81.3 / 82.8 s).

Also on 80 GB: CA0's leaf hand-over and transition now fit (20.76 GB of
kept factors staged), and it dies one level up with the same ghost
mismatch as at 384^3 ("near block shape (hi): block 300 x 357, box
147528, neighbor 133486 c 356 (raw)"), so that is CA0's blocker at
every size now; CA3 exhausts even the 57.5 GB heap in the owner pass
(3-6 GB requests with 34-39 GB in use).  Color ran end to end: 190.7 s,
levels 6-2 in 14.1 / 7.0 / 3.1 / 1.3 / 0.9 s (14.1 against 16.2 on 40
GB: level 6 packs nothing there), solve 1.85 s, acc_forward 1.43e-2,
host RSS 30.3 GB.

Left open, measured: the kept X_NR packed finer where it is cheap.  At
384^3 on 64 ranks, the kept X_NR of levels <= 5 at 0.002 (lib_dg) against
0.2 (lib_dgf), both with the fine fill store, read acc_mvp 5.1e-4 against
1.97e-3 and acc_forward 7.5e-3 against 1.36e-2; at 768^3 the levels above
the leaf and the one under it hold 4.6 GB of the 30 GB kept, so it would
cost about 4 GB of host memory a rank at 8x.  Also open: whether the
leaf's pressure packs (parents at 0.2, owner groups at 0.05) inflate level
6 - 72.02 at 768^3 on 8x against 65.6 for the same box size at 384^3 on 1x.

#### CA3 at 384^3 on 40 GB end to end (2026-10-11, night)

Three changes, after which CA3 at 8x runs the factorization and the solve:

1. *The payloads installed at each wave boundary* (owner_install_arrived,
   called from eliminate_wave; install_owner_messages is the former
   colour-end install, now shared).  A wave boundary installs the payloads
   of the current colour that have arrived; the colour's end installs the
   rest.  The GPU sketch filter keeps the colour-end meaning: a box counts
   as eliminated when its component's colour is done, or when this rank
   eliminated it (colour_now, done_here in factorize_owner_level_gpu).  With
   this the leaf no longer holds a colour's ~18 GB of received payloads, and
   the factorization completed (90.9 s).  CA3 at 192^3 is unchanged in
   accuracy: logdet -89253589.36 / .37 against .359 before, acc_forward 7.4e-3.
2. *The owned ghosts' host copies freed after the transition*
   (FactorizationDriver::transition).  An owner level's owned ghosts kept
   their X_NR (temp2) and X_RR_full on the host as fill sources the
   transition may restore.  Nothing read them afterwards, yet they stayed
   through every later level: 6.2 GB a rank at the 384^3 leaf.  With ranks
   4-7 at 57-59 GB, the second node's job had no room for the level-5
   parents' packed host copies, so the parents stayed raw on the card and
   the transition's heap fragmented (4.3 GB free, the largest range 157 MB,
   a 235 MB request).  This only showed under CUDA_LAUNCH_BLOCKING; without
   it the run got through.  With the copies freed, level 6's ghosts hold
   0.04 GB at level 5's end and the RSS there drops from 53.9 to 46.5 GB.
3. *The streamed solve over the waves that read blocks*.  The solve then
   died with an illegal access, reported by the next sync ("solve
   messages").  The leaf's kept factors are streamed (18.75 GB through a
   10.17 GB ring).  CA3's solve waves are the level's global elimination
   sequence (ca_home_order), so a rank has waves with none of its boxes.
   Such a wave has no X_NR segments, and the decompression's prefetch of
   "the next segment" then reached a wave whose blocks were not in the ring
   yet.  It read stale bytes; under CUDA_LAUNCH_BLOCKING: "H2 GPU solve:
   decompressing a wave's X_NR failed".  The kept-block prefetch
   (prepare_wave_kept) now goes to the next wave that reads blocks, and the
   ring plan (plan_kept_streaming) keeps consecutive reading waves apart,
   not just w and w + 1.  Color's and CA0's waves are never empty, so they
   are unchanged: Color's logdet is bitwise the earlier runs'.

Same job, 2 nodes, 8 ranks, 384^3 (lib_det22):

| | total | levels 6 / 5 / 4 / 3 / 2 (s) | logdet | acc_mvp | acc_forward | acc_backward | solve |
|---|---|---|---|---|---|---|---|
| Color | 81.4 s | 49.0 / 12.3 / 5.1 / 2.6 / 0.4 | -792567588.01 | 1.38e-3 | 1.57e-2 | 1.2e-15 | 2.39 s |
| CA3 | 85.8 s | 50.2 / 13.7 / 5.1 / 3.0 / 0.4 | -792567667.82 | 1.23e-3 | 1.54e-2 | 8.2e-16 | 2.49 s |

CA3's logdet moves by ~4 in the 9th digit from run to run (-663.4, -665.4,
-667.8): payloads now install in arrival order.  CA0 in an earlier job of
the day: 82.8 s, leaf 46.6 s.  So at 384^3 CA3 runs but is not faster
than Color.  Its leaf spends 45.0 s in the owner schedule's waves and 3.1
s in colour transports, and sends 25.5 GB a rank.

#### Color, CA0 and CA3 at 768^3 on 64 ranks, 40 and 80 GB (2026-10-11, late)

CA3 at 768^3 on 64 ranks (8x the boxes a rank, CA at the leaf, level 7)
died in four ways after d29804d.  Each fix is in level_eliminator.hpp.

1. *The generators' original X_NR in one block a wave* (`d_xnr`, 4.27 GB,
   with 1.7-3.0 GB the largest free range).  Now a block per production
   chunk, freed after the chunk's generator copies.  Allocated this way,
   the source and temp2 groups have to follow the chunks exactly.  On a CA
   level, though, the wave's boxes were sorted by LU size class *after*
   being assigned to chunks.  A box whose class moved it to an earlier
   chunk would then read a temp2 block not yet allocated.  This never
   happened: at 216-point leaf boxes nearly every r is in the 256 class.
   The groups and chunks are now formed after the sort.  Side effect: with
   less held at the pressure decisions, Color's leaf packs a few blocks
   differently.  It is no longer bitwise the earlier runs' (logdet
   -792567588.01 -> -583.07 at 384^3, acc_forward 1.57e-2 -> 1.51e-2), as
   is CA0 (-541.11 -> -543.11).
2. *Owner payloads, 10.9-14.8 GB a rank on the main heap*.  A payload
   lived about three waves: issued at the sender's next wave boundary
   (once its copies had run), received at the receiver's, installed at the
   one after.  The transport (owner_progress, owner_install_arrived) now
   also runs after each production chunk and at the wave's end.  At a
   chunk boundary the near store is consistent (the chunk's packs have
   run), and `meta_` is rebuilt by each chunk.  Components of one colour
   are not adjacent, so an install never touches a block of this colour's
   waves.  Tried first and dropped: keeping the payloads on the main heap
   only within the hold limit, and in pinned host memory beyond it.  That
   killed the host at 384^3 and at 768^3.  At 384^3 the node has no room
   (ranks at 52-59 GB).  And the pinned pool reuses an idle chunk only if
   it is at least the request's size, so payloads of varying size pin new
   chunks while the old ones sit idle.
3. *The host, in the interior waves* (RSS 63 GB a rank, the node's share
   is 60).  The owned ghosts' temp2 was copied to the host as a fill
   source for the transition, which reads the packed temp2^T first anyway.
   With the packed fill store (packed_restore_) the copy is no longer
   made (host_drops_ghost_xnr); X_RR_full still is, since no kept LU
   rebuilds it.  CA3's pageable host memory grew by 13 GB through the
   green colour before; the leaf's RSS peak is now 55.3 GB against
   Color's 54.0.  At 384^3 the leaf also got faster (49.8 -> 47.8 s:
   less copied off the card).
4. BPACK_TRACE=wave now prints the owner payloads held at each wave's end
   (sends and receives).  The largest at 768^3: 14.2 GB on 40 GB (6.9 +
   7.3), 17.9 GB on 80 GB.

384^3 on 8 ranks with all four (lib_det26): CA3 82.0 s (leaf 47.8), logdet
-792567657.24, acc_forward 1.57e-2, acc_backward 6.7e-16; 192^3 logdet
-89253589.360, acc_forward 7.5e-3.

768^3 on 64 ranks (16 nodes), lib_det26, each card size in one job.  Times
in s; the levels are the shortest rank's elimination (the longest is
within 0.15 s at the leaf).

| card | variant | total | levels 7 / 6 / 5 / 4 / 3 / 2 | solve | logdet | acc_forward | acc_backward | leaf RSS peak |
|---|---|---|---|---|---|---|---|---|
| 40 GB | Color | 100.0 | 52.9 / 15.7 / 7.4 / 5.5 / 2.7 / 0.9 | 2.99 | -6968677311.1 | 1.65e-2 | 7.9e-16 | 54.0 GB |
| 40 GB | CA0 | 106.5 | 54.7 / 19.2 / 7.9 / 5.5 / 2.7 / 0.9 | 3.61 | -6968676793.7 | 1.72e-2 | 6.7e-16 | 59.0 GB |
| 40 GB | CA3 | 105.6 | 53.9 / 17.0 / 7.5 / 6.1 / 3.6 / 0.9 | 3.23 | -6968678017.3 | 1.81e-2 | 7.6e-16 | 55.3 GB |
| 80 GB | Color | 90.1 | 50.0 / 14.4 / 7.0 / 3.0 / 1.3 / 0.9 | 1.94 | -6968673482.4 | 1.47e-2 | 8.3e-16 | 42.9 GB |
| 80 GB | CA0 | 99.3 | 51.6 / 20.9 / 7.1 / 3.0 / 1.3 / 0.9 | 2.14 | -6968673036.6 | 1.45e-2 | 8.1e-16 | 47.6 GB |
| 80 GB | CA3 | 90.7 | 50.5 / 13.7 / 7.0 / 2.8 / 1.3 / 0.9 | 2.00 | -6968674648.0 | 1.47e-2 | 7.4e-16 | 43.1 GB |

Neither CA variant beats Color at the leaf at 8x on 64 ranks.  CA3's leaf
is within 1-2% of Color's on both cards (53.9 vs 52.9, 50.5 vs 50.0);
CA0's is 3% slower.  The leaf's exchange wait drops from 10.4 s (40 GB) /
13.6 s (80 GB) for Color to 3.1 / 2.2 s for CA3, but the owner schedule's
extra compute takes that back.  Above the leaf, CA3 is slower than Color
on 40 GB (levels 6, 4, 3: +1.3, +0.7, +0.9 s) and level with or faster
than it on 80 GB.  CA0's level 6 is 3.5-6.5 s slower than Color's
on both.  The logdets differ by ~1e-7 between variants and by 5e-7
between card sizes (what packs differs); for scale, the variation run to
run on 40 GB was Color 100.0-104.8 s, CA0 106.5-113.0 s.

#### Why CA0 stopped beating Color on 80 GB, and what was left (2026-10-12)

The "CA0's leaf 49.0 s against Color's 52.1 s" on 80 GB at 768^3 (11 Oct,
noon) had two things in it that were not CA0's speed:

1. That CA0 diverged.  Its ghost copies were packed by each rank's own
   pressure, and it stopped above the leaf ("near block shape (hi)").  The
   fix (blocks on more than one rank packed by a level-start rule) cost ~3
   s of the leaf's owner pass on 80 GB, where every write of a shared block
   was packed.
2. A one-time MPI cost.  GPU-aware MPI meets a buffer's whole allocation
   on a pair's first transfer from it, and a level's larger messages land
   in the main heap (one allocation of most of the card) once the 2.65 GB
   exchange arena is full.  A Color leaf paid this in its first four
   exchanges (up to 1.1 s of payload wait each at 384^3; Color's leaf
   exchange 6.56 -> 4.39 s once warmed).  A CA0 leaf exchanges nothing, so
   the level above it paid instead: +3.4 s at 384^3, +6.5 s at 768^3 on
   80 GB.  CA0 does no communication between its waves (checked: no
   transport, no MPI in the CA loop); its leaf only moved the cost.

Changes, all in one commit:

- The exchange warm-up also sends 64 KB through the main heap to every
  peer.  At 768^3 it takes 3.4 s (40 GB) / 6.7 s (80 GB), outside the
  factor time but inside the "total time" line, so for one factorization
  the cost moves rather than disappears.  The heap's capacity is unchanged.
- The shared-block rule ends with the ghosts' groups (replicas_done at the
  interior group, the same wave on every rank: no ghost is eliminated
  after it, so no copy on another rank is read again).  On 40 GB this
  changes nothing measurable (pressure packs the interior anyway: 504 GB
  packed either way).
- The level-start estimate counts the boundary groups' boxes (the rank's
  boundary boxes and its ghosts) at their full-size near set.  Half that
  count let 38 of 64 ranks fall back in the green group on 80 GB (37.0 GB
  estimated against 40.9 GB of room; the raw peak was 41.5 GB at the wave
  ends), and the level above stopped.  At the full count both cards pack.
- Replica check, 192^3 CA0: 0 of 9488 / 3536 ghost copies differ, at the
  default heap (shared blocks left to pressure) and at a 35% heap (packed;
  that run then exhausts the heap at level 4, as it did before).

768^3 on 64 ranks, each card in one job, Color and CA0 twice each (factor
time excludes the warm-up; total includes it):

| card | variant | factor | total | warm-up exch. | leaf | level 6 | solve |
|---|---|---|---|---|---|---|---|
| 40 GB | Color | 105 / 100 s | 108.4 / 103.5 | 3.4 / 3.6 | 54.2 / 51.8 | 17.0 / 16.6 | 3.03 / 2.99 |
| 40 GB | CA0 | 106 / 106 s | 110.0 / 109.1 | 3.6 / 3.4 | 56.0 / 55.9 | 17.1 / 17.0 | 3.59 / 3.59 |
| 40 GB | CA3 | 106 s | 109.2 | 3.4 | 54.1 | 17.3 | 3.23 |
| 80 GB | Color | 86.1 / 85.1 s | 92.7 / 91.8 | 6.7 / 7.1 | 44.8 / 44.6 | 14.5 / 14.5 | 1.92 / 1.92 |
| 80 GB | CA0 | 93.2 / 93.3 s | 99.7 / 99.7 | 6.9 / 6.8 | 51.3 / 50.9 | 14.2 / 14.7 | 2.13 / 2.13 |

CA0's level 6 now matches Color's.  Its leaf does not: the leaf ends with
the slowest rank, an interior one with 41,176 boxes (its ghost layer, +26%)
against Color's 32,768.  Its waves take 49.3 s of host time against
the slowest Color rank's 37.4 s (80 GB), and 53.6 s against 44.3 s on
40 GB, more than Color's 8-9 s of leaf exchange.  At 8x the boxes a rank on
64 ranks the replication costs more than the communication it saves, on
both cards.

#### The heap's host time, and an owner group's raw blocks in one block (2026-10-12)

The box region's "buffers" timer was not restarted per production chunk, so
from the second chunk on it also counted the previous chunk's launches,
compressions and keeps: 6-8.5 s a leaf, of which 0.26-0.56 s is buffers
(fixed in 39b9fa6, which also times the heap under BPACK_TRACE=wave).
Measured at the 384^3 leaf on 8 ranks (the same 32,768 boxes a rank as
768^3 on 64): the heap does 3.8M allocations and 3.9M frees, 5.3 s of
host time, on 40 GB cards (2.1M and 1.7M, 5.0 s, on 80 GB).  Most are the
owner pass's: each group puts its packed near blocks back raw, a heap
block each (~105 KB), and under pressure packs them again and frees them
(2.85M blocks a leaf on 40 GB, 0.66M on 80 GB).  The GPU idles 5.4-5.8 s
inside the leaf's waves.

Fix 1: an owner group allocates one heap block for the raw blocks its pack
will take back, foreseen by the pack's own rule (pressure from what the
group adds; on a replicated level's ghost groups the shared-block rule),
and cuts them out with a bump pointer.  The pack's decision stays where it
was, after the GEMMs; a block it leaves raw after all is copied to a block
of its own before the group's goes back.  Sizes round as the heap's do, so
the heap's count at every decision is unchanged and the results are
bitwise the same.

384^3 on 8 ranks, alternating builds on the same nodes (Color, two runs
each):

| card | | heap allocations / frees | heap host time | leaf | factor time |
|---|---|---|---|---|---|
| 40 GB | before | 3.81M / 3.91M | 5.26 s | 47.35 / 47.17 s | 79.7 / 79.3 s |
| 40 GB | fix 1 | 1.28M / 1.38M | 2.31 s | 45.15 / 45.12 s | 77.5 / 77.3 s |
| 80 GB | before | 2.06M / 1.68M | 4.87 s | 39.10 / 39.10 s | 68.7 / 68.6 s |
| 80 GB | fix 1 | 1.44M / 1.06M | 4.03 s | 38.94 / 38.67 s | 68.6 / 68.4 s |

On 40 GB 716 groups took 2.53M blocks this way, none moved out; on 80 GB
(where 409 of 776 groups stay raw) 369 groups took 0.63M blocks, 11,892
moved out.  Color, CA0 at 384^3 and CA0 at 192^3 bitwise unchanged; the
192^3 replica check: 0 of 9488 / 3536 ghost copies differ.  CA0's leaf
took 3.0M blocks this way, CA3's 2.6M.

What is left on 80 GB is the transient search (19.8M range look-ups a
leaf: each step of the size-class walk looks the range up again in the
address map) and the per-call bookkeeping (five or six node containers):
fix 2, not done.

#### CA0's boundary phase at 768^3 on 64 ranks, measured by rank and wave (2026-10-12)

BPACK_TRACE=wave now prints a line per wave on every rank ([wave-row]: boxes
and ghosts, host phases, the owner pass's host parts, device times, the
shared-block writes it packed); h2_c1_work/wave_rows.py sums them by phase.
One job, 80 GB cards, 768^3 on 64 ranks (lib_det32); "CA0 np" is a
measurement-only build in which the blocks on more than one rank are never
packed (no level rule, no fallback).

The leaf's slowest rank (35, interior: 5,768 own boundary boxes and 8,408
ghosts) and a corner rank (0: 3,826 ghosts):

| rank 35 | boundary phase | of which waiting | owner GEMMs (dev) | owner host mat + pack | interior phase |
|---|---|---|---|---|---|
| Color (8 waves) | 5.12 s | +4.75 s between waves | 1.42 s | 0.07 s | 31.6 s |
| CA0 (24 waves) | 18.64 s | 0 | 9.91 s | 5.30 s | 30.1 s |
| CA0 np | 12.30 s | 0 | 3.43 s | 0.14 s | 31.4 s |

| rank 0 | boundary phase | of which waiting | interior phase |
|---|---|---|---|
| Color | 3.95 s | +5.67 s | 29.6 s |
| CA0 | 10.50 s | 0 | 28.9 s |
| CA0 np | 6.76 s | 0 | 29.8 s |

So the slowest rank's 18.6 s is: its own boundary boxes ~5 s (as Color's),
the ghosts ~7.2 s (0.80-0.86 ms each, as an own box), and packing the
blocks on more than one rank ~6.3 s (1.66M packed writes, 276 GB raw).  The
wave count costs nothing measurable: a fit over 1,536 boundary waves gives
no per-wave constant (wall = 0.865 ms per own box + 0.804 ms per ghost,
CA0 np).  Leaves: Color 44.1 s, CA0 50.6 s, CA0 np 45.1 s; factor time
84.8 / 93.1 / 87.5 s.  CA0 np completed every level (no rank packs a
shared block, so the copies cannot differ; quick verification passes) and
held at most 53.2 GB at a wave end, against 40.0 GB packed and a 43.3 GB
hold limit.  Without the packing CA0 is 2.9 s ahead of Color on a corner
rank's boundary phase and 2.4 s behind on an interior one; the leaf ends
with the interior ranks.

#### FLOP rates of CA0's colour groups (2026-10-12)

The trace rows now carry each wave's floating-point work by device phase
(the sketch's GEMMs and stored-block products, the ID's pivoted QR at its
full size - an upper bound -, the box region's GEMMs, LU and solves, the
owner GEMMs; kernel values, gathers and adds count none);
h2_c1_work/flop_groups.py sums them by colour group.  768^3 on 64 ranks, 80
GB, one job; the slowest rank (interior), GFLOP/s per GPU-busy second and
per wall second:

| run | group | waves (boxes) | GFLOP | busy / wall | rate busy | rate wall | owner pass | sketch |
|---|---|---|---|---|---|---|---|---|
| CA0 | blue/orange | 8 x 128 | 5,324 | 1.25 / 1.36 s | 4,257 | 3,906 | 1,610 | 17,004 |
| CA0 | purple | 8 x 468 | 19,690 | 4.46 / 4.90 s | 4,413 | 4,021 | 1,663 | 17,773 |
| CA0 | green | 8 x 1,176 | 49,877 | 11.18 / 12.33 s | 4,462 | 4,045 | 1,739 | 17,981 |
| CA0 | interior | 8 x 3,375 | 135,364 | 25.25 / 30.17 s | 5,362 | 4,486 | 2,223 | 16,935 |
| CA0 np | blue/orange | 8 x 128 | 5,324 | 0.69 / 0.80 s | 7,667 | 6,643 | 5,017 | 17,019 |
| CA0 np | purple | 8 x 468 | 19,690 | 2.54 / 2.99 s | 7,743 | 6,587 | 4,923 | 17,794 |
| CA0 np | green | 8 x 1,176 | 49,878 | 6.80 / 8.43 s | 7,339 | 5,915 | 4,879 | 17,244 |
| CA0 np | interior | 8 x 3,375 | 135,365 | 26.55 / 31.30 s | 5,098 | 4,325 | 2,093 | 16,923 |
| Color | boundary | 8 x 721 | 28,388 | 3.95 / 4.76 s | 7,184 | 5,966 | 5,245 | 18,573 |
| Color | interior | 8 x 3,375 | 135,009 | 25.34 / 31.87 s | 5,327 | 4,236 | 2,297 | 16,967 |

The interior is not the efficient phase.  Every box costs ~5 GFLOP (ghosts
5.2-5.3, interior 5.0), even the 128-box waves keep the GPU as busy as the
3,375-box ones (idle 14-19% of the wall in every group), and with no
packing the boundary groups run at 7.3-7.7 TFLOP/s busy against the
interior's 5.1-5.4.  The difference is the owner pass's interval: ~5
TFLOP/s in the boundary groups, 2.1-2.4 in the interior, where the heap is
past the hold limit and the near blocks are put back and packed again
inside it (the same in Color), plus ~5 s a rank of kernels with no
floating-point work (kernel values, gathers, near-block adds).  The shared
packing takes CA0's boundary owner pass from ~5 down to 1.6-1.7 TFLOP/s.
The sketch's GEMMs run at 17-18.6 TFLOP/s (FP64 tensor cores, peak 19.5),
the box region's at 6-7, the ID's QR at ~1.1.

#### Owner-pass updates read nowhere on the rank (2026-10-12): ~1%, not 1-2 s

The idea: on a replicated level a rank updates near blocks that only the
other ranks' copies need.  The owner pass already skipped pairs of two
eliminated ghosts.  Now also an eliminated ghost's own Schur block (its
parent is another rank's; build_parent reads only pairs with a box of this
rank) and its pairs with boxes beyond the halo (never eliminated here).
Counted at the 384^3 leaf on rank 0 (8 ranks, 40 GB): pairs of two
eliminated ghosts 59 GFLOP, the Schur blocks 16, the pairs beyond the halo
394, against 38,141 GFLOP made - 1.2% in all, 1.1% newly skipped.  CA0's
leaf 45.40 / 45.38 s before, 45.58 / 45.52 s after (factor time 78.1-78.2 s
in all four).  192^3 replica check (CA at the leaf; and at levels 4 and 5):
0 differ, logdet bitwise unchanged; Color bitwise unchanged.  The estimate
of 1-2 s that proposed it assumed the dead updates were a large share;
they are not - once two ghosts are both eliminated, few later
eliminations touch their pair.

#### What the transient actually is (2026-10-08): the sketch's [P; W], 5.41 GB

With transient blocks tagged, the 192^3 leaf reads, identically on both paths:

| transient at the peak | GB |
|---|---|
| **sketch [P; W]** | **5.41** |
| sketch fill slots | 0.54 |
| codec staging | 0.46 |
| sketch lists | 0.29 |
| sketch vectors | 0.13 |
| wave result | 0.14 |
| wave sources | 0.08 |
| total | 7.05 (accounted 7.05) |

So the single largest item in the whole 15.40 GB peak is neither the near
blocks (3.24 resident) nor temp2^T (2.34): it is the sketch's per-box
`[P; W]`.  The transient peak is reached in the *sketch*, at the head of a
wave, before the solves have produced anything - which is why freeing 5.14 GB
of wave outputs moved it by 0.61.

`P` is `fill_rank x n` and `W` is `d x fill_rank`, where
`fill_rank = sum of r over the box's fill sources`, and **every box of the
wave has its pair allocated at once**: one loop allocates them all
(level_eliminator.hpp, the sketch's plan loop) and one frees them all after
the sketch.  At the leaf that is 256 boxes x about 20 MB.

**Done (2026-10-08), and it was free.**  The sum is over sources, so a group's
terms are accumulated into Y and the scratch reused: for each fill group, build
that group's P_E, sketch its W_E, then `Y -= W_g P_g` with beta = 1 before the
next group overwrites the block.  The groups already existed (a group's
temp2^T slots are only in the fill scratch while its batches run), so the only
new thing is the per-group offsets, and one block sized for the widest group
replaces a stack per box.

| 192^3 leaf | peak | transient | Factor |
|---|---|---|---|
| Color before | 15.40 | 7.05 | 9.48 s |
| Color after | **13.34** | 5.95 | 9.39 s |
| CA0 before | 17.47 | 7.04 | 10.6 s |
| CA0 after | **15.40** | 5.95 | 10.4 s |

-2.06 GB on both paths, **bit-identical logdet** on both at both sizes, and
marginally faster.  I expected the summation order to move the result by
rounding and it did not: the groups are swept in the order the stacked product
summed over k, so the terms are added in the same order.  That is this
layout's property and not a guarantee, so the gate for any change here stays
logdet.

`sketch P;W` is now below the reporting threshold, and the transient reads
wave work 1.09, **wave outputs 3.36**, source window 0.54, codec staging 0.52,
wave result 0.14, sketch vectors 0.13, wave sources 0.08.  **The peak has
moved into the window where the raw outputs are live**, which is what stage C
shortens - so stage C, worth 0.62 GB when it landed, is now load-bearing, and
the outputs are the next target again at 3.36 GB.

#### X_RR_full is rebuilt, not stored (2026-10-08)

The user's observation: X_RR_full is completely reconstructible.  It is the
symmetrized X_RR as it stood before getrf factored it in place
(level_eliminator.hpp, the gather beside the symmetrize), and a kept box's LU
is resident for the solve, so `X_RR_full = P L U` of what is already there.
Rebuilding costs rounding - 1e-16 relative, against the 1e-3 the fill sources
are held to - where compressing it would have cost 1e-3.

A box gives its stored copy back as soon as its factors are kept, and the two
readers rebuild what they need: the sketch into a window a fill group at a
time (one rebuild a (group, source), not one per reader), and the transition
into a block a source.  A remote box keeps the stored copy, since its LU is on
its owner's rank, and so does a level that keeps no factors.

| 192^3 leaf | peak | resident | transient | Factor |
|---|---|---|---|---|
| Color before | 13.34 | 7.38 | 5.95 | 9.39 s |
| Color after | **12.84** | **6.11** | 6.73 | 9.42 s |
| CA0 before | 15.40 | 9.45 | 5.95 | 10.4 s |
| CA0 after | **14.57** | **7.85** | 6.72 | 10.4 s |

`fill X_RR_full` 0.90 -> 0.30 GB (what remains is the current wave's, given
back at the keep), and the resident falls 1.27 GB on Color and 1.60 on CA0.
The peak falls less - 0.50 and 0.83 - because the peak moment moves and the
transient takes some of it back.  logdet moves by 1.7e-16 relative at most.

**The first version cost 2.1 s a run, +22%.**  I wrote the product as a kernel
with a dot product a thread and stride-ldlu loads - a hand-rolled GEMM, and a
bad one.  Splitting L and U out with an elementwise kernel, multiplying with
magma's batched GEMM and undoing the pivots afterwards made it free (+0.03 s).
The lesson is narrow but worth keeping: a product of two r x r matrices is a
GEMM, and there is a batched one to hand.

**The fix had a shape this codebase has used twice.**  A box's sketch is
independent of the other boxes', so the wave's boxes can be swept in groups
bounded by a byte budget, each group's `[P; W]` freed before the next - which
is step 3a's owner groups and the fill-source groups again.  Nothing about
the arithmetic changes, so the gate is bit-identical logdet on both paths.
With a 512 MB budget, 5.41 -> about 0.5 GB, so the transient 7.05 -> about
2.2 and the peak 15.40 -> about 10.5 GB on Color and 17.47 -> 12.6 on CA0.
The cost is more, smaller launches inside the 1.39 s sketch.


The solve factor the code calls X_NR is `temp2` = -X_NR X_RR^{-1}: the host's
Box::X_NR is downloaded from `wb.off_temp2` and the kept device factor is
compressed from the same buffer (`srcs.push_back(source + wb->off_temp2)`).
`wb.xnr` = A_NR - A_NS T is an intermediate, surviving only as
box->lazy_original_x_nr where has_remote.

So after a wave, temp2 is resident twice:

- `persist_temp2t` / `t2t_packed` = temp2^T (r x nt), the fill source later
  waves' sketches read - "fill temp2^T packed", 1.38 GB at the 192^3 leaf;
- `kb.xnr` = temp2 (nt x r), the kept solve factor - "kept factors", 1.08 GB
  at the leaf's peak with C4 on, 2.57 GB by the level's end.

Two compressed copies of one matrix, one of them transposed, and
launch_transpose already makes the second from the first. Deriving one on
demand would free about 1.1-1.4 GB of the leaf's resident set.

Before treating that as a saving: (a) do the two sets cover the same boxes -
the fill source exists only where gpu_sketch_ and lazy Schur apply, the kept
factor covers every box the level keeps; (b) the lifetimes differ - the fill
source goes at the transition, the kept factor must reach the solve - so the
kept factor would have to be the survivor, with the sketch path paying a
decompress plus transpose. Queued behind the near-field work, which is
7x larger.

### Step 2: the parent blocks, 74% of what is left

Measured and de-risked, not started.  They are **not** a separate store: the
next level adopts them as its own st.schur and edges_ (begin() does
`st.schur = it->second`).  I had assumed they were accumulated into across
waves, which would have made this the hard case; they are not.  Step 5
allocates a *fresh* block per edge and replaces the old entry, and st.schur is
only read during a wave (the three gathers for X_RR, X_SR and A_SS).  So they
are written once and read about once, like the fill sources:

| level | Schur: largest wave / reads per block | edges: largest wave / per block |
|---|---|---|
| 5 | 0.13 GB / 1.01 | 2.74 GB / 1.32 |
| 4 | 0.04 GB / 1.00 | 0.57 GB / 1.21 |

So the level that adopts 9.47 GB reads **0.61 GB** of it in its largest wave
and 6.63 GB over the level - a 0.6 GB scratch and 0.018 s of decompression -
and the transition's own live set is one chunk (0.82 GB, 0 revisits).

The plumbing is the work, not the design.  DeviceMatrixT is a plain
{ptr, rows, cols}, so the compressed form needs a parallel store in
DeviceLevelBlocksT, threaded through build_parent (compress per chunk),
begin() (adopt packed), the wave reads (decompress per wave),
download_level_blocks (the host path) and the between-level compaction (which
relocates blocks and must carry the packed ones).

**The ratio for these blocks is still unmeasured.**  They are Schur
complements rather than the smooth one-hop couplings, and temp2^T already came
in at half X_NR's ratio, so 9.47 GB may not go to 0.5 GB.  The round trip
reports it on the first run.

### Step 2 piece 2: freeing the raw blocks is where it stands (2026-10-06)

Piece 1 is in (4b01433): build_parent packs each chunk's blocks into
DeviceLevelBlocksT::schur_packed / edges_packed, one allocation per group, and
the raw blocks stay so results are bit-identical.  7.2x at the leaf over 7099
blocks.

Piece 2 - freeing the raw blocks once packed and putting them back at the
adopting level - **works, and buys almost no peak** (2026-10-06).  Three
configurations at 96^3 over 8 ranks, all exit 0:

| run | logdet | acc_forward | acc_backward | Factor |
|---|---|---|---|---|
| Color, compress 0 | -9928892.0403490998 | 3.46e-03 | 1.24e-15 | 1.16 s |
| Color, compress 1 | -9928907.2811287567 | 9.33e-03 | 9.62e-16 | 1.29 s |
| CA0, compress 1 | -9928906.4595639911 | 9.26e-03 | 7.29e-16 | 1.49 s |

logdet moves by 1.5e-6 relative and acc_forward by 2.7x.  That is a real cost
and it is expected: unlike X_NR, which only the solve reads, and unlike the
fill sources, which are an update, a parent block *is* what the next level
factorizes, so its quantization error enters the factorization.  At tol 1e-3
acc_forward stays the same order it already was.

**But the peak does not move, and 192^3 says so exactly.**  Color, 8 ranks,
compress 1, with piece 2: the leaf peak is **12.78 GB**, the same number
piece 1 recorded, and `Factor time` 9.20 s against piece 1's 8.56 s and 7.71
s uncompressed - piece 2 costs 0.6 s for no memory.  logdet
-89253766.641603738, 2.0e-6 relative from uncompressed, in line with 96^3.

The leaf breakdown is the whole story:

    level 5 resident at the peak: 12.78 GB | fill X_RR_full 0.73 | fill
    temp2^T packed 1.38 | level Schur/edges 1.14 | parent blocks 9.47 |
    near 0.02 | points 0.04
    level 5 parent blocks compress: 6.80 -> 0.94 GB (7.2x) over 7099 blocks

**That 9.47 GB is not the parent blocks, and reading the label as though it
were cost two runs.**  ResidentKind::Transition prints as "parent blocks" but
four different allocations shared it: build_parent's parent blocks,
begin()'s upload of the host's blocks, begin()'s restore of the adopted
blocks, and the owner pass's new near-field blocks.  Level 5 is the *leaf*:
it adopts nothing, and its own parent blocks are 6.80 GB raw / 0.94 GB
packed, so its 9.47 GB could never have been parent blocks.  Piece 2
producing a byte-identical 12.78 GB peak was the signal, and it was read as
"the peak sits elsewhere in the level" rather than "the label is wrong".

Splitting the kind into Uploaded / Adopted / OwnerNear says exactly where the
memory is (192^3, 8 ranks, Color, compress 1, piece 3):

| level | peak | the dominant item | the rest |
|---|---|---|---|
| 5 (leaf) | 12.78 GB | **owner near blocks 9.47** | Schur/edges 1.14, fill 2.11 |
| 4 | 7.74 GB | **adopted blocks 6.29** | fill 1.30, kept 0.06 |
| 3 | 3.52 GB | adopted 1.76 | kept 0.85, fill 0.89 |
| 2 | 2.25 GB | kept 1.05 | adopted 0.54, fill 0.65 |

So:

- **the global peak is the leaf's, and 74% of it is the owner pass's
  near-field blocks** (level_eliminator.hpp, the join that allocates
  `edges_[t.key]` for each new owner target).  Step 2 never touched them.
  This is the same 9.55 GB measured earlier as "the part that is not fill
  sources", now located at its allocation site;
- **piece 4 is in, and it is worth 1.5 GB at level 4, not the 6 GB
  projected** (see below); the 12.78 GB global peak is untouched.

Which means the lever on peak device memory - the thing that decides how many
boxes fit on a rank - is the near-field blocks, not the parent blocks.  They
are harder than anything in step 2: the waves *accumulate* into them, so
holding them packed means decompress-update-recompress per wave, which is the
trade of time for memory already anticipated.  Piece 4 is still worth doing,
but it should be understood as flattening levels 4-2, not as raising the
number of boxes per rank.

**The 96^3 breakdown says the same thing and says why.**  At level 4,
compression on: peak 1.91 GB of which *parent blocks 0.73 GB*, against a
compress report of 0.77 GB raw in all.  So at the moment of the peak the
whole raw set is still allocated, because build_parent allocates every parent
block up front (the two loops over `pairs_of`) and only then runs the chunk
loop that packs and frees them.  Piece 2's window - between a chunk's pack
and the next level's begin() - contains nothing else large, and begin()'s
restore puts the whole raw set back at once.

So the saving the step was measured for needs two more pieces, and the
measurement says which one is worth doing first:

- **piece 3: allocate a parent block when its chunk is planned, not up
  front.**  A chunk's blocks are final (0 revisits, measured), so the
  transition's live set becomes one chunk's raw (0.82 GB at 192^3) plus the
  packed accumulation, in place of 9.47 GB.  This is the leaf level's peak -
  the leaf adopts nothing, so its 12.78 GB peak is its own transition - which
  makes it the largest single item left.  The work: `chunk_blocks` and the
  `if (M.ptr != nullptr)` guards dedupe by pointer today and must dedupe by
  key, since ptr is null until the chunk starts; shapes stay where they are.
- **piece 4: the adopting level reads the packed form per wave** instead of
  begin() restoring all of it - 0.61 GB of scratch at level 4 against the
  6.63 GB it reads over the level.  This is the original step 2 design.

### Piece 4, done (2026-10-06): 6.29 -> 3.85 GB at level 4, and not at the peak

begin() keeps the packed store (packed_schur_, packed_edges_, strides,
packed_owned_) and eliminate_wave materializes what its wave reads, in one
batch, before the planner reads it.  sketch_wave does not touch these blocks,
so the hook is in one place.  Verified at 96^3 Color and CA0 and at 192^3
Color: logdet -9928907.2811287567, -9928906.4595639911 and
-89253766.641603738, each exactly the pieces 1-3 value; acc_backward 9.6e-16;
`0 left packed` on every level; materialization 0.026 s for 7099 blocks in 12
batches.  `Factor time` 8.63 s at 192^3, against 8.84 piece 3 and 9.20 piece
2, so it pays for itself in time as well.

| level | piece 3 | piece 4 | adopted set |
|---|---|---|---|
| 5 (leaf) | 12.78 GB | **12.78 GB** | adopts nothing |
| 4 | 7.74 | **6.24** | 6.29 -> 3.85 |
| 3 | 3.52 | 3.62 | 1.76 -> 1.60 |
| 2 | 2.25 | 2.31 | 0.54 -> 0.54 |

3.85 GB rather than the 0.9 GB the read union suggested, because the Schur
blocks are now materialized up front (below) and a block a wave materializes
but does not retire stays for the level.  Levels 3 and 2 come out marginally
worse, the packed form and the up-front Schur blocks costing a little more
than the per-wave saving at that size.  **The global peak does not move**: it
is the leaf's, and the leaf adopts nothing.

#### What it cost to find: four readers and one class of writer

A null ptr means "there is no block here" everywhere in the level, and every
reader answers that from the kernel or by making a fresh block.  An adopted
block left packed is therefore read as the original A, with no error - so
every site that acts on a null ptr now calls assert_not_packed() and names
itself.  That is what found the readers, in four runs: the owner pass's join,
replace_edge in the remote-generator pass, ensure_edge's inner loop, and
build_parent's diagonal child-Schur gather.  Each is now pre-materialized in
one batch (run_fill_batch plans and synchronizes twice per call, so a call per
block is not an option), with the generator pass's three conditions hoisted
into lambdas so the pass and its pre-pass cannot drift apart.

The last one was not a reader at all.  logdet stayed wrong at -9813044.20
with every guard silent and 50 of 863 blocks never materialized; a run that
materialized everything through the piece-4 path reproduced the exact
uncompressed logdet, which cleared the codec and the store; and a dump of the
leftovers said all 50 were **Schur** blocks, with a pointer in their slot.  A
box's Schur block has several writers that reinstall it - ensure_schur builds
it from the kernel when it finds no block, the wave replaces it with the
skeleton block, the exchange installs one - and each drops a packed one.
That is a class of writer, not a site, so the adopted Schur blocks are
materialized in begin() instead: 64 of those 863 blocks, about 7% of the
bytes, so the near blocks keep the saving and no writer has a null pointer to
misread.

### Piece 4 as built: the near blocks lazily, the Schur blocks not

begin() keeps schur_packed / edges_packed / strides / packed_owned in the
eliminator, and a block goes back in the wave that reads it
(materialize_adopted, three passes and one run_fill_batch, as the piece 2
restore did).  Correct at 96^3 for Color and CA0: logdet
-9928907.2811287567 and -9928906.4595639911, both the pieces 1-3 values,
acc_forward 9.33e-03, acc_backward 9.62e-16, and 0 blocks left packed.
863 blocks and 0.77 GB go back in 12 batches in 0.007 s.

**What it is worth, at 192^3 over 8 ranks (Color, compress 1):**

| level | peak, piece 3 | peak, piece 4 | adopted blocks | materialized |
|---|---|---|---|---|
| 5 (leaf) | 12.78 GB | **12.78 GB** | - (adopts nothing) | - |
| 4 | 7.74 GB | **6.24 GB** | 6.29 -> 3.85 GB | 6.80 GB, 12 batches, 0.024 s |
| 3 | 3.52 GB | 3.62 GB | 1.76 -> 1.60 GB | 1.94 GB, 0.009 s |
| 2 | 2.25 GB | 2.31 GB | 0.54 GB | 0.57 GB, 0.002 s |

logdet -89253766.641603723, the piece 3 value exactly; acc_backward
1.35e-15; `Factor time` 8.90 s against piece 3's 8.84 s, so the
decompression is free at this size.

**1.50 GB off level 4, and the global peak does not move** - it is the
leaf's, and 9.47 GB of its 12.78 GB is the owner pass's near blocks.

**The projection of "level 4 to about 1.5 GB" was wrong.**  The adopted set
goes to 3.85 GB resident, not to the 0.61 GB per-wave read union: that
measurement covered the *wave gathers* only, while the Schur blocks are now
eager, the remote-generator pass materializes the union of all its
generators at once, and a block materialized by a wave that then does not
read it is never retired.  Closing that gap is possible but it would be
work spent below the peak.

**The near blocks are lazy; the Schur blocks are materialized in begin().**
A near block has one reader and one writer: the wave gathers from it and
replaces it.  A box's *Schur* slot has several writers that each treat a null
pointer as "there is no block here" and fill a fresh one from the kernel -
ensure_schur, the wave's post-elimination landing, the exchange's install -
so a packed Schur block is dropped and its Schur-complement contents are
replaced by kernel values.  That was logdet out by 1.2e-2 with no crash and
no guard: 50 of 863 blocks at level 3, every one a Schur block of an
ordinary local writable box.  Schur blocks are 64 of those 863 and about 7%
of the bytes, so making them eager keeps the saving.

**What the readers cost to find: five runs, and the guards earned their
keep.**  Outside the wave materialization a null pointer already means "no
block", and the reader answers from the kernel or makes a fresh one, so a
packed block there is a wrong answer and not an error.  Every site that acts
on a null pointer therefore calls assert_not_packed() and names itself.  It
found the owner pass's join and replace_edge() in the remote-generator pass -
both needed a batch of their own, since run_fill_batch() plans and
synchronizes twice per call and one call per block is the per-call-overhead
shape that has cost this work three times.  What it did *not* find was the
Schur lifecycle, because ensure_schur was the one null-pointer site left
unguarded; and one gate (`if (!packed_edges_.empty())`) switched the
generator pre-pass off once the edges were gone, stranding the Schur blocks
it was also responsible for.

Two process notes, both of which cost a run:

- *a bisect beats a hypothesis.*  Materializing everything in begin() through
  the piece 4 path gave the right logdet and proved the packed store and the
  decompression were sound, which moved the search to the consumers;
- *check the build before reading the number.*  One gate run used a library
  linked two seconds before the source it was meant to test, and its
  unchanged result was read as "the fix does nothing".

### The piece 2 bug, and the two it was hiding

Two attempts ended in

    what():  H2 GPU heap: free of an unknown pointer

during the transition itself.  What was tried: collecting the chunk's blocks
and freeing them right after the gather that packs them, first through the
DeviceMatrix* the Pair carries, then by key through out->schur / out->edges
with a null check so a repeat is a no-op.  Both still hit it, so the pointer is
either not heap-owned or has already been freed by another path before the
chunk ends.

**It is a double free, and four hypothesis-bisects have not found it.**  The
heap now says which kind of bad pointer it got:

    free of an unknown pointer 0x... (offset 420509696 of 33398682368);
    inside the arena; nearest block below starts at 418640128, is 1713408
    bytes, so the pointer is past it

Inside the arena and in free space, so a block is being freed twice - not an
interior pointer into a group allocation, and not a foreign one.  The offset
is the same every run.  What has been ruled out:

- *the free in build_parent is clean.*  With a DeviceHeap::is_allocated()
  query it reports, at the Color leaf, `raw freed 7099, null 0, missing 0,
  stale 0 (foreign 0)` - every block still held, freed once, and the entry
  nulled afterwards;
- *it is not the packed_owned free in begin().*  Disabling it (leaking them
  instead) gives the identical error at the identical offset;
- *it is not the restore in begin() either.*  The run that looked like a
  restore fault was free-without-restore, where the next level reads blocks
  that are gone - the illegal access in the fill transposes was that, not a
  bad restore.

**The heap tag ended it in one run, and the double free was a mask.**  Each
of the four bisects above cost a 25 minute run to rule out one guess.  A
twenty line tag in DeviceHeap behind BPACK_CHECK=heap - offset to the last
place that allocated and freed it, resolved with dladdr - named all three
sites at once:

    H2 GPU heap: double free at offset 420509696
      allocated by  LevelEliminator::begin()+15751
      last freed by LevelEliminator::release_level_data()+2085
      freed again by LevelEliminator::release_level_data()+3633

objdump on those offsets (no line info in the build, but the loop shapes are
unmistakable) says what the three are:

- *begin()+15751* is the restore's own alloc_resident: `mov $0x5,%edx` is
  ResidentKind::Transition and the size comes from a PackedBlock's rows and
  cols;
- *release_level_data()+2085* is a stride-0x120 loop over BoxState freeing
  +0x28 and zeroing 16 bytes there - `heap_.free(st.schur.ptr); st.schur =
  DeviceMatrix{};`
- *release_level_data()+3633* is the `if (adopt_)` tail: two unordered_map
  walks (`mov 0x10(%rbx)`, then `mov 0x48(%rbx)`, nodes linked at +0x0)
  followed by nulling adopt_, i.e. `for (auto& kv : adopt_->schur)
  heap_.free(...)`.

So one block had two owners, `st.schur` and `adopt_->schur`.  begin() creates
that alias deliberately - `st.schur = it->second` copies the pointer and
leaves the map's copy - and unmakes it four lines later with
`adopt_.reset()`.  **The tail can therefore only run if begin() threw
between the two.**  make_level_eliminator holds the eliminator in a
unique_ptr and calls begin() itself, so a throw there destroys the object,
the destructor calls release_level_data(), the second free throws *while
unwinding*, and std::terminate prints that message in place of the real one.

Which is why four bisects found nothing: there is no double free to find.
The bug is whatever begin() threw, and it has never been printed.

Three changes, all of them right independently of what that error turns out
to be:

1. **one owner per adopted block** - `it->second.ptr = nullptr` after
   `st.schur = it->second`.  The tail still frees what an owner level did
   not adopt, which is what it is for, and a throw in begin() now reports
   itself;
2. **the restore reads a raw slot, it never creates one** - find() instead of
   operator[], throwing a message that names the box or the box pair.
   operator[] would insert an empty 0x0 block for a key the transition did
   not write, which surfaces much later as a mis-shaped block;
3. **edges-only levels are restored too** - the guard was
   `!schur_packed.empty()`, so a level that packed near blocks but no Schur
   block would have had its raw edges freed and never put back.

The heap's double-free message now also carries the block's resident kind and
size, so a recurrence does not need another objdump.

#### Step 2 piece 4: the adopting level reads the packed form (2026-10-06)

In: begin() keeps the packed store (packed_schur_, packed_edges_,
packed_strides_, packed_owned_) instead of restoring it, and each wave
materializes what it is about to read, in one batch, through
materialize_adopted().  The wave then retires the block as it always did.

Correct at 96^3 over 8 ranks and at 192^3, Color and CA0, with logdet
identical to pieces 1-3 to the last digit and 0 blocks left packed at every
level:

| run | logdet | acc_backward | Factor |
|---|---|---|---|
| 96^3 Color | -9928907.2811287567 | 9.62e-16 | 1.44 s (1.29 pieces 1-3, 1.16 off) |
| 96^3 CA0 | -9928906.4595639911 | 7.29e-16 | 1.47 s |
| 192^3 Color | -89253766.641603738 | - | 8.63 s (8.59 pieces 1-3, 7.71 off) |

Materializing costs 0.026 s for 7099 blocks in 12 batches at the 192^3 leaf,
so the codec is not the price.  **The memory is less than projected:**

| level | pieces 1-3 | piece 4 | adopted, raw -> materialized + packed |
|---|---|---|---|
| 5 (leaf) | 12.78 GB | 12.78 GB | adopts nothing |
| 4 | 7.74 GB | **6.24 GB** | 6.29 -> 3.85 + 0.94 |
| 3 | 3.52 GB | 3.62 GB | 1.76 -> 1.60 + 0.25 |
| 2 | 2.25 GB | 2.31 GB | 0.54 -> 0.54 + 0.07 |

The projection assumed the resident set would be the per-wave read union
(0.61 GB at level 4).  It is not: **a materialized block is only retired when
the wave that read it replaces it**, and the ones no wave replaces stay, so
the level accumulates materialized blocks as it goes.  Add the packed store,
which must stay resident for the whole level, and levels 3 and 2 come out
marginally *worse* than holding the raw blocks.  Level 4 still gains 1.5 GB
because its adopted set is large relative to its packed form.

#### What it cost to get right, and the lesson

Five gate runs and four wrong fixes.  Outside the wave materialization a null
ptr already means "there is no block here", and every reader answers that from
the kernel or by making a fresh one - so a block left packed is not a crash,
it is the original A in place of the updated block.  The guards
(assert_not_packed, naming the site) caught three readers the wave hook
missed: the owner pass's join, replace_edge in the remote-generator pass, and
ensure_edge's inner loop.  Two more were silent:

- build_parent reads a *child's Schur* block for a diagonal pair and falls
  back to eval() when the pointer is null - guarded now, and the transition
  materializes whatever the waves left;
- and the one that mattered: **a box's Schur block has several writers that
  reinstall it** (ensure_schur from the kernel, the wave's skeleton block
  after elimination, the exchange), each treating a null pointer as "no
  block".  50 of 863 blocks at level 3 of 96^3 were dropped this way, every
  one a Schur block, for a 1.2e-2 error in logdet.  Chasing the writers one
  run at a time was not converging, so the Schur blocks are now materialized
  in begin() - 64 of those 863 blocks, about 7% of the bytes, so the near
  blocks keep the saving and no writer has a null pointer to misread.

What found it was not reasoning about the readers.  It was (a) a counter for
blocks *left packed* at level end, which said 50 immediately, (b) a bisect
that materialized everything through the piece-4 path and reproduced the
uncompressed logdet exactly, clearing the codec and the store, and (c) a dump
of what the leftovers were, which said "all Schur" in one run.  Each of the
four fixes before those three diagnostics was a guess, and each cost a run.

One false alarm worth recording: a run at 4.94 s against 1.29 s looked like a
3.8x regression and was BPACK_CHECK=fillsrc, left on to see the dump - it
turns on the fill-source read-back.  The tell was that the leaf slowed down
too, and the leaf adopts nothing.

### A gap found while reading, not yet the failure

device_reduction.hpp moves parent blocks to their new owner's rank when the
parent level runs on fewer ranks, and it works on the raw blocks only:

- it sends from `blocks->schur.at(pm).ptr` and `blocks->edges.at(key).ptr`,
  which piece 2 has freed;
- it *erases* the entries it hands away, which leaves a schur_packed or
  edges_packed entry with no raw slot;
- it *emplaces* the blocks it receives, which have no packed form (harmless -
  they stay raw);
- and its drop path frees every block of a Blocks it discards, which with
  piece 2 leaks that Blocks' packed_owned.

So a transition whose parents are reduced to other ranks must not pack, or
must unpack first.  Skipping the pack there costs nothing worth having: the
reduction happens only where the rank count shrinks, which is the coarse
levels - level 2 of this test packs 8 blocks and 0.02 GB, against 863 blocks
and 0.77 GB at the leaf.  **This is not the current failure**: the run that
fails reports `0 rounds, sent 0.00 GB in 0 blocks` at level 2, so no block
was moved, and the rank that aborts is rank 0, the one that goes on to
level 1 and runs the restore.

### The older note, kept because its measurement still stands

**The diagnostic has now been run, and the free is not the problem.**  With a
DeviceHeap::is_allocated() query (owns() only says the address is inside the
arena, which does not tell a foreign pointer from a freed one) the free loop
reports, at the Color leaf:

    raw freed 7099, null 0, missing 0, stale 0 (foreign 0)

Every block the chunk packed was still held by the heap and was freed exactly
once.  So looking up by key and nulling the entry is correct, and the earlier
"free of an unknown pointer" came from somewhere downstream of it.

With the free clean the failure moves to the **next level**, in the
fill-source sub-group transposes:

    CUDA error in fill transposes: an illegal memory access was encountered

at level 4, after level 5's transition reported normally (6 chunks, blocks
kept on the device).  An illegal access there is a bad pointer in a
TransposeItem, not bad data, so the suspect is the restore in begin() that
puts the packed blocks back - it allocates each block and decompresses into
it, and alloc_resident can call the reclaimer part way through, which frees
released fill sources while h_src/h_dst already hold addresses.

Next, in order:

- test the halves apart: free without restoring (the next level then has no
  blocks, so only the transition is exercised) and restore without freeing
  (force it past the null-pointer guard), to say which side faults;
- if it is the restore, allocate every block first, then build h_src/h_dst,
  then decompress - never interleave allocation with address collection;
- check whether the between-level compaction runs in these configurations at
  all.  It relocates blocks and rewrites the pointers it knows about;
  packed_owned and PackedBlock::packed are not in its list, so if it does run
  they are stale afterwards.  Nothing in the Color runs so far has proved it
  fires here.

Older hypotheses, now unlikely given the above:

- a parent block's pointer is replaced during the transition (out->edges[key]
  assigned twice for one key), so the chunk's record and the map disagree;
- something else already frees it - the retired list in build_parent is the
  candidate, and the Pair's block may not always be a fresh allocation;
- the free is correct but runs while the gather still reads it, and the heap
  then sees the second free from the normal teardown.

The diagnostic that settles it is cheap and was not run: at the free, print the
pointer, whether DeviceHeap::owns() accepts it, and its key, and compare
against the allocation log for that key.  Worth doing before another attempt -
two hypotheses have already been wrong.

Note the peak this is aimed at: the leaf's 9.47 GB of parent blocks is resident
*during the transition*, which is where the leaf's peak is, so freeing has to
happen per chunk to move it.  Packing at the end of build_parent would not.

### After step 2

The transient (7.02 GB) becomes the largest single term, at about 60% of the
leaf's peak.  It is the wave's per-box work blocks, all live at once
(level_eliminator.hpp, `work_blocks`), so it is bounded by tiling rather than
compression - but tiling it splits step 5's accumulations across sub-batches
and reorders sums into shared near blocks, which the CA replica check relies
on being identical.  That question has to be settled before the tiling, not
after.

Two small remainders with known routes: X_RR_full (0.73 GB) by reconstructing
it from its own LU - exact, no tolerance - once C4 lets the leaf keep its
factors, and the packed temp2^T (1.38 GB) only by FP32 storage, which the user
has deferred.

### Scaling, honestly

Tiling makes the working memory O(1) in problem size; the compressed stores
stay proportional to it.  So 8x more boxes a rank looks reachable and 64x does
not - at 64x the compressed data alone exceeds the card, and that needs host
offload with streaming, which is a different mechanism.

## Backlog from the factor compression (H2_compress_factors)

X_NR is 92% of a kept block and compresses 16x at tol_comp/5 (measured at
192^3 on 8 ranks), taking the device's factors from 27.76 to 2.96 GB per
rank. Three things are left:

- C1 (done, 2026-10-05). The solve was 4.7x slower with compression on
  (0.264 -> 1.239 s at 192^3 on 8 ranks, 40 GB nodes), past the 4x the user
  accepted. Overlapping the decompression with the sweeps brought it to
  0.954 s, **3.6x**, with the same results to the last digit (logdet
  -89253588.194065422 either way, acc_backward 1.3e-15).

  The first diagnosis in this entry was wrong and is kept here as a warning.
  It claimed the codec was not the cost and blamed per-call overhead and the
  box-table upload. Timing the calls on the stream (CUDA events, now in
  DeviceSolveTimes::xnr, printed by BPACK_TRACE=phase under "GPU solve") says
  otherwise:

  | | off | on, per wave | on, overlapped |
  |---|---|---|---|
  | GPU solve | 0.257 s | 1.231 s | 0.903 s |
  | forward / backward | 0.105 / 0.134 | 0.592 / 0.619 | 0.430 / 0.454 |
  | MPI | 0.069 | 0.393 | 0.356 |
  | X_NR decompression | - | 0.666 s, 19.3 GB/s | 0.669 s, 19.2 GB/s |

  The decompression is 0.666 s of the 0.973 s added, and the MPI inflation
  (+0.324 s, ranks waiting for peers that are still decompressing) is the
  rest - the whole slowdown, accounted for. **bitcomp's lossy decompression
  runs at 19 GB/s here**, not the 400-600 GB/s an isolated benchmark
  suggested, and that is the floor: the 80 calls of a solve average 160 MB
  each, so batching waves into bigger windows had nothing to gain. (Small
  batches do cost: CA0, whose waves are ~52 MB, gets 5.4 GB/s.)

  What was done instead: each level's scratch holds kXnrBuffers = 2 waves and
  a wave decompresses on a stream of its own while the previous wave's
  kernels and messages run (prepare_wave_xnr, which now takes the sweep's
  direction so it knows which wave comes next; events gate each buffer
  against its last reader). That hid half the decompression. The second
  buffer is taken only while the heap keeps kXnrSecondBufferHeadroom free
  afterwards, so a tight level falls back to one buffer and the old strictly
  per-wave behaviour - without that, a level could gain +0.7 GB (947 MiB at
  level 4 of this case, against 473 MiB) where there is no room for it.
  The box-table patch kernel (launch_patch_xnr) also went in; it removes
  123 MB of H2D traffic per sweep but measured as nothing (1.2384 ->
  1.2380 s), which is the other half of the correction above.

  **Not irreducible after all (2026-10-05, same day).** The 19 GB/s was the
  block count, not the codec. Lossy Bitcomp decompresses a batch in parallel
  across its blocks, and a wave was one block per box: 27-37 of them. A sweep
  over the dumped X_NR of a real level-4 wave (37 boxes, 436 MB) says what
  that costs, at the tolerance the H2 path uses:

  | chunk | blocks | ratio | compress GB/s | decompress GB/s | staging MB |
  |---|---|---|---|---|---|
  | whole box | 37 | 19.04 | 60.5 | 63.0 | 436.3 |
  | 4 MB | 127 | 19.04 | 169.2 | 189.5 | 436.9 |
  | 1 MB | 456 | 19.02 | 362.0 | 447.0 | 439.5 |
  | 0.5 MB | 891 | 19.00 | 411.2 | 548.1 | 443.3 |
  | 0.25 MB | 1759 | 18.97 | 462.3 | 612.4 | 443.6 |
  | 0.125 MB | 3503 | 18.90 | 468.3 | 641.1 | 457.8 |

  The ratio barely moves and max_err does not move at all, so the fix is to
  cut each box's X_NR into chunks (kXnrChunkBytes, 0.5 MB) and let the batch
  be chunks rather than boxes - the same bytes, the same scratch, just more
  blocks. A chunk is a contiguous range of the box's elements, so the buffer
  comes back as its chunks laid end to end; the delta stays per box, so the
  tolerance a box is held to does not depend on how it was cut. Each chunk
  gets a 256-byte start in the kept block (bitcomp reads its compressed bytes
  with wide loads and faults on a 4-byte start - that was the one bug), and
  the padded lengths go in KeptSolveLevel::chunk_strides, which is how the
  solve walks to a chunk.

  Measured at 192^3 on 8 ranks, 40 GB nodes:

  | | off | per box | in 0.5 MB chunks |
  |---|---|---|---|
  | X_NR decompression | - | 0.669 s, 19.2 GB/s | **0.035 s, 365 GB/s** |
  | time_solve | 0.264 s | 0.954 s (3.6x) | **0.284 s (1.07x)** |
  | GPU solve | 0.257 s | 0.903 s | 0.273 s |
  | of which MPI | 0.069 s | 0.356 s | 0.074 s |
  | Factor time | 7.68 s | ~8.2 s | **7.38 s** |
  | blocks per level-4 wave | - | 27 / 37 / 37 | 153 / 625 / 965 |
  | ratio, levels 4/3/2 | - | 18.7 / 19.5 / 17.6 | 18.6 / 19.4 / 17.5 |
  | acc_backward | 1.20e-15 | 1.30e-15 | 1.30e-15 |

  So compression now costs 7% of the solve instead of 260%, and the
  factorization is *faster* with it on than off - the gather copies 0.97 GB
  into the kept blocks instead of 7.03 GB, which more than pays for the
  compression. CA0 gains the most in rate (its waves were 52 MB): 5.4 ->
  168 GB/s, time_solve 1.450 -> 0.946 s, the remainder being its host CA
  level. That also settles C3 (no measurable per-wave cost left) and takes
  away the speed argument for FP32 X_NR; FP32 remains a memory argument.

  The two-buffer overlap below is now mostly redundant - 0.035 s has little
  left to hide - but it costs nothing when the second buffer fits and falls
  back cleanly when it does not, so it stays.
- C2. CA3 keeps 5.71 GB of 19.90 GB (3.5x) where Color and CA0 reach
  7.4-7.7x, because land_full_factors - the component-owner path that
  stores the factors of boxes landing from other ranks - lays them out
  uncompressed. Its per-box composition is the same as keep_solve_factors'
  (92% X_NR, 16x), so compressing it should bring CA3 in line. Its boxes
  arrive one at a time rather than in waves, so they need collecting into
  a batch before compressing.
- C3 (resolved by the chunking, 2026-10-05). keep_solve_factors synchronizes
  the stream once per wave to keep the compressor's staging alive through the
  gather, which cost about 10% of the factorization (7.49 -> 8.22 s at
  tol_comp). With the compression at 400+ GB/s instead of 60 the sync has
  almost nothing to wait for: the factorization is now 7.38 s with
  compression against 7.68 s without. No double buffering needed.
- C4 (found while measuring C1, 2026-10-05). **The leaf level never keeps
  its solve factors, so the compression misses 73% of X_NR.** At 192^3 on 8
  ranks (40 GB nodes) the leaf's X_NR is 17.45 GB per rank of the 23.87 GB
  in the whole factorization; levels 4, 3 and 2 hold 5.13, 1.15 and 0.14 GB
  and compress 18.7x, 19.5x and 17.6x, taking the kept factors from 7.03 to
  0.97 GB. The leaf's go to the host and come back **uncompressed** at the
  first solve (19.81 GB per rank, 2.15 s of it host packing), so the solve
  holds about 20.8 GB of factors on the device whether compression is on or
  not. BPACK_TRACE=wave prints this per level ("X_NR read in place").

  decide_solve_keep rejects the leaf because solve_keep_bound is both
  uncompressed and loose - it uses each box's n where the factors have r
  columns (levels 4-2 come in at 53-68% of their bound), and at the leaf it
  is tens of GB against a budget of 0.8 x capacity less what the level's own
  working set already holds (26.18 GB resident of a 33.40 GB heap).

  Dividing the bound's X_NR term by 4 when compression is on does flip the
  decision, and that was tried and reverted: the leaf then fails the keep a
  second way ("no free range of the heap held a wave's factors" -
  largest_free < one wave's block), and **the rescue path after that failure
  is broken** - level 4 then reads NaN/Inf in its ID input on every rank.
  rescue_skipped_xnr(true) does not restore what skip_xnr_ told the waves not
  to copy to the host. That is a latent bug in the existing keep-failure
  path, now reachable; it needs fixing before any attempt to keep the leaf.

  **Why it does not fit, measured** (BPACK_TRACE=wave now prints this at the
  failure). During the leaf's waves the heap holds 31.84 of 33.40 GB, so
  1.56 GB is free, of which 1.43 GB is in one range - the heap is full, not
  fragmented, and the free space is the leaf's own working set's leftovers
  (26.18 GB resident + 7.03 GB transient at the peak). Then:

  1. A big leaf wave's factors are 2.70 GB uncompressed.
  2. WaveXnrCompressor::compress stages **the whole wave at
     bitcompMaxBuflen** (>= the input size), so it asks for 2.70 GB and
     try_alloc fails against 1.56 GB free. The wave falls back to plain.
  3. The plain 2.70 GB block then fails the contiguity check against
     1.43 GB, and the level drops everything it had kept (0.87 GB).

  So the binding constraint is the compressor's own staging, not the stored
  factors: compressed, a leaf wave is about 145 MB and would fit with room
  to spare.

  **Note what the C1 chunking did and did not do.** Cutting a box's X_NR
  into kXnrChunkBytes blocks gave the batch its parallelism, but the staging
  is still the whole wave - it is the sum of bitcompMaxBuflen over the
  chunks, and that bound is at least the input size whatever the chunk size
  (the sweep measured 443 MB of staging for 436 MB of data at 0.5 MB
  chunks). So C4's blocker is exactly as it was. What the chunking did change
  is that the chunks now exist, so bounding the staging is a matter of
  grouping them rather than of cutting anything new.

  Two ways out, the first much cheaper:

  - **Compress a wave in groups bounded by a fixed staging buffer.** A group
    is a run of the wave's chunk list whose bitcompMaxBuflen sum fits the
    buffer; the staging becomes O(buffer) instead of O(wave). 512 MB fits in
    the 1.56 GB free today and holds about 1024 chunks of 0.5 MB, which is
    past the point where the throughput curve flattens, so it costs nothing
    in rate. The work is in keep_solve_factors: it compresses the whole wave
    and then gathers from the staging with one launch_copy_spans, so each
    group has to be gathered out before the buffer is reused.
  - **Free working set first** (step 4, compressing the fill sources). The
    leaf's 26.18 GB resident is what leaves only 1.56 GB, so step 4 and C4
    are the same constraint seen from two sides.

  **Tried, 2026-10-05: the first two are done, the third does not work yet.**

  1. *The rescue bug is fixed.* rescue_skipped_xnr read `ntot * r` elements
     from a kept X_NR that holds `xnr_bytes` of compressed chunks, running
     off the end of the block and handing the host nonsense - the same bug
     that materialize_host_factors already had fixed in the solve, missed in
     this second site. It now brings the compressed bytes over and
     decompresses them, so a mid-level keep failure no longer poisons the
     next level. Verified: with the leaf attempting a keep, Color runs to
     completion where it used to give NaN/Inf in level 4's ID input.
  2. *The staging is bounded.* keep_solve_factors now keeps a wave in groups
     of boxes whose X_NR comes to at most kXnrStagingBudget (512 MB), each
     group compressed and gathered into a block of its own
     (keep_solve_factor_group). At the leaf the request that used to be
     2.70 GB becomes 559 MB. At 192^3 the groups change nothing - a level-4
     wave is at most 473 MiB, so it is one group - which is the point: it
     costs nothing now and is what C4 needs later. Color and CA0 both verified
     unchanged (time_solve 0.282 s, 365 GB/s; CA0 168 GB/s).
  3. *The compression-aware keep bound does not work yet, and is left off*
     (kXnrBoundDivisor = 1). Dividing the bound's X_NR term by 4 does flip
     the leaf's verdict, and the leaf then keeps 1.13 GB of its ~2.4 GB
     before failing again - this time on fragmentation, 2.92 GB free but
     7.7 MB in the largest range, because the compressor cannot stage its
     group and falls back to plain. It also segfaults CA0 at CA_level 4,
     where level 4 attempts a keep with the heap at 33.19 of 33.40 GB; that
     is a second unhandled mid-level keep-failure path and it is **not
     diagnosed**.

  The reason 3 cannot be fixed by choosing a better divisor: decide_solve_keep
  runs at level start and measures the heap as it is *then*, before the level
  has allocated the 26.18 GB it holds at its peak. For the leaf the budget is
  therefore over-optimistic however the bound is scaled - the level says yes,
  then runs out of room half way through its own waves. Making the leaf keep
  its factors needs the working set reduced (step 4) or a keep decision that
  estimates the level's own working set. The groups and the rescue fix are
  both prerequisites that are now in place.

  The cost of succeeding is now small. At the chunked rate (365 GB/s) the
  leaf's 17.45 GB over two sweeps is about 0.10 s, not the 1.8 s it would
  have been at 19 GB/s, so C4 is no longer a memory-for-time trade - it is
  close to free, and worth doing for the ~19 GB it frees at solve time.
