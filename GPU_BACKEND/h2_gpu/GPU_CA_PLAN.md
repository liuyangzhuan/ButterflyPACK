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

## Backlog from the factor compression (H2_compress_factors)

X_NR is 92% of a kept block and compresses 16x at tol_comp/5 (measured at
192^3 on 8 ranks), taking the device's factors from 27.76 to 2.96 GB per
rank. Three things are left:

- C1. The solve is 5.3x slower with compression on (0.251 -> 1.337 s at
  192^3 on 8 ranks), past the 4x the user accepted. The codec is not the
  cost: about 24 GB of X_NR per sweep over six passes is 144 GB in 3.4 s,
  40 GB/s against the 400-600 GB/s a wave-sized batch reaches on its own.
  The cost is per-call: a wave is one bitcompBatchUncompress (about 300 per
  solve, and the coarse levels' waves are far below the batch size that
  fills the GPU - 9.8 GB/s at 8 blocks), and each wave re-uploads the whole
  SolveBox array (410 KB) to move its pointers. MPI time triples
  (0.069 -> 0.392 s), so the decompression is serializing against the
  exchanges rather than overlapping. Fix: decompress a window of waves in
  one batch, and patch only the wave's box entries with a small kernel
  instead of re-uploading the table.
- C2. CA3 keeps 5.71 GB of 19.90 GB (3.5x) where Color and CA0 reach
  7.4-7.7x, because land_full_factors - the component-owner path that
  stores the factors of boxes landing from other ranks - lays them out
  uncompressed. Its per-box composition is the same as keep_solve_factors'
  (92% X_NR, 16x), so compressing it should bring CA3 in line. Its boxes
  arrive one at a time rather than in waves, so they need collecting into
  a batch before compressing.
- C3. keep_solve_factors synchronizes the stream once per wave to keep the
  compressor's staging alive through the gather, which costs about 10% of
  the factorization (7.49 -> 8.22 s at tol_comp). Double buffer the staging
  across waves instead.
