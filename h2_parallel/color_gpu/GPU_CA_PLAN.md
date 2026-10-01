# GPU CA (communication-avoiding) factorization: plan and milestones

Status: M0 done; M1, M2 and M4b implemented and validated, M4b with
per-level solve residency and the solve's halo in device memory (results
under each). Keep this list current as milestones finish, so deferred
improvements are not lost.

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

Files: no new folder. Device code goes into `color_gpu/` (a CA mode of the
level eliminator, a new file for the CA glue: ghost upload, CA schedule, CA
sketch semantics; a CA entry in the shared GPU driver); `color_CA/` and
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
3. Device ghost-block exchange (new, `color_gpu/ca_halo.hpp`), after the
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

**M2 as built (2026-09-29).** `color_gpu/ca_halo.hpp`, the transition's
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
`color_gpu/device_solve.hpp`: a CA level's `DeviceSolveLevel` holds the
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
  `refresh_CA_assisting_skeletons` (color_CA/serialization.hpp) sends only
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

### M5. CA3 (component-owner) on the GPU
Evaluate after M1-M3 numbers. Unique owner per boundary component, the
corner -> edge -> face DAG, FULL/COMPACT/SKELETON routing, deterministic
replay (R1/R2/R3). Risk: small component batches underfill the GPU.


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

## Backlog from the GPU Color work (not CA-specific)

- B1. Level-start wait: the first transport of each leaf level (0.4-0.7 s at
  64 ranks) waits for the slowest rank's level setup; trim the setup it
  waits on.
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
