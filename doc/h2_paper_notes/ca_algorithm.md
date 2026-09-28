# Communication-Avoiding and Component-Owner H2 Factorization

Status: paper-development notes describing the ButterflyPACK integration.

## 1. Candidate Paper Thesis

The communication-avoiding (CA) ordering groups process-boundary boxes by
geometric components rather than global parity.  The replicated CA baseline
exchanges immutable level-entry state once and then redundantly eliminates a
halo, trading extra work and memory for fewer dependent communication rounds.

The improved algorithm introduces a second point in this tradeoff:

- **mode 0, replicated CA**: one level-entry dependency, replicated boundary
  elimination;
- **mode 3, component-owner CA**: every boundary component is eliminated
  once by a selected owner, with an asynchronous corner-to-edge-to-face DAG;
- **Color**: every box is eliminated once by its fixed home owner, with
  `2^d` parity dependencies.

In 3D, component ownership replaces the eight-round parity chain with three
boundary-color dependencies, while eliminating replicated factorization.  It
pays for owned-component entry fetches, routed eliminated state, and runtime
bookkeeping.  The paper should present this as a work/latency/volume design
space, not as a universal ordering of the three methods.

## 2. Geometric CA Ordering

The process decomposition partitions a level into rectangular rank bricks.
Boundary boxes are classified by the process-grid feature they touch:

- **blue**, including the orange thickening: corner clusters;
- **purple** in 3D: edge runs;
- **green**: face patches in 3D and edge runs in 2D;
- **interior**: boxes whose update neighborhood is rank local.

The elimination order is

```text
blue(+orange) -> purple -> green -> interior        (3D)
blue(+orange) -> green -> interior                  (2D)
```

Boxes inside a component may be adjacent and therefore require internal
waves.  Different components of the same geometric color are independent.
That distinction is central: owning individual boxes does not remove CA
replication safely, while owning an entire connected component does.

The full boundary stratification requires sufficient local geometry.  The
integrated threshold is `7^3=343` boxes per active rank in 3D and `5^2=25` in
2D.  A requested CA level below this threshold automatically runs Color.
Single-rank levels are allowed, and all 1D levels use Color.

## 3. Mode 0: Replicated CA

### 3.1 Level-entry invariant

Before elimination, each active rank imports the entry state of the halo
needed to reproduce every boundary component touching its brick.  The entry
state contains:

- box geometry, point/index sets, and neighbor lists;
- the box Schur block;
- modified one-hop interaction blocks.

Strict two-hop blocks are not communicated.  At level entry they are either
kernel-fresh or, in lazy mode, recoverable from retained generators.  Fills
created during the level remain inside the replicated local computation and
are assembled into parent near interactions at transition.

After the gather, each rank independently eliminates its local boxes plus
the required ghost boundary boxes in the same CA order.  There is no
cross-rank dependency between boundary groups inside that level because each
rank has a sufficiently deep replica.

### 3.2 Computation cost

Let `H_{r,l}` be the ghost boxes eliminated by rank `r` at level `l`, and let
`C_B` be the ID plus elimination cost of box `B`.  Then

```math
F_{r,l}^{CA0}
 \simeq \sum_{B\in Home_{r,l}\cup H_{r,l}} C_B.
```

The global replication factor is

```math
\rho_l^{rep}
 = \frac{\sum_r |Home_{r,l}\cup H_{r,l}|}
        {|Boxes_l|}.
```

`rho_l^rep` approaches one when rank bricks are large compared with the halo
depth, but can approach two or more for small bricks.  Interior process-grid
ranks replicate more faces, edges, and corners than domain-corner ranks, so
replication also creates geometric load imbalance.

### 3.3 Memory cost

A homogeneous leading model is

```math
M_{r,l}^{CA0}
 \simeq M_{home,l}
 + \sigma\sum_{B\in H_{r,l}}
   \left(k_B^2+\sum_{A\in N_1(B)}n_Bn_A\right)
 + M_{ghost\ factors,l}+M_{scratch,l}.
```

Symmetric pair pruning reduces bytes on the wire, but the receiver restores
the reciprocal view, so the live replicated halo still contains the state
needed by its local elimination.  With lazy far updates, dense two-hop blocks
are absent and generator storage replaces them as described in the Color
note.

### 3.4 Communication cost

Ignoring metadata, the monolithic halo volume is

```math
V_l^{CA0}
 \simeq \sigma\sum_r\sum_{B\in RequestedGhosts_{r,l}}
 \left(k_B^2+\sum_{A\in N_1(B)}n_Bn_A\right).
```

This counts duplicates when one boundary component is needed by several
sharers.  Symmetric-pair pruning sends only one orientation of a pair when
possible and reconstructs the other by transpose.

The logical dependency distance is one level-entry exchange rather than
`2^d` elimination rounds.  A simple exposed-time model is

```math
T_{comm,l}^{CA0}\approx \chi_l\alpha+\beta V_l^{CA0}+T_{pack/install,l},
```

before overlap.  The staged protocol below uses more physical messages but
can hide their later stages behind earlier elimination.

## 4. Staged Halo Exchange

`H2_CA_staged_halo=2` divides the level-entry transfer into four priority
stages in 3D:

| Stage | Payload |
|---:|---|
| 0 | shells: geometry, points, index sets, neighbor lists; no matrices |
| 1 | blue/orange Schur and near blocks |
| 2 | purple Schur and near blocks |
| 3 | green Schur and near blocks |

All receives are posted up front.  Senders pack and initiate stages in this
order.  Elimination waits only for the stage required by the next group, so
purple and green data can move while blue and purple compute.  If lazy Schur
updates are disabled, the implementation waits for all block stages before
elimination, making staging primarily a transport reorganization.

Updates targeting an unarrived ghost are diverted into pending state and
merged when that stage is installed.  For symmetric matrices, staged pruning
is enabled automatically: one pair orientation is sent and the reciprocal
view is reconstructed only after its dependency arrives.

Staging does not fundamentally reduce bytes.  Its goals are:

1. prioritize small shell metadata and early colors;
2. overlap later payloads with useful elimination;
3. shorten the exposed communication tail;
4. parallelize packing and unpacking.

Its additional startup term is roughly `O(4 chi_l alpha)` in 3D, but only the
waits not hidden behind earlier work lie on the critical path.

## 5. Lazy Schur and Streamed ID on CA Levels

`H2_use_sketch=2` is supported on CA levels.  It streams the ID target and
applies lazy far contributions in sketch space exactly as described in
[color_algorithm.md](color_algorithm.md).

The CA path supports effective lazy modes:

- `H2_lazy_schur=0`: eager local far storage;
- `H2_lazy_schur=1`: no dense far storage; regenerate from generators;
- requested mode 2: automatically treated as mode 1.

Generated-near communication is a Color transport optimization and is not
implemented for CA.  CA's level-entry near blocks encode prior local
elimination history; the receiver cannot reconstruct that immutable entry
state from the kernel alone without receiving a larger child-level history.

For replicated CA, lazy far is mainly a memory optimization because two-hop
blocks were never part of the halo payload.  It changes local storage from

```math
\frac{\sigma}{2}(b_l+h_l)z_2n_l^2
```

to generator storage of approximate size

```math
\sum_{E\ local+ghost}\sigma
\left(r_E^2+r_Ek_E+r_E\sum_{A\in N_1(E)}n_A\right).
```

## 6. Why Component Ownership Is Necessary

A CA boundary component can straddle several ranks and contains adjacent
boxes.  Assigning its boxes independently would allow two owners to eliminate
neighboring boxes without seeing each other's updates.  Replicated CA avoids
this by giving every sharer a complete component copy.

Mode 3 instead assigns the entire connected component to one eligible rank.
The owner fetches remote entry boxes of that component, eliminates all its
boxes, and publishes the resulting state.  Consequently:

- every boundary box is factored once;
- dependencies inside a component are owner local;
- only dependencies between geometric colors cross owners;
- ownership can be selected to balance work among component sharers.

## 7. Component DAG and Owner Assignment

### 7.1 Static graph

The graph is a pure function of the process grid and CA classification:

```text
corner components -> edge components -> face components -> rank interiors
```

In 2D, green edge components depend on their endpoint corner components.  In
3D, a purple edge depends on its corner predecessors, and a green face
depends on bounding edge/corner predecessors.  Each rank interior has a join
dependency on the boundary components touching its brick.

The graph includes:

- component boxes and geometric color;
- sharer/eligible-owner ranks;
- predecessor and successor lists;
- recipient ranks at one- and two-hop distance;
- per-rank interior join nodes.

Every rank builds the same graph without MPI.  An assignment hash is checked
so divergent geometry fails before communication can deadlock.

### 7.2 Assignment

The implemented owner policy is deterministic largest-processing-time (LPT)
assignment over eligible sharers, followed by deterministic relief and
smoothing passes.  Coarse levels balance each color separately; fine levels
balance cumulative ownership across colors so a rank light in an early color
can own more later work.

The present cost is primarily box count, not measured factorization cost.
This is an explicit limitation: two boxes with different numerical ranks can
have substantially different costs.

## 8. Mode-3 Dataflow Execution

Each rank maintains counters for owned components, expected payloads, and its
interior join.  A component becomes runnable when all required lower-color
state and auxiliary two-hop skeletons are installed.

The high-level event loop is:

1. post entry-state and component-payload receives;
2. eliminate ready owned components;
3. serialize and push completed component state;
4. install arrived shared components;
5. replay deferred updates in canonical order;
6. poll communication at wave boundaries;
7. start interior work after the local boundary join is satisfied.

The runtime does not require a global barrier between components.  A late
component delays only its graph successors.  The critical path is the longest
weighted corner-to-edge-to-face-to-interior chain, not the sum of the slowest
rank in each global phase.

### 8.1 Batching and component-aware waves

Running one small component at a time would starve an OpenMP team.  Runnable
components of one color are merged into a batch.  Inside each component, a
greedy graph coloring produces the fewest practical nonadjacent waves:

- a face sheet typically needs about four classes;
- an edge run typically needs two;
- the general 3D upper target is eight.

At fine levels, batches are arrival driven and can be preempted at a wave
boundary by newly runnable lower-color work.  At a coarse level, currently
identified by fewer than eight local boundary boxes per OpenMP worker, the
runtime uses two maximally full batches per color: owned components, then
installed shared copies.

### 8.2 Serial oracle

`H2_CA_owner_serial=1` executes the same component graph in a serialized
schedule.  It is a correctness oracle for the asynchronous runtime, not a
performance mode.  Fast and serial mode should be bitwise identical for the
same graph and numerical chunking.

## 9. Routing and Current Payload Policy

Mode 3 routes each eliminated box according to what the recipient will read:

| Kind | Recipient | Leading contents |
|---|---|---|
| `FULL` | the box's home rank | factors and all state needed by transition and solve |
| `COMPACT` | ranks with a one-hop pair of interest | generator/factor state and required near blocks |
| `SKELETON` | ranks with a two-hop interest only | index sets needed by lazy far regeneration |

The home rank must receive `FULL` state when a remote owner eliminates its
box because parent construction and later solves are organized by home
ownership.

### 9.1 Current COMPACT tradeoff

The current implementation sends full `X_NR` in COMPACT state and recomputes
`temp2=-X_NR X_RR^{-1}` on the receiver.  It also sends partner-modified near
blocks needed by that recipient; kernel-only blocks can be rebuilt locally.

An earlier row-diet design regenerated many `X_NR` rows from kernel or partner
copies and reduced bytes substantially.  It also made installs depend on the
receiver's pair copies being at the exact source-entry wave.  The resulting
wave-by-wave install rounds cost more at 64 ranks than the saved bytes.
Shipping full `X_NR` permits every copy in a batch to install in one parallel
round, after which its update passes run wave by wave.

This is a useful paper result: minimum bytes did not minimize time because it
increased latency and reduced batch width.

### 9.2 Need-based routing

A completed box is not broadcast globally:

- `FULL` is sent only to home when necessary;
- `COMPACT` is sent only to ranks holding relevant one-hop pairs;
- `SKELETON` is sent only to ranks needing a canonical two-hop view;
- recipient-specific near blocks are pruned.

Send/receive buffers come from a capacity-retaining pool and serialization
writes directly into final offsets.  This avoids product-buffer copies and
repeated allocation/page-fault churn.

## 10. Deterministic Deferred Replay

Arrival order cannot determine floating-point accumulation order.  Mode 3
therefore separates forming an update from committing it to a shared pair.
The runtime uses three replay regions:

- **R1**: sources in the active batch update pairs having an endpoint in that
  batch, in wave order;
- **R2**: before a batch starts, available earlier sources are applied to its
  endpoint pairs above each pair's watermark;
- **R3**: before interior starts, remaining sources above pair watermarks are
  applied to all pairs of interest.

The canonical source key is `(color,wave,Morton)`.  Pair watermarks prevent a
source from being applied twice when two component copies become live at
different times.

Determinism also depends on fixed numerical calls:

- 512-row right-side solve chunks;
- 256-output-column GEMM chunks;
- canonical source lists for lazy far regeneration;
- fixed per-target slots for solve contributions.

Thus message timing changes readiness and overlap, but not arithmetic order.

## 11. Mode-3 Computation Cost

Ignoring runtime overhead, unique component ownership gives

```math
F_l^{CA3,elim}\simeq \sum_{B\in Boxes_l}C_B,
```

instead of the replicated sum over local plus halo copies.  Additional work
is

```math
F_l^{CA3,extra}
 = F_{install}+F_{near\ rebuild}+F_{replay}
 + F_{pack/unpack}+F_{scheduling}.
```

Some receiver GEMMs are relocated work rather than new work: mode 0 performs
them while eliminating replicas, whereas mode 3 rebuilds the receiving copy
from the unique owner's state.  The genuine overheads are serialization,
copy installation, graph scheduling, canonical replay bookkeeping, and any
loss of BLAS/OpenMP efficiency in thinner component waves.

Mode 3 is favored when replicated box work dominates those overheads.  Mode 0
can remain competitive at fine levels where halo replication is a small
fraction and broad replicated waves use the node efficiently.

## 12. Mode-3 Communication Cost

For component `C` with owner `o(C)`, define:

- `Fetch(C)`: entry state of boxes in `C` not already held by `o(C)`;
- `Full(E)`: home payload for source box `E` when owner differs from home;
- `Compact(E,r)`: one-hop payload required by recipient `r`;
- `Skel(E,r)`: two-hop index payload required by recipient `r`;
- `Rec_1(E)` and `Rec_2(E)`: the corresponding one- and two-hop recipient
  rank sets.

Then

```math
V_l^{CA3}
 \simeq \sum_C |Fetch(C)|
 + \sum_E\left[
     I(o(E)\ne home(E))|Full(E)|
     +\sum_{r\in Rec_1(E)}|Compact(E,r)|
     +\sum_{r\in Rec_2(E)}|Skel(E,r)|
   \right].
```

This expression is preferable to a single asymptotic bound because component
sharer counts and recipient overlap drive the real volume.  Homogeneously,
the dense part of a COMPACT payload is dominated by

```math
\Theta\left(\sigma[z_1n_lr_l+r_l^2+\text{retained near blocks}]
\right).
```

Compared with replicated CA, owned-entry fetch removes duplicate elimination
state but adds routed post-elimination state.  Compared with Color mode 2,
mode 3 generally moves more bytes because it must fetch component entry state
and restore home factors.  Historical measurements suggest roughly
`1.5-2x` Color bytes as a useful hypothesis, not a universal bound.

The latency structure is different from both baselines.  In 3D there are
three boundary-color dependencies after entry fetch rather than eight parity
dependencies.  Because execution is asynchronous, these are DAG depths, not
global barriers.  A schematic critical-path model is

```math
T_l^{CA3}
 \approx \max_{\pi\in DAG_l}
 \sum_{v\in\pi} C_v
 + \sum_{e\in\pi}(\alpha+\beta V_e),
```

plus the interior and transition.  Counting the initial entry fetch gives
approximately `d+1` dependent communication stages; the exposed cost can be
smaller when transfers overlap owner work.

## 13. Mode-3 Memory Cost

The main live categories are

```math
M_{r,l}^{CA3} = M_{home,l}
 + M_{owned\ fetch,l}
 + M_{installed,l}
 + M_{generator,l}
 + M_{deferred,l}
 + M_{buffers,l}
 + M_{scratch,l}.
```

- `M_owned fetch`: remote entry boxes for components owned here; releasable
  after elimination and payload completion.
- `M_installed`: shared component copies needed for local pair replay.
- `M_generator`: `X_NR/temp2`, `X_RR_full`, and metadata retained for lazy
  far reconstruction and transition.
- `M_deferred`: formed or pending near updates waiting for canonical replay.
- `M_buffers`: pooled sends, receives, serialization offsets, and in-flight
  payloads.
- `M_scratch`: streamed ID, right-side solve, and replay workspaces.

Mode 3 removes full replicated halo factorization state, but it can still
have a higher transient peak than Color because fetched entry boxes, installed
copies, routed payloads, and pool capacity overlap.  The current full-`X_NR`
policy raises bytes and live payload memory in exchange for fewer install
rounds.

Tracked live memory and RSS must be reported separately:

- pooled capacity is reserved and reusable but may not be logically live;
- tcmalloc with release rate zero retains freed pages in RSS;
- parent construction/reduction can overlap child boxes and payload buffers;
- `MaxRSS` records a historical peak and cannot reveal later release.

## 14. Comparing Color, CA0, and CA3

| Property | Color, lazy mode 2 | CA0, staged mode 2 | CA3, staged mode 2 |
|---|---|---|---|
| Boundary elimination | unique fixed box owner | replicated halo | unique component owner |
| In-level dependency depth, 3D | 8 parity colors | none after entry halo | 3 component colors |
| Entry transfer | assisting/neighbor state as needed | full replicated halo | owned-component entry fetch |
| Update representation | generators for far and remote near | explicit near, generators for far | routed eliminated state, explicit near, generators for far |
| Dense far storage | no | no with lazy mode 1 | no |
| Replicated factor work | none | potentially large | none |
| Load balance | fixed by home partition | geometric halo imbalance | assignable by component |
| Main latency risk | eight communication tails | initial halo tail | longest component DAG chain |
| Main volume risk | generator fanout | halo duplication | fetch plus FULL/COMPACT/SKELETON routing |
| Main memory risk | near blocks and generators | replicated halo | payload/install/buffer overlap |

The expected regime map is:

- large rank bricks and fine levels: CA0 can win because replication is cheap
  and one staged entry transfer replaces many Color rounds;
- small bricks/coarser distributed levels: CA0 replication and imbalance grow,
  favoring CA3 or Color;
- latency-dominated strong scaling: CA3's component depth can beat Color's
  parity depth if its extra bytes remain hidden;
- bandwidth/memory-constrained cases: Color may win because it avoids owned
  entry fetch and large COMPACT state;
- thin components with high-rank boxes: CA3 can lose node efficiency unless
  batching and fixed-chunk splits recover parallelism.

## 15. Hybrid Level Selection

`CA_level=L` requests CA on eligible levels `l>=L`; smaller level numbers are
closer to the root.  Therefore:

- `CA_level=10000`: all practical levels use Color;
- small `CA_level`: all eligible fine/deep levels request CA;
- intermediate values: coarse levels use Color and deeper levels use CA.

The automatic geometry fallback may still replace a requested CA level with
Color.  This gives the desired hybrid policy: use CA where rank bricks are
large enough to amortize replication or component fetch, and Color where CA
geometry is absent or its payload overhead is unfavorable.

For mode 3, a future model-based selector should compare estimated

```math
T_l^{Color},\quad T_l^{CA0},\quad T_l^{CA3}
```

using measured rank, boundary/halo counts, peer volume, and prior-level box
cost rather than a fixed level alone.

## 16. Owner-Based Solve and Multiply

### 16.1 Solve

Replicated CA solve applies factors on ghost boxes, which is especially costly
because triangular-solve sweeps are memory-bandwidth bound.  Mode 3 records
component ownership and uses a distributed vector-only sweep:

- forward substitution routes additive `X_NR x_R` contributions to holders
  of later components;
- backward substitution routes finalized neighbor values in reverse
  elimination order;
- factors remain at or are restored to home ranks by the FULL payload;
- contributions land in fixed per-target slots and fold in canonical order;
- one aggregated message is sent per graph-derived neighbor, not an
  all-to-all factor gather.

The communicated object shrinks from matrices/factors to vectors.  For
multiple right-hand sides the current implementation is not yet a fully
blocked owner solve, so this remains an optimization opportunity.

### 16.2 Multiply

The mode-3 multiply follows the owner graph:

- the upward/`W` sweep pulls stepped box state to owners of higher-color
  neighbors;
- the downward/`V` sweep pushes contributions to the target holder in reverse
  order;
- only vectors move; factors are not gathered.

This preserves component ownership beyond factorization and avoids paying a
replicated-halo penalty on every subsequent operator application.

## 17. Recommended Paper Experiments

### 17.1 Fair baselines

All three methods should share:

- the same right-side solve and diagonal solve kernels;
- fixed 512-row and 256-column chunks;
- streamed sketching;
- lazy far regeneration;
- symmetric edge storage;
- sequential BLAS and identical allocator policy.

Recommended configurations:

| Label | Configuration |
|---|---|
| Color | `CA_level=10000`, sketch 2, lazy 2 |
| CA0-legacy | CA levels, owner 0, staged halo 0, sketch 2, lazy 1 |
| CA0-staged | CA levels, owner 0, staged halo 2, sketch 2, lazy 1 |
| CA3-serial | owner 3, staged halo 2, serial oracle 1, sketch 2, lazy 1 |
| CA3-fast | owner 3, staged halo 2, serial oracle 0, sketch 2, lazy 1 |

### 17.2 Per-level measurements

- unique and replicated boxes eliminated;
- component count, size distribution, owner load, and wave width;
- critical-path component timestamps and idle waits;
- entry-fetch, FULL, COMPACT, and SKELETON bytes separately;
- message counts, wait tails, pack/install/replay time;
- streamed-ID, box, right-solve, near-update, transition, and reduction time;
- live home, ghost/fetch, installed, generator, deferred, buffer, scratch, and
  parent memory;
- RSS, page faults, and allocator free/unmapped counters;
- solve/multiply vector bytes and time;
- accuracy, ranks, skeleton counts, and fast-versus-serial checksums.

### 17.3 Scaling questions

1. At what boxes/rank does CA0 replication become more expensive than its
   latency advantage?
2. At what scale does CA3's three-stage DAG beat Color's eight-stage chain?
3. How much extra volume does the current full-`X_NR` COMPACT policy create,
   and is it recovered by one-round installs?
4. Does component-count balancing correlate with factor time when numerical
   ranks are anisotropic?
5. Which levels should the hybrid assign to Color, CA0, or CA3?
6. How do real Laplace and complex Helmholtz/VIE change the balance through
   scalar size, rank, and Bunch-Kaufman cost?

## 18. Limitations and Open Questions

- Mode 3 currently requires the symmetric non-Hermitian path, streamed
  sketching, and lazy far updates.
- CA does not implement generated-near mode 2.
- The owner assignment balances box counts rather than measured rank-weighted
  work or communication volume.
- A large face component is indivisible; component granularity can limit load
  balance.
- Full `X_NR` improves install latency but dominates COMPACT payload bytes.
- Entry fetch, payload, and parent construction can overlap in memory even
  when their logical lifetimes are short.
- Progress overlap benefits from `MPICH_ASYNC_PROGRESS=1`, but correctness
  must not depend on asynchronous progress.
- Transition/process-reduction overlap is not yet the full cross-level
  pipeline suggested by the dataflow graph.
- A blocked owner solve for many right-hand sides remains future work.
- Nonsymmetric matrices need paired left/right generator and payload rules.

## 19. Suggested Paper Structure

1. Replicated communication-avoiding skeletonization and its work/halo cost.
2. Staged and pruned immutable halo exchange.
3. Why box ownership fails and component ownership is sufficient.
4. Component graph, deterministic owner assignment, and payload routing.
5. Asynchronous execution, batching, component waves, and canonical replay.
6. Work, memory, volume, and critical-path models for CA0/CA3/Color.
7. Owner-based solve and multiply.
8. Hybrid level selection.
9. Ablation, scaling, accuracy, and memory experiments.
10. Regime map, limitations, and nonsymmetric extensions.
