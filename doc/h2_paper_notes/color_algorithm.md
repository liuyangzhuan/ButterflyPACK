# Optimized Distributed Color H2 Factorization

Status: paper-development notes describing the ButterflyPACK integration.

## 1. Candidate Paper Thesis

The Color algorithm gives every box a unique MPI owner and uses a `2^d`-color
ordering to make same-wave boxes independent.  It avoids redundant
factorization but traditionally pays for dense Schur-product storage,
pair-product communication, and `2^d` dependent communication rounds per
level.  The recent implementation attacks the first two costs without
changing the color dependency depth:

1. BLAS-3 right-side solves and deterministic task splitting improve the
   underfilled coarse-level kernels.
2. Streamed sketching removes the materialized ID target.
3. Lazy Schur updates retain low-rank elimination generators instead of dense
   two-hop products.
4. Generated-near transport communicates each source representation rather
   than every endpoint-pair product.

The strongest defensible claim is therefore not that Color becomes
communication avoiding.  It is that the same unique-owner algorithm can move
from dense pair state to streamed/generator state, reducing peak memory and
bandwidth while retaining its original synchronization structure.

## 2. Baseline Color Schedule

Color a box by the parity of its grid coordinates, equivalently by the low
`d` Morton bits.  There are

```math
q=2^d
```

colors, and no two one-hop neighbors have the same color.  At each distributed
level, ButterflyPACK currently forms:

1. `q` boundary waves, one for each parity color;
2. `q` interior sub-waves, also parity colored.

After a boundary wave, factor updates needed by remote boxes are exchanged
and installed before a dependent later color starts.  Interior work is local,
but the final boundary payload must be drained before the first interior
sub-wave.  Thus the algorithm has `q` elimination-dependent communication
steps per level, plus the separate parent transition/process redistribution.

For an eliminated box `E`, all update endpoints lie in `N_1(E)`.  Two such
endpoints can be separated by two box widths, so eliminating `E` modifies
both one-hop and strict two-hop interactions.  In 3D, one source has at most
26 neighbor endpoints and `binom(26,2)=325` endpoint pairs before geometric,
boundary, ownership, and symmetry pruning.

## 3. Per-Box Algebra

Let the transformed box block be ordered as skeleton `S`, redundant `R`, and
neighbor rows `N`.  The implementation factors the redundant pivot block
`X_RR` and forms

```math
T_1 = -X_{SR}X_{RR}^{-1},
\qquad
T_2 = -X_{NR}X_{RR}^{-1}.
```

`T_1` is retained for the triangular solve.  Endpoint row slots of `T_2` are
the Schur generators used to update neighboring boxes.  With the original
unfactored redundant block `D_E=X_RR_full`, a pair update can be regenerated
as

```math
\Delta^{(E)}_{A,C}
 = -T_{2,E}[A]D_ET_{2,E}[C]^T.
```

The eager and lazy paths differ only in when this product is formed and how
long its factors are retained.

## 4. BLAS-Level and Shared-Memory Optimizations

### 4.1 Solve tall matrices from the right

The expensive operands in `T_1` and `T_2` are tall matrices multiplied by
`X_RR^{-1}` on the right.  The old path transposed a tall matrix, called a
left-side LAPACK solve, and transposed it back.  The new path operates in its
native layout.

For Cholesky, with `X_RR=LL^T`, it applies two right-side `TRSM`s:

```math
B X_{RR}^{-1} = (B L^{-T})L^{-1}.
```

For LU, it applies right-side triangular solves with `U` and `L`, followed by
column swaps dictated by the LAPACK pivots.  For Bunch-Kaufman, it converts
the factor representation with `syconv`, applies the pivot permutations,
uses right-side triangular solves, and handles the 1-by-1 and 2-by-2 diagonal
pivots directly.  This also covers complex-symmetric factors without
conjugating the transpose.

Benefits:

- no full tall-matrix transpose buffers;
- no pair of strided transpose passes;
- Level-3 kernels over many right-hand sides;
- natural row blocking for OpenMP tasks.

The independent diagonal solve in the H2 solve phase uses the blocked
`dsytrs2/zsytrs2` path instead of the older Level-2 solve.

### 4.2 Fixed numerical chunks

Changing an OpenBLAS call's dimensions can change its internal blocking and
last-bit rounding.  A naive `K`-way split therefore made results depend on
wave size, thread count, and arrival timing.  The integrated implementation
uses fixed numerical chunks:

- right-side solves: 512 matrix rows per solve call;
- split GEMMs: 256 output columns per GEMM call.

`H2_GEMM_split` limits how many OpenMP tasks execute the fixed chunks, but it
does not change chunk boundaries.  The arithmetic call sequence for each
output chunk is consequently independent of the worker count.

### 4.3 Split only compute-heavy work

If a wave contains fewer boxes than OpenMP workers, the runtime gives a box
up to

```math
K_B = \min\left(K_{max},\left\lfloor t_l/|\mathcal W|\right\rfloor\right)
```

tasks.  Splitting is retained for kernel columns, right-side solves, and the
large GEMMs used by lazy regeneration and streamed sketches.  Sparse sketch
AXPYs, row gathers, and transposes remain on the box's owning thread; tests
showed that splitting those memory-bound loops increased cache traffic and
task barriers.

The outer box loops use dynamic scheduling because numerical rank and
neighbor state make box costs nonuniform.

### 4.4 In-place accumulation

The deferred owner pass accumulates Schur updates directly into their target
with GEMM `beta=1`.  It does not allocate a dense product followed by a
separate add.  This removes one product-sized temporary and one memory pass.

### 4.5 Sequential BLAS requirement

These task splits assume one BLAS thread per call.  OpenMP provides the
parallelism across boxes and fixed chunks.  Running a threaded BLAS inside
these tasks causes nested oversubscription and defeats the deterministic call
shape model.

## 5. Streamed Sketching for the ID

### 5.1 Materialized formulation

Let the ID training matrix for box `B` be

```math
W_B \in \mathbb C^{m_B\times n_B}.
```

Its columns correspond to the active indices of `B`; its rows concatenate
selected interactions with the two-hop environment and optional training
rows.  The randomized path applies a sparse embedding

```math
S_B\in\mathbb R^{s_B\times m_B},
\qquad Y_B=S_BW_B,
```

then runs pivoted QR/ID on `Y_B`.  `H2_use_sketch=1` still materializes all of
`W_B` before forming `Y_B`.

### 5.2 Streaming identity

Partition the rows into blocks,

```math
W_B = [W_1^T, W_2^T,\ldots,W_J^T]^T,
\qquad
S_B = [S_1,S_2,\ldots,S_J].
```

Then

```math
Y_B = \sum_{j=1}^{J} S_jW_j.
```

`H2_use_sketch=2` evaluates or gathers one `W_j`, immediately accumulates it
into `Y_B`, and reuses the block buffer.  The sketch distribution is the same
as in the materialized sparse-sketch path; only the evaluation order and
storage lifetime change.

### 5.3 Applying lazy fill in sketch space

For a row endpoint `A`, target box `B`, and eliminated common-neighbor source
`E`, write

```math
W_{A,B}^{(E)}=-G_{E,A}D_EG_{E,B}^T.
```

Define

```math
P_{E,B}=D_EG_{E,B}^T,
\qquad
Q_E=\sum_A S_A G_{E,A}.
```

Associativity gives

```math
\sum_A S_AW_{A,B}^{(E)}=-Q_EP_{E,B}.
```

The code therefore never regenerates a dense far block for ID.  It first
streams kernel/stored rows into the sketch while recording the sparse draws,
then, source by source in canonical order, forms `Q_E` and applies one GEMM
against `P_{E,B}`.  At any instant only one source's row slice, `Q_E`, and
`P_{E,B}` need be live.

If the sparse embedding has `zeta` signed nonzeros per input row, the lazy
fill work for source `E` is approximately

```math
O\left(\zeta m_E r_E + s_B r_E n_B\right),
```

rather than forming every dense pair contribution before sketching it.

### 5.4 Workspace memory

The materialized target needs, per concurrent box worker,

```math
M_{ID,mat}(B)=\Theta(\sigma m_B n_B).
```

The streamed path needs approximately

```math
M_{ID,str}(B)=O\left(\sigma(s_Bn_B+a_Bn_B+s_Br_E+r_En_B+a_Br_E)
                    + m_B\zeta\right).
```

The final term stores compact sketch indices/signs, not scalars.  The main
distinction is that scalar storage no longer grows as `m_B n_B`.  Rank-level
peak scratch multiplies these per-worker quantities by the number of boxes
processed concurrently, not necessarily by all OpenMP workers.

### 5.5 Current applicability

Streamed sketching is active on both Color and CA levels for the symmetric
non-Hermitian path above level 1.  Adaptive proxy selection
(`H2_ID_proxy=2`) currently uses the materialized route because its selected
row set is constructed adaptively before sketching.

## 6. Lazy Schur Modes

### 6.1 Mode 0: eager products

`H2_lazy_schur=0` forms dense Schur products when `E` is eliminated.

- Near products update persistent one-hop blocks.
- Far products update persistent strict two-hop blocks.
- Products needed remotely are buffered and exchanged after color waves.
- The ID reads stored far blocks.
- Parent construction consumes stored child interactions.

This minimizes regeneration but has the largest block storage and product
payloads.

### 6.2 Mode 1: lazy far, explicit near

`H2_lazy_schur=1` keeps one-hop behavior eager but does not store two-hop
modified blocks.  It retains `T_2` endpoint slots and `D_E`, and ships remote
generators where another rank must reconstruct a far interaction.

- Near products are still accumulated and transported explicitly.
- Far products are regenerated in sketch space for ID.
- Far contributions needed by the parent are regenerated during transition;
  cached `P_{E,C}=D_EG_{E,C}^T` blocks are reused across sibling rows.
- Remote generator boxes persist until the level no longer needs them.

This is primarily a memory optimization, with a partial communication
benefit from removing far pair products.

### 6.3 Mode 2: lazy far, generated near

`H2_lazy_schur=2` retains the mode-1 far behavior and also suppresses dense
near products whose remote endpoint would otherwise receive them.  A source
ships generators once to the required recipients.  Receivers form diagonal,
near, and replacement updates in canonical source order and accumulate them
into their persistent near blocks.

This changes the wire representation from one dense matrix per endpoint pair
to one low-rank source representation per recipient.  It does **not** remove
persistent one-hop blocks: those blocks encode accumulated level history and
are required by later eliminations and the parent transition.

Mode 2 is currently implemented only on Color levels.  CA levels reduce it to
mode 1.

## 7. Per-Level Communication Cost

Let `E_l^x` be boundary sources whose effects cross a rank boundary.  Let
`P_E^near` and `P_E^far` be the remotely relevant near and far endpoint-pair
sets.  Ignoring headers and aggregation, the leading byte volumes are

```math
V_l^{(0)} \simeq
  \sigma\sum_{E\in E_l^x}
  \sum_{(A,C)\in P_E^{near}\cup P_E^{far}} n_A n_C,
```

```math
V_l^{(1)} \lesssim
  \sigma\sum_{E\in E_l^x}
  \sum_{(A,C)\in P_E^{near}} n_A n_C
  + \sum_{E\in E_l^x} f_E g_E,
```

```math
V_l^{(2)} \lesssim
  \sum_{E\in E_l^x} f_E g_E.
```

The inequalities emphasize that payload routing can send only endpoint slots
needed by a recipient.  Under homogeneous sizes, one eager interior source
has an upper bound of `c_d n_l^2` product scalars, while one unpruned generator
has about

```math
z_1n_lr_l+r_l^2+r_lk_l
```

scalars.  In 3D, with `r_l/n_l=rho`, the ideal source-size ratio before
fanout is roughly

```math
\frac{26\rho+\rho^2+\rho(1-\rho)}{325}.
```

This is a geometry-constant reduction even when `r_l=Theta(n_l)`, and it
becomes rank-dependent when `r_l<<n_l`.

All three modes retain the same Color dependency depth.  With payloads
aggregated by peer, a simple communication model is

```math
T_{comm,l}^{(j)}
 \approx q\,\chi_l\alpha + \beta V_l^{(j)} + T_{apply,l}^{(j)}.
```

Lazy modes reduce `V_l`; they do not remove the `q` round latency.  A late
rank can therefore expose a tail at every color boundary.

## 8. Per-Level Memory Cost

For symmetric edge storage and homogeneous box sizes, persistent one-hop
interaction storage is approximately

```math
M_{near,l}\simeq
  \frac{\sigma}{2} b_l z_1 n_l^2.
```

Mode 0 additionally stores strict two-hop interactions,

```math
M_{far,l}^{(0)}\simeq
  \frac{\sigma}{2} b_l z_2 n_l^2.
```

In an interior 3D region, the far/near ratio from geometry alone is
`98/26=3.77`.  This explains why dense two-hop blocks dominated earlier peak
breakdowns even after symmetric edge storage.

Modes 1 and 2 replace that term with retained generators,

```math
M_{gen,l}\simeq
  \sum_{E\ eliminated}\sigma
  \left(r_E^2+r_Ek_E+r_E\sum_{A\in N_1(E)}n_A\right),
```

including remote generator boxes needed later in the level.

The transient communication state has the following leading behavior:

| Mode | Dominant pending state before a transport/application point |
|---|---|
| 0 | dense near and far pair products, up to `O(sigma |wave_bdry| c_d n_l^2)` |
| 1 | dense near products plus generator sends/receives |
| 2 | generator sends/receives; near products are formed directly into persistent targets |

Thus a useful peak model is

```math
M_{peak,l}^{(j)} = M_{factors,l}+M_{near,l}
 + [j=0]M_{far,l}^{(0)}+[j>0]M_{gen,l}
 + M_{ID,l}^{(j)}+M_{pending,l}^{(j)}+M_{transition,l}.
```

`M_transition,l` matters because parent boxes, serialization buffers, and
child-level state can overlap.  Allocator-retained free pages affect RSS but
must not be counted as algorithmically live memory.

## 9. Computation Tradeoffs

Mode 0 forms each required dense product once and reads it later.  Modes 1/2
reduce writes, network bytes, and dense storage, but add:

- generator-row gathers;
- source-wise sketch-space GEMMs;
- transition-time regeneration;
- receiver-side near formation in mode 2.

The mode-2 near GEMMs are mostly relocated work: eager Color forms them on the
source/owner before transport, while generated-near Color forms them on the
recipient.  Performance depends on whether avoided memory traffic and bytes
outweigh poorer receiver occupancy or repeated regeneration.

The expected regime map is:

- fine levels: many boxes fill the OpenMP team; bandwidth and memory capacity
  favor streaming and lazy modes;
- coarse levels: few boxes make fixed-chunk task splitting essential;
- strong scaling: mode 2 reduces bytes, but the unchanged eight-round 3D
  latency can dominate;
- high numerical rank: generator sizes approach dense-neighbor sizes, but
  still avoid the endpoint-pair geometry factor;
- very low rank: generator transport and storage gain both geometrically and
  algebraically.

## 10. Accuracy and Reproducibility

1. Streamed and materialized sketching represent the same random projection
   in exact arithmetic.
2. Lazy and eager Schur updates represent the same block in exact arithmetic.
3. Lazy regrouping changes floating-point association, so eager and lazy
   factors need not be bitwise equal.
4. Sources are replayed in canonical `(color,wave,Morton)` order.
5. Fixed 256-column and 512-row chunks make the result independent of the
   number of task workers consuming those chunks.
6. Numerical validation should report `acc_mvp`, forward error, backward
   error, rank distributions, and skeleton counts.  Residual alone cannot
   distinguish approximation error from conditioning.

## 11. Recommended Paper Experiments

### 11.1 Core ablation

Use identical matrix, tree, process grid, tolerance, BLAS, and allocator:

| Case | Sketch | Lazy mode | Purpose |
|---|---:|---:|---|
| C0 | 1 | 0 | materialized-sketch eager reference |
| C1 | 2 | 0 | isolate streamed ID |
| C2 | 2 | 1 | isolate lazy far |
| C3 | 2 | 2 | add generated near |

Optionally include `H2_use_sketch=0` to quantify sketching accuracy and cost,
not as the scalable baseline.

### 11.2 BLAS/task ablation

- old transpose/left-solve versus right-side solve;
- unsplit versus `H2_GEMM_split={2,4,8,16}`;
- report fine-level and coarse-level box times separately;
- verify fixed-thread and dynamic-thread runs have matching checksums.

### 11.3 Required per-level measurements

- total elimination, box, deferred-owner, exchange, transition, and process
  reduction times;
- sent/received payload bytes and message count by color;
- live near, far, generator, pending, scratch, and transition bytes;
- process RSS and minor faults, clearly separated from tracked live memory;
- average/max `n_l`, `k_l`, `r_l`, actual endpoint-pair count, and generator
  fanout `f_E`;
- numerical errors and factorization success.

### 11.4 Scaling axes

- 2D and 3D;
- Laplace and oscillatory Helmholtz/VIE kernels;
- weak scaling at fixed boxes/rank;
- strong scaling into the latency-dominated regime;
- at least two tolerances to vary numerical rank.

## 12. Limitations and Open Questions

- Color still has `2^d` dependent communication steps per level.
- Generated-near mode is symmetric-only in the current integration.
- Adaptive proxy selection is not streamed.
- Lazy transition regeneration can become visible if parent construction has
  too little parallelism or cache reuse.
- Generator fanout can duplicate `D_E` and common metadata across recipients;
  measured recipient-specific bytes are needed for a tight model.
- The current allocator policy improves page reuse but may increase retained
  RSS; paper memory plots must show both live bytes and RSS.
- Future nonsymmetric support needs two interaction orientations and separate
  left/right generators, changing the one-half symmetric-storage factors.

## 13. Suggested Paper Structure

1. Recursive skeletonization and distributed parity coloring.
2. Bottlenecks of eager Color: tall solves, ID target, far storage, pair
   communication.
3. BLAS-3 and deterministic task kernels.
4. Streamed randomized ID, including the sketch-space lazy-fill derivation.
5. Lazy Schur modes and generator transport.
6. Per-level memory, bandwidth, and latency model.
7. Implementation and reproducibility.
8. Accuracy, ablation, and scaling experiments.
9. Regime map and limitations.
