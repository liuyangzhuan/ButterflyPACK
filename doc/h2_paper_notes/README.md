# H2 Color and CA Paper Notes

Status: working technical notes for paper development, 2026-09-15.

These notes record the algorithmic understanding used while integrating the
updated `Communication-Avoiding-FMM_New` implementation into ButterflyPACK.
They are not a claim that every bound below is tight.  In particular, the
cost models expose geometry constants and recipient counts that should be
measured rather than hidden inside big-O notation.

## Documents

- [color_algorithm.md](color_algorithm.md): the optimized parity-color
  algorithm, BLAS and task-level kernels, streamed sketching, and lazy Schur
  modes 0/1/2.
- [ca_algorithm.md](ca_algorithm.md): replicated CA (mode 0), staged halo
  exchange, and the asynchronous component-owner engine (mode 3).

## Common Algorithmic Setting

At level `l`, each tree box `B` owns an active index set.  An interpolative
decomposition (ID) splits it into skeleton indices `S_B` and redundant indices
`R_B`.  After applying the interpolation transform, the redundant block is
factored and eliminated.  If `E` is the eliminated source box, its update to
two boxes `A,C` in its one-hop neighborhood has the generator form

```math
\Delta^{(E)}_{A,C}
  = -G_{E,A} D_E G_{E,C}^{T}.
```

Here `D_E` is the retained redundant pivot block (`X_RR_full` in the code)
and `G_{E,A}` is the row slot for endpoint `A` in the finalized `temp2` /
`X_NR` generator.  The transpose is non-conjugating because the currently
supported path is symmetric or complex-symmetric, not Hermitian.

Both algorithms perform the same algebra.  They differ in:

1. the legal elimination ordering;
2. which rank eliminates a boundary box;
3. whether a Schur product is materialized or retained as generators;
4. when and in what representation cross-rank state is communicated.

## Notation

Quantities may vary by box.  A subscript `B` or `E` denotes the exact local
quantity; a subscript `l` denotes a representative level average.

| Symbol | Meaning |
|---|---|
| `d` | spatial dimension |
| `q = 2^d` | number of parity colors |
| `p_l` | active MPI ranks at level `l` |
| `b_l` | home boxes on one active rank |
| `b_l^bdry` | home boundary boxes on one active rank |
| `h_l` | halo/ghost boxes imported by a replicated CA rank |
| `n_B` | active degrees of freedom entering box `B` |
| `k_B` | skeleton count after ID |
| `r_B = n_B-k_B` | redundant count |
| `m_B` | row count in the ID training matrix `W_B` |
| `s_B` | row count of the randomized ID sketch |
| `a_B` | largest row block streamed into the ID |
| `z_1 = 3^d-1` | maximum number of strict one-hop neighbors |
| `z_2 = 5^d-3^d` | maximum number of strict two-hop neighbors |
| `c_d = binom(z_1,2)` | upper bound on endpoint pairs updated by one interior source |
| `chi_l` | number of communicating peer ranks at the level |
| `f_E` | number of recipients of source `E`'s generator payload |
| `sigma` | bytes per scalar, normally 8 real or 16 complex |
| `alpha`, `beta` | message startup time and inverse network bandwidth |
| `t_l` | OpenMP workers assigned to an MPI rank at level `l` |

For 3D, `q=8`, `z_1=26`, and `z_2=98`.  Domain boundaries and already
eliminated neighbors reduce these maxima substantially.

## Cost-Model Conventions

The notes keep four costs separate:

- **Persistent level memory**: matrices and factors needed until transition or
  solve.
- **Transient peak memory**: ID workspaces, pending messages, product buffers,
  serialization buffers, and parent-construction overlap.
- **Communication volume**: payload bytes crossing MPI rank boundaries.
- **Latency depth**: dependent communication stages on the critical path.

For a set of source boxes `E`, define the leading generator size

```math
g_E = \sigma\left(r_E^2 + r_E k_E
                  + r_E\sum_{A\in N_1(E)} n_A\right).
```

The terms are `X_RR_full`, `X_RS_entry`, and the endpoint row slots in
`temp2/X_NR`.  Headers, index vectors, pivots, and recipient-specific pruning
are lower-order terms in the dense regimes considered here.  Under the
homogeneous approximation,

```math
g_l = \Theta\left(\sigma(z_1 n_l r_l+r_l^2+r_l k_l)\right).
```

For a source `E`, materializing every Schur pair has leading size

```math
p_E = \sigma\sum_{(A,C)\in P(E)} n_A n_C,
```

where `P(E)` is the subset of endpoint pairs that must be stored or sent.
The interior upper bound is `Theta(sigma c_d n_l^2)`.  These formulas should
be summed over actual boxes for paper plots; replacing all sizes by averages
is useful only for interpretation.

## ButterflyPACK Option Map

| Option | Integrated behavior |
|---|---|
| `CA_level=L` | levels `l>=L` request CA; `10000` gives full Color for practical tree depths |
| `H2_use_sketch=0` | materialize `W_B` and run ID on the full target |
| `H2_use_sketch=1` | materialize `W_B`, then form the randomized sketch |
| `H2_use_sketch=2` | stream row blocks directly into the randomized sketch |
| `H2_lazy_schur=0` | eager near and far Schur products |
| `H2_lazy_schur=1` | retain generators and regenerate far products; near products remain explicit |
| `H2_lazy_schur=2` | Color only: lazy far plus near updates communicated as generators |
| `H2_GEMM_split=K` | maximum task split for one compute-heavy box kernel; `0/1` disables splitting |
| `H2_XRR_factor=0` | factor each `X_RR` with Bunch-Kaufman (default) |
| `H2_XRR_factor=1` | factor each `X_RR` with LU and partial pivoting; symmetric half storage is unchanged |
| `H2_CA_staged_halo=0` | monolithic CA level-entry gather |
| `H2_CA_staged_halo=2` | shell/blue/purple/green staged CA gather with overlap |
| `H2_CA_owner_component=0` | replicated CA boundary elimination |
| `H2_CA_owner_component=3` | asynchronous component-owner elimination |
| `H2_CA_owner_serial=1` | serialize mode-3 scheduling as a correctness oracle |

Important interactions in the current integration:

1. `H2_use_sketch=2` is supported on both Color and CA levels.  It is enabled
   for levels above level 1 when the matrix is symmetric non-Hermitian.
   Adaptive proxy mode (`H2_ID_proxy=2`) then selects its rows against the
   streamed sketch; with `H2_use_sketch=0` or `1` it uses the materialized
   target.
2. A nonzero `H2_lazy_schur` requires `H2_use_sketch=2`.
3. On a CA level, requested lazy mode 2 is intentionally reduced to mode 1;
   generated-near transport is not implemented in the CA path.
4. Mode-3 CA additionally requires streamed sketching, lazy far updates, and
   the symmetric non-Hermitian factorization path.
5. Distributed CA is used only when a rank has enough local geometry: at
   least `7^3=343` boxes in 3D or `5^2=25` in 2D.  Otherwise that level
   automatically runs Color.  One-dimensional problems always run Color.

## Execution Assumptions

- BLAS is sequential.  Parallelism is controlled by the H2 OpenMP wave and
  fixed-chunk task layers.  A threaded vendor BLAS would oversubscribe these
  teams and can invalidate the intended scheduling/reproducibility model.
- Symmetric edge storage keeps one dense orientation and reconstructs the
  reverse view by transpose.  The cost tables include the corresponding
  factor of approximately one half where stated.  A future nonsymmetric path
  must store and transport both orientations.
- Dynamic threading may vary `t_l` with the active communicator.  Fixed
  numerical chunks make the arithmetic independent of how many tasks consume
  those chunks.
- The current run scripts preload tcmalloc and use
  `TCMALLOC_RELEASE_RATE=0`.  This reduces page release/refault overhead but
  retains free spans in RSS.  It is an allocator policy, not an algorithmic
  memory reduction.

## Claims and Evidence

The intended paper hierarchy is:

1. **Algebraic claims**, such as streamed sketch equivalence and generator
   regeneration, should be proved from the formulas.
2. **Complexity claims** should use the symbols above and report actual
   recipient/pair counts alongside asymptotic expressions.
3. **Performance claims** should come only from new ButterflyPACK runs with
   identical kernels, tolerances, process grids, allocator policy, and BLAS.
4. **Bitwise reproducibility** means independence from arrival order and
   thread count within one implementation.  It does not imply bitwise
   equality between eager and lazy algebraic reorderings.

## Implementation Source Map

- Public options and defaults:
  `SRC/BPACK_defs.f90`, `SRC/BPACK_utilities.f90`.
- Per-level algorithm selection and Color driver:
  `h2_parallel/structured/butterfly_factorization.hpp`.
- Elimination, streamed ID, lazy regeneration, right-side solves, and fixed
  chunks: `h2_parallel/core/factorization.hpp`.
- Runtime options and tree geometry: `h2_parallel/core/tree.hpp` and
  `tree_impl.hpp`.
- Staged CA gather: `h2_parallel/core/staged_halo.hpp`.
- Component graph: `h2_parallel/core/dataflow.hpp`.
- Component-owner execution: `h2_parallel/core/owner_schedule.hpp` and
  `owner_exchange.hpp`.
- Owner solve and multiply: `h2_parallel/core/owner_solve.hpp` and
  `owner_mul.hpp`.

The upstream development record is in
`Communication-Avoiding-FMM_New/concise-algorithm-ref/`.  Those files contain
historical experiments and superseded designs; these notes describe the
behavior integrated into ButterflyPACK.
