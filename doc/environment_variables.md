# Environment variables

ButterflyPACK's H2 code (`h2_parallel/`, format 7) and the GPU backends of
H2 (option `H2_use_gpu`) and HODLR (format 1, option `HODLR_use_gpu`) read
six environment variables. Algorithmic choices are solver options; the
environment variables cover what is set once per process (device memory, MPI,
CPU cores) and debugging aids (checks and traces).

| Variable | Purpose |
|---|---|
| [`BPACK_GPU_HEAP_FRACTION`](#bpack_gpu_heap_fraction) | Size of each rank's device memory pool |
| [`BPACK_GPU_EXCHANGE_MB`](#bpack_gpu_exchange_mb) | Size of the device arena for MPI messages |
| [`BPACK_GPU_AWARE_MPI`](#bpack_gpu_aware_mpi) | Whether MPI may send from and receive into device memory |
| [`BPACK_MAX_CPUS_PER_NODE`](#bpack_max_cpus_per_node) | Cores per node shared by the ranks' threads (H2 factorization) |
| [`BPACK_CHECK`](#bpack_check) | Checks of GPU results against the CPU, and comparison modes |
| [`BPACK_TRACE`](#bpack_trace) | Extra output for performance analysis and debugging |

At the start of an H2 compression or factorization and of a HODLR GPU
construction, rank 0 prints once per process the checks and traces that are
on (`BPACK: BPACK_CHECK: solve, replica`), any value of the two lists it does
not know (ignored), and any variable of an earlier version that is set (no
longer read; see the [last section](#variables-of-earlier-versions)).

## Resources

These are per process: with several matrices in one process (a Gaussian
process and its derivative matrices, or an H2 and a HODLR matrix), they are
shared by all of them.

### BPACK_GPU_HEAP_FRACTION

Default `0.85`. The device memory pool of each rank, as a fraction of the GPU
memory free when the pool is reserved (the first GPU work of the process). It
must lie in (0, 1). Ranks of a node that use the same GPU (found by its
device UUID) split it: each takes an equal part of the fraction of the memory
they found free together.

Everything the H2 and HODLR backends keep on the device lives in the pool.
Within it, the H2 factorization keeps a level's solve factors on the device
when an upper bound of their size fits under 0.8 of the pool (the others are
uploaded at the first solve), and the H2 compression keeps the blocks of the
device matvec while the pool stays under 0.9 of its size. With several H2
matrices in one process, each keeps its device data; the least recently used
matrices' data is released when the others hold more than half the pool.

### BPACK_GPU_EXCHANGE_MB

Default: min(2560, GPU memory in MB / 16), divided among the ranks that share
the GPU. A device arena of this many MB for the messages MPI sends from and
receives into device memory, reserved before the pool. It exists only when
the run has more than one rank and MPI is GPU-aware (`BPACK_GPU_AWARE_MPI`).
MPI registers it once; device buffers outside it make the transfers much
slower. A message that does not fit goes through pinned host buffers instead.

### BPACK_GPU_AWARE_MPI

`1`: MPI is GPU-aware, so the H2 and HODLR backends send and receive device
memory directly (through the exchange arena). `0`: every message is staged
through pinned host buffers. Unset: the value of `MPICH_GPU_SUPPORT_ENABLED`
(Cray MPICH; set it to `1` on Perlmutter, with the `craype-accel-nvidia80`
module loaded at build time).

### BPACK_MAX_CPUS_PER_NODE

Default: unset (off). Dynamic threading of the H2 factorization: the ranks of
a node share at most this many cores, pinned per rank, and a level on which
only some of a node's ranks are active gives them the idle ranks' cores. Off,
every rank keeps `OMP_NUM_THREADS` threads.

## BPACK_CHECK

A comma-separated list of the values below, for example
`BPACK_CHECK=solve,replica`. Each check does extra work, mostly a CPU
reference run of the same step, so timings taken with a check on are not
meaningful. The checks print on rank 0 unless noted.

| Value | Applies to | What it does | Expected result |
|---|---|---|---|
| `solve` | H2 GPU (`H2_use_gpu`), factorization | Every solve and every multiply with the factors on the GPU is run again on the host. Prints `GPU solve check: \|x_gpu - x_host\| / \|x_host\|` and `GPU multiply check: \|y_gpu - y_host\| / \|y_host\|` for each call. | Round-off: about 1e-16 to 1e-13. |
| `matvec` | H2 GPU, compression only (`precon 2`, the H2 matvec without a factorization) | Every matvec of the compressed matrix on the GPU is run again on the host. Prints `GPU matvec check: \|y_gpu - y_host\| / \|y_host\|`. | Round-off. |
| `kernel` | H2 GPU and HODLR GPU with a registered GPU evaluator of the entries ([gpu_kernels.md](gpu_kernels.md)) | Compares the evaluator's entries with the application's (its CPU entry callback) on up to 32 x 32 entries per block. H2: at the start of every GPU level, three blocks per rank: a box with itself (the singular self terms), with a neighbor, and with the last local box; prints `[gpu] level L kernel check (rank R, form): ...`, the form being `entry evaluator`, `entry evaluator (source, NVRTC)`, `block evaluator` or `entry-list evaluator`. HODLR: before the GPU construction, a dense leaf with itself and its rows with the last columns of the matrix, the largest over the ranks; prints `HODLR GPU kernel check: ...`. Each value is the largest difference relative to the block's largest entry. | Round-off; a large value means the evaluator does not compute the application's kernel (parameters, ids off by one, rows and columns swapped, coordinates read without `BPACK_GPU_COORDINATES`). |
| `replica` | H2 GPU, replicated CA levels (`CA_level`, `H2_CA_owner_component 0`) | After each CA level on the GPU, compares every ghost copy of a box (each rank eliminates its ghosts again) with its owner's copy, bitwise, by hashes of the skeleton, T, LU(X_RR), pivots and the other factors. Prints `[gpu] level L CA replica check: N ghost copies, M differ from their owners' (bitwise)` and, per rank, the first copies that differ and in which factors. Also keeps the host copies of the factors, which are compared. | `0 differ`. Any difference breaks the CA solve's consistency. |
| `parity` | H2 CPU, component-owner CA (`H2_CA_owner_component 3`) | Not a check of results but a comparison mode: the component-owner schedule uses the raw Morton parity of the boxes as its waves (the order of the replicated mode 0) instead of its partition with the fewest waves, so a mode-3 run can be compared box by box with a mode-0 run. Must be set alike on every rank. | Mode 3 then matches mode 0's elimination order; it runs slower than the default mode 3. |
| `hodlr` | HODLR GPU (`HODLR_use_gpu`) | Runs the CPU routine next to each GPU step and prints the difference: each construction level and the dense leaves; the upload of the factors; their round trip back to the host (compared, then taken, so the CPU checks that follow use the fetched copies); the entries of the construction check, extracted again with the host factors overwritten by a sentinel so that every row must come from the GPUs; every multiply; the factorization (log-determinant); and the solves. It runs the CPU construction and factorization as well, and keeps host copies of the factors built on the GPU (which are otherwise left on the GPU only). | Round-off, except the construction levels where the GPU and the CPU draw different random numbers: there the difference is of the order of the tolerance, and the log says so (use `hodlr-exact` for round-off there too). |
| `hodlr-transpose` | HODLR GPU | Everything `hodlr` does, plus the transposed products of each multiply. | As `hodlr`. |
| `hodlr-exact` | HODLR GPU | Not a check by itself but a comparison mode: the GPU construction splits a block shared by several ranks exactly as the CPU does (the block's own group of ranks, the same random numbers, one piece per rank whatever `HODLR_gpu_pieces`) instead of spreading the pieces over the node's group. With `hodlr`, the construction levels then compare at round-off; it can also be used to compare a GPU run with a separate CPU run. | Same ranks as the CPU; a slower construction of the symmetric HODLR. |

## BPACK_TRACE

A comma-separated list of the values below, for example
`BPACK_TRACE=phase,wave`. Traces add output and some synchronization; they
do not change results.

| Value | Applies to | Output |
|---|---|---|
| `wave` | H2 GPU, factorization | Per GPU level (Color or CA), a host timeline of its waves printed by the level's first rank at the level end: the level start, each wave's host phases (sketch plan and launch, sketch wait, ID store, elimination plan and launch, owner pass, host store) and device times, the time between waves, and the level finish; also the spread of the ranks' wave totals and the wave tables of the slowest rank. |
| `exchange` | H2 GPU, Color levels on several ranks | One `[exchange]` line per rank, level and wave: the time of the message sizes and of the payloads, the GB sent and received, and the messages that did not fit in the exchange arena. |
| `phase` | H2, factorization | Per level, the elimination's phase breakdown (`[phase] level L rank-root` and `max-ranks` lines: waves, box, candidate, owner, mirror, finish and staged times, per-thread sketch, ID, factor and near-field times) and, on the GPU, the device box path's statistics. |
| `bk` | H2 CPU, Bunch-Kaufman (LDL^T) factorizations | Diagnostics of the Bunch-Kaufman factorizations (symmetry \|\|A - A^T\|\|_F / \|\|A\|\|_F, diagonal range, NaN or Inf), with matrix dumps when one fails. |
| `memory` | H2, factorization | Host memory by phase: `H2_MEM_DIAG` lines per rank with the resident set size, its high-water mark, and the bytes of the factorization's data structures (local and halo metadata, factors, near-field, far-field and deferred blocks, H2 blocks, assisting boxes, solve data, pending updates, scratch, staged parent data, communication). |
| `id:<prefix>` | H2 | One line per compressed box describing its ID target, appended by each rank to `<prefix>.rank<r>`. |
| `hodlr-qr` | HODLR GPU | Each block and chunk of the BACA recompression QR of the construction (sizes, ranks, pointers). |
| `sync` | HODLR GPU | Synchronizes the device after each step of the BACA recompression, so that an asynchronous device fault is reported by the kernel that caused it. |

## Variables of earlier versions

These are no longer read; rank 0 names any of them that is set.

| Earlier variable | Now |
|---|---|
| `H2_GPU_HEAP_FRACTION` | `BPACK_GPU_HEAP_FRACTION` |
| `H2_GPU_HEAP_GB` | Removed: ranks that share a GPU split its memory by themselves. |
| `H2_GPU_EXCHANGE_MB`, `HODLR_GPU_EXCHANGE_MB` | `BPACK_GPU_EXCHANGE_MB` (one arena for H2 and HODLR) |
| `H2_GPU_AWARE_MPI`, `H2_GPU_DEVICE_EXCHANGE=0` | `BPACK_GPU_AWARE_MPI` (`1`, `0`) |
| `FMM_MAX_CPUS_PER_NODE` | `BPACK_MAX_CPUS_PER_NODE` |
| `H2_GPU_SOLVE_CHECK`, `H2_GPU_MATVEC_CHECK`, `H2_GPU_KERNEL_CHECK`, `H2_CA_REPLICA_CHECK` | `BPACK_CHECK=solve`, `matvec`, `kernel`, `replica` |
| `H2_OWNER_PARITY_WAVES=1` | `BPACK_CHECK=parity` |
| `HODLR_GPU_CHECK=1`, `=2` | `BPACK_CHECK=hodlr`, `hodlr-transpose` |
| `HODLR_GPU_SPLIT=0` | `BPACK_CHECK=hodlr-exact` |
| `H2_GPU_WAVE_TRACE`, `H2_GPU_EXCHANGE_TRACE`, `H2_PHASE_REPORT`, `H2_BK_DIAGNOSTICS`, `FMM_MEMORY_DIAGNOSTICS`, `H2_ID_TRACE=<prefix>` | `BPACK_TRACE=wave`, `exchange`, `phase`, `bk`, `memory`, `id:<prefix>` |
| `HODLR_GPU_DEBUG=1`, `=2` | `BPACK_TRACE=hodlr-qr`, `hodlr-qr,sync` |
| `H2_GPU_SOLVE`, `H2_GPU_SKETCH`, `H2_GPU_MATVEC`, `H2_GPU_MATVEC_SYMMETRIC`, `H2_GPU_CA_SOLVE_HALO`, `H2_GPU_CA_DEVICE_HALO`, `H2_GPU_KEEP_OPERATORS`, `H2_GPU_SOLVE_KEEP`, `HODLR_GPU_KEEP_FACTORS`, `HODLR_GPU_DEFER_HOST` | Removed; always on (the HODLR factors' host copies are still filled with `BPACK_CHECK=hodlr` or `ErrFillFull=1`). |
| `H2_GPU_SOLVE_DIRECT`, `H2_GPU_MATVEC_DIRECT` | Removed; messages go through device memory whenever MPI is GPU-aware. |
| `H2_GPU_SOLVE_KEEP_FRACTION`, `H2_GPU_MATVEC_KEEP_FRACTION`, `H2_GPU_COPY_THREADS` | Removed; fixed (0.8 and 0.9 of the pool; the copy threads by the rank's OpenMP threads). |
| `H2_GPU_WARMUP` | Removed: the GPU warm-ups always run, and their time is reported apart from the factor, compression and construction times. H2: the batched kernels and a first message with each rank the levels exchange with, before the first level of the process's first factorization, and the GPU evaluator before its first use (`[gpu] warm-up: kernels, exchanges, evaluator`). HODLR: the kernels of the factorization, the routines and cuSOLVER handles of the GPU construction and the GPU evaluator, before the construction (`HODLR GPU warm-up: kernels, evaluator`). |
