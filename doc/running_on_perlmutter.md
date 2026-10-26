# Running on Perlmutter: preloads, OpenMP binding and BLAS

## Do not export an `LD_PRELOAD` that contains an OpenMP library into a wrapper

A launch like

```bash
export OMP_PROC_BIND=spread OMP_PLACES=cores OMP_NUM_THREADS=16
export LD_PRELOAD="$TCMALLOC/libtcmalloc.so:$BUILD/libbutterflypack.so.5.0.0"   # wrong
srun -n 8 -c 32 --cpu-bind=cores /usr/bin/time ./cvie3d_t_so ...
```

runs **every OpenMP thread of each rank on one core**. `libbutterflypack` links
libgomp, so the preload loads libgomp into `/usr/bin/time` (or into `bash` in a
`bash -c "exec ..."` launcher). With `OMP_PROC_BIND` set, libgomp's constructor
pins that wrapper to its first place: one core, or one hardware thread with
`OMP_PLACES=threads`. The application inherits the wrapper's CPU mask, its own
libgomp builds the place list from that mask, and `omp_get_num_places()` is 1.
The run is correct and only slow. Nothing in the output says so.

Measured with a 10-line test program (16 cores per rank):

| launch | `omp_get_num_places()` |
|---|---|
| `srun ... ./app` | 16 |
| `srun ... /usr/bin/time ./app` | 16 |
| `LD_PRELOAD=...libbutterflypack... srun ... /usr/bin/time ./app` | **1** |
| `srun ... /usr/bin/time env LD_PRELOAD=... ./app` | 16 |

What it cost (2026-10, tensor VIE with H2 block-diagonal preconditioner, CPU
nodes): a 640^3 run took 1112 s instead of 291 s. Each H2 block took 51 s
instead of 14 s, and the tensor matvec 17.8 s instead of 5.3 s. Several days of
CPU timings were wrong by 3-5x before this was found. Under that placement an
OpenMP-threaded OpenBLAS looked 25x *slower*, because its 16 spin-waiting
threads shared one core. With correct placement it is faster.

**Do this instead:** keep the preload out of the environment and hand it to the
application, after any wrapper:

```bash
PRELOAD="$TCMALLOC/libtcmalloc.so:$BUILD/libbutterflypack.so.5.0.0"
srun -n 8 -c 32 --cpu-bind=cores /usr/bin/time env LD_PRELOAD="$PRELOAD" ./cvie3d_t_so ...
```

Put profilers and debuggers before `env` (`perf record ... env LD_PRELOAD=... ./app`;
in gdb, `set exec-wrapper env LD_PRELOAD=...`). The same applies to any
preloaded library that links an OpenMP runtime, such as an OpenMP build of
OpenBLAS. A preload of tcmalloc and the sequential OpenBLAS alone does not link
libgomp and is safe, as in `example_scripts/runit_3dlaplace_perlmutter_h2.sh`.

**Check when a CPU timing looks off:** print `omp_get_num_places()` in the
application's environment (or from a 10-line program launched the same way). It
must equal the cores per rank. `example_scripts/runit_3dvie_tensor_bdiag_perlmutter_cpu.sh`
launches correctly.

## BLAS for the H2 code

- Build against the sequential OpenBLAS (CA0 needs it).
- The Color H2 factorization and solve also run correctly on an OpenMP-threaded
  OpenBLAS. Calls inside the code's parallel regions run on one thread and the
  rest on all threads. With correct placement, a 40^3 tensor-VIE H2 block was 17%
  faster with it (13.8 -> 11.5 s), with the same results.
- Preload the threaded library at run time over the sequential one (same soname
  `libopenblas.so.0`); there is no need to rebuild. `OPENBLAS_NUM_THREADS` does
  not cap the OpenMP build of OpenBLAS 0.3.30; it follows `omp_get_max_threads()`.
- A region guarded by `#pragma omp parallel if (n > 1)` is still a region when
  n = 1: an inactive one, inside which `omp_in_parallel()` is false. A threaded
  BLAS called there opens a nested team, and libgomp starts fresh threads for
  every nested team. The H2 code uses `fmm::omp_parallel_if`
  (`h2_parallel/core/omp_region.hpp`) for such regions; a single item runs with no
  region.
