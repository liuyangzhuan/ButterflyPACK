# Partial H2 factorization (`H2_partial`)

The H2 code (format 7) can factor part of a matrix: eliminate a leading set of
points and keep the Schur complement on the others, as a multifrontal solver
does with a front. The object then serves the front's partial solves and the
Schur complement's entries and products to its parent, without ever forming
the Schur complement as a dense matrix.

This page describes the options, the C calls, what runs on the GPU, and the
test driver. The design notes and measurements are in
[`h2_parallel/unstructured/PARTIAL_FACTORIZATION_PLAN.md`](../h2_parallel/unstructured/PARTIAL_FACTORIZATION_PLAN.md).

## What it computes

Let A be the N x N matrix with its points in input order (the global input
ids, 0-based), and n11 = `H2_partial`. The points with id < n11 form F11, the
others F22 (N22 = N - n11 points):

```
A = [ F11  F12 ]      CB = F22 - F21 F11^-1 F12   (N22 x N22)
    [ F21  F22 ]
```

`c_bpack_factor` eliminates the F11 points with the RS-S factorization of a
full factorization. The F22 points are compressed level by level with the
others but never eliminated; what the eliminations do to them is kept as
records of the F22 boxes, and the F22 points left at the top form a dense root
block. The Schur complement CB is these records and the root block: its
entries and its products with vectors are computed from them on request. Its
accuracy is that of the factorization (`tol_comp`).

## Options

| Option | Value | |
|---|---|---|
| `H2_partial` | n11, in [0, N] | The number of leading input points that form F11. 0: off. |
| `format` | 7 | The H2 format. |
| `H2_unstructured` | 1 | |
| `sym` | 0 | A nonsymmetric (LU) factorization, also for a symmetric matrix. |
| `H2_use_sketch`, `H2_lazy_schur` | 2, 2 | |
| `H2_XRR_factor` | 1 | LU with partial pivoting of the X_RR blocks. |
| `precon` | 1 or 3 | A factorization (not the compression-only `precon 2`). |
| `H2_ID_proxy` | 0, 1 or 2 | Optional: far-field rows of the IDs beyond their neighbourhood, as for a full factorization: 1 a proxy surface (`H2_ID_proxy_points` points), 2 adaptive row sampling (`BACA_Batch` rows a node). Both the F11 IDs and the F22 parts' IDs take them, on the CPU and the GPU; mode 2 selects its rows against the columns each ID is over (the F11 columns, or the F22 ones). Modes 1 and 2 need the points' global order on every rank: an H2 set up with `c_bpack_construct_init` (a rank's rows are its range of the tree order), not with `c_bpack_construct_init_distributed64`. |
| `CA_level` | at least the number of levels | Color levels only (as for every `H2_unstructured 1` factorization). |
| `H2_compress_factors` | 0 or 1 | Optional. 1: the F22 records are kept compressed, lossy with an error bound of 0.02 `tol_comp` (CPU and GPU builds). The F11 factors of the partial solves stay uncompressed. |
| `H2_partial_extract` | `auto` (0), `device` (1), `host` (2) | Optional, GPU builds: where the Schur complement's entries are computed ([below](#entries-h2_partial_extract)). |
| `H2_use_gpu` | 0 or 1 | Optional: the factorization, the partial solves and the Schur complement's products and entries on the GPU. |
| `elem_extract` | 0, 1 or 2 | As for a full factorization. With 1 (a block callback all ranks call together, as STRUMPACK's fronts give their entries) the F22 step asks for its blocks in rounds too: on the CPU a box that lacks one runs again after the round, on the GPU the device fetches it on the spot. |

The H2 setup (`c_bpack_construct_init` and its variants) checks these and
throws with the option's name when one does not hold. They are set like the
other options, with `c_bpack_set_I_option(&option, "H2_partial", n11)` or on
the command line (`--H2_partial n11`, `--H2_partial_extract auto|device|host`).

## Calling sequence

```cpp
// format 7, H2_unstructured 1, sym 0, H2_use_sketch 2, H2_lazy_schur 2, H2_XRR_factor 1
d_c_bpack_set_I_option(&option, "H2_partial", n11);
// ... construct the H2 as usual (c_bpack_construct_init[_distributed64], c_bpack_construct_element_compute[...]) ...
d_c_bpack_factor(&bmat, &option, &stats, &ptree, &msh);  // eliminates F11 only

int N22 = N - n11;
// the front's partial solves: b, x local rows (as for c_bpack_solve); out2, x2 N22 x nrhs on every rank
d_c_bpack_partial_forward(out2, b, &nloc, &N22, &nrhs, &bmat, &option, &stats, &ptree);  // out2 = b2 - F21 F11^-1 b1
// ... the parent solves for x2 ...
d_c_bpack_partial_backward(x, x2, &nloc, &N22, &nrhs, &bmat, &option, &stats, &ptree);   // x = [F11^-1 (b1 - F12 x2); x2]

// the Schur complement
d_c_bpack_partial_schur_mult(y2, x2, &N22, &ncol, &bmat, &option);           // y2 = CB x2
d_c_bpack_partial_schur_extract(blk, rows, &nr, cols, &nc, &bmat, &option);  // blk = CB(rows, cols), 0-based
// or c_bpack_extractelement on the same object: blocks of CB, 1-based F22 indices

d_c_bpack_partial_schur_release(&bmat, &option);  // when the parent has read CB; the partial solves stay
```

The prefix selects the data type as for the rest of the C interface (`d_` real
double, `z_` complex double). [`EXAMPLE/Laplace3D_H2Test_Driver.cpp`](../EXAMPLE/Laplace3D_H2Test_Driver.cpp)
(`partial_checks`, `partial_extract_bench`, `partial_release_check`) calls
each of them.

## The calls

All are collective over the matrix's communicator. Where a call takes F22
vectors, indices or columns, every rank passes the same ones. The F22 index of
an F22 point is its input id minus n11.

| Call | Computes | Arguments |
|---|---|---|
| `c_bpack_partial_forward(out2, b, Nloc, N22, Nrhs, ...)` | out2 = b2 - F21 F11^-1 b1 | `b`: Nloc x Nrhs, the rank's local rows as for `c_bpack_solve` on the same handle (a `distributed64` handle: the rank's input points, Nloc its input-local size). `out2`: N22 x Nrhs, column-major, the same on every rank. |
| `c_bpack_partial_backward(x, x2, Nloc, N22, Nrhs, ...)` | x = [F11^-1 (b1 - F12 x2); x2] | `x2`: N22 x Nrhs, the same on every rank. `x`: the local rows, laid out as `b`; its F22 rows are x2's. b1 is the last forward call's b, so a backward call needs a forward call with the same Nrhs before it. |
| `c_bpack_partial_schur_mult(y2, x2, N22, Ncol, ...)` | y2 = CB x2 | `x2`, `y2`: N22 x Ncol, column-major, the same on every rank. |
| `c_bpack_partial_schur_extract(out, rows, nrows, cols, ncols, ...)` | out = CB(rows, cols) | 0-based F22 indices; `out` nrows x ncols, column-major, the same on every rank. |
| `c_bpack_extractelement`, `c_bpack_extractelement_bridge` | blocks of CB | The usual arguments, with CB as an N22 x N22 matrix: 1-based F22 indices in `allrows` and `allcols`, and a 1 x 1 process grid for each intersection (its `pmaps` row: 1, 1, the rank that receives it). Every intersection of a call is computed in one pass, so a parent should ask for all the blocks it reads in one call. |
| `c_bpack_partial_schur_release(...)` | frees CB | The records and the root block, nearly all of the object's memory (host and device). The partial solves keep working; Schur complement calls after it throw. |
| `c_bpack_logdet` | log det F11 | |

`c_bpack_solve` and `c_bpack_inv_mult` throw on a partial object (it is
applied by the forward and backward calls), and the partial calls throw on an
object without a partial factorization.

## Memory

The records of the F22 boxes and the root block are most of the object's
memory. On a 64^3 grid with the first 4 planes in F11 (n11 = 16384, one GPU,
`tol_comp` 1e-6) the records took 31.7 GB, or 5.9 GB with
`H2_compress_factors 1`, and the root block 2.1 GB. Release them with
`c_bpack_partial_schur_release` as soon as the parent has read the Schur
complement. The factorization prints the records' size per level on rank 0
(`verbosity` 1), here with `H2_compress_factors 1` on the GPU:

```
  F22 phase: 512 boxes, 156886 points sparsified, 85639 kept; records 3.57 GB (17.59 GB uncompressed), 3.57 GB of them kept on the device only
```

## GPU

With `H2_use_gpu 1` the factorization runs on the GPU as a full one does, and
the object's calls use the GPU when their data fit in the device memory pool
(`BPACK_GPU_HEAP_FRACTION`, [environment_variables.md](environment_variables.md)).
Each falls back to the host on its own when they do not.

**The partial solves.** The forward and backward sweeps run on the device
when the device solve does (the factors kept on the device or uploaded at the
first solve, as for `c_bpack_solve`). Rank 0 prints
`GPU partial forward sweeps: ...` and `GPU partial backward sweeps: ...`.

**Where the records live.** The factorization keeps the F22 records on the
device, without a host copy, while they fit under 0.2 of the pool; the others
go to the host (the `kept on the device only` part of the F22 phase line). At
the first Schur complement call that runs on the device, the ranks decide
together whether all the records, the root block and the working space fit
in the free pool. If they do, the host records are uploaded and everything
stays there until the release:

```
GPU partial CB: the records (compressed) and the root block on the device, 4.69 GB per rank at most (records 0.37, root 2.13), set up in 0.25 s (upload 0.19 s on rank 0); 73 columns a pass
```

(the records counted there are those uploaded; the ones the factorization
kept are read in place).

If they do not, the records kept on the device are copied to the host and
freed, and the products and entries run on the host:

```
GPU partial CB: host product (not enough device memory on some rank; records 24.98 GB, root block 2.13 GB, need 27.44 GB against 24.74 GB on rank 0)
```

With several H2 matrices in one process, the least recently used ones'
device data is released when a later build needs the room
([environment_variables.md](environment_variables.md#bpack_gpu_heap_fraction));
a partial object's device-only records are copied to the host first.

**Products.** `c_bpack_partial_schur_mult` runs on the device when the records
are there, in passes of as many columns as the working space holds.

### Entries (`H2_partial_extract`)

| Value | The entries of a call are computed |
|---|---|
| `auto` (0, default) | on the device when the records are there and either the factorization kept some of them (they have no host copy) or the call asks for at least 65536 entries in all; else on the host. |
| `device` (1) | on the device whenever the records are there (uploaded at the first call if needed). |
| `host` (2) | on the host. Records kept on the device are first copied to the host (once). |

The threshold keeps small calls, where the device's fixed cost per call
would dominate, on the host; it has not yet been calibrated against a
multifrontal solver's calls. The device gains the most with
`H2_compress_factors 1`, as it decompresses the records it reads in batches:
in the 64^3 case above, parent-like rounds of 2.4 M entries took 0.099 s on
the device against 0.43 s on the host.

On the device, the pair terms are added into a device copy of the call's
result, which comes to the host once at the end (when it fits in the working
space; else each group's terms come to the host and the host adds them in).
A point asked for at several positions of a block gets its sums at its first
position, copied to the others on the host. In STRUMPACK's H2-only fronts
(Poisson 200^3, 16 ranks, the same calls timed on both under
`BPACK_CHECK=solve`) this took the device from 287 s to 80 s summed over the
calls, against 162 s on the host.

### Messages between ranks

On several ranks, the entries' blocks that other ranks read go from device
memory to device memory when MPI is GPU-aware (`BPACK_GPU_AWARE_MPI`, through
the exchange arena of `BPACK_GPU_EXCHANGE_MB`); only their index lists go
through the host. Without GPU-aware MPI their values are staged through host
memory. In the 64^3 case above on 4 ranks and 4 GPUs, the parent-like rounds
took 0.062 s (0.172 s staged through the host, 0.268 s before the values left
the host messages).

The products keep their vectors on the host (only the records' and the root
block's products run on the device), so their messages are host MPI. The
partial solves on the device use the device exchanges of the full solve.

### Checks

`BPACK_CHECK=solve` runs every partial solve, Schur complement product and
Schur complement entry call that ran on the device again on the host, and
prints on rank 0:

```
GPU partial forward check: |gpu - host| / |host| = 5.476e-16
GPU partial backward check: |gpu - host| / |host| = 3.334e-16
GPU partial CB product check: |y_gpu - y_host| / |y_host| = 2.653e-15
GPU partial CB entries check: |gpu - host| / |host| = 1.209e-16 (2438 entries in 4 blocks; device 0.0029 s, host 0.0052 s)
```

Round-off is expected (about 1e-17 to 1e-14), also with compressed records:
both sides decompress the same bytes.

## Test driver

`EXAMPLE/Laplace3D_H2Test_Driver.cpp` (`claplace3d_h2`) factors a 3D Laplace
kernel on an n^3 grid; with `--H2_partial n11` the first n11 points (the
grid's first planes) form F11. It then calls the partial solves, the product,
the entries both ways (comparing them) and the release, and prints their
times.

| Option | |
|---|---|
| `--H2_partial <n11>` | The partial factorization. With `--distributed64 1` a rank's vector rows are its range of the input order, with `--distributed64 0` (which `--H2_ID_proxy 1` and `2` need) its range of the tree order. |
| `--H2_partial_extract <auto\|device\|host>` | As the option. |
| `--partial_check 1` | Compares the partial solves, the product and blocks of the Schur complement with a dense LU on rank 0 (N <= 40000). |
| `--extract_ranks <P>` | The partial H2 takes its entries through the collective block callback (`elem_extract 1`) from a compression-only H2 of the kernel on the first P ranks, as from a STRUMPACK front (`--distributed64 0`). The checks then also carry the source's compression error. |
| `--partial_extract_bench <R>` | R rounds of Schur complement blocks as a parent front reads them through `c_bpack_extractelement` (64 x 64 near-field blocks and 16 x 4096 panels), timed. |

```
srun -n 4 ./claplace3d_h2 --grid-size 20 --tol-comp 1e-6 --Nmin_leaf 64 --elem_extract 2 \
  --H2_unstructured 1 --sym 0 --CA_level 10000 --H2_XRR_factor 1 --h2_use_sketch 2 --h2_lazy_schur 2 \
  --distributed64 1 --H2_use_gpu 1 --H2_partial 1200 --partial_check 1
```

## Limitations

- An MPI build (format 7 needs it).
- `sym 0` only: there is no symmetric (LDL^T) partial factorization.
- Color levels only; no CA levels.
- Complex-valued (`z_`) fronts compile but have not been tested.
- On several ranks, the Schur complement's products exchange their vectors
  through host memory ([above](#messages-between-ranks)).
