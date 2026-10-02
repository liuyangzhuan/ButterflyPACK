# GPU evaluators: the matrix's entries on the GPU

The GPU backends of ButterflyPACK (H2, format 7, with the option `H2_use_gpu`;
HODLR, format 1, with `HODLR_use_gpu`) compute the entries of the matrix on
the device. They need the device counterpart of the application's CPU entry
callback: a **GPU evaluator**, which the application writes for its kernel
and registers with the matrix. This page describes the three forms an
evaluator can take, how to write and register each from C++, C, Fortran and
Python, and how to check one.

Without an evaluator the backends still run, with less on the GPU: the H2
GPU factorization runs only its Schur-update pass on the device, and the
HODLR GPU backend evaluates the entries on the host.

The interface is in `GPU_INTERFACE/` (installed with the library headers):

| Header | For |
|---|---|
| `bpack_gpu.h` | The C interface: the flags, the block and entry-list types, the registrations. Included by `dBPACK_wrapper.h` and `zBPACK_wrapper.h`. |
| `bpack_gpu_entry.cuh` | Entry evaluators from C++: `bpack::gpu::set_entry_evaluator`. Include it in a `.cu` file compiled by `nvcc`. |
| `bpack_gpu_kernels.cuh` | The library's evaluation kernels, which the two routes above instantiate for the application's evaluator. Not used directly. |

`EXAMPLE/gpu/` has one example of each form for the C++ drivers, and
`EXAMPLE/user_block_funcs_george_gpu.py` the two Python routes.

## Conventions

All forms compute the same thing: the entry K(i, j) of the matrix, for the
points with **0-based global ids** i (row) and j (column), in the
application's ordering. These are the 1-based ids of the CPU callbacks minus
one; ButterflyPACK's internal (tree) order never appears.

The **coordinates** of a point are those the application gave to
`c_bpack_construct_init` (`dim` values per point). An evaluator that reads
them says so with `BPACK_GPU_COORDINATES`; without the flag the library may
pass none.

The value type is that of the library the matrix belongs to: `double` for
the `d_` library, double complex for the `z_` library (`bpack::gpu::dcomplex`
in device code: two doubles `re`, `im`, the layout of C's `double _Complex`,
C++'s `std::complex<double>` and Fortran's `complex(8)`). Only these two
libraries have GPU backends; the `s_` and `c_` libraries accept the
registrations and ignore them.

An evaluator is registered after `c_bpack_construct_init` (which creates the
matrix) and before the construction (`c_bpack_construct_element_compute`).
A matrix has one evaluator; registering again replaces it.

### Flags

The registrations take an or-ed combination of:

| Flag | Value | Meaning |
|---|---|---|
| `BPACK_GPU_SYMMETRIC` | 1 | K(i, j) == K(j, i). **Required by the H2 format**, which factors symmetric matrices: H2 refuses an evaluator without it (and says so) and runs as without one. HODLR does not need it. |
| `BPACK_GPU_COORDINATES` | 2 | The evaluator reads the points' coordinates. |
| `BPACK_GPU_HOST_IDS` | 4 | Block evaluators only: the evaluator reads the ids of each block on the host (`row_ids_host`, `col_ids_host`). The library then waits for their copy to the host before each batch, which costs a synchronization per batch. |

## Which form?

| Form | What the application writes | From | Speed |
|---|---|---|---|
| **Entry evaluator** | A device function returning one entry K(i, j). | C++ (CUDA), compiled by `nvcc` with the application; or its CUDA source text, compiled by the library at run time with NVRTC (any language, e.g. Python). | The fastest: the library inlines it in its own kernels, including the fused sketch of the H2 factorization, so the entries are never stored. |
| **Block evaluator** | A function that fills a batch of blocks with its own launches. | C, C++, Fortran (any language that can launch GPU work). | For kernels whose entries share work within a block (a quadrature shared by neighboring basis functions, a table lookup per row). The library gathers each block's ids (and coordinates) on the device. |
| **Entry-list evaluator** | A function that fills the values of a list of (row, column) pairs. | Any; from Python, CuPy array code. | A block evaluator in flat form, for vectorized array code. The library flattens the blocks into id lists and scatters the values back. |

Start with an entry evaluator. Use a block evaluator when computing the
entries one at a time repeats expensive work, and an entry-list evaluator
when the kernel is easiest to write as array code.

With the block and entry-list forms, the H2 factorization's sketch first
evaluates its kernel rows into scratch and then sketches them, where an
entry evaluator is fused into the sketch.

## Entry evaluators from C++

Write a trivially copyable functor with a `value_type` and a `__device__`
call operator, in a `.cu` file:

```cpp
#include "bpack_gpu_entry.cuh"

// K(i, j) = scale / |x - y|, and `diagonal` when i == j
struct LaplaceEntry {
    using value_type = double;     // or bpack::gpu::dcomplex for the z_ library
    double scale;
    double diagonal;

    __device__ double operator()(const double* x, int64_t i, const double* y, int64_t j) const {
        if (i == j) return diagonal;
        const double dx = x[0] - y[0], dy = x[1] - y[1], dz = x[2] - y[2];
        return scale / sqrt(dx * dx + dy * dy + dz * dz);
    }
};

void register_gpu_entries(void** bmat, double scale, double diagonal) {
    bpack::gpu::set_entry_evaluator(bmat, LaplaceEntry{scale, diagonal},
                                    BPACK_GPU_SYMMETRIC | BPACK_GPU_COORDINATES);
}
```

`x` and `y` are the coordinates of points `i` and `j` (`x[0 .. dim - 1]`).
`bmat` is the matrix handle of `c_bpack_construct_init` (`F2Cptr*`).

- **Data.** The functor is copied by the library at registration and passed
  by value to its kernels: it holds parameters and device pointers (at most
  4000 bytes), never host pointers. Tables the evaluator reads (coefficients
  per point, a mesh) are device arrays the application allocates and keeps
  until it deletes the matrix (`EXAMPLE/gpu/vie3d_gpu.cu` reads a
  coefficient per column point).
- **What it is called with.** Any ids of the matrix, in any order, from many
  threads at once, and both (i, j) and (j, i) for symmetric kernels. It must
  not synchronize or allocate.
- **Building.** Compile the `.cu` file with `nvcc` (C++17) with
  `GPU_INTERFACE/` on the include path, for the GPU architecture the library
  targets, and link the shared CUDA runtime (`libcudart.so`, as the library
  does; with CMake, the target property `CUDA_RUNTIME_LIBRARY Shared`).
  `EXAMPLE/CMakeLists.txt` (`bpack_example_gpu`) does this for the examples.
- **Versions.** The call instantiates the library's kernels in the
  application, so the application must be rebuilt with each new version of
  the library. The headers carry `BPACK_GPU_INTERFACE_VERSION`; a library
  of another version refuses the registration with a message that says so.

Fortran or C applications can keep the functor and the registration in a
small `.cu` file and call it through a C function, as the C++ drivers do
(`EXAMPLE/gpu/gpu_evaluators.h`).

## Entry evaluators from source text (NVRTC)

Any language can give the entry evaluator as CUDA source text, which the
library compiles with NVRTC on each rank at its first use (once per
process) and inlines in its kernels like a C++ one. The text defines

```c
__device__ double bpack_entry(const double* x, long long i, const double* y, long long j,
                              const double* params, int dim);
```

(for the `z_` library returning `bpack_dcomplex`, a struct of two doubles
`re` and `im` with the usual arithmetic operators, constructed as
`bpack_dcomplex(re, im)`), and may define helper `__device__` functions.
`params` is a device copy of the parameters given at registration (any
number of doubles: constants, and tables such as one value per point by
global id); `dim` is the number of coordinates per point. The CUDA math
functions (`exp`, `sqrt`, `sincos`, ...) are available; headers are not.

```c
// C
const char* source =
    "__device__ double bpack_entry(const double* x, long long i, const double* y, long long j,\n"
    "                              const double* params, int dim) {\n"
    "    double q = 0.0;\n"
    "    for (int d = 0; d < dim; ++d) { double t = x[d] - y[d]; q += t * t; }\n"
    "    return params[0] * exp(-0.5 * q / params[1]) + (i == j ? params[2] : 0.0);\n"
    "}\n";
const double params[3] = {amplitude, length2, noise};
const int nparams = 3, flags = BPACK_GPU_SYMMETRIC | BPACK_GPU_COORDINATES;
d_c_bpack_set_gpu_entry_source(&bmat, source, params, &nparams, &flags);
```

A source that does not compile stops the run at the first GPU construction
with NVRTC's messages (the line numbers count from the first line of the
source).

## Block evaluators

A block evaluator fills a batch of blocks:

```c
void my_blocks(const bpack_gpu_block_batch* batch, void* user);
int flags = BPACK_GPU_SYMMETRIC | BPACK_GPU_HOST_IDS;
z_c_bpack_set_gpu_block_evaluator(&bmat, &my_blocks, my_data, &flags);
```

Block `b` of the batch is `batch->blocks[b]` (`bpack_gpu.h`):

| Field | |
|---|---|
| `m`, `n`, `ld` | Its size and leading dimension. |
| `out` | Device memory of the matrix's type: `out[r + c * ld] = K(row_ids[r], col_ids[c])`, `r < m`, `c < n` (column-major). |
| `row_ids`, `col_ids` | The ids of its rows and columns, on the device. |
| `row_ids_host`, `col_ids_host` | The same ids on the host, with `BPACK_GPU_HOST_IDS` (else null). |
| `row_coords`, `col_coords` | The coordinates of its rows and columns on the device, `dim` per point (point `r` at `row_coords[r * dim]`), with `BPACK_GPU_COORDINATES` (else null). |

and the batch carries `count`, the same array in device memory
(`blocks_device`, for a kernel that walks the blocks), the largest `m` and
`n`, `dim`, the CUDA stream, and an allocator.

- **Stream.** Launch everything on `batch->stream` (a `cudaStream_t`) and
  return without waiting: the library orders its next work after it. Do not
  synchronize the device.
- **Scratch.** Take device scratch from `batch->allocate(batch->allocator,
  bytes)` and give it back with `batch->release(batch->allocator, ptr)`
  before returning: the library holds most of the device memory in its own
  pool, and memory released here is reused by later launches on the same
  stream, so it can be released right after the launches that use it. The
  allocator returns NULL when the memory is short: split the batch and try
  again with fewer blocks (as `EXAMPLE/gpu/emsurf_gpu.cu` does).
- **Host data.** Host buffers may be copied to the device with
  `cudaMemcpyAsync` on the stream; pageable host memory can be freed once the
  call returns.
- **Batches.** A batch holds many blocks of different sizes (the library
  bounds the total of their ids); the same block may come again later.

`EXAMPLE/gpu/emsurf_gpu.cu` evaluates the EFIE of RWG edges of a triangle
mesh this way: an entry is a sum over the 4 pairs of triangles of its two
edges, and a triangle pair serves up to 9 edge pairs, so the evaluator finds
each block's distinct triangles on the host (`BPACK_GPU_HOST_IDS`), runs the
quadrature of every triangle pair once into scratch, and combines the pair
sums into the entries.

## Entry-list evaluators

An entry-list evaluator fills `values[e] = K(rows[e], cols[e])` for
`e < count`; the three arrays are on the device, and the work goes on
`stream`:

```c
void my_list(int64_t count, const int64_t* rows, const int64_t* cols, void* values, void* stream, void* user);
int flags = BPACK_GPU_SYMMETRIC;
d_c_bpack_set_gpu_list_evaluator(&bmat, &my_list, my_data, &flags);
```

A call can hold tens of millions of entries.

## From Python

`Py_BPACK_worker.py` registers an evaluator given in the payload of
`bpack_factor` (`Test_python_master.py --gpu entry|list` shows both):

```python
# an entry evaluator as CUDA source text (NVRTC), the faster route
payload["gpu_entry"] = {"source": SOURCE, "params": params, "flags": 3}

# or an entry-list evaluator in CuPy: a function of block_func_module
payload["gpu_list"] = {"func_name": "compute_entries_gpu", "flags": 1}
```

`SOURCE` and `params` are as in [the source-text route](#entry-evaluators-from-source-text-nvrtc).
The workers' options select the GPU backend as for any driver (H2:
`--format 7 --H2_use_gpu 1` and the options of [H2 specifics](#h2-specifics);
HODLR: `--format 1 --HODLR_use_gpu 1` and a BACA construction,
`--reclr_leaf 4` or `5`). The worker registers the evaluator between the two
steps of the initialization (`py_bpack_init`, `py_bpack_compute`), then warms
it up there (`GPU evaluator warm-up: ...`), so its first-use costs count in
the initialization rather than in the factorization.
The entry-list function is called on each rank of the workers as

```python
def compute_entries_gpu(rows, cols, values, meta):
    # rows, cols: CuPy int64 arrays (0-based global ids); values: the CuPy
    # array to fill (float64 or complex128); meta: the payload's
    ...
```

within the library's stream (`cupy.cuda.ExternalStream`), so ordinary CuPy
code is ordered correctly. Its temporaries come from CuPy's memory pool,
outside the library's: process large calls in chunks, and leave the GPU
memory CuPy needs by lowering `BPACK_GPU_HEAP_FRACTION`
([environment variables](environment_variables.md)). An exception in the
function aborts the run. For the best speed, write the same kernel as
source text instead (`EXAMPLE/user_block_funcs_george_gpu.py` has both for a
Gaussian-process kernel).

## H2 specifics

- The evaluator must be symmetric (`BPACK_GPU_SYMMETRIC`).
- The factorization calls it on the levels its GPU box path takes, which
  needs `h2_use_sketch 2`, `H2_XRR_factor 1`, and `h2_lazy_schur 2` on levels
  shared by several ranks (`h2_lazy_schur 1` or 2 on CA levels). With
  `verbosity 1` each such level prints `GPU box path: on`, or `off` and why.
  Other levels run on the host, with only the Schur-update pass on the GPU,
  and do not call the evaluator.
- The interpolative decompositions sample the kernel against training
  points. With training rows from the whole tree (`H2_ID_radius` above 2,
  `H2_ID_proxy 1`, or the adaptive rows of `H2_ID_proxy 2`, which the GPU
  selects with `h2_use_sketch 2`), an evaluator that reads coordinates needs
  the coordinates of every point, which `H2_ID_proxy 1` and, with
  `H2_use_gpu`, `H2_ID_proxy 2` keep; otherwise those levels sample on the
  host.

## First use

An evaluator's first call can be slow: the NVRTC compile of a source-text
evaluator, the loading of the kernels a C++ evaluator instantiates, CuPy
compiling the kernels of an entry-list function (it keeps them in
`~/.cupy/kernel_cache` for later runs), the application's own setup. The
library therefore calls each evaluator once, on a block of a few of the
rank's points and its sketch, before the first level that uses it, and
reports that time apart from the factorization, compression and
construction times:

- H2: `[gpu] warm-up: kernels X s, exchanges Y s, evaluator Z s, before the levels (not in the factor time)`
  (and `[gpu] warm-up: evaluator Z s` before a compression).
- HODLR: `HODLR GPU warm-up: kernels X s, evaluator Z s (not in the construction time)`.

An H2 application can also warm its evaluator up right after registering it,
with `c_bpack_gpu_warm_up` (`bpack_gpu.h`; collective over the matrix's
ranks): the cost then falls there, on the GPU the factorization will use,
and the factorization's warm-up reports `evaluator 0.00 s`. The Python
worker does so.

## Checking an evaluator

`BPACK_CHECK=kernel` ([environment variables](environment_variables.md))
compares the evaluator with the application's CPU entries and prints the
largest difference relative to the largest entry:

- H2: at the start of every GPU level, on up to 32 x 32 entries of three
  blocks: a box with itself (the singular self terms), with a neighbor, and
  with the last local box (`[gpu] level L kernel check (rank R, form): ...`,
  the form being `entry evaluator`, `entry evaluator (source, NVRTC)`,
  `block evaluator` or `entry-list evaluator`).
- HODLR: before the GPU construction, on up to 32 x 32 entries of a dense
  leaf with itself and of its rows with the last columns of the matrix
  (`HODLR GPU kernel check: ...`).

Expect round-off. A large difference means the evaluator does not compute
the application's kernel: wrong parameters, ids off by one, rows and columns
swapped, or coordinates read without `BPACK_GPU_COORDINATES`.
`BPACK_CHECK=hodlr` compares every block of a HODLR GPU construction with
the CPU's.
