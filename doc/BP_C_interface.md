# C interface for standalone BP matrices

The `c_bp_*` functions follow the existing `c_bf_*` calling convention. Precision
preprocessing exports `d_c_bp_*`, `z_c_bp_*`, `s_c_bp_*`, and `c_c_bp_*` for double
real, double complex, single real, and single complex respectively.

`c_bp_construct_init` has the same arguments as `c_bf_construct_init`. It consumes
**existing row and column meshes** and delegates to `BP_Construct_Init_from_mshrc`.
It does not expose the coordinate/permutation signature of the Fortran
`BP_Construct_Init`; prepare those meshes first, as in the BF examples.

A typical double-precision call sequence is:

```cpp
F2Cptr bp = nullptr, bp_mesh = nullptr, bp_kernel = nullptr, bp_stats = nullptr;
d_c_bpack_createstats(&bp_stats);
d_c_bpack_set_I_option(&option, "cpp", 1);
d_c_bpack_set_I_option(&option, "format", 2); // H-matrix representation

d_c_bp_construct_init(&M, &N, &Mloc, &Nloc, nnsr, nnsc, &row_mesh, &col_mesh,
    &bp, &option, &bp_stats, &bp_mesh, &bp_kernel, &ptree,
    distance_callback, nearfar_callback, context);
d_c_bp_construct_element_compute(&bp, &option, &bp_stats, &bp_mesh, &bp_kernel,
    &ptree, entry_callback, block_entry_callback, context);

char trans = 'N';
d_c_bp_mult(&trans, x, y, &Nloc, &Mloc, &nrhs, &bp, &option, &bp_stats, &ptree);

d_c_bp_delete(&bp);
d_c_bpack_deletekernelquant(&bp_kernel);
d_c_bpack_deletemesh(&bp_mesh);
d_c_bpack_deletestats(&bp_stats);
```

## Conventions and limitations

- **Handles:** BP, BF, and BPACK handles are not interchangeable. BP handles also
  own persistent callback-pointer storage. Only `c_bp_delete` releases them;
  deleting a null BP handle is harmless. The output mesh/kernel are separate
  allocations. Input meshes, options, statistics, process tree, and application
  context remain caller-owned. Do not invoke callbacks through the BP kernel
  after deleting the BP handle.
- **Formats:** the standalone BP backend supports `format=1` (one BF/LR block),
  `2` (H-matrix), `3` (butterfly-plus), and `5` (BLR). H2 and tensor formats are not
  supported here. Keep `format` unchanged between initialization and compression.
- **Callbacks:** set `cpp=1`. Use `elem_extract=0` for the scalar entry callback or
  `elem_extract=1` for the block-extraction callback, with the same block layout
  contract as BF. Unused callbacks may be null. Keep the application context live
  for all calls that use it. With `nogeo=2`, supply distance and near/far callbacks;
  the near/far callback receives cluster IDs, not element IDs.
- **Ordering:** entry and distance callbacks receive signed, **1-based reordered**
  indices: positive indices identify the row mesh, negative indices identify the
  column mesh. Either argument can identify the row. Convert each axis through
  its own permutation if the application stores entries in original order.
  `c_bf_new2old_row/col` can still be used on the corresponding input meshes.
- **Neighbors:** for `nogeo=3/4` with `knn>0`, supply both arrays. They are laid out
  as C `[M][knn]` and `[N][knn]`, equivalent to Fortran `(knn,M)` and `(knn,N)`.
  Values are 1-based reordered indices on the opposite axis, with zero denoting
  a missing neighbor. These arrays may be null when unused. With `nogeo=0/4`,
  the input meshes must still contain coordinates when BP initialization runs.
- **Multiplication:** calls are collective over the process tree. Buffers contain
  local, reordered vectors in column-major order. For `N`, use input size `Nloc`
  and output size `Mloc`; swap them for `T` and `C`. Input and output must not
  overlap. `C` uses conjugation around the backend's transpose operation and
  leaves the input unchanged. As with `c_bf_mult`, the result is divided by
  `option%scale_factor`; keep that value consistent with construction.
- **Statistics:** allocate fresh statistics for a BP construction instead of
  reusing the statistics populated while building its row/column meshes.

## Regression example

`EXAMPLE/InterfaceBPTest.cpp` uses `M=512`, `N=768`, `Nmin_leaf=32`, and three
right-hand sides. It requires nonzero compressed storage, a reduced rank, and
total matrix storage below the dense equivalent for every case, and prints
those measurements. It checks products against a dense reference with normal/transpose/conjugate-transpose
operations, all four supported formats, and deletion. It also exercises the
rectangular nearest-neighbor path and `forwardN15flag=2`.

After configuring an MPI build, explicitly build the desired test targets:

```sh
cmake --build build --target dbp_interface_test zbp_interface_test sbp_interface_test cbp_interface_test
srun -n 2 build/EXAMPLE/dbp_interface_test
srun -n 2 build/EXAMPLE/zbp_interface_test
```

The test targets are excluded from the default build. Their definitions select
single/complex precision automatically; the same source can also be compiled
manually with `BP_TEST_SINGLE` and/or `BP_TEST_COMPLEX`.
