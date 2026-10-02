"""GPU evaluators of the Gaussian-process kernel of user_block_funcs_george.py,
for H2_use_gpu / HODLR_use_gpu (Test_python_master.py --gpu entry|list;
doc/gpu_kernels.md):

    K(i, j) = amplitude exp(-sum_d (x_d - y_d)^2 / (2 metric_d)) + yerr_i^2 [i == j]

the kernel george evaluates on the host (a constant times ExpSquaredKernel(metric)).
meta: the payload's, with "coordinates", "yerr", and the kernel's "amplitude"
and "metric" (one value per dimension).

Two routes from Python:
  gpu_entry(meta): an entry evaluator as CUDA source text, which the library
      compiles with NVRTC and inlines in its kernels (the fastest);
  compute_entries_gpu(rows, cols, values, meta): an entry-list evaluator in
      CuPy, called with the ids of a batch of entries.
"""
import numpy as np

from user_block_funcs_george import compute_block  # noqa: F401 (the entries on the host)

BPACK_GPU_SYMMETRIC = 1
BPACK_GPU_COORDINATES = 2

# params: amplitude, the inverse metric (dim values), then yerr^2 per point
# (by 0-based global id)
ENTRY_SOURCE = r"""
__device__ double bpack_entry(const double* x, long long i, const double* y, long long j,
                              const double* params, int dim) {
    double q = 0.0;
    for (int d = 0; d < dim; ++d) {
        const double t = x[d] - y[d];
        q += params[1 + d] * t * t;
    }
    double v = params[0] * exp(-0.5 * q);
    if (i == j) v += params[1 + dim + i];
    return v;
}
"""


def gpu_entry(meta):
    """The payload's "gpu_entry": the source text and its parameters."""
    params = np.concatenate((
        [float(meta["amplitude"])],
        1.0 / np.asarray(meta["metric"], dtype=np.float64),
        np.asarray(meta["yerr"], dtype=np.float64) ** 2,
    ))
    return {
        "source": ENTRY_SOURCE,
        "params": params,
        "flags": BPACK_GPU_SYMMETRIC | BPACK_GPU_COORDINATES,
    }


_device = {}  # the points and parameters on the device, per process


def compute_entries_gpu(rows, cols, values, meta):
    """values[e] = K(rows[e], cols[e]): CuPy arrays (int64 ids, float64
    values) on the device, run on the library's stream.  A batch can hold
    tens of millions of entries: they are done in chunks, to bound CuPy's
    temporaries (taken from CuPy's memory pool, outside the library's)."""
    import cupy as cp

    if not _device:
        _device["xyz"] = cp.asarray(meta["coordinates"], dtype=cp.float64)
        _device["inverse_metric"] = cp.asarray(1.0 / np.asarray(meta["metric"], dtype=np.float64))
        _device["noise"] = cp.asarray(np.asarray(meta["yerr"], dtype=np.float64) ** 2)
    xyz = _device["xyz"]
    amplitude = float(meta["amplitude"])
    chunk = 1 << 22
    for start in range(0, rows.size, chunk):
        r = rows[start:start + chunk]
        c = cols[start:start + chunk]
        diff = xyz[r] - xyz[c]
        q = (diff * diff) @ _device["inverse_metric"]
        v = amplitude * cp.exp(-0.5 * q)
        values[start:start + chunk] = v + cp.where(r == c, _device["noise"][r], 0.0)
