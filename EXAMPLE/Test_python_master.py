import numpy as np
import os
import ctypes
import time
import sys
import pickle
import argparse
from dPy_BPACK_wrapper import *


def positive_int(value):
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError(f"expected a positive integer, got {value}")
    return parsed


parser = argparse.ArgumentParser(
    description="Run the ButterflyPACK Python file-interface example."
)
parser.add_argument(
    "--Npo", "--npo", "--np", dest="npo", type=positive_int, default=1000000,
    help="number of points (default: 1000000)",
)
parser.add_argument(
    "--Ndim", "--ndim", dest="ndim", type=positive_int, default=3,
    help="coordinate dimension (default: 3)",
)
parser.add_argument(
    "--gpu", choices=["none", "entry", "list"], default="none",
    help="GPU evaluator of the george kernel for H2_use_gpu / HODLR_use_gpu "
         "(user_block_funcs_george_gpu.py): CUDA source text compiled with NVRTC "
         "(entry) or a CuPy entry-list function (list); default none",
)
args = parser.parse_args()


####################################################################################################
####################################################################################################
####################### create the matrix
seed=12345
rng = np.random.default_rng(seed=seed)
nrhs = 1
verbosity=True
Npo = args.npo
Ndim = args.ndim
coordinates = rng.random((Npo, Ndim)).astype(np.float64)




############################################## the following tests a 1/r kernel
# meta = {"coordinates": coordinates}
# payload = {
#     "block_func_filepath": os.path.abspath(__file__),
#     "block_func_module": "user_block_funcs_1_r",
#     "block_func_name": "compute_block",
#     "meta": meta
# }


############################################### the following tests a RBF kernel in george
import george
input_dim=Ndim
intialguess=[5e-2, 1] + [1]*input_dim
#### Note that intialguess contains theta, but george needs theta^2
K = george.kernels.ExpSquaredKernel(metric=np.array(intialguess[2:]), ndim=input_dim)
amplitude = intialguess[1]
K *= amplitude
err=np.sqrt(intialguess[0])
meta = {
    "coordinates": coordinates,
    "kernel": K,
    "yerr": np.repeat(err, Npo).astype(np.float64),
    # the kernel's parameters, for its GPU evaluators
    "amplitude": amplitude,
    "metric": np.array(intialguess[2:], dtype=np.float64)
}
payload = {
    "block_func_filepath": os.path.abspath(__file__),
    "block_func_module": "user_block_funcs_george",
    "block_func_name": "compute_block",
    "meta": meta
}
if args.gpu == "entry":
    from user_block_funcs_george_gpu import gpu_entry
    payload["gpu_entry"] = gpu_entry(meta)
elif args.gpu == "list":
    payload["block_func_module"] = "user_block_funcs_george_gpu"
    payload["gpu_list"] = {"func_name": "compute_entries_gpu", "flags": 1}  # BPACK_GPU_SYMMETRIC





####################################################################################################
####################################################################################################
####################### call the APIs


bpack_factor(payload, verbosity)
sign,logd = bpack_logdet(verbosity)
xb = rng.random((Npo,nrhs)).astype(np.float64,order="F")
y=bpack_solve(xb, verbosity)
xb = rng.random((Npo,nrhs)).astype(np.float64,order="F")
y=bpack_mult(xb, "N", verbosity)
bpack_free(verbosity)
bpack_terminate(verbosity)
