import numpy as np
import os
import ctypes
import time
import sys
import pickle
from dPy_BPACK_wrapper import *



####################################################################################################
####################################################################################################
####################### create the matrix
seed=12345
rng = np.random.default_rng(seed=seed)
nrhs = 1
verbosity=True
Npo = 1000
Ndim = 3
coordinates = rng.random((Npo, Ndim)).astype(np.float64)




############################################## the following compresses a 1/r kernel with channel_id 0
channel_id_0 = 0
meta = {"coordinates": coordinates}
payload = {
    "block_func_filepath": os.path.abspath(__file__),
    "block_func_module": "user_block_funcs_1_r",
    "block_func_name": "compute_block",
    "meta": meta
}
bpack_factor(payload, verbosity, fid=channel_id_0)




############################################### the following compresses a RBF kernel in george with channel_id 1
import george
channel_id_1 = 1
input_dim=Ndim
intialguess=[5e-6, 1] + [1]*input_dim
#### Note that intialguess contains theta, but george needs theta^2
K = george.kernels.ExpSquaredKernel(metric=np.array(intialguess[2:]), ndim=input_dim)
amplitude = intialguess[1]
K *= amplitude
err=np.sqrt(intialguess[0])
meta = {
    "coordinates": coordinates,
    "kernel": K,
    "yerr": np.repeat(err, Npo).astype(np.float64)
}
payload = {
    "block_func_filepath": os.path.abspath(__file__),
    "block_func_module": "user_block_funcs_george",
    "block_func_name": "compute_block",
    "meta": meta
}
bpack_factor(payload, verbosity, fid=channel_id_1)




####################################################################################################
####################################################################################################
####################### call the APIs when both instances are active

sign,logd = bpack_logdet(verbosity, fid=channel_id_0)
xb = rng.random((Npo,nrhs)).astype(np.float64,order="F")
y=bpack_solve(xb, verbosity, fid=channel_id_0)
xb = rng.random((Npo,nrhs)).astype(np.float64,order="F")
y=bpack_mult(xb, "N", verbosity, fid=channel_id_0)


sign,logd = bpack_logdet(verbosity, fid=channel_id_1)
xb = rng.random((Npo,nrhs)).astype(np.float64,order="F")
y=bpack_solve(xb, verbosity, fid=channel_id_1)
xb = rng.random((Npo,nrhs)).astype(np.float64,order="F")
y=bpack_mult(xb, "N", verbosity, fid=channel_id_1)

bpack_free(verbosity, fid=channel_id_0)
bpack_free(verbosity, fid=channel_id_1)


