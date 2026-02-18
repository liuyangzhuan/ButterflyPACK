import numpy as np
import os
import ctypes
import time
import sys
import dPy_BPACK_wrapper


try:
    import mpi4py
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    bcast = comm.bcast
    if(rank==0):
        print('mpi4py version: ', mpi4py.__version__)
        print('MPI count:', size)
except ImportError:
    comm=None
    rank=0
    size=1
    def bcast(obj, root=0):
        return obj




# ####################################################################################################
# ####################################################################################################
# ####################### create the matrix
# from user_block_funcs_1_r import * # this is the file that defines the compute_block function
# seed=12345
# rng = np.random.default_rng(seed=seed)
# nrhs = 1
# Npo = 1000
# Ndim = 3
# coordinates = rng.random((Npo, Ndim)).astype(np.float64)
# coordinates = bcast(coordinates, root=0)
# meta = {"coordinates": coordinates}





############################################### the following tests a RBF kernel in george
from user_block_funcs_george import * # this is the file that defines the compute_block function
import george
seed=12345
rng = np.random.default_rng(seed=seed)
nrhs = 1
Npo = 1000
Ndim = 3
coordinates = rng.random((Npo, Ndim)).astype(np.float64)
coordinates = bcast(coordinates, root=0)
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





####################################################################################################
####################################################################################################
####################### handle options
argv=sys.argv
if(len(argv)==1): # options are not passed via command line, set them manually here. If they are not set here, default values are used
    argv.extend(['-option'])
    argv.extend(['--tol_comp', '1e-6'])

argv.append(compute_block)
argv.append(meta)

argc = len(argv)
if(rank==0):
    print('BPACK options: ',argv[1:])
argv_ctypes = (ctypes.c_void_p * argc)(
    *[ctypes.c_void_p(id(a)) for a in argv]
)



####################################################################################################
####################################################################################################
####################### call the APIs
start = time.time()
sp = dPy_BPACK_wrapper.d_py_bpack_load()
####################### initialization
pyobj = ctypes.c_void_p()
maxrank = ctypes.c_int(0)

sp.d_py_bpack_init_compute(
    Npo,
    Ndim,
    coordinates.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    ctypes.byref(pyobj),            # void **pyobj
    ctypes.byref(maxrank),
    argc,                           # int argc
    argv_ctypes                     # char *argv[]
)

end = time.time()
if(rank==0):
    print(f"Time spent in py_bpack_init_compute: {end - start} seconds. Maxrank: {maxrank.value}")



####################### factor
start = time.time()
sp.d_py_bpack_factor(
    ctypes.byref(pyobj),            # void **pyobj
)
end = time.time()
if(rank==0):
    print(f"Time spent in d_py_bpack_factor: {end - start} seconds")


####################### solve
start = time.time()
xb = rng.random((Npo*nrhs)).astype(np.float64) # d_py_bpack_solve will broadcast xb on rank 0 to all ranks

# nrhs=2
# positive_values = np.arange(1, Npo + 1)  # Array [1, 2, ..., Npo]
# negative_values = -np.arange(1, Npo + 1)  # Array [-1, -2, ..., -Npo]
# xb = np.vstack((positive_values, negative_values)).astype(np.float64, copy=False)
# print(xb)


sp.d_py_bpack_solve(
    ctypes.byref(pyobj),            # void **pyobj
    nrhs,                           # int nrhs
    xb.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),  # double *nzval
)
end = time.time()
if(rank==0):
    print(f"Time spent in d_py_bpack_solve: {end - start} seconds")


####################### mult
start = time.time()
trans="N"
xb = rng.random((Npo*nrhs)).astype(np.float64) # d_py_bpack_mult will broadcast xb on rank 0 to all ranks
sp.d_py_bpack_mult(
    ctypes.byref(pyobj),            # void **pyobj
    nrhs,                           # int nrhs
    trans.encode("ascii"),
    xb.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),  # double *nzval
)
end = time.time()
if(rank==0):
    print(f"Time spent in d_py_bpack_mult: {end - start} seconds")


####################### log-determinant
start = time.time()
sign = ctypes.c_double(1)
logdet = ctypes.c_double(0.0)
sp.d_py_bpack_logdet(
    ctypes.byref(pyobj),            # void **pyobj
    ctypes.byref(sign),                           # int nrhs
    ctypes.byref(logdet),  # double *nzval
)

if(rank==0):
    print("bpack logdet:",int(sign.value),logdet.value)
    rows=np.arange(Npo)
    cols=np.arange(Npo)
    fullmat = compute_block(rows,cols,meta)
    sign, logdet = np.linalg.slogdet(fullmat)
    print("numpy logdet:",int(sign),logdet)
end = time.time()
if(rank==0):
    print(f"Time spent in d_py_bpack_logdet: {end - start} seconds")


####################### free stuff
start = time.time()
sp.d_py_bpack_free(ctypes.byref(pyobj))
end = time.time()
if(rank==0):
    print(f"Time spent in d_py_bpack_free: {end - start} seconds")


####################### terminate bpack
start = time.time()
sp.d_py_bpack_terminate()
end = time.time()
if(rank==0):
    print(f"Time spent in d_py_bpack_terminate: {end - start} seconds")
