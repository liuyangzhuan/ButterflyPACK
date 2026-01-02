import numpy as np
import os
import ctypes
import time
import sys
import mpi4py
from mpi4py import MPI
import Py_BPACK_wrapper
import pickle
import importlib


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if(rank==0):
    print('mpi4py version: ', mpi4py.__version__)
    print('MPI count:', size)

####################### handle options
argv=sys.argv
if(len(argv)==1): # options are not passed via command line, set them manually here. If they are not set here, default values are used
    argv.extend(['-option'])
    argv.extend(['--tol_comp', '1e-6'])

######################## define the files used to communicate between masters and workers
CONTROL_FILE=os.getenv("CONTROL_FILE", "control.txt")
DATA_FILE=os.getenv("DATA_FILE", "data.bin")
RESULT_FILE=os.getenv("RESULT_FILE", "result.bin")
MAX_ID_FILE = int(os.getenv("MAX_ID_FILE", "10"))
poll_interval = 0.001
pyobjs = [None] * (MAX_ID_FILE + 1)
VALID_FLAGS = {"init", "factor", "solve", "mult", "logdet", "free", "terminate"}

# Ensure the file exists; if not, wait a moment and try again.
while True:
    flag = None
    fid = None
    if rank == 0:
        while True:
            for i in range(MAX_ID_FILE + 1):
                fname = f"{CONTROL_FILE}.{i}"
                if os.path.exists(fname):
                    with open(fname, "r") as f:
                        content = f.read().strip()
                    if content in VALID_FLAGS:
                            flag = content
                            fid = i
                            break
            if flag is not None:
                    break
            time.sleep(poll_interval)
    flag = comm.bcast(flag, root=0)
    if(flag=="init"):
        #####  read payload by rank 0 and broadcast
        payload=None
        if rank == 0:
            with open(f"{DATA_FILE}.{fid}", "rb") as f:
                payload = pickle.load(f)
        payload = comm.bcast(payload, root=0)
        # Resolve user function
        mod = importlib.import_module(payload["block_func_module"])
        compute_block = getattr(mod, payload["block_func_name"])
        meta = payload["meta"]

        argv_tmp = argv.copy()
        argv_tmp.append(compute_block)
        argv_tmp.append(meta)

        argc = len(argv_tmp)
        if(rank==0):
            threshold=np.get_printoptions()["threshold"]
            np.set_printoptions(threshold=50)
            print('BPACK options: ',argv_tmp[1:])
            np.set_printoptions(threshold=threshold)
        argv_ctypes = (ctypes.c_void_p * argc)(
            *[ctypes.c_void_p(id(a)) for a in argv_tmp]
        )

        sp = Py_BPACK_wrapper.py_bpack_load()
        ####################### initialization
        pyobjs[fid] = ctypes.c_void_p()

        sp.py_bpack_init_compute(
            meta["coordinates"].shape[0],
            meta["coordinates"].shape[1],
            meta["coordinates"].ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.byref(pyobjs[fid]),            # void **pyobj
            argc,                           # int argc
            argv_ctypes                     # char *argv[]
        )

    elif(flag=="factor"):
        ####################### factor
        sp.py_bpack_factor(
            ctypes.byref(pyobjs[fid]),            # void **pyobj
        )
    elif(flag=="solve"):
        ####################### solve
        #####  read in the RHS by rank 0
        nrhs=-1
        xb = np.random.rand(1).astype(np_dt)
        if rank == 0:
            with open(f"{DATA_FILE}.{fid}", "rb") as f:
                xb,nrhs = pickle.load(f)
            xb = np.ascontiguousarray(xb, dtype=np_dt)

        sp.py_bpack_solve(
            ctypes.byref(pyobjs[fid]),            # void **pyobj
            nrhs,                           # int nrhs
            xb.ctypes.data_as(ctypes.POINTER(ctypes_dt)),
        )
        if rank == 0:
            with open(f"{RESULT_FILE}.{fid}", "wb") as f:
                pickle.dump(xb, f)
    elif(flag=="mult"):
        ####################### mult
        #####  read in the RHS by rank 0
        nrhs=-1
        xb = np.random.rand(1).astype(np_dt)
        if rank == 0:
            with open(f"{DATA_FILE}.{fid}", "rb") as f:
                xb,nrhs,trans = pickle.load(f)
            xb = np.ascontiguousarray(xb, dtype=np_dt)

        sp.py_bpack_mult(
            ctypes.byref(pyobjs[fid]),            # void **pyobj
            nrhs,                           # int nrhs
            trans.encode("ascii"),
            xb.ctypes.data_as(ctypes.POINTER(ctypes_dt)),
        )
        if rank == 0:
            with open(f"{RESULT_FILE}.{fid}", "wb") as f:
                pickle.dump(xb, f)

    elif(flag=="logdet"):
        ####################### log-determinant
        sign = ctypes_dt(1)
        logdet = ctypes_rdt(0.0)
        sp.py_bpack_logdet(
            ctypes.byref(pyobjs[fid]),            # void **pyobj
            ctypes.byref(sign),                           # int nrhs
            ctypes.byref(logdet),  # double *nzval
        )
        if rank == 0:
            with open(f"{RESULT_FILE}.{fid}", "wb") as f:
                pickle.dump((sign.value, logdet.value),f)

        if(rank==0):
            print("bpack logdet:",sign.value,logdet.value)
            # sign, logdet = np.linalg.slogdet(m.toarray())
            # print("numpy logdet:",int(sign),logdet)

    elif(flag=="free"):
        ####################### free stuff
        sp.py_bpack_free(ctypes.byref(pyobjs[fid]))
    elif(flag=="terminate"):
        sp.py_bpack_terminate()
        break

    if(rank==0):
        ####################### signal the master that the work (init, factor, solve, logdet, free) has been completed
        with open(f"{CONTROL_FILE}.{fid}", "w") as f:
            if(flag=="free"):
                f.write("clean")
            else:
                f.write("done")

if(rank==0):
    ####################### signal the master that the work (terminate) has been completed and the MPI communicator will be released
    with open(f"{CONTROL_FILE}.{0}", "w") as f:
        f.write("finish")