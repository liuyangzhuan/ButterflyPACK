import os
import ctypes
import sys
from sys import platform
import time
import pickle
import numpy as np

def py_bpack_setup(sp):
    # Define the function signatures as shown in your original code

    sp.py_bpack_init_compute.restype = None
    sp.py_bpack_init_compute.argtypes = [
        ctypes.c_int, ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p)
    ]
    sp.py_bpack_factor.restype = None
    sp.py_bpack_factor.argtypes = [ctypes.POINTER(ctypes.c_void_p)]
    sp.py_bpack_solve.restype = None
    sp.py_bpack_solve.argtypes = [ctypes.POINTER(ctypes.c_void_p), ctypes.c_int, ctypes.POINTER(ctypes_dt)]
    sp.py_bpack_mult.restype = None
    sp.py_bpack_mult.argtypes = [ctypes.POINTER(ctypes.c_void_p), ctypes.c_int, ctypes.c_char_p, ctypes.POINTER(ctypes_dt)]
    sp.py_bpack_logdet.restype = None
    sp.py_bpack_logdet.argtypes = [ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes_dt), ctypes.POINTER(ctypes_rdt)]
    sp.py_bpack_free.restype = None
    sp.py_bpack_free.argtypes = [ctypes.POINTER(ctypes.c_void_p)]
    sp.py_bpack_terminate.restype = None
    sp.py_bpack_terminate.argtypes = None

def py_bpack_load():
    # Check platform and set library extension
    if platform == "linux" or platform == "linux2":
        pos = '.so'
    elif platform == "darwin":
        pos = '.dylib'
    elif platform == "win32":
        raise Exception("Windows is not yet supported")

    DLLFOUND = False
    INSTALLDIR = os.getenv('BPACK_PYTHON_LIB_PATH')

    DLL = os.path.abspath(__file__ + "/../../") + '/butterflypack_python' + pos
    if os.path.exists(DLL):
        DLLFOUND = True
    elif INSTALLDIR is not None:
        DLL = os.path.join(INSTALLDIR, 'butterflypack_python' + pos)
        if os.path.exists(DLL):
            DLLFOUND = True
    else:
        DLL = os.path.join('./butterflypack_python' + pos)
        if os.path.exists(DLL):
            DLLFOUND = True
    if DLLFOUND:
        sp = ctypes.cdll.LoadLibrary(DLL)
        py_bpack_setup(sp)
        return sp
    else:
        raise Exception("Cannot find the butterflypack_python library. Try to set the BPACK_PYTHON_LIB_PATH environment variable correctly.")




###################################################################################################
###########  define the APIs

def wait_for_flag(expected_flag, control_file, poll_interval=0.001):
    """Poll the control file until its content equals the expected flag."""
    while True:
        if os.path.exists(control_file):
            with open(control_file, "r") as f:
                flag = f.read().strip()
            if flag == expected_flag:
                return True
        time.sleep(poll_interval)


####################### initialization and factorization
def bpack_factor(payload, verbosity=False, nofactor=False, fid=0):

    MAX_ID_FILE = int(os.getenv("MAX_ID_FILE", "10"))
    if fid > MAX_ID_FILE:
        raise ValueError(f"fid={fid} exceeds MAX_ID_FILE={MAX_ID_FILE}. Please increase MAX_ID_FILE. ")

    start = time.time()
    CONTROL_FILE=os.getenv("CONTROL_FILE", "control.txt")

    # The following if test makes sure bpack cleans up the factorization is there is one
    if os.path.exists(f"{CONTROL_FILE}.{fid}"):
        with open(f"{CONTROL_FILE}.{fid}", "r") as f:
            flag = f.read().strip()
        if flag != "clean":
            bpack_free(verbosity,fid=fid)

    DATA_FILE=os.getenv("DATA_FILE", "data.bin")
    with open(f"{DATA_FILE}.{fid}", "wb") as f:
        pickle.dump(payload, f)
    with open(f"{CONTROL_FILE}.{fid}", "w") as f:
        f.write("init")
    wait_for_flag("done", f"{CONTROL_FILE}.{fid}")
    end = time.time()
    if verbosity==True:
        print(f"ID {fid}: Time spent in py_bpack_init_compute: {end - start} seconds")

    if nofactor==False:
        start = time.time()
        with open(f"{CONTROL_FILE}.{fid}", "w") as f:
            f.write("factor")
        wait_for_flag("done", f"{CONTROL_FILE}.{fid}")
        end = time.time()
        if verbosity==True:
            print(f"ID {fid}: Time spent in py_bpack_factor: {end - start} seconds")


####################### solve
def bpack_solve(vec, verbosity=False,fid=0):
    vec = np.asarray(vec, dtype=np_dt)
    orig_shape = vec.shape
    if vec.ndim == 1:
        vec = vec.reshape(-1, 1)
    nrhs=vec.shape[-1]
    start = time.time()
    CONTROL_FILE=os.getenv("CONTROL_FILE", "control.txt")
    DATA_FILE=os.getenv("DATA_FILE", "data.bin")
    RESULT_FILE=os.getenv("RESULT_FILE", "result.bin")
    with open(f"{DATA_FILE}.{fid}", "wb") as f:
        pickle.dump((vec,nrhs), f)
    with open(f"{CONTROL_FILE}.{fid}", "w") as f:
        f.write("solve")
    wait_for_flag("done", f"{CONTROL_FILE}.{fid}")
    with open(f"{RESULT_FILE}.{fid}", "rb") as f:
        vec_out = pickle.load(f)
    end = time.time()
    if verbosity==True:
        print(f"ID {fid}: Time spent in py_bpack_solve: {end - start} seconds")
    vec_out = vec_out.reshape(orig_shape)
    return vec_out



####################### mult
def bpack_mult(vec, trans, verbosity=False,fid=0):
    vec = np.asarray(vec, dtype=np_dt)
    orig_shape = vec.shape
    if vec.ndim == 1:
        vec = vec.reshape(-1, 1)
    nrhs=vec.shape[-1]
    start = time.time()
    CONTROL_FILE=os.getenv("CONTROL_FILE", "control.txt")
    DATA_FILE=os.getenv("DATA_FILE", "data.bin")
    RESULT_FILE=os.getenv("RESULT_FILE", "result.bin")
    with open(f"{DATA_FILE}.{fid}", "wb") as f:
        pickle.dump((vec,nrhs,trans), f)
    with open(f"{CONTROL_FILE}.{fid}", "w") as f:
        f.write("mult")
    wait_for_flag("done", f"{CONTROL_FILE}.{fid}")
    with open(f"{RESULT_FILE}.{fid}", "rb") as f:
        vec_out = pickle.load(f)
    end = time.time()
    if verbosity==True:
        print(f"ID {fid}: Time spent in py_bpack_mult: {end - start} seconds")
    vec_out = vec_out.reshape(orig_shape)
    return vec_out



####################### log-determinant
def bpack_logdet(verbosity=False,fid=0):
    start = time.time()
    CONTROL_FILE=os.getenv("CONTROL_FILE", "control.txt")
    RESULT_FILE=os.getenv("RESULT_FILE", "result.bin")
    with open(f"{CONTROL_FILE}.{fid}", "w") as f:
        f.write("logdet")
    wait_for_flag("done", f"{CONTROL_FILE}.{fid}")
    with open(f"{RESULT_FILE}.{fid}", "rb") as f:
        sign,log_det = pickle.load(f)
    end = time.time()
    if verbosity==True:
        print(f"ID {fid}: Time spent in py_bpack_logdet: {end - start} seconds")
    return sign,log_det

####################### free stuff
def bpack_free(verbosity=False,fid=0):
    start = time.time()
    CONTROL_FILE=os.getenv("CONTROL_FILE", "control.txt")
    with open(f"{CONTROL_FILE}.{fid}", "w") as f:
        f.write("free")
    wait_for_flag("clean", f"{CONTROL_FILE}.{fid}")
    end = time.time()
    if verbosity==True:
        print(f"ID {fid}: Time spent in py_bpack_free: {end - start} seconds")


####################### terminate all workers if no more bpack calls are needed, only f"{CONTROL_FILE}.{0}" is written the signal
def bpack_terminate(verbosity=False):
    start = time.time()
    CONTROL_FILE=os.getenv("CONTROL_FILE", "control.txt")
    with open(f"{CONTROL_FILE}.{0}", "w") as f:
        f.write("terminate")
    end = time.time()
    if verbosity==True:
        print(f"ID {0}: Time spent in py_bpack_terminate: {end - start} seconds")



