#!/bin/bash
# Build the bf_gpu harness against the installed GPU build of ButterflyPACK.
# The .cu holds the exchange kernels and needs nvcc; the rest is plain C++.
set -e
cd "$(dirname "$0")/.."
module load PrgEnv-gnu cudatoolkit >/dev/null 2>&1 || true
ROOT=$PWD
MAGMA=/global/cfs/cdirs/m2957/lib/magma_v2.10.0
MPIINC=$(CC --cray-print-opts=cflags 2>/dev/null | tr " " "\\n" | grep "^-I" | tr "\\n" " ")
INC="-I $ROOT/GPU_BACKEND -I $ROOT/SRC -I $ROOT/build_gpu/include -I $MAGMA/include $MPIINC"
OBJ=$ROOT/bf_gpu_test/obj
mkdir -p "$OBJ"
nvcc -std=c++17 -O2 -c -DH2_HAVE_GPU -arch=sm_80 $INC \
     -o "$OBJ/bf_exchange.o" GPU_BACKEND/bf_gpu/bf_exchange.cu
CC -std=c++17 -O2 -DH2_HAVE_GPU $INC \
   -o bf_gpu_test/bf_layout_check bf_gpu_test/bf_layout_check.cpp "$OBJ/bf_exchange.o" \
   -L"$ROOT/build_gpu/lib64" -lbutterflypack -L"$MAGMA/lib" -lmagma -lcudart -lcublas \
   -Wl,-rpath,"$ROOT/build_gpu/lib64" -Wl,-rpath,"$MAGMA/lib"
echo "HARNESS BUILD OK"
