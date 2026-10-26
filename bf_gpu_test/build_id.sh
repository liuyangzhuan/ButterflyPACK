#!/bin/bash
# Build the batched-ID unit check.  bf_id.cu needs nvcc; it calls
# fmm::gpu::launch_qrcp, which lives in the installed library.
set -e
cd /global/cfs/cdirs/m2957/liuyangz/my_research/ButterflyPACK_H2_HODLR_towards_STRUMPACK
module load PrgEnv-gnu cudatoolkit >/dev/null 2>&1 || true
ROOT=$PWD
MAGMA=/global/cfs/cdirs/m2957/lib/magma_v2.10.0
MPIINC=$(CC --cray-print-opts=cflags 2>/dev/null | tr ' ' '\n' | grep '^-I' | tr '\n' ' ')
INC="-I $ROOT/GPU_BACKEND -I $ROOT/SRC -I $ROOT/build_gpu/include -I $MAGMA/include $MPIINC"
OBJ=$ROOT/bf_gpu_test/obj; mkdir -p "$OBJ"
nvcc -std=c++17 -O2 -c -DH2_HAVE_GPU -arch=sm_80 $INC -o "$OBJ/bf_id.o" GPU_BACKEND/bf_gpu/bf_id.cu
CC -std=c++17 -O2 -DH2_HAVE_GPU $INC -o bf_gpu_test/bf_id_check bf_gpu_test/bf_id_check.cpp "$OBJ/bf_id.o" \
   -L"$ROOT/build_gpu/lib64" -lbutterflypack -L"$MAGMA/lib" -lmagma -lcudart -lcublas \
   -Wl,-rpath,"$ROOT/build_gpu/lib64" -Wl,-rpath,"$MAGMA/lib"
echo "ID CHECK BUILD OK"
