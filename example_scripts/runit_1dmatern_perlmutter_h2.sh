#!/bin/bash
#SBATCH --job-name=matern1d_h2
#SBATCH -A m2957
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:30:00
#SBATCH --output=./matern1d_h2_%j.log
#SBATCH --error=./matern1d_h2_%j.err

module load PrgEnv-gnu cray-fftw
module unload cray-libsci

export LD_PRELOAD=/global/cfs/cdirs/m2957/lib/lib/PrgEnv-gnu/OpenBLAS_sequential/build/install/lib/libopenblas.so.0

export OMP_NUM_THREADS=32
export OMP_DYNAMIC=FALSE
export OMP_MAX_ACTIVE_LEVELS=1

export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export ARMPL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

NMPI=4
NNODES=1
GRID_SIZE=1048576
TOL=1e-8
NMIN_LEAF=64
REDUCTION_THRESHOLD=2
CA_LEVEL=10000
PRECON=1
NRHS=1
ELEM_EXTRACT=2
VERBOSITY=1
LENGTH_SCALE=0.1
NUGGET=1e-3

srun --export=ALL -N "${NNODES}" -n "${NMPI}" --cpu-bind=none \
  /usr/bin/time -f "MaxRSS=%M KB" \
  ./EXAMPLE/cmatern1d_h2 \
  --grid-size "${GRID_SIZE}" \
  --tol-comp "${TOL}" \
  --Nmin_leaf "${NMIN_LEAF}" \
  --reduction_threshold "${REDUCTION_THRESHOLD}" \
  --CA_level "${CA_LEVEL}" \
  --length-scale "${LENGTH_SCALE}" \
  --nugget "${NUGGET}" \
  --precon "${PRECON}" \
  --nrhs "${NRHS}" \
  --kernel matern52 \
  --number-type real \
  --dimension 1 \
  --num-proxy 0 \
  --format 7 \
  --sym 1 \
  --elem_extract "${ELEM_EXTRACT}" \
  --verbosity "${VERBOSITY}" \
  2>&1 | tee matern1d_h2_${GRID_SIZE}_1_calv_${CA_LEVEL}_rt${REDUCTION_THRESHOLD}_N${NNODES}_nmpi${NMPI}.log
