#!/bin/bash
#SBATCH --job-name=matern2d_h2
#SBATCH -A m2957
#SBATCH --constraint=cpu
#SBATCH --output=./matern2d_h2_%j.log
#SBATCH --error=./matern2d_h2_%j.err
#SBATCH --qos=premium
#SBATCH --nodes=4
#SBATCH --time=05:00:00

module load PrgEnv-gnu cray-fftw
module unload cray-libsci

export LD_PRELOAD=/global/cfs/cdirs/m2957/lib/lib/PrgEnv-gnu/OpenBLAS_sequential/build/install/lib/libopenblas.so.0

# export FMM_MAX_CPUS_PER_NODE=128
export OMP_NUM_THREADS=8
# export FMM_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export OMP_MAX_ACTIVE_LEVELS=1

export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export ARMPL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

NMPI=256
NNODES=16
GRID_SIZE=16384
TOL=1e-8
NMIN_LEAF=64
REDUCTION_THRESHOLD=4
CA_LEVEL=0                 # 0: full CA where eligible; 10000: full color
ELEM_EXTRACT=2
VERBOSITY=1
LENGTH_SCALE=0.1
NUGGET=1e-3
CG_TOL=1e-13
CG_MAX_ITERATIONS=100

srun --export=ALL -N "${NNODES}" -n "${NMPI}" --cpu-bind=none \
  /usr/bin/time -f "MaxRSS=%M KB" \
  ./EXAMPLE/cmatern2d_h2 \
  --grid-size "${GRID_SIZE}" \
  --tol-comp "${TOL}" \
  --Nmin_leaf "${NMIN_LEAF}" \
  --reduction_threshold "${REDUCTION_THRESHOLD}" \
  --CA_level "${CA_LEVEL}" \
  --length-scale "${LENGTH_SCALE}" \
  --nugget "${NUGGET}" \
  --cg-tol "${CG_TOL}" \
  --cg-max-iterations "${CG_MAX_ITERATIONS}" \
  --precon 3 \
  --iter_solver 4 \
  --elem_extract "${ELEM_EXTRACT}" \
  --verbosity "${VERBOSITY}" \
  2>&1 | tee matern2d_h2_${GRID_SIZE}_2_calv_${CA_LEVEL}_rt${REDUCTION_THRESHOLD}_tol${TOL}.log

CA_LEVEL=10000

srun --export=ALL -N "${NNODES}" -n "${NMPI}" --cpu-bind=none \
  /usr/bin/time -f "MaxRSS=%M KB" \
  ./EXAMPLE/cmatern2d_h2 \
  --grid-size "${GRID_SIZE}" \
  --tol-comp "${TOL}" \
  --Nmin_leaf "${NMIN_LEAF}" \
  --reduction_threshold "${REDUCTION_THRESHOLD}" \
  --CA_level "${CA_LEVEL}" \
  --length-scale "${LENGTH_SCALE}" \
  --nugget "${NUGGET}" \
  --cg-tol "${CG_TOL}" \
  --cg-max-iterations "${CG_MAX_ITERATIONS}" \
  --precon 3 \
  --iter_solver 4 \
  --elem_extract "${ELEM_EXTRACT}" \
  --verbosity "${VERBOSITY}" \
  2>&1 | tee matern2d_h2_${GRID_SIZE}_2_calv_${CA_LEVEL}_rt${REDUCTION_THRESHOLD}_tol${TOL}.log
