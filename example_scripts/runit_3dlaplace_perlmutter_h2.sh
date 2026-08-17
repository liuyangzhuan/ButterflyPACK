#!/bin/bash
#SBATCH --job-name=laplace3d_h2
#SBATCH -A m2957
#SBATCH --constraint=cpu
#SBATCH --output=./laplace3d_h2_%j.log
#SBATCH --error=./laplace3d_h2_%j.err
#SBATCH --qos=premium
#SBATCH --nodes=16
#SBATCH --time=01:00:00


module load PrgEnv-gnu cray-fftw
module unload cray-libsci

export LD_PRELOAD=/global/cfs/cdirs/m2957/lib/lib/PrgEnv-gnu/OpenBLAS_sequential/build/install/lib/libopenblas.so.0

export FMM_MAX_CPUS_PER_NODE=128
export OMP_NUM_THREADS=1
# export FMM_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export OMP_MAX_ACTIVE_LEVELS=1

export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export ARMPL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

NMPI=8
NNODES=4
GRID_SIZE=192
TOL=1e-3
NMIN_LEAF=216
REDUCTION_THRESHOLD=8
CA_LEVEL=2                # 10000: full color; 0,1,2: full CA
ELEM_EXTRACT=0
VERBOSITY=1

srun --export=ALL -N "${NNODES}" -n "${NMPI}" --cpu-bind=none \
  /usr/bin/time -f "MaxRSS=%M KB" \
  ./EXAMPLE/claplace3d_h2 \
  --grid-size "${GRID_SIZE}" \
  --tol-comp "${TOL}" \
  --Nmin_leaf "${NMIN_LEAF}" \
  --reduction_threshold "${REDUCTION_THRESHOLD}" \
  --CA_level "${CA_LEVEL}" \
  --elem_extract "${ELEM_EXTRACT}" \
  --verbosity "${VERBOSITY}" |& tee laplace3d_h2_${GRID_SIZE}_3_calv_${CA_LEVEL}_rt${REDUCTION_THRESHOLD}.log


CA_LEVEL=10000               # 10000: full color; 0,1,2: full CA 

  srun --export=ALL -N "${NNODES}" -n "${NMPI}" --cpu-bind=none \
  /usr/bin/time -f "MaxRSS=%M KB" \
  ./EXAMPLE/claplace3d_h2 \
  --grid-size "${GRID_SIZE}" \
  --tol-comp "${TOL}" \
  --Nmin_leaf "${NMIN_LEAF}" \
  --reduction_threshold "${REDUCTION_THRESHOLD}" \
  --CA_level "${CA_LEVEL}" \
  --elem_extract "${ELEM_EXTRACT}" \
  --verbosity "${VERBOSITY}" |& tee laplace3d_h2_${GRID_SIZE}_3_calv_${CA_LEVEL}_rt${REDUCTION_THRESHOLD}.log

