#!/bin/bash
#SBATCH --job-name=mpi_omp_job
#SBATCH -A m2957
#SBATCH --constraint=cpu
#SBATCH --output=./vie_h2%j.log
#SBATCH --error=./vie_h2_%j.err
#SBATCH --qos=regular
#SBATCH --nodes=512
#SBATCH --ntasks-per-node=1        # MPI ranks per node
#SBATCH --time=02:59:00
#SBATCH --hint=nomultithread       # use physical cores only, no hyperthreading

module load PrgEnv-gnu cray-fftw
module unload cray-libsci
export LD_PRELOAD=/global/cfs/cdirs/m2957/lib/lib/PrgEnv-gnu/OpenBLAS/build/install/lib/libopenblas.so.0


# ── thread settings ──────────────────────────────────────────────
##########################################
# setting FMM_MAX_CPUS_PER_NODE will turn on dynamic threading in H2
##########################################
export FMM_MAX_CPUS_PER_NODE=128
##########################################
##########################################
NTH=1
export OMP_NUM_THREADS=$NTH                      # start small; the code raises this per level
export FMM_NUM_THREADS=$OMP_NUM_THREADS       # keep helper fallback/debug output aligned
export FMM_PIN_THREADS=1                      # disable helper-side pinning; rely on OpenMP binding
# export OMP_PROC_BIND=close                     # threads stay near their rank's cores
# export OMP_PLACES=cores                        # each "place" is a physical core
export OMP_DYNAMIC=FALSE                       # don't let runtime change thread count
export OMP_MAX_ACTIVE_LEVELS=1
THREADS_PER_RANK=`expr $NTH \* 2`	

# ── BLAS settings (prevent thread explosion) ─────────────────────
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1                       # if using MKL instead of OpenBLAS

nmpi=8

NTH=1							 
export OMP_NUM_THREADS=$NTH

h=0.05
omega=1.4
tols2s=1e-6
format=7

srun --export=ALL -n ${nmpi} -N 1 --cpu-bind=none /usr/bin/time -f "MaxRSS=%M KB" ./EXAMPLE/cvie3d_h2  --ivelo 9 --scaleGreen 0 --omega $omega --h $h --x0max 1.0 --y0max 1.0 --z0max 1.0 --L 1.0 --H  1.0 --W 1.0  --vs 1 --shape 4 --tol_comp 1e-8 --tol_comp_s2s $tols2s --format_s2s $format --sym 1 --reduction_threshold 64 --Nmin_leaf 0 --precon 1 --verbosity 1 --lrlevel 0

