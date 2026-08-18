#!/bin/bash
#SBATCH --job-name=h2
#SBATCH -A m2957
#SBATCH --constraint=cpu
#SBATCH --output=./vie_h2%j.log
#SBATCH --error=./vie_h2_%j.err
#SBATCH --qos=premium
#SBATCH --nodes=512
#SBATCH --ntasks-per-node=1        # MPI ranks per node
#SBATCH --time=00:30:00
#SBATCH --hint=nomultithread       # use physical cores only, no hyperthreading

module load PrgEnv-gnu cray-fftw
module unload cray-libsci
export LD_PRELOAD=/global/cfs/cdirs/m2957/lib/lib/PrgEnv-gnu/OpenBLAS_sequential/build/install/lib/libopenblas.so.0


# ── thread settings ──────────────────────────────────────────────
##########################################
# setting FMM_MAX_CPUS_PER_NODE will turn on dynamic threading in H2
##########################################
# export FMM_MAX_CPUS_PER_NODE=128
##########################################
##########################################


NTH=128
export OMP_NUM_THREADS=$NTH                      # start small; the code raises this per level
# export FMM_NUM_THREADS=$OMP_NUM_THREADS       # keep helper fallback/debug output aligned
# export FMM_PIN_THREADS=1                      # disable helper-side pinning; rely on OpenMP binding
# export OMP_PROC_BIND=close                     # threads stay near their rank's cores
# export OMP_PLACES=cores                        # each "place" is a physical core
export OMP_DYNAMIC=FALSE                       # don't let runtime change thread count
export OMP_MAX_ACTIVE_LEVELS=1
# THREADS_PER_RANK=`expr $NTH \* 2`	

# ── BLAS settings (prevent thread explosion) ─────────────────────
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export ARMPL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

tols2s=1e-3
NMIN_LEAF=216
format=7
NMPI=512
NNODES=512
REDUCTION_THRESHOLD=8
ELEM_EXTRACT=2
CA_LEVEL=10000                # 10000: full color; 0,1,2: full CA


# h=0.025
# omega=1.4
# srun --export=ALL -n ${nmpi} -N 4 --cpu-bind=none /usr/bin/time -f "MaxRSS=%M KB" ./EXAMPLE/cvie3d_h2  --ivelo 9 --scaleGreen 0 --omega $omega --h $h --x0max 1.0 --y0max 1.0 --z0max 1.0 --L 1.0 --H  1.0 --W 1.0  --vs 1 --shape 4 --tol_comp 1e-8 --tol_comp_s2s $tols2s --format_s2s $format --sym 1 --reduction_threshold 8 --CA_level "$CA_LEVEL" --Nmin_leaf 0 --precon 1 --verbosity 1 --lrlevel 0

# nmpi=1
# h=0.0104166666666667
# omega=4.71238898038469
# srun --export=ALL -n ${nmpi} -N 1 --cpu-bind=none /usr/bin/time -f "MaxRSS=%M KB" ./EXAMPLE/cvie3d_h2  --ivelo 9 --scaleGreen 0 --omega $omega --h $h --x0max 1.0 --y0max 1.0 --z0max 1.0 --L 1.0 --H  1.0 --W 1.0  --vs 1 --shape 4 --tol_comp $tols2s --tol_comp_s2s $tols2s --format_s2s $format --sym 1 --reduction_threshold 8 --CA_level "$CA_LEVEL" --Nmin_leaf 216 --precon 1 --elem_extract "$ELEM_EXTRACT" --verbosity 1 --lrlevel 0 | tee vie3d_h2_96_3.log



# h=0.0104166666666667
# n=96
# CA_LEVEL=2
# omega=4.71238898038469


# h=0.005208333333333333
# n=192
# CA_LEVEL=5            
# omega=4.71238898038469


# h=0.00260416666
# n=384
# CA_LEVEL=6            
# omega=4.71238898038469


h=0.00130208333
n=768
CA_LEVEL=7            
omega=4.71238898038469

# srun --export=ALL -n ${NMPI} -N ${NNODES} --cpu-bind=none /usr/bin/time -f "MaxRSS=%M KB" ./EXAMPLE/cvie3d_h2  --ivelo 9 --scaleGreen 0 --omega $omega --h $h --x0max 1.0 --y0max 1.0 --z0max 1.0 --L 1.0 --H  1.0 --W 1.0  --vs 1 --shape 4 --tol_comp $tols2s --tol_comp_s2s $tols2s --format_s2s $format --sym 1 --reduction_threshold ${REDUCTION_THRESHOLD} --CA_level "$CA_LEVEL" --Nmin_leaf ${NMIN_LEAF} --precon 1 --elem_extract "$ELEM_EXTRACT" --verbosity 1 --lrlevel 0  |& tee vie3d_h2_${n}_3_calv_${CA_LEVEL}_rt${REDUCTION_THRESHOLD}_N${NNODES}_nmpi${NMPI}.log

CA_LEVEL=10000           
srun --export=ALL -n ${NMPI} -N ${NNODES} --cpu-bind=none /usr/bin/time -f "MaxRSS=%M KB" ./EXAMPLE/cvie3d_h2  --ivelo 9 --scaleGreen 0 --omega $omega --h $h --x0max 1.0 --y0max 1.0 --z0max 1.0 --L 1.0 --H  1.0 --W 1.0  --vs 1 --shape 4 --tol_comp $tols2s --tol_comp_s2s $tols2s --format_s2s $format --sym 1 --reduction_threshold ${REDUCTION_THRESHOLD} --CA_level "$CA_LEVEL" --Nmin_leaf ${NMIN_LEAF} --precon 1 --elem_extract "$ELEM_EXTRACT" --verbosity 1 --lrlevel 0  |& tee vie3d_h2_${n}_3_calv_${CA_LEVEL}_rt${REDUCTION_THRESHOLD}_N${NNODES}_nmpi${NMPI}.log

