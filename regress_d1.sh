#!/bin/bash -l
#SBATCH --account=m2957
#SBATCH -q debug
#SBATCH -C cpu
#SBATCH -N 4
#SBATCH -t 0:28:00
#SBATCH -J mdd1
ulimit -s unlimited
module load PrgEnv-gnu
module unload cray-libsci
module load cray-fftw
module load python
export LD_LIBRARY_PATH=/global/cfs/cdirs/m2957/liuyangz/my_software/scalapack-2.2.0/build/install/lib/:$LD_LIBRARY_PATH
export BPACK_BF_MVP_BATCH=1
export OMP_NUM_THREADS=8
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
ROOT=/global/cfs/cdirs/m2957/liuyangz/my_research/ButterflyPACK_gcc_libsci_perlmutter_MixedPrecision_SIE_2026
EXE=${EXE:-$ROOT/build_diag/EXAMPLE/cvie3d_t_so}
TAG=${TAG:-d1}
export BPACK_MD_SLOT_BUDGET_MB=${BUDGET_MB:-0}
echo "BPACK_MD_SLOT_BUDGET_MB=$BPACK_MD_SLOT_BUDGET_MB"

# 320^3 (h=0.0025), omega scaled from the 1280^3 production case to keep PPW:
# 804.24/4 = 201.06. 64 ranks on 4 nodes so the trans_invariant=0 branch, which
# allocates boundary_map over ALL row groups, has memory headroom.
COMMON="--ivelo 9 --scaleGreen 1 --omega 201.06192983 --htensor_mvp_level_batch 1 \
--h 0.0025 --x0max 0.8 --y0max 0.8 --z0max 0.8 --bdiag_precon 0 --L 0.8 --H 0.8 --W 0.8 \
--vs 1 --shape 4 --knn 0 --lrlevel 5 --nbundle 8 --n_iter 2 --format 6 --iter_solver 1 \
--baca_batch 8 --knn_near_para 20.0 --elem_extract 0 --near_para 2.1d0 --use_zfp 1 \
--use_qtt 0 --xyzsort 1 --use_fft_circulant 1 --nmin_leaf 32 --sample_para_outer 0.8 \
--sample_para 0.8 --rmax 4 --rankrate 2.0 --itermax 1 --tol_comp 1e-2 --tol_rand 1e-5 \
--tol_Rdetect 0.3e-5 --fastsample_tensor 2 --verbosity 2 --bdiag_lrlevel 0 \
--bdiag_tol_rand 1e-2 --tol_itersol 1e-2 --precon 2"

for ti in 2; do
  case $ti in
    0) pad=0 ;;
    2) pad=1 ;;
  esac
  echo "###### trans_invariant=${ti} pad_s2s=${pad} 320^3 64ranks ######"
  date
  srun -N 4 -n 64 -c 16 --cpu_bind=cores $EXE $COMMON \
     --trans_invariant ${ti} --pad_s2s ${pad} 2>&1 | tee $ROOT/reg_${TAG}_ti${ti}.out
  date
  echo "###### END ti=${ti} rc=$? ######"
done
