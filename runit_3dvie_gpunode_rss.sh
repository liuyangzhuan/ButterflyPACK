#!/bin/bash -l

#SBATCH --account=m4872_g
#SBATCH -q premium
#SBATCH -C gpu
#SBATCH -N 64
#SBATCH --ntasks-per-node=4
#SBATCH -c 32
#SBATCH --gpus-per-task=1
#SBATCH -t 0:12:00
#SBATCH -J vierss
#SBATCH --mail-user=liuyangzhuan@lbl.gov
ulimit -s unlimited

module load PrgEnv-gnu
module load cmake
module unload cray-libsci
module load cray-fftw
module load python

ROOT=/global/cfs/cdirs/m2957/liuyangz/my_research/ButterflyPACK_gcc_libsci_perlmutter_MixedPrecision_SIE_2026

export LD_LIBRARY_PATH=/global/cfs/cdirs/m2957/liuyangz/my_software/scalapack-2.2.0/build/install/lib/:$LD_LIBRARY_PATH
export BPACK_BF_MVP_BATCH=1   # batched BF MVP (OpenBLAS)
export OMP_NUM_THREADS=16     # 4 ranks x 16 threads = 64 cores = full GPU node
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

nmpi=256

# ---- problem/solver parameters: as in runit_3dvie_perlmutter_nobdiag.sh ----
bdiag_precon=0
vs=1
shape=4
scaleGreen=1
L=0.8; H=0.8; W=0.8
tol_rand=1e-5
tol_Rdetect=0.3e-5
tol_itersol=1e-2
nmin_leaf_t=32
xyzsort=1
BACA_Batch=8
near_para=2.1d0
precon=2
n_iter=3            # a few iterations is enough for per-iteration timing
verbosity=2
knn_near_para=20.0
knn=0
fastsample_tensor=2
iter_solver=1
nbundle=8
bdiag_lrlevel=0
bdiag_tol_rand=1e-2

omega=804.24
ivelo=9
tol_comp=1e-2
h=0.000625
x0max=0.8; y0max=0.8; z0max=0.8

elem_extract=0
sample_para_outer=0.8
sample_para=0.8

# HTENSOR with trans_invariant
format=6
pad_s2s=1
trans_invariant=2
rmax=4
LRlevel=5
rankrate=2.0
itermax=1
use_zfp=1
use_qtt=0
use_fftw=1

run_one () {
  local tag=$1
  local exedir=$2
  local lbatch=$3
  cd ${exedir}
  echo "###### RUN tag=${tag} exe=${exedir}/EXAMPLE/cvie3d_t_so OMP=${OMP_NUM_THREADS} level_batch=${lbatch} ######"
  date
  srun -N 64 -n ${nmpi} -c 32 --cpu_bind=cores \
    ./EXAMPLE/cvie3d_t_so --ivelo ${ivelo} --pad_s2s ${pad_s2s} --trans_invariant ${trans_invariant} \
    --scaleGreen ${scaleGreen} --omega ${omega} --htensor_mvp_level_batch ${lbatch} --h ${h} \
    --x0max $x0max --y0max $y0max --z0max $z0max --bdiag_level ${bdiag_level} --bdiag_precon ${bdiag_precon} \
    --L ${L} --H ${H} --W ${W} --nshape ${nshape} --smin_ivelo11 ${slowness_min} --smax_ivelo11 ${slowness_max} \
    --vs ${vs} --shape ${shape} --knn ${knn} --lrlevel ${LRlevel} --nbundle ${nbundle} --n_iter ${n_iter} \
    --format ${format} --iter_solver ${iter_solver} --baca_batch ${BACA_Batch} --knn_near_para ${knn_near_para} \
    --elem_extract ${elem_extract} --near_para ${near_para} --use_zfp ${use_zfp} --use_qtt ${use_qtt} \
    --xyzsort ${xyzsort} --use_fft_circulant ${use_fftw} --nmin_leaf ${nmin_leaf_t} \
    --sample_para_outer ${sample_para_outer} --sample_para ${sample_para} --rmax $rmax --rankrate ${rankrate} \
    --itermax ${itermax} --tol_comp ${tol_comp} --tol_rand ${tol_rand} --tol_Rdetect ${tol_Rdetect} \
    --fastsample_tensor ${fastsample_tensor} --verbosity ${verbosity} --bdiag_lrlevel ${bdiag_lrlevel} \
    --bdiag_tol_rand ${bdiag_tol_rand} --tol_itersol ${tol_itersol} --precon ${precon} \
    2>&1 | tee ${ROOT}/cmp_${tag}.out
  date
  echo "###### END tag=${tag} ######"
}

# all fixes; vary the number of levels redistributed per exchange.
# level_batch=1 gives 12 sync points per matvec, 6 gives 2.
# per-level load balance + memory high-water mark, level_batch=1 so each
# MVPPROF/MDBAL level entry corresponds to exactly one Bplus level
run_one rss $ROOT/build_diag 1
