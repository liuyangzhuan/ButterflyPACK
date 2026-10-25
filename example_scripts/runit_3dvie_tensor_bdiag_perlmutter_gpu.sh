#!/bin/bash
#SBATCH -A m4293_g
#SBATCH -C gpu
#SBATCH -q premium
#SBATCH -N 2
#SBATCH --gpus-per-node=4
#SBATCH -t 0:30:00
#SBATCH -J vie_bdiag_gpu
#SBATCH -o vie_bdiag_gpu_%j.out
# Tensor VIE (cvie3d_t_so: format 6, trans_invariant 2, FFT circulant) on NS^3 points at 6.25 points per wavelength
# of the background slowness (h = 0.00125, omega 402.12), solved by GMRES with a block-diagonal preconditioner on
# Perlmutter GPU nodes: 4 ranks a node, 16 cores and one GPU each; 2^BDIAG_LEVEL blocks (NS^3 / 2^BDIAG_LEVEL points
# each: NS=640, BDIAG_LEVEL=12 gives 4096 blocks of 40^3).  Runs the configurations in CASES, in order:
#   h2 = H2 blocks on the GPU (bdiag_format 7, structured Color), hodlr = HODLR blocks (bdiag_lrlevel 0; on the GPU,
#   HODLR_GPU=0: on the CPU), none = no preconditioner.
#
#   sbatch example_scripts/runit_3dvie_tensor_bdiag_perlmutter_gpu.sh             # 160^3, 64 blocks of 40^3, 2 nodes
# (the incident field vanishes inside cubes smaller than about 0.1, i.e. NS <= 80: GMRES returns at once)
#   NS=640 BDIAG_LEVEL=12 N_ITER=200 sbatch -N 64 -t 2:00:00 example_scripts/runit_3dvie_tensor_bdiag_perlmutter_gpu.sh
#
# Memory (measured, 640^3): 40^3 blocks at 1e-2 keep 2.1 GB of factors each, 16 a 40 GB GPU on 64 nodes; at 1e-3
# 2.9 GB, some blocks a GPU are solved on the host; 80^3 blocks at 1e-3 keep 42.5 GB each and need one 80 GB GPU
# each (-C gpu&hbm80g, NS=640 BDIAG_LEVEL=9 on 128 nodes).
#
# Environment: NS, BDIAG_LEVEL, CASES, BDIAG_TOL (bdiag_tol_rand), N_ITER, TOL_ITERSOL, ITER_SOLVER (2 GMRES, 1 TFQMR),
# GMRES_RESTART, IVELO (9: linear velocity profile), IVELO9_GRAD_SCALE (factor on its gradient; passed only when set),
# H2_COMPRESS (H2_compress_factors), BDIAG_OVERLAP (restricted additive Schwarz), HODLR_GPU, BDIAG_KNN, EXTRA_ARGS,
# CHECK (BPACK_CHECK), HEAP_FRACTION (BPACK_GPU_HEAP_FRACTION), BUDGET_MB (BPACK_MD_SLOT_BUDGET_MB), TAG,
# OUTDIR (logs and field slices; default the submit directory), FROZEN (a directory with copies of cvie3d_t_so and
# libbutterflypack.so.5.0.0 to run instead of the build), REPO.
# Build: build_gpu (enable_h2_gpu ON, sequential OpenBLAS).
set -u
module load PrgEnv-gnu cray-fftw cudatoolkit craype-accel-nvidia80 >/dev/null 2>&1
module unload cray-libsci >/dev/null 2>&1
ulimit -s unlimited

if [[ -n "${REPO:-}" ]]; then
  repo=$(cd "${REPO}" && pwd)
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -x "${SLURM_SUBMIT_DIR}/build_gpu/EXAMPLE/cvie3d_t_so" ]]; then
  # Slurm copies the batch script into /var/spool/slurmd, so BASH_SOURCE no
  # longer identifies the checkout when the script is submitted with sbatch.
  repo=$(cd "${SLURM_SUBMIT_DIR}" && pwd)
else
  repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fi
if [[ -n "${FROZEN:-}" ]]; then
  EXE=$FROZEN/cvie3d_t_so; LIB=$FROZEN/libbutterflypack.so.5.0.0
else
  EXE=$repo/build_gpu/EXAMPLE/cvie3d_t_so; LIB=$repo/build_gpu/libbutterflypack.so.5.0.0
fi
OUTDIR=${OUTDIR:-${SLURM_SUBMIT_DIR:-$PWD}}
cd $OUTDIR

TCMALLOC_ROOT=/global/cfs/cdirs/m2957/lib/lib/PrgEnv-gnu/gperftools-2.18.1
OPENBLAS_LIBRARY=/global/cfs/cdirs/m2957/lib/lib/PrgEnv-gnu/OpenBLAS_sequential/build/install/lib/libopenblas.so.0
export LD_LIBRARY_PATH="$TCMALLOC_ROOT/lib:/global/cfs/cdirs/m2957/liuyangz/my_software/scalapack-2.2.0/build/install/lib/:${LD_LIBRARY_PATH:-}"
export TCMALLOC_RELEASE_RATE=1
# The preload goes to the application only (env after /usr/bin/time, inside bash -c): libbutterflypack pulls in
# libgomp, which under OMP_PROC_BIND pins any wrapper that loads it (bash, time) to its first place - one hardware
# thread with OMP_PLACES=threads - and the application then inherits that mask: all its OpenMP threads, and so the
# tensor construction and matvec on the CPU, on one hardware thread.  See doc/running_on_perlmutter.md.
PRELOAD="$TCMALLOC_ROOT/lib/libtcmalloc.so:$OPENBLAS_LIBRARY:$LIB"
export OMP_NUM_THREADS=16 OMP_PLACES=threads OMP_PROC_BIND=spread OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1 GOTO_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export BPACK_MD_SLOT_BUDGET_MB=${BUDGET_MB:-0}
[ -n "${CHECK:-}" ] && export BPACK_CHECK=${CHECK}
[ -n "${HEAP_FRACTION:-}" ] && export BPACK_GPU_HEAP_FRACTION=${HEAP_FRACTION}
NN=${SLURM_JOB_NUM_NODES}
NR=$((NN * 4))
NS=${NS:-160}
BDIAG_LEVEL=${BDIAG_LEVEL:-6}
CASES=${CASES:-h2}
N_ITER=${N_ITER:-30}
TOL_ITERSOL=${TOL_ITERSOL:-1e-2}
BDIAG_TOL=${BDIAG_TOL:-1e-2}
TAG=${TAG:-gpu_n${NS}}
H=0.00125
X0=$(python3 -c "print(${NS}*${H})")
OMEGA=402.1238597   # 6.25 points per wavelength of the background slowness 2 at h = 0.00125
COMMON="--ivelo ${IVELO:-9} ${IVELO9_GRAD_SCALE:+--ivelo9_grad_scale ${IVELO9_GRAD_SCALE}} --pad_s2s 1 --trans_invariant 2 --scaleGreen 1 \
 --omega ${OMEGA} --htensor_mvp_level_batch 1 --h ${H} \
 --x0max ${X0} --y0max ${X0} --z0max ${X0} --L ${X0} --H ${X0} --W ${X0} --vs 1 --shape 4 --knn 0 --lrlevel 5 --nbundle 8 \
 --n_iter ${N_ITER} --format 6 --iter_solver ${ITER_SOLVER:-2} --gmres_restart ${GMRES_RESTART:-200} --baca_batch 8 --knn_near_para 20.0 --elem_extract 0 --near_para 2.1d0 \
 --use_zfp 1 --use_qtt 0 --xyzsort 1 --use_fft_circulant 1 --nmin_leaf 32 --sample_para_outer 0.8 --sample_para 0.8 \
 --rmax 4 --rankrate 2.0 --itermax 1 --tol_comp 1e-2 --tol_rand 1e-5 --tol_Rdetect 0.3e-5 --fastsample_tensor 2 \
 --verbosity 1 --bdiag_level ${BDIAG_LEVEL} --bdiag_lrlevel 0 --bdiag_tol_rand ${BDIAG_TOL} --tol_itersol ${TOL_ITERSOL} --precon 2 ${EXTRA_ARGS:-}"
H2OPTS="--H2_use_gpu 2 --H2_XRR_factor 1 --h2_use_sketch 2 --h2_lazy_schur 2 --CA_level 10000 --H2_compress_factors ${H2_COMPRESS:-0}"
for c in $CASES; do
  case $c in
    none)  EXTRA="--bdiag_precon 0" ;;
    hodlr) EXTRA="--bdiag_precon 1 --bdiag_format 1 --HODLR_use_gpu ${HODLR_GPU:-2} --bdiag_knn ${BDIAG_KNN:-0}" ;;
    h2)    EXTRA="--bdiag_precon 1 --bdiag_format 7 --bdiag_overlap ${BDIAG_OVERLAP:-0} ${H2OPTS}" ;;
    *)     echo "unknown case $c (h2, hodlr, none)"; continue ;;
  esac
  log=vie_${TAG}_${c}.log
  echo "=== ${c} (GPU): ${NR} ranks on ${NN} nodes, ${NS}^3, bdiag_level ${BDIAG_LEVEL}, $EXE -> ${OUTDIR}/${log}  $(date +%H:%M:%S)"
  srun -N ${NN} -n ${NR} --ntasks-per-node=4 -c 32 --cpu-bind=cores --gpus-per-node=4 \
    bash -c "export CUDA_VISIBLE_DEVICES=\$SLURM_LOCALID; exec /usr/bin/time -f 'maxrss_kb=%M wall=%e' env LD_PRELOAD=$PRELOAD $EXE ${COMMON} ${EXTRA}" > ${log} 2>&1
  echo "    exit $? $(date +%H:%M:%S)"
done
