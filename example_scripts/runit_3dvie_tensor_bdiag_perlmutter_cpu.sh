#!/bin/bash
#SBATCH -A m4293
#SBATCH -C cpu
#SBATCH -q debug
#SBATCH -N 1
#SBATCH -t 0:30:00
#SBATCH -J vie_bdiag_cpu
#SBATCH -o vie_bdiag_cpu_%j.out
# Tensor VIE (cvie3d_t_so: format 6, trans_invariant 2, FFT circulant) on NS^3 points at 6.25 points per wavelength
# of the background slowness (h = 0.00125, omega 402.12), solved by GMRES with a block-diagonal preconditioner of H2
# blocks (bdiag_format 7, structured Color, H2_use_gpu 0) on Perlmutter CPU nodes (128 cores, 512 GB):
# RPN ranks a node with 128/RPN cores each; 2^BDIAG_LEVEL blocks (NS=640, BDIAG_LEVEL=12: 4096 blocks of 40^3).
#
#   sbatch -N 1 example_scripts/runit_3dvie_tensor_bdiag_perlmutter_cpu.sh                  # 160^3, 64 blocks
#   NS=640 BDIAG_LEVEL=12 BDIAG_TOL=1e-2 N_ITER=3 GMRES_RESTART=30 \
#     sbatch -q regular -N 32 -t 0:30:00 example_scripts/runit_3dvie_tensor_bdiag_perlmutter_cpu.sh
#
# Measured (2026-10-24, the line above): 291 s, tensor construction 18 s, H2 block setup 177 s (16 blocks of 40^3
# a rank), matvec 5.3 s, preconditioner application 2.5 s, 47.8 GB a rank.
#
# Environment: NS, BDIAG_LEVEL, BDIAG_TOL (H2 block tolerance, bdiag_tol_rand), N_ITER, TOL_ITERSOL, ITER_SOLVER
# (2 GMRES, 1 TFQMR), GMRES_RESTART, IVELO (9: linear velocity profile, 11: random shapes, see EXTRA_ARGS),
# IVELO9_GRAD_SCALE (factor on the ivelo 9 velocity gradient), CASE (h2 | none: no preconditioner), RPN, BLAS
# (omp: the OpenMP-threaded OpenBLAS preloaded, 17% faster H2 blocks | seq: the build's sequential one),
# EXTRA_ARGS, EXTRA_PRELOAD, WRAP (a launcher before the application, e.g. a perf wrapper), TAG, OUTDIR (logs and
# field slices; default the submit directory), FROZEN (a directory with copies of cvie3d_t_so and
# libbutterflypack.so.5.0.0 to run instead of the build), REPO.
# Build: a CPU-only build in build_cpu (enable_parsec OFF, enable_h2_gpu OFF, TPL_BLAS_LIBRARIES and
# TPL_LAPACK_LIBRARIES the sequential OpenBLAS in /global/cfs/cdirs/m2957/lib/lib/PrgEnv-gnu/OpenBLAS_sequential);
# its TPLs need the cudatoolkit module at run time (CPU nodes have libcuda).
set -u
module load cpu >/dev/null 2>&1
module load PrgEnv-gnu cray-fftw cudatoolkit >/dev/null 2>&1
module unload cray-libsci craype-accel-nvidia80 >/dev/null 2>&1
ulimit -s unlimited

if [[ -n "${REPO:-}" ]]; then
  repo=$(cd "${REPO}" && pwd)
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -x "${SLURM_SUBMIT_DIR}/build_cpu/EXAMPLE/cvie3d_t_so" ]]; then
  # Slurm copies the batch script into /var/spool/slurmd, so BASH_SOURCE no
  # longer identifies the checkout when the script is submitted with sbatch.
  repo=$(cd "${SLURM_SUBMIT_DIR}" && pwd)
else
  repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fi
if [[ -n "${FROZEN:-}" ]]; then
  EXE=$FROZEN/cvie3d_t_so; LIB=$FROZEN/libbutterflypack.so.5.0.0
else
  EXE=$repo/build_cpu/EXAMPLE/cvie3d_t_so; LIB=$repo/build_cpu/libbutterflypack.so.5.0.0
fi
OUTDIR=${OUTDIR:-${SLURM_SUBMIT_DIR:-$PWD}}
cd $OUTDIR

TCMALLOC_ROOT=/global/cfs/cdirs/m2957/lib/lib/PrgEnv-gnu/gperftools-2.18.1
export LD_LIBRARY_PATH="$TCMALLOC_ROOT/lib:${LD_LIBRARY_PATH:-}"
export TCMALLOC_RELEASE_RATE=1
# The preload goes to the application only (env after /usr/bin/time): libbutterflypack pulls in libgomp, which
# under OMP_PROC_BIND pins any wrapper that loads it (bash, time) to its first place, and the application then
# inherits that one-core mask and builds a one-place list - all its OpenMP threads on one core.
OMP_OPENBLAS=/global/cfs/cdirs/m2957/lib/lib/PrgEnv-gnu/OpenBLAS/build/install/lib/libopenblas.so.0
BLAS=${BLAS:-omp}
BLAS_PRELOAD=""; [ "$BLAS" = omp ] && BLAS_PRELOAD=$OMP_OPENBLAS
PRELOAD="$TCMALLOC_ROOT/lib/libtcmalloc.so:${EXTRA_PRELOAD:+$EXTRA_PRELOAD:}${BLAS_PRELOAD:+$BLAS_PRELOAD:}$LIB"

RPN=${RPN:-8}
TPR=$((128 / RPN))
export OMP_NUM_THREADS=$TPR OMP_PLACES=cores OMP_PROC_BIND=spread OMP_MAX_ACTIVE_LEVELS=1
# the threaded OpenBLAS runs on TPR threads outside the code's parallel regions and on one inside them (the
# sequential one ignores this)
export OPENBLAS_NUM_THREADS=${BLAS_THREADS:-$TPR} GOTO_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export MPICH_GPU_SUPPORT_ENABLED=0
NN=${SLURM_JOB_NUM_NODES}
NR=$((NN * RPN))
NS=${NS:-160}
BDIAG_LEVEL=${BDIAG_LEVEL:-6}
N_ITER=${N_ITER:-30}
TOL_ITERSOL=${TOL_ITERSOL:-1e-2}
BDIAG_TOL=${BDIAG_TOL:-1e-3}
TAG=${TAG:-cpu_n${NS}}
H=0.00125
X0=$(python3 -c "print(${NS}*${H})")
OMEGA=402.1238597   # 6.25 points per wavelength of the background slowness 2 at h = 0.00125
COMMON="--ivelo ${IVELO:-9} --ivelo9_grad_scale ${IVELO9_GRAD_SCALE:-1} --pad_s2s 1 --trans_invariant 2 --scaleGreen 1 \
 --omega ${OMEGA} --htensor_mvp_level_batch 1 --h ${H} \
 --x0max ${X0} --y0max ${X0} --z0max ${X0} --L ${X0} --H ${X0} --W ${X0} --vs 1 --shape 4 --knn 0 --lrlevel 5 --nbundle 8 \
 --n_iter ${N_ITER} --format 6 --iter_solver ${ITER_SOLVER:-2} --gmres_restart ${GMRES_RESTART:-200} --baca_batch 8 --knn_near_para 20.0 --elem_extract 0 --near_para 2.1d0 \
 --use_zfp 1 --use_qtt 0 --xyzsort 1 --use_fft_circulant 1 --nmin_leaf 32 --sample_para_outer 0.8 --sample_para 0.8 \
 --rmax 4 --rankrate 2.0 --itermax 1 --tol_comp 1e-2 --tol_rand 1e-5 --tol_Rdetect 0.3e-5 --fastsample_tensor 2 \
 --verbosity 1 --bdiag_level ${BDIAG_LEVEL} --bdiag_lrlevel 0 --bdiag_tol_rand ${BDIAG_TOL} --tol_itersol ${TOL_ITERSOL} --precon 2 ${EXTRA_ARGS:-}"
H2OPTS="--H2_use_gpu 0 --H2_XRR_factor 1 --h2_use_sketch 2 --h2_lazy_schur 2 --CA_level 10000"
CASE=${CASE:-h2}   # h2: H2 diagonal blocks on the CPU; none: no preconditioner
if [ "$CASE" = none ]; then EXTRA="--bdiag_precon 0"; else EXTRA="--bdiag_precon 1 --bdiag_format 7 ${H2OPTS}"; fi
log=vie_${TAG}_${CASE}.log
echo "=== ${CASE} (CPU): ${NR} ranks x ${TPR} threads on ${NN} nodes, ${NS}^3, bdiag_level ${BDIAG_LEVEL}, BLAS ${BLAS}, $EXE -> ${OUTDIR}/${log}  $(date +%H:%M:%S)"
srun -N ${NN} -n ${NR} --ntasks-per-node=${RPN} -c $((2 * TPR)) --cpu-bind=cores \
  /usr/bin/time -f 'maxrss_kb=%M wall=%e' ${WRAP:-} env LD_PRELOAD="$PRELOAD" $EXE ${COMMON} ${EXTRA} > ${log} 2>&1
echo "    exit $? $(date +%H:%M:%S)"
