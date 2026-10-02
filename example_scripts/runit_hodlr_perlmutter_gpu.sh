#!/bin/bash
#SBATCH --job-name=hodlr_gpu
#SBATCH -A m2957
#SBATCH --constraint=gpu
#SBATCH --qos=regular
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --time=00:30:00
#SBATCH --output=./hodlr_gpu_%j.log

# HODLR (format 1, LRlevel 0) with the GPU backend (HODLR_use_gpu) on
# Perlmutter GPU nodes, with the GPU build of
# run_cmake_build_gnu_perlmutter_openblas_sequential_h2gpu.sh (in ../build_gpu).
# The construction, the factorization, the multiply and the solve run on the
# GPUs, for sym=1 (symmetric HODLR) and sym=0, over any number of nodes: one
# MPI rank per A100 by default (4 per node, 16 cores each); with
# RANKS_PER_NODE=8 two ranks share each GPU and split its memory.  See
# ../GPU_BACKEND/hodlr_gpu/README.md for the options and environment switches.
#
# Overrides: NODES, RANKS_PER_NODE, CASE (laplace, vie, efie, cfie), SYM (0 or
# 1; default 1 for laplace, 0 otherwise), USE_GPU (1: FP64 GEMMs, 2: FP64 tensor-core GEMMs), PIECES (the option
# HODLR_gpu_pieces, default 4), CHECK (1: BPACK_CHECK=hodlr, compares every
# GPU step with the CPU; 2: hodlr-transpose, also the transposed products;
# slow; doc/environment_variables.md), GRID_SIZE, VIE_H,
# SCALE_GREEN, MESH (sphere mesh stem), WAVELENGTH, JOB_ID (run inside an
# existing allocation), REPO.
#
# Examples (inside an allocation of 4 nodes):
#   JOB_ID=<id> NODES=4 CASE=laplace GRID_SIZE=64 bash runit_hodlr_perlmutter_gpu.sh     # 262k unknowns
#   JOB_ID=<id> NODES=4 CASE=cfie MESH=<dir>/sphere_128000 bash runit_hodlr_perlmutter_gpu.sh
#   JOB_ID=<id> NODES=4 CASE=vie SCALE_GREEN=0 VIE_H=0.025 bash runit_hodlr_perlmutter_gpu.sh

set -uo pipefail

if [[ -n "${REPO:-}" ]]; then
  repo=$(cd "${REPO}" && pwd)
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -d "${SLURM_SUBMIT_DIR}/build_gpu/EXAMPLE" ]]; then
  repo=$(cd "${SLURM_SUBMIT_DIR}" && pwd)
else
  repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fi
bin=${repo}/build_gpu/EXAMPLE

nodes=${NODES:-${SLURM_JOB_NUM_NODES:-1}}
ranks_per_node=${RANKS_PER_NODE:-4}
ranks=$((nodes * ranks_per_node))
case_name=${CASE:-laplace}
sym=${SYM:-}
if [[ -z "${sym}" ]]; then  # (the validated defaults: symmetric Laplace, unsymmetric VIE and surface IEs)
  sym=0
  [[ "${case_name}" == laplace ]] && sym=1
fi
use_gpu=${USE_GPU:-2}
job_id=${JOB_ID:-}

module load PrgEnv-gnu cray-fftw cudatoolkit craype-accel-nvidia80
module unload cray-libsci >/dev/null 2>&1 || true

export TCMALLOC_ROOT=/global/cfs/cdirs/m2957/lib/lib/PrgEnv-gnu/gperftools-2.18.1
export OPENBLAS_LIBRARY=/global/cfs/cdirs/m2957/lib/lib/PrgEnv-gnu/OpenBLAS_sequential/build/install/lib/libopenblas.so.0
export LD_LIBRARY_PATH="${TCMALLOC_ROOT}/lib:${LD_LIBRARY_PATH:-}"
export LD_PRELOAD="${TCMALLOC_ROOT}/lib/libtcmalloc.so:${OPENBLAS_LIBRARY}"
export TCMALLOC_RELEASE_RATE=1
export OMP_NUM_THREADS=$((64 / ranks_per_node))
export OMP_DYNAMIC=FALSE
export OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export MPICH_GPU_SUPPORT_ENABLED=1  # (CUDA-aware MPI: the backend's messages go between the GPUs)
if [[ "${CHECK:-0}" == 1 ]]; then export BPACK_CHECK=hodlr; fi
if [[ "${CHECK:-0}" == 2 ]]; then export BPACK_CHECK=hodlr-transpose; fi

gpu_opts=(--HODLR_use_gpu "${use_gpu}")
if [[ -n "${PIECES:-}" ]]; then gpu_opts+=(--HODLR_gpu_pieces "${PIECES}"); fi

srun_args=(-N "${nodes}" -n "${ranks}" --ntasks-per-node="${ranks_per_node}"
  -c $((128 / ranks_per_node)) --cpu-bind=cores --gpus-per-node=4)
if [[ -n "${job_id}" ]]; then
  srun_args=(--jobid="${job_id}" "${srun_args[@]}")
fi

case "${case_name}" in
  laplace)
    grid_size=${GRID_SIZE:-32}
    cmd=("${bin}/claplace3d_h2" --grid-size "${grid_size}" --tol-comp 1e-4 --Nmin_leaf 64
         --format 1 --sym "${sym}" --elem_extract 2 --precon 1 --nrhs 2 --verbosity 0 "${gpu_opts[@]}")
    ;;
  vie)
    # --scaleGreen 1: a symmetric kernel; 0: column-scaled, not symmetric (keep SYM=0)
    cmd=("${bin}/cvie3d_h2" --ivelo 9 --scaleGreen "${SCALE_GREEN:-1}" --omega 1.4 --h "${VIE_H:-0.05}"
         --x0max 1.0 --y0max 1.0 --z0max 1.0 --L 1.0 --H 1.0 --W 1.0 --vs 1 --shape 4
         --tol_comp 1e-6 --tol_comp_s2s 1e-6 --format_s2s 1 --sym "${sym}" --Nmin_leaf 64 --precon 1
         --elem_extract 2 --verbosity 0 --lrlevel 0 "${gpu_opts[@]}")
    ;;
  efie|cfie)
    # EFIE (cfie_alpha 1, a symmetric kernel) or CFIE (cfie_alpha 0.5, not symmetric: keep SYM=0)
    alpha=1
    [[ "${case_name}" == cfie ]] && alpha=0.5
    mesh=${MESH:-${repo}/EXAMPLE/EM3D_DATA/sphere_2300}
    cmd=("${bin}/cie3d" -quant --data_dir "${mesh}" --wavelength "${WAVELENGTH:-2}" --scaling 1 --cfie_alpha "${alpha}"
         --rcs_static 2 --rcs_nsample 16
         -option --format 1 --sym "${sym}" --lrlevel 0 --tol_comp 1e-6 --tol_rand 1e-6 --nmin_leaf 128
         --xyzsort 1 --reclr_leaf 5 --baca_batch 16 --knn 10 --errsol 1 --verbosity 0 "${gpu_opts[@]}")
    ;;
  *)
    echo "unknown CASE ${case_name} (laplace, vie, efie, cfie)" >&2
    exit 2
    ;;
esac

# (each rank on its node's GPU of its local rank)
srun "${srun_args[@]}" bash -c 'export CUDA_VISIBLE_DEVICES=$((SLURM_LOCALID % 4)); exec "$@"' hodlr "${cmd[@]}" \
  |& tee "hodlr_gpu_${case_name}_sym${sym}_gpu${use_gpu}_N${nodes}_nmpi${ranks}.log"
