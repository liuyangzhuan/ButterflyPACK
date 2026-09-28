#!/bin/bash
#SBATCH --job-name=laplace3d_h2_gpu
#SBATCH -A m2957
#SBATCH --constraint=gpu
#SBATCH --qos=regular
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --time=00:30:00
#SBATCH --output=./laplace3d_h2_gpu_%j.log

# H2 Color factorization of the 3D Laplace kernel on Perlmutter GPU nodes,
# with the GPU build of run_cmake_build_gnu_perlmutter_openblas_sequential_h2gpu.sh
# (in ../build_gpu).  One MPI rank per A100, four ranks per node sharing its 64
# cores.  The rank count must be a power of 8 (1, 8, 64, 512).
#
# Overrides: NODES, RANKS_PER_NODE, GRID_SIZE, USE_GPU (1: FP64 GEMMs, 2: FP64
# tensor-core GEMMs), JOB_ID (run inside an existing allocation), REPO.

set -uo pipefail

if [[ -n "${REPO:-}" ]]; then
  repo=$(cd "${REPO}" && pwd)
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -x "${SLURM_SUBMIT_DIR}/build_gpu/EXAMPLE/claplace3d_h2" ]]; then
  # Slurm copies the batch script into /var/spool/slurmd, so BASH_SOURCE no
  # longer identifies the checkout when the script is submitted with sbatch.
  repo=$(cd "${SLURM_SUBMIT_DIR}" && pwd)
else
  repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fi
exe=${repo}/build_gpu/EXAMPLE/claplace3d_h2

nodes=${NODES:-16}
ranks_per_node=${RANKS_PER_NODE:-4}
ranks=$((nodes * ranks_per_node))
grid_size=${GRID_SIZE:-192}
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
export MPICH_ASYNC_PROGRESS=1
# CUDA-aware MPI: the GPU backend sends its messages from device memory
export MPICH_GPU_SUPPORT_ENABLED=1

if [[ ! -x "${exe}" ]]; then
  echo "Missing executable: ${exe}" >&2
  exit 2
fi

srun_args=(-N "${nodes}" -n "${ranks}" --ntasks-per-node="${ranks_per_node}"
  -c $((128 / ranks_per_node)) --cpu-bind=cores --gpus-per-node=4)
if [[ -n "${job_id}" ]]; then
  srun_args=(--jobid="${job_id}" "${srun_args[@]}")
fi

# Each rank sees only its node-local GPU (shared round-robin when a node
# runs more ranks than GPUs).  Cray MPICH sets up its GPU
# transport on the current device in MPI_Init: with all four GPUs visible that
# is device 0 for every rank, the H2 backend then moves ranks 1-3 to their own
# GPUs, and MPICH disables CUDA IPC between the GPUs of a node ("This process
# is not using the same device it did during gtl_init"), which cost 12-22% of
# the factorization time at 8 and 64 ranks.  --H2_XRR_factor 1 (LU of X_RR)
# is required by the GPU box path; without it only the owner pass runs on the
# GPU.
srun "${srun_args[@]}" bash -c 'export CUDA_VISIBLE_DEVICES=$((SLURM_LOCALID % ${SLURM_GPUS_ON_NODE:-4})); exec "$@"' h2 \
  "${exe}" \
  --grid-size "${grid_size}" --tol-comp 1e-3 --Nmin_leaf 216 --reduction_threshold 8 --sym 1 \
  --CA_level 10000 --elem_extract 2 --verbosity 1 --distributed64 1 \
  --h2_use_sketch 2 --h2_lazy_schur 2 --h2_gemm_split 16 \
  --H2_CA_staged_halo 2 --H2_CA_owner_component 0 \
  --H2_XRR_factor 1 --H2_use_gpu "${use_gpu}" \
  |& tee "laplace3d_h2_gpu_${grid_size}_N${nodes}_nmpi${ranks}_gpu${use_gpu}.log"
