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

# H2 factorization of the 3D Laplace kernel on Perlmutter GPU nodes, with the
# Color algorithm on every level (METHOD=color, the default) or the replicated
# communication-avoiding algorithm CA0 on the finest levels (METHOD=ca0), with
# the GPU build of run_cmake_build_gnu_perlmutter_openblas_sequential_h2gpu.sh
# (in ../build_gpu).  One MPI rank per A100, four ranks per node sharing its 64
# cores.  The rank count must be a power of 8 (1, 8, 64, 512).
#
# Overrides: NODES, RANKS_PER_NODE, GRID_SIZE, NMIN_LEAF, USE_GPU (1: FP64
# GEMMs, 2: FP64 tensor-core GEMMs), METHOD (color, ca0), CA_LEVEL, JOB_ID (run
# inside an existing allocation), REPO.
#
# CA levels.  CA_LEVEL is the first (coarsest) level that uses CA: levels
# >= CA_LEVEL use CA0, levels above it Color (10000: Color everywhere).  The
# leaf level follows from the grid and Nmin_leaf: the leaf boxes have edge
# d = round(Nmin_leaf^(1/3)) points or a little more, and the leaf level is
# log2(grid / d) rounded down (Nmin_leaf 216: grid 192 -> leaf 5, 384 -> 6,
# 768 -> 7).  METHOD=ca0 sets CA_LEVEL = leaf - 1 (CA on the two finest
# levels); CA_LEVEL=<leaf> gives CA on the leaf only.  A CA level needs at
# least 343 (7^3) boxes per rank, otherwise it falls back to Color by itself.
# On the GPU only CA0 is supported (--H2_CA_owner_component 0); CA3 (3) runs
# on the CPU.
#
# Memory.  --constraint=gpu can give 40 GB or 80 GB A100 nodes; for comparable
# timings use --constraint="gpu&hbm40g" or "gpu&hbm80g".  With Nmin_leaf 216,
# 884K unknowns per rank (e.g. 384^3 on 64 ranks) needs the 80 GB nodes;
# ~600K or less fits in 40 GB.  If a CA0 run fails with "H2 GPU heap
# exhausted" after its CA leaf (the leaf's solve factors kept on the device
# fragment the heap), run it with H2_GPU_SOLVE_KEEP=0.

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
nmin_leaf=${NMIN_LEAF:-216}
use_gpu=${USE_GPU:-2}
method=${METHOD:-color}
job_id=${JOB_ID:-}

# the leaf level, as calc_num_levels (h2_parallel/butterfly_init.hpp) sets it
edge=$(awk -v n="${nmin_leaf}" 'BEGIN { printf "%d", n^(1/3) + 0.5 }')
levels=0
while (( edge <= grid_size )); do
  edge=$((edge * 2))
  levels=$((levels + 1))
done
leaf_level=$((levels - 1))
case "${method}" in
  color) ca_level=${CA_LEVEL:-10000} ;;
  ca0)   ca_level=${CA_LEVEL:-$((leaf_level - 1))} ;;
  *) echo "METHOD must be color or ca0" >&2; exit 2 ;;
esac
echo "grid ${grid_size}, Nmin_leaf ${nmin_leaf}: leaf level ${leaf_level}; METHOD=${method}, CA_level ${ca_level}"

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
# load the batched GPU kernels before the first level, so first-use module
# loading does not land inside a level's time
export H2_GPU_WARMUP=1

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
# GPU.  --H2_CA_staged_halo 0: a GPU CA level gathers its halo before the
# level (the staged overlap of mode 2 is a CPU feature).
srun "${srun_args[@]}" bash -c 'export CUDA_VISIBLE_DEVICES=$((SLURM_LOCALID % ${SLURM_GPUS_ON_NODE:-4})); exec "$@"' h2 \
  "${exe}" \
  --grid-size "${grid_size}" --tol-comp 1e-3 --Nmin_leaf "${nmin_leaf}" --reduction_threshold 8 --sym 1 \
  --CA_level "${ca_level}" --elem_extract 2 --verbosity 1 --distributed64 1 \
  --h2_use_sketch 2 --h2_lazy_schur 2 --h2_gemm_split 16 \
  --H2_CA_staged_halo 0 --H2_CA_owner_component 0 \
  --H2_XRR_factor 1 --H2_use_gpu "${use_gpu}" \
  |& tee "laplace3d_h2_gpu_${grid_size}_N${nodes}_nmpi${ranks}_gpu${use_gpu}_${method}_calv${ca_level}.log"
