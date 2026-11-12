#!/bin/bash
#SBATCH --job-name=laplace3d_h2_partial_gpu
#SBATCH -A m2957
#SBATCH --constraint=gpu
#SBATCH --qos=regular
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --time=00:15:00
#SBATCH --output=./laplace3d_h2_partial_gpu_%j.log

# Partial H2 factorization of the 3D Laplace kernel on Perlmutter GPU nodes
# (--H2_partial: the first H2_PARTIAL grid points eliminated, the Schur
# complement on the rest served through c_bpack_partial_* and
# c_bpack_extractelement), then --partial_extract_bench rounds of Schur
# complement blocks read as a STRUMPACK parent front reads a child's update
# block, on every rank.  Built by
# run_cmake_build_gnu_perlmutter_openblas_sequential_h2gpu.sh (../build_gpu).
#
# NCCL.  Each block's sum over the ranks that hold the factorization runs on
# the device through NCCL (one ncclReduce a landing rank) when the build has
# NCCL (TPL_NCCL_*) and the group can use it: always on one node, across nodes
# only with NCCL's network plugin, which the nccl module brings (NCCL_NET).
# NCCL=0 leaves the module out: a group across nodes then sums on the host
# (MPI_Reduce_scatter after downloading every rank's whole contribution).  A
# STRUMPACK H2 front run (Poisson 400^3, 64 A100 on 16 nodes) took 52.5 s
# without the device sums and 49.3 s with them; this script's defaults (48^3,
# 2 nodes) 33.4 ms a bench round with NCCL=0, 22.8 ms with NCCL=1.  The first
# sum on a group makes its communicator (0.2-0.6 s).
#
# Overrides: NODES, RANKS_PER_NODE, GRID_SIZE, TOL, H2_PARTIAL (default: all
# but the last quarter of the grid's planes, at least 16 - the bench needs 16
# whole F22 planes), ROUNDS, USE_GPU, NCCL (1: load the nccl module), CHECK (1: BPACK_CHECK=solve,
# every device sum also against the host's, "sum check" lines), TRACE (1:
# BPACK_TRACE=extract, a line per extraction call with its sum's time), JOB_ID
# (run inside an existing allocation), REPO.  Environment variables:
# doc/environment_variables.md.

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

nodes=${NODES:-2}
ranks_per_node=${RANKS_PER_NODE:-4}
ranks=$((nodes * ranks_per_node))
grid_size=${GRID_SIZE:-48}
tol=${TOL:-1e-4}
f22_planes=$((grid_size / 4 > 16 ? grid_size / 4 : 16))
h2_partial=${H2_PARTIAL:-$((grid_size * grid_size * (grid_size - f22_planes)))}
rounds=${ROUNDS:-20}
use_gpu=${USE_GPU:-1}
nccl=${NCCL:-1}
job_id=${JOB_ID:-}

module load PrgEnv-gnu cray-fftw cudatoolkit craype-accel-nvidia80
module unload cray-libsci >/dev/null 2>&1 || true
if [[ "${nccl}" == 1 ]]; then
  # NCCL_NET (its Slingshot plugin) and the plugin's path; it also sets
  # FI_CXI_DISABLE_HOST_REGISTER and FI_CXI_RDZV_THRESHOLD, which reach Cray
  # MPICH too (no slowdown measured)
  module load nccl/2.29.2-cu13
fi

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
[[ "${CHECK:-0}" == 1 ]] && export BPACK_CHECK=solve
[[ "${TRACE:-0}" == 1 ]] && export BPACK_TRACE=extract

if [[ ! -x "${exe}" ]]; then
  echo "Missing executable: ${exe}" >&2
  exit 2
fi

srun_args=(-N "${nodes}" -n "${ranks}" --ntasks-per-node="${ranks_per_node}"
  -c $((128 / ranks_per_node)) --cpu-bind=cores --gpus-per-node=4)
if [[ -n "${job_id}" ]]; then
  srun_args=(--jobid="${job_id}" "${srun_args[@]}")
fi

# Each rank sees only its node-local GPU (see runit_3dlaplace_perlmutter_h2_gpu.sh).
# --H2_partial needs the unstructured tree (--H2_unstructured 1) and sym 0;
# --H2_XRR_factor 1 is required by the GPU box path.
srun "${srun_args[@]}" bash -c 'export CUDA_VISIBLE_DEVICES=$((SLURM_LOCALID % ${SLURM_GPUS_ON_NODE:-4})); exec "$@"' h2 \
  "${exe}" \
  --grid-size "${grid_size}" --tol-comp "${tol}" --Nmin_leaf 64 --reduction_threshold 8 --sym 0 --skew 0.5 \
  --CA_level 10000 --elem_extract 2 --verbosity 1 --distributed64 1 \
  --h2_use_sketch 2 --h2_lazy_schur 2 --H2_unstructured 1 --H2_XRR_factor 1 --H2_use_gpu "${use_gpu}" \
  --H2_partial "${h2_partial}" --partial_extract_bench "${rounds}" \
  |& tee "laplace3d_h2_partial_gpu_${grid_size}_N${nodes}_nmpi${ranks}_gpu${use_gpu}_nccl${nccl}.log"
