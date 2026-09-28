#!/bin/bash
#SBATCH --account=m2957
#SBATCH --qos=regular
#SBATCH --constraint=gpu
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --time=00:20:00
#SBATCH --job-name=emsurf_h2_gpu
#SBATCH --output=slurm-%j.out

# EMSURF EFIE on a sphere mesh (cie3d, H2 format 7 on the unstructured Color
# backend) on Perlmutter GPU nodes, with the GPU build of
# run_cmake_build_gnu_perlmutter_openblas_sequential_h2gpu.sh (in ../build_gpu).
# One MPI rank per A100, four ranks per node sharing its 64 cores.  The rank
# count must be a power of 8 (1, 8, 64, 512).
#
# Overrides: NODES, RANKS_PER_NODE, MESH (sphere size, e.g. 9000, 128000,
# 512000), MESH_DIR (mesh prefix, <prefix>_node.inp and <prefix>_elem.inp),
# RESULT_DIR, JOB_ID (run inside an existing allocation), REPO.

set -uo pipefail

if [[ -n "${REPO:-}" ]]; then
  repo=$(cd "${REPO}" && pwd)
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -x "${SLURM_SUBMIT_DIR}/build_gpu/EXAMPLE/cie3d" ]]; then
  # Slurm copies the batch script into /var/spool/slurmd, so BASH_SOURCE no
  # longer identifies the checkout when the script is submitted with sbatch.
  repo=$(cd "${SLURM_SUBMIT_DIR}" && pwd)
else
  repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fi
exe=${repo}/build_gpu/EXAMPLE/cie3d

nodes=${NODES:-2}
ranks_per_node=${RANKS_PER_NODE:-4}
ranks=$((nodes * ranks_per_node))
mesh=${MESH:-128000}
mesh_dir=${MESH_DIR:-${repo}/EXAMPLE/EM3D_DATA/preprocessor_3dmesh/sphere_${mesh}}
result_dir=${RESULT_DIR:-${repo}/build_gpu/emsurf_h2_sphere${mesh}_gpu_nmpi${ranks}}
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
if [[ ! -f "${mesh_dir}_node.inp" || ! -f "${mesh_dir}_elem.inp" ]]; then
  echo "Missing sphere mesh ${mesh_dir}_{node,elem}.inp (set MESH_DIR)" >&2
  exit 2
fi

srun_args=(-N "${nodes}" -n "${ranks}" --ntasks-per-node="${ranks_per_node}"
  -c $((128 / ranks_per_node)) --cpu-bind=cores --gpus-per-node=4)
if [[ -n "${job_id}" ]]; then
  srun_args=(--jobid="${job_id}" "${srun_args[@]}")
fi

mkdir -p "${result_dir}"
cd "${result_dir}" || exit 2
# Each rank sees only its node-local GPU.  Cray MPICH sets up its GPU
# transport on the current device in MPI_Init: with all four GPUs visible that
# is device 0 for every rank, the H2 backend then moves ranks 1-3 to their own
# GPUs, and MPICH disables CUDA IPC between the GPUs of a node ("This process
# is not using the same device it did during gtl_init"), which cost 12-22% of
# the factorization time at 8 and 64 ranks.  --H2_XRR_factor 1 (LU of X_RR)
# is required by the GPU box path; --CFIE_alpha 1 (EFIE) selects the device
# kernel of the EMSURF entries.
srun "${srun_args[@]}" bash -c 'export CUDA_VISIBLE_DEVICES=${SLURM_LOCALID}; exec "$@"' h2 \
  "${exe}" \
  --data_dir "${mesh_dir}" --wavelength 2 --model 1 --CFIE_alpha 1 --scaling 1 \
  --rcs_static 2 --rcs_nsample 1 --format 7 --sym 1 --elem_extract 2 \
  --tol_comp 1e-4 --tol_rand 1e-4 --tol_Rdetect 1e-4 --tol_itersol 1e-8 \
  --precon 1 --nmin_leaf 8 --baca_batch 64 --LR_BLAS 2 \
  --H2_unstructured 1 --H2_ID_proxy 1 --H2_ID_radius 2 --h2_lazy_schur 2 \
  --H2_CA_level 10000 --H2_GEMM_split 16 --h2_use_sketch 2 \
  --H2_XRR_factor 1 --H2_use_gpu 1 --verbosity 1 \
  > run.log 2>&1
status=$?
grep -E "total time:|data exchange communication time:|H2_CheckError|GPU solve:|Error on rank" run.log
exit "${status}"
