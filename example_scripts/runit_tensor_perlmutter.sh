#!/bin/bash -l

#SBATCH -q premium
#SBATCH -N 2
#SBATCH -t 10:00:00
#SBATCH -J paralleltest
#SBATCH --mail-user=liuyangzhuan@lbl.gov
#SBATCH -C haswell

module load PrgEnv-gnu
NTH=32
CORES_PER_NODE=128
THREADS_PER_RANK=`expr $NTH \* 2`								 

export OMP_NUM_THREADS=$NTH
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
  


if [[ $(uname -s) == 'Darwin' ]]; then
    export GPTUNEROOT=/Users/liuyangzhuan/Desktop/GPTune/
    export MPIRUN="$GPTUNEROOT/openmpi-4.1.5/bin/mpirun"
else
    export MPIRUN=mpirun
fi



tol=1e-6
export OMP_NUM_THREADS=$NTH

############## 3D seperated cubes
# srun -n 4 -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 3 --wavelen 0.015625 -option --nmin_leaf 8 --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 1.0 --sample_para_outer 1.0 --fastsample_tensor 2 | tee a.out_tensor_3d_green
# srun -n 4 -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben -quant --tst 3 --wavelen 0.015625 -option --nmin_leaf 8 --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 2.0 --sample_para_outer 2.0  | tee a.out_matrix_3d_green



############## 2D parallel plates
# srun -n 4 -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 2 --wavelen 0.0156 -option --nmin_leaf 8 --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 0.8 --sample_para_outer 0.8 --fastsample_tensor 2 | tee a.out_tensor_2d_green
# # srun -n 4 -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 2 --wavelen 0.0156 -option --nmin_leaf 8 --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 1.0 --sample_para_outer 1.0 --fastsample_tensor 1 | tee a.out_tensor_2d_green
# srun -n 4 -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 2 --wavelen 0.0156 -option --nmin_leaf 8 --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 0.8 --sample_para_outer 0.8 | tee a.out_tensor_2d_green
# srun -n 4 -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben -quant --tst 2 --wavelen 0.0156 -option --nmin_leaf 64  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 | tee a.out_matrix_2d_green


############## DFT
tol=1e-3
# srun -n 4 -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben -quant --tst 4 --ndim_FIO 4 --N_FIO 8 -option --nmin_leaf 64  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 | tee a.out_matrix_DFT
srun -N 4 -n 8 -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 4 --ndim_FIO 4 --N_FIO 16 -option --nmin_leaf 8  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 --fastsample_tensor 2 | tee a.out_tensor_DFT




