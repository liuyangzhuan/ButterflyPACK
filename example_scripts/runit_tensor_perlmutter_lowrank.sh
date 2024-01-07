#!/bin/bash -l

#SBATCH -q premium
#SBATCH -N 2
#SBATCH -t 10:00:00
#SBATCH -J paralleltest
#SBATCH --mail-user=liuyangzhuan@lbl.gov
#SBATCH -C haswell
module load PrgEnv-gnu
CORES_PER_NODE=128
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
  


if [[ $(uname -s) == 'Darwin' ]]; then
    export GPTUNEROOT=/Users/liuyangzhuan/Desktop/GPTune/
    export MPIRUN="$GPTUNEROOT/openmpi-4.1.5/bin/mpirun"
else
    export MPIRUN=mpirun
fi






# ############## 3D seperated cubes
NTH=1
THREADS_PER_RANK=`expr $NTH \* 2`								 
export OMP_NUM_THREADS=$NTH
tol=1e-2
nmpi=1
wavelen=0.03125
zdist=1.0
ppw=4.0
nmin_leaf_t=16
nmin_leaf_m=4096
lrlevel=0
srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 3 --wavelen ${wavelen} --zdist ${zdist} --ppw ${ppw} -option --xyzsort 1 --nmin_leaf ${nmin_leaf_t} --lrlevel ${lrlevel} --verbosity 1 --tol_comp $tol --sample_para 0.8 --sample_para_outer 0.8 --fastsample_tensor 2 | tee a.out_tensor_3d_green_wavelen${wavelen}_zdist${zdist}_tol${tol}_mpi${nmpi}_omp${NTH}_nmin_leaf_t${nmin_leaf_t}_ppw${ppw}_lrlevel${lrlevel}
# srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben -quant --tst 3 --wavelen ${wavelen} --zdist ${zdist} --ppw ${ppw} -option --xyzsort 1 --nmin_leaf ${nmin_leaf_m} --lrlevel ${lrlevel} --verbosity 1 --tol_comp $tol --sample_para 2.0 --sample_para_outer 2.0  | tee a.out_matrix_3d_green_wavelen${wavelen}_zdist${zdist}_tol${tol}_mpi${nmpi}_omp${NTH}_nmin_leaf_m${nmin_leaf_m}_ppw${ppw}_lrlevel${lrlevel}



# ############## 2D parallel plates
# NTH=8
# THREADS_PER_RANK=`expr $NTH \* 2`								 
# export OMP_NUM_THREADS=$NTH
# tol=1e-6
# nmpi=64
# # wavelen=0.00012207031 #0.00024414062 #0.00048828125 #0.0009765625 #0.001953125 #0.00390625 #0.0078125 #0.015625
# zdist=1.0
# ppw=4.0
# nmin_leaf_t=16
# nmin_leaf_m=256
# lrlevel=0
# # ppw=2.0
# # nmin_leaf_t=8
# # nmin_leaf_m=64
# for wavelen in 0.015625 #0.0078125 # 0.00390625 #
# do
# # srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 2 --wavelen ${wavelen} --zdist ${zdist} --ppw ${ppw} -option --nmin_leaf ${nmin_leaf_t} --xyzsort 1 --lrlevel ${lrlevel} --verbosity 1 --tol_comp $tol --sample_para 0.8 --sample_para_outer 0.8 --fastsample_tensor 2 | tee a.out_tensor_2d_green_wavelen${wavelen}_zdist${zdist}_tol${tol}_mpi${nmpi}_omp${NTH}_nmin_leaf_t${nmin_leaf_t}_ppw${ppw}_lrlevel${lrlevel}
# srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben -quant --tst 2 --wavelen ${wavelen} --zdist ${zdist} --ppw ${ppw} -option --nmin_leaf ${nmin_leaf_m}  --xyzsort 1 --lrlevel ${lrlevel} --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 | tee a.out_matrix_2d_green_wavelen${wavelen}_zdist${zdist}_tol${tol}_mpi${nmpi}_omp${NTH}_nmin_leaf_m${nmin_leaf_m}_ppw${ppw}_lrlevel${lrlevel}
# done 
# # srun -n 4 -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 2 --wavelen 0.0156 --zdist 1.0 --ppw 2.0 -option --nmin_leaf 8 --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 0.8 --sample_para_outer 0.8 --fastsample_tensor 2 | tee a.out_tensor_2d_green_wavelen${wavelen}_zdist${zdist}_tol${tol}_mpi${nmpi}_omp${NTH}
# # # srun -n 4 -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 2 --wavelen 0.0156 --zdist 1.0 --ppw 2.0 -option --nmin_leaf 8 --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 1.0 --sample_para_outer 1.0 --fastsample_tensor 1 | tee a.out_tensor_2d_green_wavelen${wavelen}_zdist${zdist}_tol${tol}_mpi${nmpi}_omp${NTH}
# # srun -n 4 -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 2 --wavelen 0.0156 --zdist 1.0 --ppw 2.0 -option --nmin_leaf 8 --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 0.8 --sample_para_outer 0.8 | tee a.out_tensor_2d_green_wavelen${wavelen}_zdist${zdist}_tol${tol}_mpi${nmpi}_omp${NTH}


# ############## DFT
# NTH=16
# THREADS_PER_RANK=`expr $NTH \* 2`								 
# export OMP_NUM_THREADS=$NTH
# tol=1e-3
# nmpi=32
# ndim_FIO=4

# N_FIO=16
# # srun -n 4 -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben -quant --tst 4 --ndim_FIO ${ndim_FIO} --N_FIO ${N_FIO} -option --nmin_leaf 64  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 | tee a.out_matrix_DFT_ndim_FIO${ndim_FIO}_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}
# srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 4 --ndim_FIO ${ndim_FIO} --N_FIO ${N_FIO} -option --nmin_leaf 8  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 0.8 --sample_para_outer 0.8 --fastsample_tensor 2 | tee a.out_tensor_DFT_ndim_FIO${ndim_FIO}_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}

# N_FIO=32
# # srun -n 4 -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben -quant --tst 4 --ndim_FIO ${ndim_FIO} --N_FIO ${N_FIO} -option --nmin_leaf 64  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 | tee a.out_matrix_DFT_ndim_FIO${ndim_FIO}_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}
# srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 4 --ndim_FIO ${ndim_FIO} --N_FIO ${N_FIO} -option --nmin_leaf 8  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 0.8 --sample_para_outer 0.8 --fastsample_tensor 2 | tee a.out_tensor_DFT_ndim_FIO${ndim_FIO}_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}



