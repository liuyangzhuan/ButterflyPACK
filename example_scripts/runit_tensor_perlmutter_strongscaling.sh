#!/bin/bash -l

#SBATCH --account=mp127
#SBATCH -q regular
#SBATCH -N 16
#SBATCH --constraint=cpu
#SBATCH -t 10:00:00
#SBATCH -J tensorbf
#SBATCH --mail-user=liuyangzhuan@lbl.gov


module load PrgEnv-gnu
module unload cray-libsci

CORES_PER_NODE=128
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
  


if [[ $(uname -s) == 'Darwin' ]]; then
    export GPTUNEROOT=/Users/liuyangzhuan/Desktop/GPTune/
    export MPIRUN="$GPTUNEROOT/openmpi-4.1.5/bin/mpirun"
else
    export MPIRUN=mpirun
fi



# # ############## 3D seperated cubes
# NTH=8
# THREADS_PER_RANK=`expr $NTH \* 2`								 
# export OMP_NUM_THREADS=$NTH
# tol=1e-2
# nmpi=8
# # wavelen=0.001953125 0.00390625 0.0078125 0.015625 0.03125 0.0625
# zdist=1.0
# # zdist=0.001
# ppw=4.0
# nmin_leaf_t=16 #4
# nmin_leaf_m=64
# use_zfp=1
# for wavelen in 0.00390625
# do
# srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 3 --wavelen ${wavelen} --zdist ${zdist} --ppw ${ppw} -option --use_zfp ${use_zfp} --xyzsort 1 --nmin_leaf ${nmin_leaf_t} --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 0.8 --sample_para_outer 0.8 --fastsample_tensor 2 | tee a.out_tensor_3d_green_wavelen${wavelen}_zdist${zdist}_tol${tol}_mpi${nmpi}_omp${NTH}_nmin_leaf_t${nmin_leaf_t}_ppw${ppw}_use_zfp${use_zfp}
# # srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben -quant --tst 3 --wavelen ${wavelen} --zdist ${zdist} --ppw ${ppw} -option --xyzsort 1 --nmin_leaf ${nmin_leaf_m} --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 2.0 --knn 10 --sample_para_outer 2.0  | tee a.out_matrix_3d_green_wavelen${wavelen}_zdist${zdist}_tol${tol}_mpi${nmpi}_omp${NTH}_nmin_leaf_m${nmin_leaf_m}_ppw${ppw}
# done


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
# use_zfp=1
# # ppw=2.0
# # nmin_leaf_t=8
# # nmin_leaf_m=64
# for wavelen in 0.00012207031 0.00024414062
# do
# srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 2 --wavelen ${wavelen} --zdist ${zdist} --ppw ${ppw} -option --nmin_leaf ${nmin_leaf_t} --use_zfp ${use_zfp} --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 0.8 --sample_para_outer 0.8 --fastsample_tensor 2 | tee a.out_tensor_2d_green_wavelen${wavelen}_zdist${zdist}_tol${tol}_mpi${nmpi}_omp${NTH}_nmin_leaf_t${nmin_leaf_t}_ppw${ppw}_use_zfp${use_zfp}
# # srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben -quant --tst 2 --wavelen ${wavelen} --zdist ${zdist} --ppw ${ppw} -option --nmin_leaf ${nmin_leaf_m}  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 | tee a.out_matrix_2d_green_wavelen${wavelen}_zdist${zdist}_tol${tol}_mpi${nmpi}_omp${NTH}_nmin_leaf_m${nmin_leaf_m}_ppw${ppw}
# done 
# # srun -n 4 -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 2 --wavelen 0.0156 --zdist 1.0 --ppw 2.0 -option --nmin_leaf 8 --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 0.8 --sample_para_outer 0.8 --fastsample_tensor 2 | tee a.out_tensor_2d_green_wavelen${wavelen}_zdist${zdist}_tol${tol}_mpi${nmpi}_omp${NTH}
# # # srun -n 4 -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 2 --wavelen 0.0156 --zdist 1.0 --ppw 2.0 -option --nmin_leaf 8 --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 1.0 --sample_para_outer 1.0 --fastsample_tensor 1 | tee a.out_tensor_2d_green_wavelen${wavelen}_zdist${zdist}_tol${tol}_mpi${nmpi}_omp${NTH}
# # srun -n 4 -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 2 --wavelen 0.0156 --zdist 1.0 --ppw 2.0 -option --nmin_leaf 8 --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 0.8 --sample_para_outer 0.8 | tee a.out_tensor_2d_green_wavelen${wavelen}_zdist${zdist}_tol${tol}_mpi${nmpi}_omp${NTH}


# ############## DFT
# NTH=16
# THREADS_PER_RANK=`expr $NTH \* 2`								 
# export OMP_NUM_THREADS=$NTH
# tol=1e-8
# nmpi=32
# ndim_FIO=3
# use_zfp=1
# for N_FIO in 32
# do
# # srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben -quant --tst 4 --ndim_FIO ${ndim_FIO} --N_FIO ${N_FIO} -option --nmin_leaf 64  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 | tee a.out_matrix_DFT_ndim_FIO${ndim_FIO}_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}
# srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 4 --ndim_FIO ${ndim_FIO} --N_FIO ${N_FIO} -option --nmin_leaf 2  --xyzsort 1 --use_zfp ${use_zfp} --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 1.0 --sample_para_outer 1.0 --fastsample_tensor 2 | tee a.out_tensor_DFT_ndim_FIO${ndim_FIO}_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}_use_zfp${use_zfp}
# done


############## NUFFT
### sample_para needs to be at least 10 to have accurate compression results. Not sure how to resolve this yet. 
NTH=16
THREADS_PER_RANK=`expr $NTH \* 2`								 
export OMP_NUM_THREADS=$NTH
tol=1e-5
nmpi=128
ndim_FIO=3
use_zfp=1
for N_FIO in 256
do
# srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben -quant --tst 10 --ndim_FIO ${ndim_FIO} --N_FIO ${N_FIO} -option --nmin_leaf 64  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 | tee a.out_matrix_DFT_ndim_FIO${ndim_FIO}_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}
srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 10 --ndim_FIO ${ndim_FIO} --N_FIO ${N_FIO} -option --nmin_leaf 2  --xyzsort 0 --use_zfp ${use_zfp} --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 0.8 --sample_para_outer 0.8 --fastsample_tensor 2 | tee a.out_tensor_DFT_ndim_FIO${ndim_FIO}_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}_use_zfp${use_zfp}
done


# ############## 2D Radon for ellipse integral
# NTH=8
# THREADS_PER_RANK=`expr $NTH \* 2`								 
# export OMP_NUM_THREADS=$NTH
# tol=1e-3
# nmpi=64
# use_zfp=1
# # N_FIO=256 512 1024 2048 4096
# for N_FIO in 1024
# do
# # srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben -quant --tst 5 --N_FIO ${N_FIO} -option --nmin_leaf 16  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 | tee a.out_matrix_Radon2D_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}
# srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 5 --N_FIO ${N_FIO} -option --nmin_leaf 4  --xyzsort 1 --use_zfp ${use_zfp} --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 0.6 --sample_para_outer 0.6 --fastsample_tensor 0 | tee a.out_tensor_Radon2D_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}_use_zfp${use_zfp}_fastsample_tensor0
# # srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 5 --N_FIO ${N_FIO} -option --nmin_leaf 4  --xyzsort 1 --use_zfp ${use_zfp} --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 0.8 --sample_para_outer 0.8 --fastsample_tensor 2 | tee a.out_tensor_Radon2D_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}_use_zfp${use_zfp}_fastsample_tensor2
# done



# ############## 3D Radon for plane integral
# # NTH=1
# # THREADS_PER_RANK=`expr $NTH \* 2`								 
# # export OMP_NUM_THREADS=$NTH
# # tol=1e-3
# # nmpi=512
# # use_zfp=1
# # # N_FIO=256 512 1024 2048 4096
# # for N_FIO in 64 128 256
# # do
# # srun -N 16 -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben -quant --tst 6 --N_FIO ${N_FIO} -option --nmin_leaf 8  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 | tee a.out_matrix_Radon3D_plane_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}
# # # srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 6 --N_FIO ${N_FIO} -option --nmin_leaf 2  --xyzsort 1 --use_zfp ${use_zfp} --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 0.6 --sample_para_outer 0.6 --fastsample_tensor 0 | tee a.out_tensor_Radon3D_plane_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}_use_zfp${use_zfp}_fastsample_tensor0
# # # srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 6 --N_FIO ${N_FIO} -option --nmin_leaf 4  --xyzsort 1 --use_zfp ${use_zfp} --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 0.8 --sample_para_outer 0.8 --fastsample_tensor 2 | tee a.out_tensor_Radon3D_plane_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}_use_zfp${use_zfp}_fastsample_tensor2
# # done

# NTH=8
# THREADS_PER_RANK=`expr $NTH \* 2`								 
# export OMP_NUM_THREADS=$NTH
# tol=1e-3
# nmpi=64
# use_zfp=1
# # N_FIO=256 512 1024 2048 4096
# for N_FIO in 64 
# do
# # srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben -quant --tst 6 --N_FIO ${N_FIO} -option --nmin_leaf 8  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 | tee a.out_matrix_Radon3D_plane_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}
# srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 6 --N_FIO ${N_FIO} -option --nmin_leaf 2  --xyzsort 1 --use_zfp ${use_zfp} --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 0.6 --sample_para_outer 0.6 --fastsample_tensor 0 | tee a.out_tensor_Radon3D_plane_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}_use_zfp${use_zfp}_fastsample_tensor0
# # srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 6 --N_FIO ${N_FIO} -option --nmin_leaf 2  --xyzsort 1 --use_zfp ${use_zfp} --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 0.8 --sample_para_outer 0.8 --fastsample_tensor 2 | tee a.out_tensor_Radon3D_plane_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}_use_zfp${use_zfp}_fastsample_tensor2
# done




# ############## 3D Radon for sphere integral
# NTH=8
# THREADS_PER_RANK=`expr $NTH \* 2`								 
# export OMP_NUM_THREADS=$NTH
# tol=1e-3
# nmpi=256
# use_zfp=1
# # N_FIO=256 512 1024 2048 4096
# for N_FIO in 64 128 256 512
# do
# # srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben -quant --tst 9 --N_FIO ${N_FIO} -option --nmin_leaf 16  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 | tee a.out_matrix_Radon3D_sphere_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}
# # srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 9 --N_FIO ${N_FIO} -option --nmin_leaf 2  --xyzsort 1 --use_zfp ${use_zfp} --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 0.6 --sample_para_outer 0.6 --fastsample_tensor 0 | tee a.out_tensor_Radon3D_sphere_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}_use_zfp${use_zfp}_fastsample_tensor0
# srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 9 --N_FIO ${N_FIO} -option --nmin_leaf 2  --xyzsort 1 --use_zfp ${use_zfp} --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 0.6 --sample_para_outer 0.6 --fastsample_tensor 2 | tee a.out_tensor_Radon3D_sphere_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}_use_zfp${use_zfp}_fastsample_tensor2
# done


# ############## 2D Radon for line integral
# NTH=8
# THREADS_PER_RANK=`expr $NTH \* 2`								 
# export OMP_NUM_THREADS=$NTH
# tol=1e-3
# nmpi=64
# use_zfp=1
# # N_FIO=256 512 1024 2048 4096
# for N_FIO in 1024
# do
# # srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben -quant --tst 7 --N_FIO ${N_FIO} -option --nmin_leaf 16  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 | tee a.out_matrix_Radon3D_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}
# srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 7 --N_FIO ${N_FIO} -option --nmin_leaf 4  --xyzsort 1 --use_zfp ${use_zfp} --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 0.6 --sample_para_outer 0.6 --fastsample_tensor 0 | tee a.out_tensor_Radon2DYing_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}_use_zfp${use_zfp}_fastsample_tensor0
# # srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 7 --N_FIO ${N_FIO} -option --nmin_leaf 4  --xyzsort 1 --use_zfp ${use_zfp} --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 0.8 --sample_para_outer 0.8 --fastsample_tensor 2 | tee a.out_tensor_Radon2DYing_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}_use_zfp${use_zfp}_fastsample_tensor2
# done

# ############## 2D hyperbolic Radon 
# NTH=8
# THREADS_PER_RANK=`expr $NTH \* 2`								 
# export OMP_NUM_THREADS=$NTH
# tol=1e-3
# nmpi=64
# use_zfp=1
# # N_FIO=256 512 1024 2048 4096
# for N_FIO in 4096
# do
# # srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben -quant --tst 8 --N_FIO ${N_FIO} -option --nmin_leaf 16  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 | tee a.out_matrix_Radon3D_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}
# srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 8 --N_FIO ${N_FIO} -option --nmin_leaf 4  --xyzsort 1 --use_zfp ${use_zfp} --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 0.6 --sample_para_outer 0.6 --fastsample_tensor 0 | tee a.out_tensor_Radon2Dhyperbolic_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}_use_zfp${use_zfp}_fastsample_tensor0
# # srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/frankben_t -quant --tst 8 --N_FIO ${N_FIO} -option --nmin_leaf 4  --xyzsort 1 --use_zfp ${use_zfp} --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 0.8 --sample_para_outer 0.8 --fastsample_tensor 2 | tee a.out_tensor_Radon2Dhyperbolic_N_FIO${N_FIO}_tol${tol}_mpi${nmpi}_omp${NTH}_use_zfp${use_zfp}_fastsample_tensor2
# done