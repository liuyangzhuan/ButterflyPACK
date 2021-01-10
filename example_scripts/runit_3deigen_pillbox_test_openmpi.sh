#!/bin/bash -l

#SBATCH -q premium
#SBATCH -N 2
#SBATCH -t 10:00:00
#SBATCH -J paralleltest
#SBATCH --mail-user=liuyangzhuan@lbl.gov
#SBATCH -C haswell


module unload cray-mpich
module swap PrgEnv-intel PrgEnv-gnu
export MKLROOT=/opt/intel/compilers_and_libraries_2019.3.199/linux/mkl
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/compilers_and_libraries_2019.3.199/linux/mkl/lib/intel64
module load openmpi/4.0.1
CCC=mpicc
CCCPP=mpicxx
FTN=mpif90


NTH=8
CORES_PER_NODE=32
THREADS_PER_RANK=`expr $NTH \* 2`								 

export EXEC=./EXAMPLE/frontaldist
export OMP_NUM_THREADS=$NTH
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
  
 
knn=100
which='LM'
nev=20
si=1
postprocess=0
sample_para=2.0d0 
lrlevel=0
mesh_normal=1
verbosity=0


# model="pillbox_4000"
# mesh_normal=-1
# # for freq in 1.826e9    
# for freq in 1.145e9 
# # for freq in 1.51e9    
# # for freq in 1.145e9 1.146e9 1.14743e9 1.148e9 1.149e9 1.826e9 1.827e9 1.82824e9 1.829e9 1.830e9 


# model="cavity_wakefield_4K_feko"
# for freq in 1136.1e6 



model="cavity_rec_17K_feko"
# for freq in 2141.3e6 1510.0e6 2321.0e6 2200.0e6
# for freq in 7495.0e6
# for freq in 1510.0e6
for freq in 2311.0e6
# for freq in 2200.0e6

# model="cavity_5cell_30K_feko"
# for freq in 9.3773e8
# for freq in 6.3377e8 6.38132e8 6.4362e8 6.4816e8 6.4996e8

do
srun -n 128 --nodes=32 -c $THREADS_PER_RANK --cpu_bind=cores ./EXAMPLE/ie3dporteigen -quant --data_dir ../EXAMPLE/EM3D_DATA/preprocessor_3dmesh/$model --freq $freq --si $si --which $which --postprocess ${postprocess} --nev ${nev} --mesh_normal ${mesh_normal} --cmmode 0 --model ${model} -option --verbosity 2 --reclr_leaf 5 --baca_batch 64 --tol_comp 1e-4 --lrlevel $lrlevel --precon 1 --verbosity ${verbosity} --xyzsort 2 --nmin_leaf 100 --near_para 0.01d0 --pat_comp 3 --format 1 --sample_para $sample_para --knn $knn | tee a.out_freq_${freq}_lrlevel_${lrlevel}
done

