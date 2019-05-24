#!/bin/bash -l

#SBATCH -q premium
#SBATCH -N 2
#SBATCH -t 10:00:00
#SBATCH -J paralleltest
#SBATCH --mail-user=liuyangzhuan@lbl.gov
#SBATCH -C haswell

module swap PrgEnv-intel/6.0.4 PrgEnv-gnu
NTH=8
CORES_PER_NODE=32
THREADS_PER_RANK=`expr $NTH \* 2`								 

export EXEC=./EXAMPLE/frontaldist
export OMP_NUM_THREADS=$NTH
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
 
 
knn=10
which='LM'
si=1
sample_para=2.0d0 
 
# 010 110 111
for freq in 1.15e9 1.83e9 2.36e9    
# for freq in 1.145e9 1.146e9 1.14743e9 1.148e9 1.149e9 1.826e9 1.827e9 1.82824e9 1.829e9 1.830e9 
do
srun -n 32 --nodes=8 -c $THREADS_PER_RANK --cpu_bind=cores ./EXAMPLE/ie3deigen -quant --data_dir ../EXAMPLE/EM3D_DATA/preprocessor_3dmesh/pillbox_4000 --freq $freq --si 1 --which LM --nev 1 --cmmode 0 -option --tol_comp 1e-4 --lrlevel 100 --precon 1 --xyzsort 2 --nmin_leaf 100 --near_para 4.01d0 --pat_comp 3 --format 1 --sample_para $sample_para --knn $knn | tee a.out_freq_$freq
done

