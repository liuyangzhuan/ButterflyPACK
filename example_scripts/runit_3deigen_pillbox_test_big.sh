#!/bin/bash -l

#SBATCH -q premium
#SBATCH -N 32
#SBATCH -t 1:00:00
#SBATCH -J paralleltest
#SBATCH --mail-user=liuyangzhuan@lbl.gov
#SBATCH -C haswell

module swap PrgEnv-intel/6.0.4 PrgEnv-gnu
NTH=1
CORES_PER_NODE=32
THREADS_PER_RANK=`expr $NTH \* 2`								 

export EXEC=./EXAMPLE/frontaldist
export OMP_NUM_THREADS=$NTH
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
  
 
knn=100
which='LM'
si=1
sample_para=2.0d0 

for lrlevel in 0 
do 
# 010 110 111
for freq in 1.06e10   
# for freq in 1.06e10   
# for freq in 1.145e9 1.146e9 1.14743e9 1.148e9 1.149e9 1.826e9 1.827e9 1.82824e9 1.829e9 1.830e9 
do
srun -n 1024 --nodes=32 -c $THREADS_PER_RANK --cpu_bind=cores ./EXAMPLE/ie3deigen -quant --data_dir ../EXAMPLE/EM3D_DATA/preprocessor_3dmesh/pillbox_100000 --freq $freq --si 1 --which LM --nev 6 --cmmode 1 -option --reclr_leaf 5 --baca_batch 64 --tol_comp 1e-4 --lrlevel $lrlevel --precon 1 --xyzsort 2 --nmin_leaf 100 --near_para 0.01d0 --pat_comp 3 --format 1 --sample_para $sample_para --knn $knn | tee a.out_freq_${freq}_lrlevel_${lrlevel}
done
done

