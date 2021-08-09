#!/bin/bash -l

#SBATCH -q premium
#SBATCH -N 2
#SBATCH -t 10:00:00
#SBATCH -J paralleltest
#SBATCH --mail-user=liuyangzhuan@lbl.gov
#SBATCH -C haswell

module swap PrgEnv-intel PrgEnv-gnu
NTH=4
CORES_PER_NODE=32
THREADS_PER_RANK=`expr $NTH \* 2`								 

export OMP_NUM_THREADS=$NTH
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
  
 
knn=100
which='LM'
si=1
sample_para=2.0d0 
lrlevel=0
noport=0
model="pillbox_4000"
# model="cavity_wakefield_4K_feko"
# model="cavity_rec_5K_feko"
# model="cavity_rec_17K_feko"

# 010 110 111
# for freq in 1.826e9    
# for freq in 765e6 760e6 755e6 750e6   
for freq in 1.15e9
# for freq in 1.51e9    
# for freq in 1.145e9 1.146e9 1.14743e9 1.148e9 1.149e9 1.826e9 1.827e9 1.82824e9 1.829e9 1.830e9 
do
srun -n 128 -c $THREADS_PER_RANK --cpu_bind=cores ./EXAMPLE/ie3dporteigen -quant --data_dir ../EXAMPLE/EM3D_DATA/preprocessor_3dmesh/$model --model $model --freq $freq --si $si --which $which --noport $noport --nev 6 --cmmode 0 --postprocess 1 -option --verbosity 2 --reclr_leaf 5 --baca_batch 64 --tol_comp 1e-7 --lrlevel $lrlevel --precon 1 --xyzsort 2 --nmin_leaf 100 --near_para 0.01d0 --pat_comp 3 --format 1 --sample_para $sample_para --knn $knn | tee a.out_freq_${freq}_lrlevel_${lrlevel}
done

