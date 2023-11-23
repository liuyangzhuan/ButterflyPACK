#!/bin/bash -l

#SBATCH -q premium
#SBATCH -N 2
#SBATCH -t 10:00:00
#SBATCH -J paralleltest
#SBATCH --mail-user=liuyangzhuan@lbl.gov
#SBATCH -C haswell

module load PrgEnv-gnu
NTH=1
CORES_PER_NODE=128
THREADS_PER_RANK=`expr $NTH \* 2`								 

export OMP_NUM_THREADS=$NTH
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
  
lrlevel=100
sample_para=2.0
knn=0
format=3
tol=1e-3
xyzsort=1
leaf=400
use_zfp=1

srun -n 64 -c $THREADS_PER_RANK --cpu_bind=cores ./EXAMPLE/full_simple -quant --tst 4 -option --verbosity 1 --forwardN15flag 0 --tol_comp $tol --format $format  --use_zfp $use_zfp --reclr_leaf 5 --baca_batch 16 --lrlevel $lrlevel --precon 1 --nogeo 0 --xyzsort $xyzsort --nmin_leaf $leaf --near_para 2.01d0 --pat_comp 3 --sample_para $sample_para --knn $knn | tee a.out_format${format}_tol${tol}_xyzsort${xyzsort}_leaf${leaf}_zfp${use_zfp}
# srun -n 4 -c $THREADS_PER_RANK --cpu_bind=cores ./EXAMPLE/full_simple -quant --tst 3 -option --verbosity 2 --tol_comp 1e-3 --format 3 --reclr_leaf 5 --baca_batch 16 --lrlevel $lrlevel --precon 1 --nogeo 1 --xyzsort 0 --nmin_leaf 100 --near_para 0.01d0 --pat_comp 3 --sample_para $sample_para --knn $knn | tee a.out


