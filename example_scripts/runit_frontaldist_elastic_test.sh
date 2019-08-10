 
#!/bin/bash -l

#SBATCH -q premium
#SBATCH -N 2
#SBATCH -t 10:00:00
#SBATCH -J paralleltest
#SBATCH --mail-user=liuyangzhuan@lbl.gov
#SBATCH -C haswell

# module swap PrgEnv-intel/6.0.4 PrgEnv-gnu
NTH=8
CORES_PER_NODE=32
THREADS_PER_RANK=`expr $NTH \* 2`								 

export EXEC=./EXAMPLE/frontaldist
export OMP_NUM_THREADS=$NTH
export OMP_PLACES=threads 
export OMP_PROC_BIND=spread


# tol=1d-4
errcheck=0
# lrcomp=5
# bACAbatch=16
LRlevel=0
leafsize=500
pat_comp=3
schulzlevel=3000
Nbundle=1
format=1

for tol in 1d-6 
do

for com_opt in 2
do
for batch in 1 
do
for blknum in 1 2 4 8 16 32
do

mpirun -n 4 ./EXAMPLE/frontaldist -quant --nunk 2500 --data_dir ~/Desktop/research/matrix/Frontal_elastic_2_5K/Frontal_elastic_2500 --explicitflag 1 -option --nmin_leaf $leafsize --lr_blk_num $blknum --tol_comp $tol --lrlevel $LRlevel --rank0 1000 --reclr_leaf $com_opt --baca_batch $batch --errfillfull $errcheck | tee elastic_2_5k.out_comp_${com_opt}_tol_${tol}_bsize_${batch}_blknum_${blknum}_errcheck_${errcheck}


# mpirun -n 16 ./EXAMPLE/frontaldist -quant --nunk 10000 --data_dir ~/Desktop/research/matrix/outHelmholtzFront100/Froot --explicitflag 1 -option --nmin_leaf $leafsize --lr_blk_num $blknum --tol_comp $tol --lrlevel $LRlevel --rank0 1000 --reclr_leaf $com_opt --baca_batch $batch --errfillfull $errcheck | tee helmholtz_10k.out_comp_${com_opt}_tol_${tol}_bsize_${batch}_blknum_${blknum}

done 
done 
done 
done 



