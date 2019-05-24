#!/bin/bash -l

#SBATCH -q premium
#SBATCH -N 128
#SBATCH -t 00:30:00
#SBATCH -J paralleltest
#SBATCH --mail-user=liuyangzhuan@lbl.gov
#SBATCH -C haswell

module swap PrgEnv-intel/6.0.4 PrgEnv-gnu
NTH=1
CORES_PER_NODE=32
THREADS_PER_RANK=`expr $NTH \* 2`								 

export EXEC=./EXAMPLE/ie2deigen
export OMP_NUM_THREADS=$NTH
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


# for nmpi in 16 18 32 50 64 98 128 200 256 512 1024
#for nmpi in  512 
for nmpi in  128
# for nmpi in  32 64
do

NODE_VAL=`expr $nmpi / $CORES_PER_NODE \* $NTH`
for pat_comp in  3 #1: from right to left 2: from left to right 3: from outter to inner
do

for precon in 1  #1: direct 2: no preconditioner 3: LU preconditioner
do

# ######## half cyclinder
#Ns=(50000 500000 5000000)
Ns=(100000)
wavelengths=(0.004)

for ((i = 0; i < ${#Ns[@]}; i++)); do
N=${Ns[i]}
wavelength=${wavelengths[i]}

blknum=1
model=10
# N=5000
# wavelength=0.08
# wavelength=0.01
tol=1d-4
errcheck=0
lrcomp=4
bACAbatch=16
LRlevel=100
xyzsort=2
leafsize=50
para=5.0
schulzlevel=100		  
Nbundle=1		  
format=1
knn=500
nev=2
which='LM'
si=0
sample_para=2.0d0

srun -n $nmpi -N $NODE_VAL -c $THREADS_PER_RANK --cpu_bind=cores $EXEC -quant --model2d $model --nunk $N --wavelength $wavelength --cmmode 0 --si $si --which $which --nev $nev -option --lr_blk_num $blknum --tol_comp $tol --errfillfull $errcheck --reclr_leaf $lrcomp --baca_batch $bACAbatch --lrlevel $LRlevel --precon $precon --xyzsort $xyzsort --nmin_leaf $leafsize --near_para $para --sample_para $sample_para --pat_comp $pat_comp --schulzlevel $schulzlevel --nbundle $Nbundle --format $format --knn $knn | tee cavity_N_${N}_w_${wavelength}_tol_${tol}_mpi_${nmpi}_nth_${OMP_NUM_THREADS}_LRlevel_${LRlevel}_precon_${precon}_sort_${xyzsort}_pat_${pat_comp}_knn_${knn}

done
done
done

done


