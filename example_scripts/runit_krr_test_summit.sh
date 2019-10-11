#!/bin/bash
# Bash script to submit many files to Cori/Edison/Queue

#BSUB -P CSC289
#BSUB -W 2:00
#BSUB -nnodes 45
#BSUB -alloc_flags gpumps
#BSUB -J superlu_gpu



EXIT_SUCCESS=0
EXIT_HOST=1
EXIT_PARAM=2


module swap xl gcc/9.1.0
# module load mercurial
module load cmake
module load essl
# module load cuda
module load netlib-lapack
module load netlib-scalapack


RANK_PER_RS=1
GPU_PER_RANK=0							 

export EXEC=./EXAMPLE/krr



FOLDER="/gpfs/alpine/csc289/scratch/liuyangz/ML/mats/checkCERR/scikit"
dim=8
DATA0=susy_d8_10k_1k/susy
hval=1.3
lval=3.11
ntrain=10000

#FOLDER="/gpfs/alpine/csc289/scratch/liuyangz/ML/mats/checkCERR"
#dim=784
#DATA0=mnist_d784_10K_1K_1K/mnist
#hval=1.3
#lval=3.11
#ntrain=10000


# for CORE_VAL in 16 18 32 50 64 98 128 200 256 512 1024
#for CORE_VAL in  512 
for CORE_VAL in  1 
# for CORE_VAL in  32 64
do

NTH=40
RS_VAL=`expr $CORE_VAL / $RANK_PER_RS`
MOD_VAL=`expr $CORE_VAL % $RANK_PER_RS`
if [[ $MOD_VAL -ne 0 ]]
then
  RS_VAL=`expr $RS_VAL + 1`
fi
OMP_NUM_THREADS=$NTH
TH_PER_RS=`expr $NTH \* $RANK_PER_RS`
GPU_PER_RS=`expr $RANK_PER_RS \* $GPU_PER_RANK`

for pat_comp in  3 #1: from right to left 2: from left to right 3: from outter to inner
do

for precon in 1  #1: direct 2: no preconditioner 3: LU preconditioner
do

blknum=1
model=7
# N=5000
# wavelength=0.08
# wavelength=0.01
tol=1d-2
errcheck=0
lrcomp=5
bACAbatch=16
LRlevel=100
xyzsort=1
leafsize=200
para=0.01
schulzlevel=100		  
Nbundle=1		  
format=1
knn=128

export OMP_NUM_THREADS=$OMP_NUM_THREADS

jsrun -n $RS_VAL -a $RANK_PER_RS -c $TH_PER_RS -g $GPU_PER_RS -b packed:$NTH $EXEC -quant --data_dir $FOLDER/${KRR_DATA}/$DATA0 --dimn $dim --ntrain $ntrain --ntest 1000 --sigma $hval --lambda $lval -option --lr_blk_num $blknum --tol_comp $tol --errfillfull $errcheck --reclr_leaf $lrcomp --baca_batch $bACAbatch --lrlevel $LRlevel --precon $precon --xyzsort $xyzsort --nmin_leaf $leafsize --near_para $para --pat_comp $pat_comp --knn ${knn} --schulzlevel $schulzlevel --nbundle $Nbundle --format $format | tee krr_tol_${tol}_mpi_${CORE_VAL}_nth_${OMP_NUM_THREADS}_LRlevel_${LRlevel}_precon_${precon}_sort_${xyzsort}_pat_${pat_comp}

done
done

done










