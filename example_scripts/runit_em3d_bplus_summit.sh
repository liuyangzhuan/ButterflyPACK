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

module load netlib-lapack/3.8.0
module load gcc/7.4.0
module load cmake/3.11.3
module load cuda/9.2.148


RANK_PER_RS=1
GPU_PER_RANK=0							 

export EXEC=./EXAMPLE/ie3d



# for CORE_VAL in 16 18 32 50 64 98 128 200 256 512 1024
#for CORE_VAL in  512 
for CORE_VAL in  4
# for CORE_VAL in  32 64
do

NTH=8
RS_VAL=`expr $CORE_VAL / $RANK_PER_RS`
MOD_VAL=`expr $CORE_VAL % $RANK_PER_RS`
if [[ $MOD_VAL -ne 0 ]]
then
  RS_VAL=`expr $RS_VAL + 1`
fi
OMP_NUM_THREADS=$NTH
TH_PER_RS=`expr $NTH \* $RANK_PER_RS`
GPU_PER_RS=`expr $RANK_PER_RS \* $GPU_PER_RANK`


precon=1
pat_comp=3
blknum=1
tol=1d-4
errcheck=0
lrcomp=2
bACAbatch=16
LRlevel=100
xyzsort=1
leafsize=200
para=0.01	  
Nbundle=1		  
format=1
nobp=0
verbosity=2
knn=50

export OMP_NUM_THREADS=$OMP_NUM_THREADS


wavelength=1.0
filename=plate_2000

jsrun -n $RS_VAL -a $RANK_PER_RS -c $TH_PER_RS -g $GPU_PER_RS -b packed:$NTH $EXEC -quant --data_dir ../EXAMPLE/EM3D_DATA/preprocessor_3dmesh/$filename  --wavelength $wavelength -option --lr_blk_num $blknum --tol_comp $tol --errfillfull $errcheck --reclr_leaf $lrcomp --baca_batch $bACAbatch --lrlevel $LRlevel --precon $precon --xyzsort $xyzsort --nmin_leaf $leafsize --near_para $para --pat_comp $pat_comp --nbundle $Nbundle --format $format --lnobp $nobp --verbosity $verbosity --knn $knn | tee ${filename}_w_${wavelength}_tol_${tol}_mpi_${nmpi}_nth_${OMP_NUM_THREADS}_LRlevel_${LRlevel}_precon_${precon}_sort_${xyzsort}_pat_${pat_comp}_lnobp_${lnobp}_knn_${knn}

done


