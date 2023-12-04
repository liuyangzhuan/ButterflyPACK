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

export EXEC=./EXAMPLE/full
export OMP_NUM_THREADS=$NTH
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


# for nmpi in 16 18 32 50 64 98 128 200 256 512 1024
#for nmpi in  512 
for nmpi in  1
# for nmpi in  32 64
do

NODE_VAL=`expr $nmpi / $CORES_PER_NODE \* $NTH`
for pat_comp in  3 #1: from right to left 2: from left to right 3: from outter to inner
do

for precon in 1  #1: direct 2: no preconditioner 3: LU preconditioner
do


# Ns=(503 1005 1005 1225 5041 5041 5041 20164 20164)
# wavelengths=(32 32 32 4 4 4 8 8 8)
# ms=(2 2 2 3 3 3 3 3 3)
# exs=(1 2 2 3 4 4 3 4 4)
# orders=(2 2 4 2 2 4 2 2 4)

Ns=(503 1005 1005 )
wavelengths=(32 32 32 )
ms=(2 2 2 )
exs=(1 2 2 )
orders=(2 2 4 )

for ((i = 0; i < ${#Ns[@]}; i++)); do
N=${Ns[i]}
wavelength=${wavelengths[i]}
m=${ms[i]}
ex=${exs[i]}
order=${orders[i]}

dir=/global/cscratch1/sd/liuyangz/my_research/Babich_Ansatz
tol=1d-4
errcheck=0
LRlevel=100
xyzsort=2
leafsize=50
para=2.0	  
format=2
knn=10
sample_para=2.0d0

srun -n $nmpi -N $NODE_VAL -c $THREADS_PER_RANK --cpu_bind=cores $EXEC -quant --tst 2 --nunk ${N} --ndim ${m} --data_dir ${dir}/N_${N}_matGO_ex_${ex}_order_${order}_w_${wavelength}pi.csv --geo_dir ${dir}/N_${N}_locationGO_ex_${ex}_order_${order}_w_${wavelength}pi.csv -option --tol_comp $tol --errfillfull $errcheck --lrlevel $LRlevel --xyzsort $xyzsort --nmin_leaf $leafsize --near_para $para --sample_para $sample_para --format $format --knn $knn | tee GO_N_${N}_w_${wavelength}_tol_${tol}_mpi_${nmpi}_nth_${OMP_NUM_THREADS}_LRlevel_${LRlevel}_sort_${xyzsort}_knn_${knn}_ex_${ex}_order_${order}

done
done
done

done


