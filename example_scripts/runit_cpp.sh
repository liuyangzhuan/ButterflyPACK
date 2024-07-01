#!/bin/bash -l

#SBATCH --account=mp127
#SBATCH -q regular
#SBATCH -N 16
#SBATCH --constraint=cpu
#SBATCH -t 10:00:00
#SBATCH -J tensorbf
#SBATCH --mail-user=liuyangzhuan@lbl.gov
ulimit -s unlimited

module load PrgEnv-gnu
module load cmake
module unload cray-libsci
module load python

export EXEC=./EXAMPLE/ctest
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

export LD_LIBRARY_PATH=/global/cfs/cdirs/m2957/liuyangz/my_software/scalapack-2.2.0/build/install/lib/:$LD_LIBRARY_PATH

THREADS_PER_RANK=`expr 2 \* $OMP_NUM_THREADS`


for nmpi in 32
do
# ######## Liza's datset MNIST 10K
# tst=1
# DATA=../EXAMPLE/KRR_DATA/susy_10Kn
# Ndim=8
# N=10000
# # DATA=/project/projectdirs/m2957/liuyangz/my_research/ML/SUSY_Origin/susy_10Kn
# # Ndim=784
# # N=10000
# ker=1
# # h=1.0 
# lambda=1.
# Nmin=200
# # tol=1e-2
# com_opt=5
# # checkerr=1
# # batch=32
# # blknum=1
# knn=0
  
# # for knn in 0 10
# for knn in 0
# do
# for batch in 16
# do
# # for tol in 1e-2 1e-6 1e-10
# for tol in 1e-2 
# do
# for h in 3.0
# do
# for blknum in 1 
# do
# for checkerr in 0
# do
# # mpirun -n 1 $EXEC $tst $DATA $Ndim $ker $h $lambda $Nmin $tol $com_opt $checkerr $batch | tee -a datset_kernel.out
# srun -n $nmpi -c 2 --cpu_bind=cores $EXEC --tst $tst --trainfile $DATA --N $N --Ndim $Ndim --ker $ker --h $h --lambda $lambda --nmin_leaf $Nmin --tol_comp $tol --reclr_leaf $com_opt --errfillfull $checkerr --baca_batch $batch --LR_BLK_NUM $blknum --knn $knn --verbosity 1  | tee mnist_datset_kernel_ker_${ker}_h_${h}_l_${lambda}_cherr_${checkerr}_bnum_${blknum}_comp_${com_opt}_tol_${tol}_mpi_${nmpi}_bsize_${batch}
# done
# done 
# done 
# done 
# done 
# done 


# tst=1
# # DATA=../EXAMPLE/KRR_DATA/susy_10Kn
# DATA=/project/projectdirs/m2957/liuyangz/my_research/ML/SUSY_Origin/susy_100Kn
# Ndim=8
# ker=1
# # h=0.1
# lambda=1.
# Nmin=200
# # tol=1e-2
# # com_opt=2
# for com_opt in 5 
# do
# # checkerr=0
# # batch=100
# for batch in 16 32 64
# do
# #for tol in 1e-2 1e-6 1e-10
# for tol in 1e-2
# do
# for h in 1.0
# do
# for blknum in 4 8 16 32 64
# do
# for checkerr in 0 
# do
# srun -n $nmpi -c 2 --cpu_bind=cores $EXEC $tst $DATA $Ndim $ker $h $lambda $Nmin $tol $com_opt $checkerr $batch $blknum | tee susy_datset_kernel_ker_${ker}_h_${h}_l_${lambda}_cherr_${checkerr}_bnum_${blknum}_comp_${com_opt}_tol_${tol}_mpi_${nmpi}_bsize_${batch}
# done
# done 
# done 
# done 
# done 										  
# done 										  


######## Gaussian cloud
tst=2
N=16384
Ndim=3
ker=9  # Wajih's paper: ker=8: Laplacian ker=9: Laplacian with low-rank update ker=10: Wave equation kernel 
knn=0
# h=0.2
lambda=10.
Nmin=256
#tol=1e-4
# com_opt=4
checkerr=0
# batch=32
# blknum=1

for com_opt in 6
do
# checkerr=0
# batch=100
for batch in 1 
do
#for tol in 1e-2 1e-6 1e-10
for tol in 1e-6
do
for h in 0.2
do
for blknum in 1 
do
for checkerr in 0 
do										  
# mpirun -n 1 $EXEC $tst $N $Ndim $ker $h $lambda $Nmin $tol $com_opt $checkerr $batch | tee -a cloud_kernel.out
srun -n $nmpi -c 2 --cpu_bind=cores $EXEC --tst $tst --N $N --Ndim $Ndim --ker $ker --h $h --lambda $lambda --nmin_leaf $Nmin --tol_comp $tol --format 2 --near_para 2.01 --rankrate 2.0 --reclr_leaf $com_opt --errfillfull $checkerr --baca_batch $batch --LR_BLK_NUM $blknum --knn $knn --verbosity 1 | tee cloud_kernel_N_${N}_ker_${ker}_h_${h}_l_${lambda}_cherr_${checkerr}_bnum_${blknum}_comp_${com_opt}_tol_${tol}_mpi_${nmpi}_bsize_${batch}
done
done 
done 
done 
done 										  
done 



# ######## product of Gaussian
# tst=3
# N=20000
# rank=800
# Nmin=200
# tol=1e-4
# # com_opt=2
# checkerr=0
# for com_opt in 5 
# do
# for batch in 8  
# do
# for tol in 1e-4
# do
# for blknum in  1
# do
# for checkerr in 0 
# do
# srun -n $nmpi -c 2 --cpu_bind=cores $EXEC $tst $N $rank $Nmin $tol $com_opt $checkerr $batch $blknum | tee randomproduct_cherr_${checkerr}_bnum_${blknum}_comp_${com_opt}_tol_${tol}_mpi_${nmpi}_bsize_${batch}_N_${N}_rank_${rank}
# done
# done 
# done 
# done 
# done 										  

done
