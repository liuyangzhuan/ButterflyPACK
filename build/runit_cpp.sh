export EXEC=./EXAMPLE/ctest
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

nmpi=2

######## Liza's datset
tst=1
DATA=../EXAMPLE/SUSY/susy_10Kn
# DATA=/project/projectdirs/m2957/liuyangz/my_research/ML/SUSY_Origin/susy_10Kn
Ndim=8
ker=1
h=0.2
lambda=1.
Nmin=200
# tol=1e-2
com_opt=4
checkerr=1
batch=256
# blknum=1

for tol in 1e-6 
do
for blknum in 1 
do
# mpirun -n 1 $EXEC $tst $DATA $Ndim $ker $h $lambda $Nmin $tol $com_opt $checkerr $batch | tee -a datset_kernel.out
srun -n $nmpi -c 2 --cpu_bind=cores $EXEC $tst $DATA $Ndim $ker $h $lambda $Nmin $tol $com_opt $checkerr $batch $blknum | tee datset_kernel_ker_${ker}_h_${h}_l_${lambda}_cherr_${checkerr}_bnum_${blknum}_comp_${com_opt}_tol_${tol}_mpi_${nmpi}_bsize_${batch}_history
done
done 


# ######## Gaussian cloud
# tst=2
# N=5000
# Ndim=50
# ker=5
# h=0.2
# lambda=10.
# Nmin=200
# tol=1e-4
# com_opt=3
# checkerr=1
# batch=32
# # blknum=1

# for tol in 1e-2
# do
# for blknum in 1 2 4 8 16 32 
# do
# # mpirun -n 1 $EXEC $tst $N $Ndim $ker $h $lambda $Nmin $tol $com_opt $checkerr $batch | tee -a cloud_kernel.out
# srun -n $nmpi -c 2 --cpu_bind=cores $EXEC $tst $N $Ndim $ker $h $lambda $Nmin $tol $com_opt $checkerr $batch $blknum | tee cloud_kernel_N_${N}_ker_${ker}_h_${h}_l_${lambda}_cherr_${checkerr}_bnum_${blknum}_comp_${com_opt}_tol_${tol}_mpi_${nmpi}
# done 
# done
  
# ######## product of Gaussian 
# tst=3
# N=5000
# rank=100
# Nmin=200
# tol=1e-4
# com_opt=4
# checkerr=0
# batch=32
# # blknum=1

# for blknum in 1 2 4 8 16 32
# do
# # mpirun -n 1 $EXEC $tst $N $rank $Nmin $tol $com_opt $checkerr $batch | tee randomproduct.out  
# srun -n $nmpi -c 2 --cpu_bind=cores $EXEC $tst $N $rank $Nmin $tol ${com_opt} $checkerr $batch $blknum | tee -a randomproduct_N_${N}_r_${rank}_cherr_${checkerr}_bnum_${blknum}_comp_${com_opt}_mpi_${nmpi} 
# done 
