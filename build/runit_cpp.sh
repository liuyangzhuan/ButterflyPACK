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
lambda=10.
Nmin=200
tol=1e-4
com_opt=3
checkerr=1
batch=100

mpirun -n $nmpi valgrind --leak-check=yes --track-origins=yes --error-limit=no $EXEC $tst $DATA $Ndim $ker $h $lambda $Nmin $tol $com_opt $checkerr $batch | tee -a datset_kernel.out 
#srun -n $nmpi -c 2 --cpu_bind=cores $EXEC $tst $DATA $Ndim $ker $h $lambda $Nmin $tol $com_opt $checkerr $batch | tee -a datset_kernel.out 



# ######## Gaussian cloud
# tst=2
# N=5000
# Ndim=8
# ker=1
# h=0.2
# lambda=10.
# Nmin=200
# tol=1e-4
# com_opt=2
# checkerr=1
# batch=100

# # mpirun -n 1 $EXEC $tst $N $Ndim $ker $h $lambda $Nmin $tol $com_opt $checkerr $batch | tee -a cloud_kernel.out
# srun -n $nmpi -c 2 --cpu_bind=cores $EXEC $tst $N $Ndim $ker $h $lambda $Nmin $tol $com_opt $checkerr $batch | tee -a cloud_kernel.out  
  
# ######## product of Gaussian 
# tst=3
# N=5000
# rank=8
# Nmin=200
# tol=1e-4
# com_opt=2
# checkerr=1
# batch=100


# # mpirun -n 1 $EXEC $tst $N $rank $Nmin $tol $com_opt $checkerr $batch | tee -a randomproduct.out  
# srun -n $nmpi -c 2 --cpu_bind=cores $EXEC $tst $N $rank $Nmin $tol $com_opt $checkerr $batch | tee -a randomproduct.out  
  
