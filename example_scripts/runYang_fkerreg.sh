for cluster in 0 2 3
do

# export EXEC=./EXAMPLE/fkerreg
# export DATA=/project/projectdirs/m2957/liuyangz/my_research/ML/SUSY_Origin/susy_10Kn
# export OMP_NUM_THREADS=1
# export OMP_PLACES=threads
# export OMP_PROC_BIND=spread
# h=0.1
# lambda=10.0
# com_opt=2

# srun -N 2 -n 2 -c 2 --cpu_bind=cores ${EXEC} ${DATA} 8 10000 1000 ${h} ${lambda} ${cluster} ${com_opt} | tee hodlr_susy_kernel10K_clus_${cluster}.out


# export EXEC=./EXAMPLE/fkerreg
# export DATA=/project/projectdirs/m2957/liuyangz/my_research/ML/GAS/gas_10Kn
# export OMP_NUM_THREADS=1
# export OMP_PLACES=threads
# export OMP_PROC_BIND=spread
# h=0.1
# lambda=10.0
# com_opt=2

# srun -N 2 -n 2 -c 2 --cpu_bind=cores ${EXEC} ${DATA} 128 10000 1000 ${h} ${lambda} ${cluster} ${com_opt} | tee hodlr_gas_kernel10K_clus_${cluster}.out


# export EXEC=./EXAMPLE/fkerreg
# export DATA=/project/projectdirs/m2957/liuyangz/my_research/ML/MNIST_Origin/mnist_10Kn
# export OMP_NUM_THREADS=1
# export OMP_PLACES=threads
# export OMP_PROC_BIND=spread
# h=3.0
# lambda=4.0
# com_opt=2
# srun -N 2 -n 2 -c 2 --cpu_bind=cores ${EXEC} ${DATA} 784 10000 1000 ${h} ${lambda} ${cluster} ${com_opt} | tee hodlr_mnist_kernel10K_clus_${cluster}.out


# export EXEC=./EXAMPLE/fkerreg
# export DATA=/project/projectdirs/m2957/liuyangz/my_research/ML/COVTYPE/covtype_10Kn
# export OMP_NUM_THREADS=1
# export OMP_PLACES=threads
# export OMP_PROC_BIND=spread
# h=0.1
# lambda=10.0
# com_opt=2
# srun -N 2 -n 2 -c 2 --cpu_bind=cores ${EXEC} ${DATA} 54 10000 1000 ${h} ${lambda} ${cluster} ${com_opt} | tee hodlr_covtype_kernel10K_clus_${cluster}.out


# export EXEC=./EXAMPLE/fkerreg
# export DATA=/project/projectdirs/m2957/liuyangz/my_research/ML/SUSY_Origin/susy_100Kn
# export OMP_NUM_THREADS=1
# export OMP_PLACES=threads
# export OMP_PROC_BIND=spread
# h=0.1
# lambda=10.0
# com_opt=2
# srun -N 32 -n 32 -c 2 --cpu_bind=cores ${EXEC} ${DATA} 8 100000 1000 ${h} ${lambda} ${cluster} ${com_opt} | tee hodlr_susy_kernel100K_clus_${cluster}.out


# export EXEC=./EXAMPLE/fkerreg
# export DATA=/project/projectdirs/m2957/liuyangz/my_research/ML/COVTYPE/covtype_500Kn
# export OMP_NUM_THREADS=1
# export OMP_PLACES=threads
# export OMP_PROC_BIND=spread
# h=0.07
# lambda=0.3
# com_opt=3
# srun -N 32 -n 128 -c 2 --cpu_bind=cores ${EXEC} ${DATA} 54 490000 1000 ${h} ${lambda} ${cluster} ${com_opt} | tee hodlr_covtype_kernel500K_clus_${cluster}.out


export EXEC=./EXAMPLE/fkerreg
export DATA=/project/projectdirs/m2957/liuyangz/my_research/ML/SUSY_Origin/susy_1Mn
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
h=0.1
lambda=10.0
com_opt=3
srun -N 32 -n 128 -c 2 --cpu_bind=cores ${EXEC} ${DATA} 8 1000000 1000 ${h} ${lambda} ${cluster} ${com_opt} | tee hodlr_susy_kernel1M_clus_${cluster}.out

done

