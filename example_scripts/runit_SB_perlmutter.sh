#!/bin/bash -l

#SBATCH --account=mp127
#SBATCH -q regular
#SBATCH -N 16
#SBATCH --constraint=cpu
#SBATCH -t 10:00:00
#SBATCH -J tensorbf
#SBATCH --mail-user=liuyangzhuan@lbl.gov


module load PrgEnv-gnu
module unload cray-libsci

CORES_PER_NODE=128
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
  


if [[ $(uname -s) == 'Darwin' ]]; then
    export GPTUNEROOT=/Users/liuyangzhuan/Desktop/GPTune/
    export MPIRUN="$GPTUNEROOT/openmpi-4.1.5/bin/mpirun"
else
    export MPIRUN=mpirun
fi


############## 2D Inverse FIO for ellipse integral with SB
NTH=1
THREADS_PER_RANK=`expr $NTH \* 2`								 
export OMP_NUM_THREADS=$NTH
ker=2
tol_comp=1e-7
tol_rand=1e-4

# ker=1
# tol_comp=1e-10
# tol_rand=1e-7

tol_itersol=1e-12
rank0=100
rankrate=2.0
mu=64.0
lambda=1.0
maxitsb=40
nmin_leaf=32
nmpi=512
precon=3
use_zfp=0
hextralevel=1
# N_FIO=256 512 1024 2048 4096
for N_FIO in 65536 #16384 #
do
srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores ../build/EXAMPLE/cifio2dsb -quant --ker ${ker} --M ${N_FIO} --N ${N_FIO} --mu ${mu} --lambda ${lambda} --maxitsb ${maxitsb} -option --use_zfp ${use_zfp} --format 2 --near_para 2.01 --tol_comp ${tol_comp} --tol_rand ${tol_rand} --errsol 0 --precon ${precon} --tol_itersol ${tol_itersol} --xyzsort 1 --rank0 ${rank0} --hextralevel ${hextralevel} --rankrate ${rankrate} --nmin_leaf ${nmin_leaf} | tee a.out_matrix_IFIO2D_SB_N_FIO${N_FIO}_ker${ker}_tol_comp${tol_comp}_tol_rand${tol_rand}_mpi${nmpi}_omp${NTH}
done
