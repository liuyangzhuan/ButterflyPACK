#!/bin/bash -l

#SBATCH -q premium
#SBATCH -N 128
#SBATCH -t 00:30:00
#SBATCH -J paralleltest
#SBATCH --mail-user=liuyangzhuan@lbl.gov
#SBATCH -C haswell

module swap PrgEnv-intel/6.0.4 PrgEnv-gnu

export EXEC=./EXAMPLE/ie2d
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
THREADS_PER_RANK=`expr 2 \* $OMP_NUM_THREADS`


# for nmpi in 16 18 32 50 64 98 128 200 256 512 1024
for nmpi in 1024
do

for pat_comp in 1 2 3 #1: from right to left 2: from left to right 3: from outter to inner
do

for precon in 1  #1: direct 2: no preconditioner 3: LU preconditioner
do

# ######## half cyclinder
#Ns=(5000 50000 500000 5000000 50000000)
#wavelengths=(0.02 0.002 0.0002 0.0002 0.00002)

Ns=(5000000)
wavelengths=(0.00005)

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
xyzsort=0
leafsize=50
para=0.01


srun -n $nmpi -N 32 -c $THREADS_PER_RANK --cpu_bind=cores $EXEC $blknum $model $N $wavelength $tol $errcheck $lrcomp $bACAbatch $LRlevel $precon $xyzsort $leafsize $para $pat_comp | tee cavity_N_${N}_w_${wavelength}_tol_${tol}_mpi_${nmpi}_nth_${OMP_NUM_THREADS}_LRlevel_${LRlevel}_precon_${precon}_sort_${xyzsort}_pat_${pat_comp}_matvec

done
done
done

done


