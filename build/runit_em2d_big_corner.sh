#!/bin/bash -l

#SBATCH -p regular
#SBATCH -N 16
#SBATCH -t 20:30:00
#SBATCH -J 50M_corner
#SBATCH --mail-user=liuyangzhuan@lbl.gov
#SBATCH -C haswell

export EXEC=./EXAMPLE/ie2d
export OMP_NUM_THREADS=8
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
THREADS_PER_RANK=`expr 2 \* $OMP_NUM_THREADS`
nmpi=16


for precon in 3  #1: direct 2: no preconditioner 3: LU preconditioner
do

# ######## corruaged corner reflector
Ns=(50000000)
wavelengths=(0.0000012)

for ((i = 0; i < ${#Ns[@]}; i++)); do
N=${Ns[i]}
wavelength=${wavelengths[i]}

blknum=1
model=9
# N=5000
# wavelength=0.08
# wavelength=0.01
tol=1d-4
errcheck=0
lrcomp=4
bACAbatch=16
LRlevel=100
xyzsort=0
leafsize=200


srun -n $nmpi -N $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC $blknum $model $N $wavelength $tol $errcheck $lrcomp $bACAbatch $LRlevel $precon $xyzsort $leafsize | tee corner_N_${N}_w_${wavelength}_tol_${tol}_mpi_${nmpi}_nth_${OMP_NUM_THREADS}_LRlevel_${LRlevel}_precon_${precon}_sort_${xyzsort}

done


done


