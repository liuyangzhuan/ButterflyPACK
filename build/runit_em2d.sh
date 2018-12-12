export EXEC=./EXAMPLE/fem2d
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

THREADS_PER_RANK=`expr 2 \* $OMP_NUM_THREADS`



for precon in 1  #1: direct 2: no preconditioner 3: LU preconditioner
do


# ######## cavity
Ns=(5000  )
wavelengths=(0.05 )

for ((i = 0; i < ${#Ns[@]}; i++)); do
N=${Ns[i]}
wavelength=${wavelengths[i]}

nmpi=1
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


srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC $blknum $model $N $wavelength $tol $errcheck $lrcomp $bACAbatch $LRlevel $precon | tee cavity_N_${N}_w_${wavelength}_tol_${tol}_mpi_${nmpi}_LRlevel_${LRlevel}_precon_$precon

done



# # ######## half cyclinder
# Ns=(5000 50000 500000)
# wavelengths=(0.02 0.002 0.0002)

# for ((i = 0; i < ${#Ns[@]}; i++)); do
# N=${Ns[i]}
# wavelength=${wavelengths[i]}

# nmpi=16
# blknum=1
# model=7
# # N=5000
# # wavelength=0.08
# # wavelength=0.01
# tol=1d-4
# errcheck=0
# lrcomp=4
# bACAbatch=16
# LRlevel=100


# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC $blknum $model $N $wavelength $tol $errcheck $lrcomp $bACAbatch $LRlevel $precon | tee hcylindar_N_${N}_w_${wavelength}_tol_${tol}_mpi_${nmpi}_LRlevel_${LRlevel}_precon_$precon

# done


done


