export EXEC=./EXAMPLE/fem2d
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread





# ######## half cyclinder
nmpi=1
# blknum=1
model=7
N=40000
wavelength=0.003
# wavelength=0.01
tol=1d-4
errcheck=0
lrcomp=4
bACAbatch=32

#for bACAbatch in 1 8 16 32
for bACAbatch in 1 
do
for blknum in 1 2 4 8 16 32
# for blknum in 1 
do
# mpirun -n 1 $EXEC $tst $DATA $Ndim $ker $h $lambda $Nmin $tol $com_opt $checkerr $batch | tee -a datset_kernel.out
srun -n $nmpi -c 2 --cpu_bind=cores $EXEC $blknum $model $N $wavelength $tol $errcheck $lrcomp $bACAbatch | tee hcylindar_N_${N}_w_${wavelength}_bnum_${blknum}_comp_${lrcomp}_mpi_${nmpi}_bsize_${bACAbatch} 
done
done


# ######## half cyclinder
# nmpi=1
# # blknum=1
# model=7
# N=40000
# # wavelength=0.003
# # wavelength=0.01
# wavelength=0.06
# tol=1d-4
# errcheck=0
# lrcomp=3
# bACAbatch=32

# for blknum in 1 2 4 8 16 32
# do
# # mpirun -n 1 $EXEC $tst $DATA $Ndim $ker $h $lambda $Nmin $tol $com_opt $checkerr $batch | tee -a datset_kernel.out
# srun -n $nmpi -c 2 --cpu_bind=cores $EXEC $blknum $model $N $wavelength $tol $errcheck $lrcomp $bACAbatch | tee hcylindar_N_${N}_w_${wavelength}_bnum_${blknum}_comp_${lrcomp}_mpi_${nmpi} 
# done



# ######## half cyclinder
# nmpi=1
# # blknum=1
# model=7
# N=40000
# # wavelength=0.003
# # wavelength=0.01
# wavelength=0.06
# tol=1d-4
# errcheck=0
# lrcomp=4
# bACAbatch=32

# for blknum in 1 2 4 8 16 32
# do
# # mpirun -n 1 $EXEC $tst $DATA $Ndim $ker $h $lambda $Nmin $tol $com_opt $checkerr $batch | tee -a datset_kernel.out
# srun -n $nmpi -c 2 --cpu_bind=cores $EXEC $blknum $model $N $wavelength $tol $errcheck $lrcomp $bACAbatch | tee hcylindar_N_${N}_w_${wavelength}_bnum_${blknum}_comp_${lrcomp}_mpi_${nmpi} 
# done


