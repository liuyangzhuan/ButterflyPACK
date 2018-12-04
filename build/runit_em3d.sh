export EXEC=./EXAMPLE/fem3d
export OMP_NUM_THREADS=8
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

THREADS_PER_RANK=`expr 2 \* $OMP_NUM_THREADS`

precon=2 #2: no preconditioner 3: LU preconditioner


for precon in 2 3
do
# ######## half sphere
# nmpi=4
# wavelength=1.0
# filename=halfsphere_2300
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon | tee $filename.out_$precon

# nmpi=4
# wavelength=0.5
# filename=halfsphere_9000
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon | tee $filename.out_$precon

nmpi=16
wavelength=0.25
filename=halfsphere_32000
srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon | tee $filename.out_$precon



# ######## sphere
# nmpi=4
# wavelength=2.0
# filename=sphere_2300
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon | tee $filename.out_$precon

# nmpi=4
# wavelength=1.0
# filename=sphere_9000
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon | tee $filename.out_$precon

nmpi=16
wavelength=0.5
filename=sphere_32000
srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon | tee $filename.out_$precon




# ######## plate
# nmpi=4
# wavelength=1.0
# filename=plate_2000
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon | tee $filename.out_$precon

# nmpi=4
# wavelength=0.5
# filename=plate_8000
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon | tee $filename.out_$precon

nmpi=16
wavelength=0.25
filename=plate_32000
srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon | tee $filename.out_$precon


# ######## corner
# nmpi=4
# wavelength=0.5
# filename=corner_2500
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon | tee $filename.out_$precon

# nmpi=4
# wavelength=0.25
# filename=corner_10000
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon | tee $filename.out_$precon

nmpi=16
wavelength=0.125
filename=corner_40000
srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon | tee $filename.out_$precon

done
