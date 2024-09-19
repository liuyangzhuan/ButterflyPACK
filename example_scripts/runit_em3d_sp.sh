# module swap PrgEnv-intel/6.0.4 PrgEnv-gnu
NTH=1
CORES_PER_NODE=32
THREADS_PER_RANK=`expr $NTH \* 2`								 

export EXENAME=ie3d_sp
export EXEC=./EXAMPLE/$EXENAME
export OMP_NUM_THREADS=$NTH
export OMP_PLACES=threads
export OMP_PROC_BIND=spread



precon=1 #1: direct 2: no preconditioner 3: LU preconditioner
# blknum=1
tol=1d-4
errcheck=0
# lrcomp=5
# bACAbatch=16
LRlevel=100
xyzsort=2
leafsize=200
para=2.01d0
sample_para=1.2d0
sample_para_outer=2.0d0
pat_comp=3
schulzlevel=3000
Nbundle=1
format=5
knn=0
less_adapt=1
verbosity=1


for com_opt in 5
do
for batch in 64 
do
for blknum in 1 
do
# sort=0 # 0:natural order 1: CKD 2: TM
for precon in 1
do
######## half sphere
# nmpi=1
# wavelength=1.4
# filename=halfsphere_1200
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon $sort | tee $filename.out_precon_$precon_sort_$sort

# nmpi=32
# wavelength=2.0
# filename=sphere_2300
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC -quant --data_dir ../EXAMPLE/EM3D_DATA/preprocessor_3dmesh/$filename --wavelength $wavelength -option --lr_blk_num $blknum --tol_comp $tol --errfillfull $errcheck --reclr_leaf $com_opt --baca_batch $batch --lrlevel $LRlevel --precon $precon --xyzsort $xyzsort --nmin_leaf $leafsize --near_para $para --sample_para $sample_para --verbosity ${verbosity} --pat_comp $pat_comp --schulzlevel $schulzlevel --nbundle $Nbundle --format $format --less_adapt ${less_adapt} --sample_para_outer $sample_para_outer --knn $knn | tee $filename.out_precon_${precon}_sort_${xyzsort}_comp_${com_opt}_tol_${tol}_bsize_${batch}_history_$EXENAME


# nmpi=32
# wavelength=4.0
# filename=sphere_9000
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC -quant --data_dir ../EXAMPLE/EM3D_DATA/preprocessor_3dmesh/$filename --wavelength $wavelength -option --lr_blk_num $blknum --tol_comp $tol --errfillfull $errcheck --reclr_leaf $com_opt --baca_batch $batch --lrlevel $LRlevel --precon $precon --xyzsort $xyzsort --nmin_leaf $leafsize --near_para $para --sample_para $sample_para --sample_para_outer $sample_para_outer --pat_comp $pat_comp --schulzlevel $schulzlevel --nbundle $Nbundle --format $format --knn $knn | tee $filename.out_precon_$precon_sort_$sort






# nmpi=4
# wavelength=0.5
# filename=halfsphere_9000
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon $sort | tee $filename.out_precon_$precon_sort_$sort

# nmpi=16
# wavelength=0.25
# filename=halfsphere_32000
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon $sort | tee $filename.out_precon_$precon_sort_$sort



# # ######## sphere
nmpi=16

# # # wavelength=0.5
# # # filename=plate_8000
# # wavelength=0.25
# # filename=halfsphere_32000
wavelength=2.0
filename=sphere_2300
# # wavelength=0.25
# # filename=plate_32000
unset OMP_PLACES
unset OMP_PROC_BIND
mpirun --allow-run-as-root -n $nmpi $EXEC -quant --data_dir ../EXAMPLE/EM3D_DATA/$filename --wavelength $wavelength -option --lr_blk_num $blknum --tol_comp $tol --errfillfull $errcheck --reclr_leaf $com_opt --baca_batch $batch --lrlevel $LRlevel --precon $precon --xyzsort $xyzsort --nmin_leaf $leafsize --near_para $para --sample_para $sample_para --verbosity ${verbosity} --pat_comp $pat_comp --schulzlevel $schulzlevel --nbundle $Nbundle --format $format --less_adapt ${less_adapt} --sample_para_outer $sample_para_outer --knn $knn | tee $filename.out_precon_${precon}_sort_${xyzsort}_comp_${com_opt}_tol_${tol}_bsize_${batch}_history_$EXENAME


# nmpi=4
# wavelength=1.0
# filename=sphere_9000
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon $sort | tee $filename.out_precon_$precon_sort_$sort

# nmpi=16
# wavelength=0.5
# filename=sphere_32000
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon $sort | tee $filename.out_precon_$precon_sort_$sort




# ######## plate
# nmpi=4
# wavelength=1.0
# filename=plate_2000
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon $sort | tee $filename.out_precon_$precon_sort_$sort

# nmpi=4
# wavelength=0.5
# filename=plate_8000
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon $sort | tee $filename.out_precon_$precon_sort_$sort

# nmpi=16
# wavelength=0.25
# filename=plate_32000
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon $sort | tee $filename.out_precon_$precon_sort_$sort


# ######## corner
# nmpi=4
# wavelength=0.5
# filename=corner_2500
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon $sort | tee $filename.out_precon_$precon_sort_$sort

# nmpi=4
# wavelength=0.25
# filename=corner_10000
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon $sort | tee $filename.out_precon_$precon_sort_$sort

# nmpi=16
# wavelength=0.125
# filename=corner_40000
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC ../EXAMPLE/preprocessor_3dmesh/$filename $wavelength $precon $sort | tee $filename.out_precon_$precon_sort_$sort

done
done
done
done
