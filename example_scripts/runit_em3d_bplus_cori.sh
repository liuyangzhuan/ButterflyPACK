module swap PrgEnv-intel/6.0.4 PrgEnv-gnu



for nmpi in  64
do

NTH=8
CORES_PER_NODE=32
THREADS_PER_RANK=`expr $NTH \* 2`								 

export EXEC=./EXAMPLE/ie3d
export OMP_NUM_THREADS=$NTH
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


precon=1
pat_comp=3
blknum=1
tol=1d-4
errcheck=0
lrcomp=2
bACAbatch=16
LRlevel=100
xyzsort=1
leafsize=200
para=0.01	  
Nbundle=1		  
format=1
verbosity=2
knn=50

for nobp in  0 100
do

export OMP_NUM_THREADS=$OMP_NUM_THREADS


wavelength=0.0625
filename=plate_512000


srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores $EXEC -quant --data_dir /global/cscratch1/sd/liuyangz/my_research/preprocessor_3dmesh/$filename --wavelength $wavelength -option --lr_blk_num $blknum --tol_comp $tol --errfillfull $errcheck --reclr_leaf $lrcomp --baca_batch $bACAbatch --lrlevel $LRlevel --precon $precon --xyzsort $xyzsort --nmin_leaf $leafsize --near_para $para --pat_comp $pat_comp --nbundle $Nbundle --format $format --knn $knn --lnobp $nobp | tee $filename.out_precon_${precon}_sort_${sort}_lnobp_${nobp}



done
done


