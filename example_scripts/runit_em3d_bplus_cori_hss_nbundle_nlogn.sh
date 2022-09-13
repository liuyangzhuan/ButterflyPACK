#!/bin/bash -l

#SBATCH -q test
#SBATCH -N 32
#SBATCH -t 04:00:00
#SBATCH -J 3Dtest
#SBATCH --mail-user=liuyangzhuan@lbl.gov
#SBATCH -C haswell

module swap PrgEnv-intel PrgEnv-gnu



for nmpi in 128 
do

NTH=4
CORES_PER_NODE=32
THREADS_PER_RANK=`expr $NTH \* 2`								 

NODE_VAL=`expr $nmpi / $CORES_PER_NODE \* $NTH`
# NODE_VAL=32

export EXEC=./EXAMPLE/ie3d
export OMP_NUM_THREADS=$NTH
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


precon=1
pat_comp=3
blknum=1
tol=1d-4
tol_rand=1d-2
errcheck=0
lrcomp=5
bACAbatch=16
LRlevel=100
xyzsort=1
leafsize=200
para=2.01	  
Nbundle=2		  
# format=1
verbosity=0
bp_cnt_lr=1
less_adapt=1


# samplepara=2
# sample_para_outer=2    # this is performance critical for nlogn algorithm
# knn=200
# forwardN15flag=0

samplepara=2
sample_para_outer=2   
knn=0
forwardN15flag=1

# samplepara=2
# sample_para_outer=2    # this is performance critical for nlogn algorithm
# knn=100
# forwardN15flag=0

for nobp in  100 
do

for format in 3
do

export OMP_NUM_THREADS=$OMP_NUM_THREADS





# wavelength=2.0
# filename=sphere_2300


# wavelength=0.25
# filename=plate_32000

wavelength=0.125
filename=plate_128000

# wavelength=0.0625
# filename=plate_512000

# wavelength=0.03125
# filename=plate_2560000

srun -n $nmpi -N $NODE_VAL -c $THREADS_PER_RANK --cpu_bind=cores $EXEC -quant --data_dir /project/projectdirs/m2957/liuyangz/my_research/ButterflyPACK_hss_factor_acc/EXAMPLE/EM3D_DATA/preprocessor_3dmesh/$filename --wavelength $wavelength -option --lr_blk_num $blknum --tol_comp $tol --tol_rand ${tol_rand} --less_adapt ${less_adapt} --errfillfull $errcheck --reclr_leaf $lrcomp --baca_batch $bACAbatch --lrlevel $LRlevel --bp_cnt_lr $bp_cnt_lr --precon $precon --xyzsort $xyzsort --forwardN15flag $forwardN15flag --nmin_leaf $leafsize --near_para $para --verbosity ${verbosity} --pat_comp $pat_comp --nbundle $Nbundle --format $format --knn $knn --sample_para $samplepara --sample_para_outer ${sample_para_outer} --lnobp $nobp | tee $filename.out_precon_${precon}_sort_${xyzsort}_lnobp_${nobp}_mpi_${nmpi}_format_${format}_less_adapt_${less_adapt}_tol_rand_${tol_rand}_nbundle_${Nbundle}_entrynlogn



done
done
done


