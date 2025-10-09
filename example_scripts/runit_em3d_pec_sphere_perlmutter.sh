#!/bin/bash -l

#SBATCH --account=m4872
#SBATCH -q premium
#SBATCH -N 4
#SBATCH --constraint=cpu
#SBATCH -t 48:00:00
#SBATCH -J ie3d_pec
#SBATCH --mail-user=liuyangzhuan@lbl.gov

module load PrgEnv-gnu
module unload matlab
module unload cray-libsci


ulimit -s unlimited



for nmpi in 128 
do

NTH=4
CORES_PER_NODE=128
THREADS_PER_RANK=`expr $NTH \* 2`								 

NODE_VAL=`expr $nmpi / $CORES_PER_NODE \* $NTH`

export EXEC=./EXAMPLE/ie3d
export OMP_NUM_THREADS=$NTH
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


#### The Default parameters
format=1
elem_extract=2
precon=1
pat_comp=3
blknum=1
tol=1d-4
tol_rand=1d-4
errcheck=0
lrcomp=5
bACAbatch=16
LRlevel=0
xyzsort=1
leafsize=200
para=0.01	  
Nbundle=1		  
verbosity=0
bp_cnt_lr=1
less_adapt=1
use_zfp=2
samplepara=2
sample_para_outer=2    # this is performance critical for nlogn algorithm
knn=0
forwardN15flag=0










########################################## HODBF
format=1
elem_extract=2
precon=3
pat_comp=3
blknum=1
tol=1d-4
tol_rand=1d-2
errcheck=0
lrcomp=5
bACAbatch=16
LRlevel=100
leafsize=200
para=0.01	  
Nbundle=2		  
bp_cnt_lr=1
less_adapt=1

# samplepara=2
# sample_para_outer=2    # this is performance critical for nlogn algorithm
# knn=200
# forwardN15flag=0


# samplepara=2
# sample_para_outer=2   
# knn=0
# forwardN15flag=1


samplepara=2
sample_para_outer=2    # this is performance critical for nlogn algorithm
knn=40
forwardN15flag=0




# ###################################### HODLR
# format=1
# precon=1
# blknum=1
# tol=1d-4
# tol_rand=1d-4
# lrcomp=5
# bACAbatch=16
# LRlevel=0
# leafsize=1000
# para=0.01	  	  
# bp_cnt_lr=1
# knn=10






# ########################################## H-LR
# format=2
# elem_extract=2
# precon=3
# pat_comp=3
# blknum=1
# tol=1d-4
# tol_rand=1d-4
# errcheck=0
# lrcomp=5
# bACAbatch=16
# LRlevel=0
# leafsize=200
# para=2.01	  		  
# bp_cnt_lr=1
# knn=0





# ########################################## H-BF
# format=2
# elem_extract=2
# precon=3
# pat_comp=3
# blknum=1
# tol=1d-4
# tol_rand=1d-2
# errcheck=0
# lrcomp=5
# bACAbatch=16
# LRlevel=100
# leafsize=200
# para=2.01	  
# Nbundle=2		  
# bp_cnt_lr=1
# less_adapt=1
# samplepara=2
# sample_para_outer=2    # this is performance critical for nlogn algorithm
# knn=0
# forwardN15flag=0





export OMP_NUM_THREADS=$OMP_NUM_THREADS





# wavelength=2.0
# filename=sphere_2300

# wavelength=1.0
# filename=sphere_9000

wavelength=0.28
filename=sphere_128000

# wavelength=0.14
# filename=sphere_512000

# wavelength=0.07
# filename=sphere_2M


# wavelength=0.25
# filename=plate_32000

# wavelength=0.125
# filename=plate_128000

# wavelength=0.0625
# filename=plate_512000

# wavelength=0.03125
# filename=plate_2560000

logfile=$filename.out_precon_${precon}_sort_${xyzsort}_mpi_${nmpi}_format_${format}_less_adapt_${less_adapt}_tol_rand_${tol_rand}_lrlevel_${LRlevel}_nbundle_${Nbundle}_forwardN15flag_${forwardN15flag}_omp${NTH}

srun -n $nmpi -N 4 -c $THREADS_PER_RANK --cpu_bind=cores /usr/bin/time -f "Time=%E MaxRSS=%M KB" $EXEC -quant --data_dir ../EXAMPLE/EM3D_DATA/preprocessor_3dmesh/$filename --wavelength $wavelength -option --lr_blk_num $blknum --use_zfp 2 --elem_extract ${elem_extract} --tol_comp $tol --tol_rand ${tol_rand} --less_adapt ${less_adapt} --errfillfull $errcheck --reclr_leaf $lrcomp --baca_batch $bACAbatch --lrlevel $LRlevel --bp_cnt_lr $bp_cnt_lr --precon $precon --xyzsort $xyzsort --forwardN15flag $forwardN15flag --nmin_leaf $leafsize --near_para $para --verbosity ${verbosity} --pat_comp $pat_comp --nbundle $Nbundle --format $format --knn $knn --sample_para $samplepara --sample_para_outer ${sample_para_outer} 2>&1 | tee $logfile


# Extract all numeric values after "MaxRSS="
values=$(grep -o 'MaxRSS=[0-9]*' "$logfile" | cut -d'=' -f2)
if [[ -z "$values" ]]; then
    echo "No MaxRSS entries found in $logfile" | tee -a "$logfile"
    exit 0
fi
# Sum all MaxRSS values (in KB)
total_kb=0
for v in $values; do
    total_kb=$((total_kb + v))
done
total_gb=$(echo "scale=6; $total_kb / 1048576" | bc -l)
printf "Total RSS: %.4f GB\n" "$total_gb" | tee -a "$logfile"




done


