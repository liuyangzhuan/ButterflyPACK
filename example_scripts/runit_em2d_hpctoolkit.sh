#!/bin/bash -l

#SBATCH -q premium
#SBATCH -N 128
#SBATCH -t 00:30:00
#SBATCH -J paralleltest
#SBATCH --mail-user=liuyangzhuan@lbl.gov
#SBATCH -C haswell

# module unload cray-mpich
module swap PrgEnv-intel PrgEnv-gnu
module unload darshan
module load hpctoolkit
module load hpcviewer
# export MKLROOT=/opt/intel/compilers_and_libraries_2019.3.199/linux/mkl
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/compilers_and_libraries_2019.3.199/linux/mkl/lib/intel64
# module load openmpi/4.0.1
# CCC=mpicc
# CCCPP=mpicxx
# FTN=mpif90


NTH=2
CORES_PER_NODE=32
THREADS_PER_RANK=`expr $NTH \* 2`								 

export APP=ie2d
export EXEC=./EXAMPLE/$APP
export OMP_NUM_THREADS=$NTH
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
  							 
											
# for nmpi in 16 18 32 50 64 98 128 200 256 512 1024
#for nmpi in  512 
for nmpi in  64 
# for nmpi in  32 64
do

NODE_VAL=`expr $nmpi / $CORES_PER_NODE \* $NTH`
for pat_comp in  3 #1: from right to left 2: from left to right 3: from outter to inner
do

for precon in 1  #1: direct 2: no preconditioner 3: LU preconditioner
do

# ######## half cyclinder
model=7
Ns=( 100000) 
wavelengths=( 0.001)
# Ns=( 10000000) 
# wavelengths=( 0.00001)
# Ns=(100000000)
# wavelengths=(0.000001)


# # ######## spiral line
# # Ns=(5000 50000 500000 )
# # wavelengths=(0.06 0.006 0.0006 )
# model=12
# #Ns=(200000 )
# #wavelengths=(0.003)
# Ns=(500000 2000000)
# wavelengths=(0.0006 0.00015)
# #Ns=(2000000 )
# #wavelengths=(0.0005 )
# # Ns=(5000000 )
# # wavelengths=(0.00006 )


# # ######## half cyclinder array
# Ns=(500000)
# wavelengths=(0.00106)
# model=13
					 

for ((i = 0; i < ${#Ns[@]}; i++)); do
N=${Ns[i]}
wavelength=${wavelengths[i]}

blknum=1
# model=7
# N=5000
# wavelength=0.08
# wavelength=0.01
tol=1d-4
tol_rand=1d-2
errcheck=0
lrcomp=5
bACAbatch=16
LRlevel=100
xyzsort=0
leafsize=200
para=0.01
schulzlevel=100		  
Nbundle=4 
format=1		  
verbosity=0
less_adapt=1

samplepara=20
knn=40
forwardN15flag=0

# samplepara=2
# knn=0
# forwardN15flag=1

rm -rf hpctoolkit-$APP-measurements*
rm -rf hpctoolkit-$APP-database*


# mpirun -n $nmpi $EXEC -quant --model2d $model --nunk $N --wavelength $wavelength -option --lr_blk_num $blknum --tol_comp $tol --tol_rand ${tol_rand}  --errfillfull $errcheck --reclr_leaf $lrcomp  --verbosity $verbosity  --forwardN15flag $forwardN15flag --sample_para $samplepara --less_adapt ${less_adapt} --baca_batch $bACAbatch --lrlevel $LRlevel --precon $precon --xyzsort $xyzsort --nmin_leaf $leafsize --near_para $para --pat_comp $pat_comp --schulzlevel $schulzlevel --nbundle $Nbundle --format $format --knn $knn | tee hcylindar_N_${N}_w_${wavelength}_tol_${tol}_mpi_${nmpi}_nth_${OMP_NUM_THREADS}_LRlevel_${LRlevel}_precon_${precon}_sort_${xyzsort}_pat_${pat_comp}_forwardN15flag_${forwardN15flag}
srun -n $nmpi -N $NODE_VAL -c $THREADS_PER_RANK --cpu_bind=cores hpcrun -e REALTIME -t $EXEC -quant --model2d $model --nunk $N --wavelength $wavelength -option --lr_blk_num $blknum --tol_comp $tol --tol_rand ${tol_rand}  --errfillfull $errcheck --reclr_leaf $lrcomp --verbosity $verbosity --forwardN15flag $forwardN15flag --sample_para $samplepara --less_adapt ${less_adapt} --baca_batch $bACAbatch --lrlevel $LRlevel --precon $precon --xyzsort $xyzsort --nmin_leaf $leafsize --near_para $para --pat_comp $pat_comp --schulzlevel $schulzlevel --nbundle $Nbundle --format $format --knn $knn | tee hcylindar_N_${N}_w_${wavelength}_tol_${tol}_mpi_${nmpi}_nth_${OMP_NUM_THREADS}_LRlevel_${LRlevel}_precon_${precon}_sort_${xyzsort}_pat_${pat_comp}_forwardN15flag_${forwardN15flag}_format_${format}



hpcstruct -j 16 $EXEC 

for entry in *; do
    entry1=$(echo $entry | grep hpctoolkit-$APP-measurements*)
    if [[ $entry1 != "" ]];then
    declare -a arr=($(echo $entry1 | grep -Eo '?[0-9]+([.][0-9]+)?'))
    num=${arr[${#arr[@]}-1]}
    fi
done


# hpcstruct -j 16 hpctoolkit-$APP-measurements-$num  # this is for GPU
hpcprof -S $APP.hpcstruct hpctoolkit-$APP-measurements-$num
hpcviewer hpctoolkit-$APP-database-$num





done
done
done

done




