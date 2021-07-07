#!/bin/bash -l

#SBATCH -q premium
#SBATCH -N 128
#SBATCH -t 00:30:00
#SBATCH -J paralleltest
#SBATCH --mail-user=liuyangzhuan@lbl.gov
#SBATCH -C haswell

module swap PrgEnv-intel PrgEnv-gnu
NTH=1
CORES_PER_NODE=32
THREADS_PER_RANK=`expr $NTH \* 2`								 

export EXEC=./EXAMPLE/ie2d
export OMP_NUM_THREADS=$NTH
export OMP_PLACES=threads
export OMP_PROC_BIND=spread



# for nmpi in 16 18 32 50 64 98 128 200 256 512 1024
for nmpi in  1024
# for nmpi in  1
# for nmpi in  32 64
do

NODE_VAL=`expr $nmpi / $CORES_PER_NODE \* $NTH`
#NODE_VAL=32
for pat_comp in  3 #1: from right to left 2: from left to right 3: from outter to inner
do

for precon in 1  #1: direct 2: no preconditioner 3: LU preconditioner
do

for format in  1  #1: hodbf/lr 2: h 3: hss-bf
do


# ######## half cyclinder
# Ns=(100000000)
# wavelengths=(0.000001)
# Ns=(1000000 ) 
# wavelengths=(0.0001 )
Ns=(10000 ) 
wavelengths=(0.01 )
model=7
geo=hcylindar
		

# # ######## spiral line
# # Ns=(5000 50000 500000 )
# # wavelengths=(0.06 0.006 0.0006 )
# model=12
# geo=spiral
# #Ns=(200000 )
# #wavelengths=(0.003)
# Ns=(500000)
# wavelengths=(0.0006)
# #Ns=(2000000 )
# #wavelengths=(0.0005 )
# # Ns=(5000000 )
# # wavelengths=(0.00006 )


# # ######## cavity
# model=10
# geo=cavity
# # Ns=(5000 50000 500000)
# # wavelengths=(0.05 0.005 0.0005)
# Ns=(500000)
# wavelengths=(0.0005)


for ((i = 0; i < ${#Ns[@]}; i++)); do
N=${Ns[i]}
wavelength=${wavelengths[i]}

less_adapt=1
blknum=1
# N=5000
# wavelength=0.08
# wavelength=0.01
tol=1d-4
tolrand=1d-4
errcheck=0
lrcomp=4
bACAbatch=64
LRlevel=100
xyzsort=2
leafsize=100
para=0.01
Nbundle=1 
# knn=20
# forwardN15flag=0
knn=0
forwardN15flag=1
format=3		  
bp_cnt_lr=1
verbosity=2
sample_para=4.0
sample_para_outer=4.0
itermax=20

# # no schulz iteration
# schulzlevel=100	
# schulzorder=7
# schulzsplitlevel=1
# schulzhardstart=0

# schulz iteration with soft or hard start
schulzlevel=0	
schulzorder=5
schulzsplitlevel=3
schulzhardstart=0



srun -n $nmpi -N $NODE_VAL -c $THREADS_PER_RANK --cpu_bind=cores $EXEC -quant --model2d $model --nunk $N --wavelength $wavelength -option --lr_blk_num $blknum --tol_comp $tol --tol_rand $tolrand --errfillfull $errcheck --less_adapt $less_adapt --reclr_leaf $lrcomp --baca_batch $bACAbatch --lrlevel $LRlevel --precon $precon --xyzsort $xyzsort --nmin_leaf $leafsize --verbosity $verbosity --near_para $para --forwardN15flag $forwardN15flag --sample_para ${sample_para} --sample_para_outer ${sample_para_outer} --pat_comp $pat_comp --schulzlevel $schulzlevel --schulzorder $schulzorder --itermax $itermax --schulzsplitlevel $schulzsplitlevel --schulzhardstart $schulzhardstart --nbundle $Nbundle --format $format --bp_cnt_lr $bp_cnt_lr --knn $knn | tee ${geo}_N_${N}_w_${wavelength}_tol_${tol}_mpi_${nmpi}_nth_${OMP_NUM_THREADS}_LRlevel_${LRlevel}_precon_${precon}_sort_${xyzsort}_pat_${pat_comp}_format_${format}_schulzlevel_${schulzlevel}_schulzorder_${schulzorder}_schulzsplitlevel_${schulzsplitlevel}_schulzhardstart_${schulzhardstart}

done
done
done
done

done



