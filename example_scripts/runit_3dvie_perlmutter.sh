#!/bin/bash -l

#SBATCH --account=mp127
#SBATCH -q regular
#SBATCH -N 16
#SBATCH --constraint=cpu
#SBATCH -t 10:00:00
#SBATCH -J tensorbf
#SBATCH --mail-user=liuyangzhuan@lbl.gov
ulimit -s unlimited

module load PrgEnv-gnu
module load cmake
module unload cray-libsci
module load python

export EXEC=./EXAMPLE/cvie2d
export OMP_NUM_THREADS=8
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

export LD_LIBRARY_PATH=/global/cfs/cdirs/m2957/liuyangz/my_software/scalapack-2.2.0/build/install/lib/:$LD_LIBRARY_PATH

THREADS_PER_RANK=`expr 2 \* $OMP_NUM_THREADS`
nmpi=64


vs=1
shape=4
scaleGreen=0
L=0.4
H=0.4
W=0.4
tol_rand=1e-2
tol_Rdetect=0.3e-2
tol_itersol=1e-12

# tol=1d-4
LRlevel=100
h0=0.005
nmin_leaf_t=4
nmin_leaf_m=16
xyzsort=1
BACA_Batch=100
near_para=2.1d0
precon=2
n_iter=3000
verbosity=1
knn_near_para=20.0
knn=0
forwardN15flag=0
rmax=1000
fastsample_tensor=2



# format=2
# LRlevel=0
# knn_near_para=20.0
# knn=0
# forwardN15flag=0
# elem_extract=2
# rmax=1000
# # sample_para=4.0   #20.0
# sample_para=2.0
# sample_para_outer=2.0

# format=1
# elem_extract=2
# knn_near_para=20.0
# knn=600
# forwardN15flag=0
# rmax=1000
# sample_para=20.0
# sample_para_outer=4.0

# # sample_para=4.0
# # sample_para_outer=2.0


# format=1
# elem_extract=0
# knn_near_para=20.0
# knn=0
# forwardN15flag=1
# rmax=1000
# sample_para=1.0
# sample_para_outer=1.0




# omega=25.132741228718345
# ivelo=9
# TNx=8
# TNy=8
# TNz=8
# tol_comp=1e-5
# h=0.02 # 0.015
# x0max=1.0
# y0max=1.0
# z0max=1.0
# xmax=1.1
# ymax=1.1
# zmax=1.1


omega=6.283185307179586
ivelo=11
TNx=8
TNy=8
TNz=8
tol_comp=1e-9
scaleGreen=1
h=0.02 # 0.015
x0max=1.0
y0max=1.0
z0max=1.0
xmax=1.1
ymax=1.1
zmax=1.1
centerx=0.5
centery=0.5
centerz=0.5
L=0.8
H=0.8
W=0.8
xoff=$(python -c "print($centerx - $L/2)")
yoff=$(python -c "print($centery - $H/2)")
zoff=$(python -c "print($centerz - $W/2)")
slowness_min=2.5
slowness_max=3
nshape=200
I=$(python -c "print(int(round($x0max / $h)+1))")
J=$(python -c "print(int(round($y0max / $h)+1))")
K=$(python -c "print(int(round($z0max / $h)+1))")
xoffint=$(python -c "print(int(round($xoff / $h)))")
yoffint=$(python -c "print(int(round($yoff / $h)))")
zoffint=$(python -c "print(int(round($zoff / $h)))")
Iobj=$(python -c "print(int(round($L / $h)+1))")
Jobj=$(python -c "print(int(round($H / $h)+1))")
Kobj=$(python -c "print(int(round($W / $h)+1))")
# python ../EXAMPLE/slowness_model3d_generator.py --shape_background ${I} ${J} ${K} --shape ${Iobj} ${Jobj} ${Kobj} --offset ${xoffint} ${yoffint} ${zoffint} --num_shapes ${nshape} --slowness_min ${slowness_min} --slowness_max ${slowness_max} --plot 0
# # echo $I $J $K $Iobj $Jobj $Kobj $xoffint $yoffint $zoffint





sample_para_outer=2.0
sample_para=2.0
elem_extract=0 # 2 is more OMP parallel, but cvie2d_t only supports 0 now. 
format=2
LRlevel=0
srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores  ./EXAMPLE/cvie3d --ivelo ${ivelo} --scaleGreen ${scaleGreen} --omega ${omega} --h ${h} --x0max $x0max --y0max $y0max --z0max $z0max --xmax $xmax --ymax $ymax --zmax $zmax --TNx ${TNx} --TNy ${TNy} --L ${L} --H ${H} --W ${W} --nshape ${nshape} --smin_ivelo11 ${slowness_min} --smax_ivelo11 ${slowness_max} --h0 ${h0} --vs ${vs} --shape ${shape} --knn ${knn} --lrlevel ${LRlevel} --n_iter ${n_iter} --format ${format} --baca_batch ${BACA_Batch} --knn_near_para ${knn_near_para} --elem_extract ${elem_extract} --near_para ${near_para} --xyzsort ${xyzsort} --nmin_leaf ${nmin_leaf_m} --sample_para_outer ${sample_para_outer} --sample_para ${sample_para} --forwardN15flag ${forwardN15flag} --reclr_leaf 2 --rmax $rmax --tol_comp ${tol_comp} --tol_rand ${tol_rand} --tol_Rdetect ${tol_Rdetect} --tol_itersol ${tol_itersol} --precon ${precon} | tee grep a.out_matrix_ivelo_${ivelo}3D_omega_${omega}_h_${h}_h0_${h0}_knn_${knn}_knn_near_para_${knn_near_para}_nmin_leaf_${nmin_leaf}_vs_${vs}_shape_${shape}_sample_para_${sample_para}_sample_para_outer_${sample_para_outer}_tol_comp${tol_comp}_tol_itersol_${tol_itersol}_precon_${precon}_format_${format}


# sample_para_outer=32.0
# sample_para=32.0
# elem_extract=0 # 2 is more OMP parallel, but cvie2d_t only supports 0 now. 
# format=3
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores  ./EXAMPLE/cvie3d --ivelo ${ivelo} --scaleGreen ${scaleGreen} --omega ${omega} --h ${h} --x0max $x0max --y0max $y0max --z0max $z0max --xmax $xmax --ymax $ymax --zmax $zmax --TNx ${TNx} --TNy ${TNy} --L ${L} --H ${H} --W ${W} --nshape ${nshape} --smin_ivelo11 ${slowness_min} --smax_ivelo11 ${slowness_max} --h0 ${h0} --vs ${vs} --shape ${shape} --knn ${knn} --lrlevel ${LRlevel} --n_iter ${n_iter} --format ${format} --baca_batch ${BACA_Batch} --knn_near_para ${knn_near_para} --elem_extract ${elem_extract} --near_para ${near_para} --xyzsort ${xyzsort} --nmin_leaf ${nmin_leaf_m} --sample_para_outer ${sample_para_outer} --sample_para ${sample_para} --forwardN15flag ${forwardN15flag} --rmax $rmax --tol_comp ${tol_comp} --tol_rand ${tol_rand} --tol_Rdetect ${tol_Rdetect} --tol_itersol ${tol_itersol} --precon ${precon} | tee grep a.out_matrix_ivelo_${ivelo}3D_omega_${omega}_h_${h}_h0_${h0}_knn_${knn}_knn_near_para_${knn_near_para}_nmin_leaf_${nmin_leaf}_vs_${vs}_shape_${shape}_sample_para_${sample_para}_sample_para_outer_${sample_para_outer}_tol_comp${tol_comp}_tol_itersol_${tol_itersol}_precon_${precon}_format_${format}


# sample_para_outer=32.0
# sample_para=16.0
# elem_extract=0 # 2 is more OMP parallel, but cvie2d_t only supports 0 now. 
# format=1
# srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores  ./EXAMPLE/cvie3d --ivelo ${ivelo} --scaleGreen ${scaleGreen} --omega ${omega} --h ${h} --x0max $x0max --y0max $y0max --z0max $z0max --xmax $xmax --ymax $ymax --zmax $zmax --TNx ${TNx} --TNy ${TNy} --L ${L} --H ${H} --W ${W} --nshape ${nshape} --smin_ivelo11 ${slowness_min} --smax_ivelo11 ${slowness_max} --h0 ${h0} --vs ${vs} --shape ${shape} --knn ${knn} --lrlevel ${LRlevel} --n_iter ${n_iter} --format ${format} --baca_batch ${BACA_Batch} --knn_near_para ${knn_near_para} --elem_extract ${elem_extract} --near_para ${near_para} --xyzsort ${xyzsort} --nmin_leaf ${nmin_leaf_m} --sample_para_outer ${sample_para_outer} --sample_para ${sample_para} --forwardN15flag ${forwardN15flag} --rmax $rmax --tol_comp ${tol_comp} --tol_rand ${tol_rand} --tol_Rdetect ${tol_Rdetect} --tol_itersol ${tol_itersol} --precon ${precon} | tee grep a.out_matrix_ivelo_${ivelo}3D_omega_${omega}_h_${h}_h0_${h0}_knn_${knn}_knn_near_para_${knn_near_para}_nmin_leaf_${nmin_leaf}_vs_${vs}_shape_${shape}_sample_para_${sample_para}_sample_para_outer_${sample_para_outer}_tol_comp${tol_comp}_tol_itersol_${tol_itersol}_precon_${precon}_format_${format}



elem_extract=0
sample_para_outer=0.8
sample_para=0.8
format=4
srun -n $nmpi -c $THREADS_PER_RANK --cpu_bind=cores  ./EXAMPLE/cvie3d_t --ivelo ${ivelo} --scaleGreen ${scaleGreen} --omega ${omega} --h ${h} --x0max $x0max --y0max $y0max --z0max $z0max --xmax $xmax --ymax $ymax --zmax $zmax --TNx ${TNx} --TNy ${TNy} --L ${L} --H ${H} --W ${W} --nshape ${nshape} --smin_ivelo11 ${slowness_min} --smax_ivelo11 ${slowness_max} --h0 ${h0} --vs ${vs} --shape ${shape} --knn ${knn} --lrlevel ${LRlevel} --n_iter ${n_iter} --format ${format} --baca_batch ${BACA_Batch} --knn_near_para ${knn_near_para} --elem_extract ${elem_extract} --near_para ${near_para} --xyzsort ${xyzsort} --nmin_leaf ${nmin_leaf_t} --sample_para_outer ${sample_para_outer} --sample_para ${sample_para} --rmax $rmax --tol_comp ${tol_comp} --tol_rand ${tol_rand} --tol_Rdetect ${tol_Rdetect} --fastsample_tensor ${fastsample_tensor} --verbosity ${verbosity} --tol_itersol ${tol_itersol} --precon ${precon} | tee grep a.out_tensor_ivelo_${ivelo}3D_omega_${omega}_h_${h}_h0_${h0}_knn_${knn}_knn_near_para_${knn_near_para}_nmin_leaf_${nmin_leaf}_vs_${vs}_shape_${shape}_sample_para_${sample_para}_sample_para_outer_${sample_para_outer}_tol_comp${tol_comp}_tol_itersol_${tol_itersol}_precon_${precon}

