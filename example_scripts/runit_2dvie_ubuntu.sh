#!/bin/bash -l


export OMP_NUM_THREADS=1

nmpi=2


vs=1
shape=4

tol_rand=1e-2
tol_Rdetect=0.3e-2

# tol=1d-4
LRlevel=100
h0=0.005
nmin_leaf_t=4
nmin_leaf_m=16
xyzsort=1
BACA_Batch=100
near_para=2.1d0


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




knn_near_para=20.0
knn=0
forwardN15flag=0
rmax=1000
sample_para=0.8
sample_para_outer=0.8
verbosity=1


# omega=25.132741228718345
# ivelo=1
# TNx=3
# TNy=3
# tol_comp=1e-8
# h=0.0125 # 0.015625 #0.025 # 0.0125 # 0.00625
# x0max=2.0
# y0max=2.0
# xmax=2.2
# ymax=2.2


# omega=31.415926535897750
# ivelo=1
# TNx=3
# TNy=3
# tol_comp=1e-8
# h=0.01 #0.000625 #0.0015625 #0.00125 #0.00125 #0.0015625
# x0max=2.0
# y0max=2.0
# xmax=2.2
# ymax=2.2


# omega=62.831853071795500
# ivelo=1
# TNx=3
# TNy=3
# tol_comp=1e-8
# h=0.005 #0.000625 #0.0015625 #0.00125 #0.00125 #0.0015625
# x0max=2.0
# y0max=2.0
# xmax=2.2
# ymax=2.2

# omega=125.663706143591
# ivelo=1
# TNx=3
# TNy=3
# tol_comp=1e-8
# h=0.0025 #0.000625 #0.0015625 #0.00125 #0.00125 #0.0015625
# x0max=2.0
# y0max=2.0
# xmax=2.2
# ymax=2.2



# omega=251.32741228718345
# ivelo=1
# TNx=3
# TNy=3
# tol_comp=1e-8
# h=0.00125 #0.000625 #0.0015625 #0.00125 #0.00125 #0.0015625
# x0max=2.0
# y0max=2.0
# xmax=2.2
# ymax=2.2


# omega=502.6548245743669
# ivelo=1
# TNx=3
# TNy=3
# tol_comp=1e-8
# h=0.000625 #0.0015625 #0.00125 #0.00125 #0.0015625
# x0max=2.0
# y0max=2.0
# xmax=2.2
# ymax=2.2

# omega=1005.309649148734
# ivelo=1
# TNx=3
# TNy=3
# tol_comp=1e-8
# h=0.0003125 #0.0015625 #0.00125 #0.00125 #0.0015625
# x0max=2.0
# y0max=2.0
# xmax=2.2
# ymax=2.2




omega=7.853981633
ivelo=9
TNx=3
TNy=3
tol_comp=1e-6
h=0.02 #0.0125 #0.005
x0max=1.0
y0max=1.0
xmax=1.2
ymax=1.2


# omega=15.70796326794896
# ivelo=9
# TNx=3
# TNy=3
# tol_comp=1e-6
# h=0.01 #0.0125 #0.005
# x0max=1.0
# y0max=1.0
# xmax=1.2
# ymax=1.2


# omega=39.269908169872402
# ivelo=9
# TNx=9
# TNy=13
# tol_comp=1e-8
# h=0.004
# x0max=1.0
# y0max=1.0
# xmax=1.2
# ymax=1.2


# omega=78.539816339744803
# ivelo=9
# TNx=9
# TNy=13
# tol_comp=1e-8
# h=0.002
# x0max=1.0
# y0max=1.0
# xmax=1.2
# ymax=1.2


# omega=157.0796326794896
# ivelo=9
# TNx=3
# TNy=3
# tol_comp=1e-8
# h=0.001
# x0max=1.0
# y0max=1.0
# xmax=1.2
# ymax=1.2

# omega=314.1592653589792
# ivelo=9
# TNx=9
# TNy=13
# tol_comp=1e-8
# h=0.0005
# x0max=1.0
# y0max=1.0
# xmax=1.2
# ymax=1.2

# omega=628.3185307179584
# ivelo=9
# TNx=9
# TNy=13
# tol_comp=1e-8
# h=0.00025
# x0max=1.0
# y0max=1.0
# xmax=1.2
# ymax=1.2



# omega=125.6637061435917
# ivelo=9
# TNx=9
# TNy=13
# tol_comp=1e-7
# h=0.00125
# x0max=1.0
# y0max=1.0
# xmax=1.2
# ymax=1.2


sample_para_outer=2.0
sample_para=2.0
elem_extract=0 # 2 is more OMP parallel, but cvie2d_t only supports 0 now. 
format=3
mpirun --allow-run-as-root -n $nmpi ./EXAMPLE/cvie2d --ivelo ${ivelo} --omega ${omega} --h ${h} --x0max $x0max --y0max $y0max --xmax $xmax --ymax $ymax --TNx ${TNx} --TNy ${TNy} --h0 ${h0} --vs ${vs} --shape ${shape} --knn ${knn} --lrlevel ${LRlevel} --format ${format} --baca_batch ${BACA_Batch} --knn_near_para ${knn_near_para} --elem_extract ${elem_extract} --near_para ${near_para} --xyzsort ${xyzsort} --nmin_leaf ${nmin_leaf_m} --sample_para_outer ${sample_para_outer} --sample_para ${sample_para} --forwardN15flag ${forwardN15flag} --rmax $rmax --tol_comp ${tol_comp} --tol_rand ${tol_rand} --tol_Rdetect ${tol_Rdetect} | tee grep a.out_matrix_ivelo_${ivelo}_TNx_${TNx}_TNy_${TNy}_omega_${omega}_h_${h}_h0_${h0}_knn_${knn}_knn_near_para_${knn_near_para}_nmin_leaf_${nmin_leaf}_vs_${vs}_shape_${shape}_sample_para_${sample_para}_sample_para_outer_${sample_para_outer}_tol_comp${tol_comp}
#elem_extract=0
# sample_para_outer=0.8
# sample_para=0.8
# format=4
# mpirun --allow-run-as-root -n $nmpi ./EXAMPLE/cvie2d_t --ivelo ${ivelo} --omega ${omega} --h ${h} --x0max $x0max --y0max $y0max --xmax $xmax --ymax $ymax --TNx ${TNx} --TNy ${TNy} --h0 ${h0} --vs ${vs} --shape ${shape} --knn ${knn} --lrlevel ${LRlevel} --format ${format} --baca_batch ${BACA_Batch} --knn_near_para ${knn_near_para} --elem_extract ${elem_extract} --near_para ${near_para} --xyzsort ${xyzsort} --nmin_leaf ${nmin_leaf_t} --sample_para_outer ${sample_para_outer} --sample_para ${sample_para} --forwardN15flag ${forwardN15flag} --rmax $rmax --tol_comp ${tol_comp} --tol_rand ${tol_rand} --tol_Rdetect ${tol_Rdetect} --fastsample_tensor 2 --verbosity ${verbosity} | tee grep a.out_tensor_ivelo_${ivelo}_TNx_${TNx}_TNy_${TNy}_omega_${omega}_h_${h}_h0_${h0}_knn_${knn}_knn_near_para_${knn_near_para}_nmin_leaf_${nmin_leaf}_vs_${vs}_shape_${shape}_sample_para_${sample_para}_sample_para_outer_${sample_para_outer}_tol_comp${tol_comp}
