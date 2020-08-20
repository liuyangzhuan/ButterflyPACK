#!/bin/bash -l

module purge
module load gcc/9.1.0
module load openmpi/gcc-9.1.0/4.0.1
module load scalapack-netlib/gcc-9.1.0/2.0.2 

NTH=1

export OMP_NUM_THREADS=$NTH

knn=100
which='LM'
si=1
nev=4
cmmode=0
sample_para=2.0d0 
lrlevel=0
tol=1d-4
scaling=1d0
# model="pillbox_4000"
model="cavity_wakefield_4K_feko"
# 010 110 111
# for freq in 1.826e9    
for freq in 765e6   
# for freq in 1.145e9 
# for freq in 1.145e9 1.146e9 1.14743e9 1.148e9 1.149e9 1.826e9 1.827e9 1.82824e9 1.829e9 1.830e9 
do
mpirun --allow-run-as-root -n 16 ./EXAMPLE/ie3deigen -quant --data_dir ../EXAMPLE/EM3D_DATA/preprocessor_3dmesh/$model --freq $freq --scaling ${scaling} --si ${si} --which ${which} --nev ${nev} --cmmode ${cmmode} -option --verbosity 2 --reclr_leaf 5 --baca_batch 64 --tol_comp 1e-4 --lrlevel $lrlevel --precon 1 --xyzsort 2 --nmin_leaf 100 --near_para 0.01d0 --pat_comp 3 --format 1 --sample_para $sample_para --knn $knn | tee a.out_freq_${freq}_lrlevel_${lrlevel}_noport
# mpirun --allow-run-as-root -n 16 ./EXAMPLE/ie3dporteigen -quant --data_dir ../EXAMPLE/EM3D_DATA/preprocessor_3dmesh/$model --freq $freq --scaling ${scaling} --si ${si} --which ${which} --nev ${nev} --cmmode ${cmmode} -option --verbosity 2 --reclr_leaf 5 --baca_batch 64 --tol_comp ${tol} --lrlevel $lrlevel --precon 1 --xyzsort 2 --nmin_leaf 100 --near_para 0.01d0 --pat_comp 3 --format 1 --sample_para $sample_para --knn $knn | tee a.out_freq_${freq}_lrlevel_${lrlevel}_port
done

