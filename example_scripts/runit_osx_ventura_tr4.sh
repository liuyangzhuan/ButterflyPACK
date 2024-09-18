if [[ $(uname -s) == 'Darwin' ]]; then
    export GPTUNEROOT=/Users/liuyangzhuan/Desktop/GPTune/
    export MPIRUN="$GPTUNEROOT/openmpi-4.1.5/bin/mpirun"
else
    export MPIRUN=mpirun
fi



tol=1e-6
export OMP_NUM_THREADS=1

# ############## 3D seperated cubes
# $MPIRUN --allow-run-as-root -n 4 ../build/EXAMPLE/frankben_t -quant --tst 3 --wavelen 0.0625 -option --nmin_leaf 8 --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 1.0 --sample_para_outer 1.0 --fastsample_tensor 2 --knn 0 | tee a.out_tensor_3d_green
# # $MPIRUN --allow-run-as-root -n 4 ../build/EXAMPLE/frankben -quant --tst 3 --wavelen 0.015625 -option --nmin_leaf 8 --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 2.0 --sample_para_outer 2.0  | tee a.out_matrix_3d_green



############## 2D parallel plates
tol=1e-6
nmpi=1
# wavelen=0.00012207031 #0.00024414062 #0.00048828125 #0.0009765625 #0.001953125 #0.00390625 #0.0078125 #0.015625
zdist=1.0
ppw=4.0
nmin_leaf_t=16
nmin_leaf_m=256
lrlevel=0
# ppw=2.0
# nmin_leaf_t=8
# nmin_leaf_m=64
for wavelen in 0.015625
do
$MPIRUN --allow-run-as-root -n ${nmpi} ../build/EXAMPLE/frankben_t -quant --tst 2 --wavelen ${wavelen} --zdist ${zdist} --ppw ${ppw} -option --nmin_leaf ${nmin_leaf_t} --xyzsort 1 --lrlevel ${lrlevel} --verbosity 1 --tol_comp $tol --sample_para 0.8 --sample_para_outer 0.8 --fastsample_tensor 2 | tee a.out_tensor_2d_green_wavelen${wavelen}_zdist${zdist}_tol${tol}_mpi${nmpi}_omp${NTH}_nmin_leaf_t${nmin_leaf_t}_ppw${ppw}_lrlevel${lrlevel}
# $MPIRUN --allow-run-as-root -n ${nmpi} ../build/EXAMPLE/frankben -quant --tst 2 --wavelen ${wavelen} --zdist ${zdist} --ppw ${ppw} -option --nmin_leaf ${nmin_leaf_m}  --xyzsort 1 --lrlevel ${lrlevel} --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 | tee a.out_matrix_2d_green_wavelen${wavelen}_zdist${zdist}_tol${tol}_mpi${nmpi}_omp${NTH}_nmin_leaf_m${nmin_leaf_m}_ppw${ppw}_lrlevel${lrlevel}
done 
# $MPIRUN --allow-run-as-root -n 4 ../build/EXAMPLE/frankben_t -quant --tst 2 --wavelen 0.0156 -option --nmin_leaf 8 --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 0.8 --sample_para_outer 0.8 | tee a.out_tensor_2d_green
# $MPIRUN --allow-run-as-root -n 4 ../build/EXAMPLE/frankben -quant --tst 2 --wavelen 0.0156 -option --nmin_leaf 64  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 | tee a.out_matrix_2d_green


# ############## DFT
# tol=1e-3
# # $MPIRUN --allow-run-as-root -n 4 ../build/EXAMPLE/frankben -quant --tst 4 --ndim_FIO 4 --N_FIO 8 -option --nmin_leaf 64  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 | tee a.out_matrix_DFT
# $MPIRUN --allow-run-as-root -n 3 ../build/EXAMPLE/frankben_t -quant --tst 4 --ndim_FIO 2 --N_FIO 16 -option --nmin_leaf 8  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 --fastsample_tensor 2 | tee a.out_tensor_DFT
