if [[ $(uname -s) == 'Darwin' ]]; then
    export GPTUNEROOT=/Users/liuyangzhuan/GPTune_testbuild/
    export MPIRUN="$GPTUNEROOT/openmpi-5.0.6/bin/mpirun"
else
    export MPIRUN=mpirun
fi



tol=1e-6
export OMP_NUM_THREADS=1


# ############## 3D seperated cubes							 
export OMP_NUM_THREADS=1
tol=1e-4
nmpi=8
# wavelen=0.001953125 0.00390625 0.0078125 0.015625 0.03125 0.0625
zdist=1.0
# zdist=0.001
ppw=4.0
nmin_leaf_t=4
nmin_leaf_m=64
use_zfp=1
use_qtt=1
for wavelen in 0.0625 
do
$MPIRUN --allow-run-as-root -n ${nmpi} ../build/EXAMPLE/frankben_t -quant --tst 3 --wavelen ${wavelen} --zdist ${zdist} --ppw ${ppw} -option --use_zfp ${use_zfp} --use_qtt ${use_qtt} --xyzsort 0 --nmin_leaf ${nmin_leaf_t} --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 0.8 --sample_para_outer 0.8 --fastsample_tensor 2 | tee a.out_tensor_3d_green_wavelen${wavelen}_zdist${zdist}_tol${tol}_mpi${nmpi}_omp${NTH}_nmin_leaf_t${nmin_leaf_t}_ppw${ppw}_use_zfp${use_zfp}_use_qtt${use_qtt}
# $MPIRUN --allow-run-as-root -n ${nmpi} ../build/EXAMPLE/frankben -quant --tst 3 --wavelen ${wavelen} --zdist ${zdist} --ppw ${ppw} -option --xyzsort 1 --nmin_leaf ${nmin_leaf_m} --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 2.0 --knn 10 --sample_para_outer 2.0  | tee a.out_matrix_3d_green_wavelen${wavelen}_zdist${zdist}_tol${tol}_mpi${nmpi}_omp${NTH}_nmin_leaf_m${nmin_leaf_m}_ppw${ppw}
done



# ############## 2D parallel plates								 
# export OMP_NUM_THREADS=1
# tol=1e-6
# nmpi=1
# # wavelen=0.00012207031 #0.00024414062 #0.00048828125 #0.0009765625 #0.001953125 #0.00390625 #0.0078125 #0.015625
# zdist=1.0
# ppw=4.0
# nmin_leaf_t=16
# nmin_leaf_m=256
# use_zfp=0
# use_qtt=1
# # ppw=2.0
# # nmin_leaf_t=8
# # nmin_leaf_m=64
# for wavelen in 0.015625
# do
# $MPIRUN --allow-run-as-root -n ${nmpi} ../build/EXAMPLE/frankben_t -quant --tst 2 --wavelen ${wavelen} --zdist ${zdist} --ppw ${ppw} -option --nmin_leaf ${nmin_leaf_t} --use_zfp ${use_zfp} --use_qtt ${use_qtt} --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 0.8 --sample_para_outer 0.8 --fastsample_tensor 2 | tee a.out_tensor_2d_green_wavelen${wavelen}_zdist${zdist}_tol${tol}_mpi${nmpi}_omp${NTH}_nmin_leaf_t${nmin_leaf_t}_ppw${ppw}_use_zfp${use_zfp}
# # $MPIRUN --allow-run-as-root -n ${nmpi} ../build/EXAMPLE/frankben -quant --tst 2 --wavelen ${wavelen} --zdist ${zdist} --ppw ${ppw} -option --nmin_leaf ${nmin_leaf_m}  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 | tee a.out_matrix_2d_green_wavelen${wavelen}_zdist${zdist}_tol${tol}_mpi${nmpi}_omp${NTH}_nmin_leaf_m${nmin_leaf_m}_ppw${ppw}
# done 










# ############## DFT
# tol=1e-3
# # $MPIRUN --allow-run-as-root -n 4 ../build/EXAMPLE/frankben -quant --tst 4 --ndim_FIO 4 --N_FIO 8 -option --nmin_leaf 64  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 | tee a.out_matrix_DFT
# $MPIRUN --allow-run-as-root -n 3 ../build/EXAMPLE/frankben_t -quant --tst 4 --ndim_FIO 2 --N_FIO 16 -option --nmin_leaf 8  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 --fastsample_tensor 2 | tee a.out_tensor_DFT








# tol=1e-4
# nmpi=8
# # wavelen=0.001953125 0.00390625 0.0078125 0.015625 0.03125 0.0625
# zdist=1.0
# # zdist=0.001
# ppw=4.0
# nmin_leaf_t=4
# nmin_leaf_m=64
# use_zfp=0
# use_qtt=1
# for wavelen in 0.0625 
# do
# $MPIRUN --allow-run-as-root -n ${nmpi} ../build/EXAMPLE/frankben_t -quant --tst 3 --wavelen ${wavelen} --zdist ${zdist} --ppw ${ppw} -option --use_zfp ${use_zfp} --use_qtt ${use_qtt} --xyzsort 1 --nmin_leaf ${nmin_leaf_t} --lrlevel 100 --verbosity 2 --tol_comp $tol --sample_para 0.8 --sample_para_outer 0.8 --fastsample_tensor 2 | tee a.out_tensor_3d_green_wavelen${wavelen}_zdist${zdist}_tol${tol}_mpi${nmpi}_omp${NTH}_nmin_leaf_t${nmin_leaf_t}_ppw${ppw}_use_zfp${use_zfp}_use_qtt${use_qtt}
# # $MPIRUN --allow-run-as-root -n ${nmpi} ../build/EXAMPLE/frankben -quant --tst 3 --wavelen ${wavelen} --zdist ${zdist} --ppw ${ppw} -option --xyzsort 1 --nmin_leaf ${nmin_leaf_m} --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 2.0 --knn 10 --sample_para_outer 2.0
# done