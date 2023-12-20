if [[ $(uname -s) == 'Darwin' ]]; then
    export GPTUNEROOT=/Users/liuyangzhuan/Desktop/GPTune/
    export MPIRUN="$GPTUNEROOT/openmpi-4.1.5/bin/mpirun"
else
    export MPIRUN=mpirun
fi



tol=1e-6
export OMP_NUM_THREADS=4

$MPIRUN --allow-run-as-root -n 4 ../build/EXAMPLE/frankben_t -quant --tst 3 --wavelen 0.015625 -option --nmin_leaf 8 --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 1.0 --sample_para_outer 1.0 --fastsample_tensor 2 | tee a.out_tensor_3d_green
$MPIRUN --allow-run-as-root -n 4 ../build/EXAMPLE/frankben -quant --tst 3 --wavelen 0.015625 -option --nmin_leaf 8 --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 2.0 --sample_para_outer 2.0  | tee a.out_matrix_3d_green



# $MPIRUN --allow-run-as-root -n 4 ../build/EXAMPLE/frankben_t -quant --tst 2 --wavelen 0.0156 -option --nmin_leaf 8 --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 0.8 --sample_para_outer 0.8 --fastsample_tensor 2 | tee a.out_tensor
# # $MPIRUN --allow-run-as-root -n 4 ../build/EXAMPLE/frankben_t -quant --tst 2 --wavelen 0.0156 -option --nmin_leaf 8 --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 1.0 --sample_para_outer 1.0 --fastsample_tensor 1 | tee a.out_tensor
# $MPIRUN --allow-run-as-root -n 4 ../build/EXAMPLE/frankben_t -quant --tst 2 --wavelen 0.0156 -option --nmin_leaf 8 --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --sample_para 0.8 --sample_para_outer 0.8 | tee a.out_tensor
# $MPIRUN --allow-run-as-root -n 4 ../build/EXAMPLE/frankben -quant --tst 2 --wavelen 0.0156 -option --nmin_leaf 64  --xyzsort 1 --lrlevel 100 --verbosity 1 --tol_comp $tol --pat_comp 3 --sample_para 2.0 --sample_para_outer 2.0 | tee a.out_matrix
