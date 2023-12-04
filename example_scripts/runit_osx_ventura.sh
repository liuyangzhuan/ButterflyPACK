export GPTUNEROOT=/Users/liuyangzhuan/Desktop/GPTune/
export MPIRUN="$GPTUNEROOT/openmpi-4.1.5/bin/mpirun"
# $MPIRUN --allow-run-as-root -n 4 ../build/EXAMPLE/frankben_t -quant --tst 3 --wavelen 0.0625 -option --nmin_leaf 10 --lrlevel 100 --verbosity 1 --sample_para 1.0 --sample_para_outer 1.0 | tee a.out_tensor
$MPIRUN --allow-run-as-root -n 4 ../build/EXAMPLE/frankben_t -quant --tst 2 --wavelen 0.03125 -option --nmin_leaf 10 --lrlevel 100 --verbosity 1 --sample_para 1.0 --sample_para_outer 1.0 | tee a.out_tensor

$MPIRUN --allow-run-as-root -n 4 ../build/EXAMPLE/frankben -quant --tst 2 --wavelen 0.03125 -option --nmin_leaf 200 --lrlevel 100 --verbosity 1 --sample_para 2.0 --sample_para_outer 2.0 | tee a.out_matrix
