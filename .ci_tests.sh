#!/bin/sh
set -e

export RED="\033[31;1m"
export BLUE="\033[34;1m"
export ROOT_DIR="$PWD"
printf "${BLUE} YL; Entered tests file:\n"

export DATA_FOLDER=$ROOT_DIR/EXAMPLE
export EXAMPLE_FOLDER=$ROOT_DIR/build/EXAMPLE
# export TEST_FOLDER=$ROOT_DIR/build/TEST

# export PATH="$PATH:$ROOT_DIR/installDir/openmpi-3.1.4/install/bin"
# export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$ROOT_DIR/installDir/openmpi-3.1.4/install/lib/" 

case "${TEST_NUMBER}" in
1) cd $DATA_FOLDER/FULLMAT_DATA
   sh file_merge.sh Smatrix.mat
   cd $ROOT_DIR/build
   mpirun --allow-run-as-root --oversubscribe "-n" "1" "$EXAMPLE_FOLDER/smat" "-quant" "--data_dir" "$DATA_FOLDER/FULLMAT_DATA" "--explicitflag" "0" "-option" "--lrlevel" "0" ;;# test sampling-based construction of HODLR with LR for an EM scattering matrix stored in file
2) mpirun --allow-run-as-root --oversubscribe "-n" "2" "$EXAMPLE_FOLDER/krr" "-quant" "--data_dir" "$DATA_FOLDER/KRR_DATA/susy_10Kn" "--dimn" "8" "--ntrain" "10000" "--ntest" "1000" "--sigma" "0.1" "--lambda" "1.0" "-option" "--xyzsort" "2" "--reclr_leaf" "2";; # test HODLR with LR for KRR
3) mpirun --allow-run-as-root --oversubscribe "-n" "4" "$EXAMPLE_FOLDER/ie2d" "-quant" "--model2d" "10" "--nunk" "5000" "--wavelength" "0.08" "-option" "--lr_blk_num" "1" "--tol_comp" "1d-4" "--errfillfull" "0" "--reclr_leaf" "4" "--baca_batch" "16" "--lrlevel" "0" "--precon" "1" "--xyzsort" "2" "--nmin_leaf" "200" "--near_para" "0.01d0" "--pat_comp" "3" "--schulzlevel" "100" "--nbundle" "1" "--format" "1"
   mpirun --allow-run-as-root --oversubscribe "-n" "4" "$EXAMPLE_FOLDER/ie2d" "-quant" "--model2d" "10" "--nunk" "5000" "--wavelength" "0.08" "-option" "--lr_blk_num" "1" "--tol_comp" "1d-4" "--errfillfull" "0" "--reclr_leaf" "4" "--baca_batch" "16" "--lrlevel" "100" "--precon" "1" "--xyzsort" "2" "--nmin_leaf" "200" "--near_para" "0.01d0" "--pat_comp" "3" "--schulzlevel" "100" "--nbundle" "1" "--format" "1";;   # test HODLR with BF and LR for 2d IE
4) mpirun --allow-run-as-root --oversubscribe "-n" "4" "$EXAMPLE_FOLDER/ie2d" "-quant" "--model2d" "10" "--nunk" "5000" "--wavelength" "0.08" "-option" "--lr_blk_num" "1" "--tol_comp" "1d-4" "--errfillfull" "0" "--reclr_leaf" "4" "--baca_batch" "16" "--lrlevel" "0" "--precon" "1" "--xyzsort" "2" "--nmin_leaf" "200" "--near_para" "2.01d0" "--pat_comp" "3" "--schulzlevel" "100" "--nbundle" "1" "--format" "2"
   mpirun --allow-run-as-root --oversubscribe "-n" "4" "$EXAMPLE_FOLDER/ie2d" "-quant" "--model2d" "10" "--nunk" "5000" "--wavelength" "0.08" "-option" "--lr_blk_num" "1" "--tol_comp" "1d-4" "--errfillfull" "0" "--reclr_leaf" "4" "--baca_batch" "16" "--lrlevel" "100" "--precon" "1" "--xyzsort" "2" "--nmin_leaf" "200" "--near_para" "2.01d0" "--pat_comp" "3" "--schulzlevel" "100" "--nbundle" "1" "--format" "2";;   # test HMAT with BF and LR for 2d IE
5) mpirun --allow-run-as-root --oversubscribe "-n" "8" "$EXAMPLE_FOLDER/ie3d" "-quant" "--data_dir" "$DATA_FOLDER/EM3D_DATA/sphere_2300" "--wavelength" "2.0" "-option" "--lr_blk_num" "1" "--tol_comp" "1e-2" "--errfillfull" "0" "--reclr_leaf" "4" "--baca_batch" "1" "--lrlevel" "0" "--precon" "1" "--xyzsort" "2" "--nmin_leaf" "100" "--near_para" "0.01d0" "--pat_comp" "3" "--format" "1"
   mpirun --allow-run-as-root --oversubscribe "-n" "8" "$EXAMPLE_FOLDER/ie3d" "-quant" "--data_dir" "$DATA_FOLDER/EM3D_DATA/sphere_2300" "--wavelength" "2.0" "-option" "--lr_blk_num" "1" "--tol_comp" "1e-2" "--errfillfull" "0" "--reclr_leaf" "4" "--baca_batch" "1" "--lrlevel" "100" "--precon" "1" "--xyzsort" "2" "--nmin_leaf" "100" "--near_para" "0.01d0" "--pat_comp" "3" "--format" "1"
   mpirun --allow-run-as-root --oversubscribe "-n" "8" "$EXAMPLE_FOLDER/ie3d_sp" "-quant" "--data_dir" "$DATA_FOLDER/EM3D_DATA/sphere_2300" "--wavelength" "2.0" "-option" "--lr_blk_num" "1" "--tol_comp" "1e-2" "--errfillfull" "0" "--reclr_leaf" "4" "--baca_batch" "1" "--lrlevel" "0" "--precon" "1" "--xyzsort" "2" "--nmin_leaf" "100" "--near_para" "0.01d0" "--pat_comp" "3" "--format" "1"
   mpirun --allow-run-as-root --oversubscribe "-n" "8" "$EXAMPLE_FOLDER/ie3d_sp" "-quant" "--data_dir" "$DATA_FOLDER/EM3D_DATA/sphere_2300" "--wavelength" "2.0" "-option" "--lr_blk_num" "1" "--tol_comp" "1e-2" "--errfillfull" "0" "--reclr_leaf" "4" "--baca_batch" "1" "--lrlevel" "100" "--precon" "1" "--xyzsort" "2" "--nmin_leaf" "100" "--near_para" "0.01d0" "--pat_comp" "3" "--format" "1";; # test HODLR with BF and LR for 3d IE (double and single precision)  
6) mpirun --allow-run-as-root --oversubscribe "-n" "8" "$EXAMPLE_FOLDER/ie3d" "-quant" "--data_dir" "$DATA_FOLDER/EM3D_DATA/sphere_2300" "--wavelength" "2.0" "-option" "--lr_blk_num" "1" "--tol_comp" "1e-2" "--errfillfull" "0" "--reclr_leaf" "4" "--baca_batch" "1" "--lrlevel" "0" "--precon" "1" "--xyzsort" "2" "--nmin_leaf" "100" "--near_para" "2.01d0" "--pat_comp" "3" "--format" "2"
   mpirun --allow-run-as-root --oversubscribe "-n" "8" "$EXAMPLE_FOLDER/ie3d" "-quant" "--data_dir" "$DATA_FOLDER/EM3D_DATA/sphere_2300" "--wavelength" "2.0" "-option" "--lr_blk_num" "1" "--tol_comp" "1e-2" "--errfillfull" "0" "--reclr_leaf" "4" "--baca_batch" "1" "--lrlevel" "100" "--precon" "1" "--xyzsort" "2" "--nmin_leaf" "100" "--near_para" "2.01d0" "--pat_comp" "3" "--format" "2";; # test HMAT with BF and LR for 3d IE   
7) mpirun --allow-run-as-root --oversubscribe "-n" "2" "$EXAMPLE_FOLDER/ctest" "--tst" "1" "--N" "10000" "--trainfile" "$DATA_FOLDER/KRR_DATA/susy_10Kn" "--Ndim" "8" "--ker" "1" "--h" "0.1" "--lambda" "10.0" "--Nmin_leaf" "200" "--tol_comp" "1e-2" "--reclr_leaf" "5" "--errfillfull" "0" "--baca_batch" "16" "--knn" "10";; # test CPP interface with HODLR-LR for KRR with SUSY dataset
8) mpirun --allow-run-as-root --oversubscribe "-n" "7" "$EXAMPLE_FOLDER/ctest" "--tst" "2" "--N" "1000" "--Ndim" "8" "--ker" "1" "--h" "0.1" "--lambda" "10.0" "--Nmin_leaf" "100" "--tol_comp" "1e-4" "--reclr_leaf" "5" "--errfillfull" "0" "--baca_batch" "16" "--knn" "10";; # test CPP interface with HODLR-LR for KRR with random points
9) mpirun --allow-run-as-root --oversubscribe "-n" "3" "$EXAMPLE_FOLDER/ctest" "--tst" "3" "--N" "5000" "--rank_rand" "20" "--Nmin_leaf" "100" "--tol_comp" "1e-4" "--reclr_leaf" "5" "--errfillfull" "0";;  # test CPP interface with HODLR-LR for product of two random matrices 
10) mpirun --allow-run-as-root --oversubscribe "-n" "5" "$EXAMPLE_FOLDER/full";; # test HODLR with LR for reading a KRR matrix from a file 
11) mpirun --allow-run-as-root --oversubscribe "-n" "2" "$EXAMPLE_FOLDER/frontal" "-quant" "--data_dir" "$DATA_FOLDER/FULLMAT_DATA/Frontal_poisson_50.mat" "--nunk" "50" "--explicitflag" "0" "-option" "--lrlevel" "0" "--tol_comp" "1d-1" "--nmin_leaf" "100" # test sampling-based construction of HODLR with LR for an 3D poisson frontal matrix stored in file
   mpirun --allow-run-as-root --oversubscribe "-n" "2" "$EXAMPLE_FOLDER/frontal" "-quant" "--data_dir" "$DATA_FOLDER/FULLMAT_DATA/Frontal_poisson_50.mat" "--nunk" "50" "--explicitflag" "0" "-option" "--lrlevel" "100" "--tol_comp" "1d-1" "--nmin_leaf" "100" # test sampling-based construction of HODLR with BF for an 3D poisson frontal matrix stored in file
   mpirun --allow-run-as-root --oversubscribe "-n" "2" "$EXAMPLE_FOLDER/frontal" "-quant" "--data_dir" "$DATA_FOLDER/FULLMAT_DATA/Frontal_poisson_50.mat" "--nunk" "50" "--explicitflag" "1" "-option" "--lrlevel" "100" "--tol_comp" "1d-1" "--nmin_leaf" "100";; # test entry-evaluation-based construction of HODLR with BF for an 3D poisson frontal matrix stored in file   
12) mpirun --allow-run-as-root --oversubscribe "-n" "4" "$EXAMPLE_FOLDER/ie2deigen" "-quant" "--model2d" "10" "--nunk" "5000" "--wavelength" "0.08" "--cmmode" "0" "--si" "0" "-option" "--lr_blk_num" "1" "--tol_comp" "1d-4" "--errfillfull" "0" "--reclr_leaf" "4" "--baca_batch" "16" "--lrlevel" "0" "--precon" "1" "--xyzsort" "2" "--nmin_leaf" "200" "--near_para" "2.01d0" "--pat_comp" "3" "--schulzlevel" "100" "--nbundle" "1" "--format" "1";; # test an eigen solver for HODLR with LR for 2d IE 
13) mpirun --allow-run-as-root --oversubscribe "-n" "8" "$EXAMPLE_FOLDER/ie3deigen" "-quant" "--data_dir" "$DATA_FOLDER/EM3D_DATA/sphere_2300" "--wavelength" "2.0" "--si" "1" "--which" "LM" "--nev" "1" "--cmmode" "0" "-option" "--tol_comp" "1e-3" "--lrlevel" "0" "--precon" "1" "--xyzsort" "2" "--nmin_leaf" "100" "--near_para" "0.01d0" "--pat_comp" "3" "--format" "1";; # test an eigen solver for HODLR with LR for 3d IE
14) cd $DATA_FOLDER/FULLMAT_DATA
   sh file_merge.sh Frontal_elastic
   cd $ROOT_DIR/build
    mpirun --allow-run-as-root --oversubscribe "-n" "4" "$EXAMPLE_FOLDER/frontaldist" "-quant" "--nunk" "2500" "--data_dir" "$DATA_FOLDER/FULLMAT_DATA/Frontal_elastic/Frontal_elastic_2500" "--explicitflag" "0" "-option" "--nmin_leaf" "100" "--tol_comp" "1e-3" "--lrlevel" "100";; # test sampling-based construction of HODLR with BF for an 3D elastic Helmholtz frontal matrix stored parallel in file
# 15) cd $DATA_FOLDER/FULLMAT_DATA
#    sh file_merge.sh FullMatKrr
#    cd $ROOT_DIR/build
#     mpirun --allow-run-as-root --oversubscribe "-n" "8" "$EXAMPLE_FOLDER/fullkrr" "-quant" "-lambda" "0.03" "-option" "--reclr_leaf" "5" "--tol_comp" "1e-7" "--xyzsort" "3";; # test full matrix kernel regression with entry-based construction of HODLR with LR for a kernel matrix stored in file, with angular gram distance-based reordering and H-BACA
15) mpirun --allow-run-as-root --oversubscribe "-n" "8" "$EXAMPLE_FOLDER/go2d" "-quant" "--nunk" "5000" "--omega" "31.4" "-option" "--tol_comp" "1e-6" "--lrlevel" "100" "--xyzsort" "2" "--nmin_leaf" "100" "--format" "1" # test the taylor-expansion-based 2D approximate Green function 
    mpirun --allow-run-as-root --oversubscribe "-n" "8" "$EXAMPLE_FOLDER/go3d" "-quant" "--nunk" "5000" "--omega" "31.4" "-option" "--tol_comp" "1e-6" "--lrlevel" "100" "--xyzsort" "2" "--nmin_leaf" "100" "--format" "1";; # test the taylor-expansion-based 3D approximate Green function 
16) mpirun --allow-run-as-root --oversubscribe "-n" "7" "$EXAMPLE_FOLDER/cfio" "--tst" "2" "--M" "1500" "--N" "1000" "--K" "2000" "--L" "1800" "--lrlevel" "100" # test CPP interface for compressing product of three FIOs 
    mpirun --allow-run-as-root --oversubscribe "-n" "7" "$EXAMPLE_FOLDER/cifio" "--M" "1000" "--N" "1000" "--ker" "1" "--errsol 1" "--format" "1" # test CPP interface for approximate inverse of a 1D FIO, using HODLR    
    mpirun --allow-run-as-root --oversubscribe "-n" "8" "$EXAMPLE_FOLDER/cifio" "--M" "1000" "--N" "1000" "--ker" "1" "--errsol 1" "--format" "2" "--near_para" "4.0";; # test CPP interface for approximate inverse of a 1D FIO, using strong admissible H matrix     
17) "$EXAMPLE_FOLDER/krr_seq" "-quant" "--data_dir" "$DATA_FOLDER/KRR_DATA/susy_10Kn" "--dimn" "8" "--ntrain" "10000" "--ntest" "1000" "--sigma" "0.1" "--lambda" "1.0" "-option" "--xyzsort" "2" "--reclr_leaf" "2";; # test HODLR with LR for KRR (non-MPI)
*) printf "${RED} ###YL: Unknown test\n" ;;
esac

