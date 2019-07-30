#!/bin/sh
set -e

export RED="\033[31;1m"
export BLUE="\033[34;1m"
printf "${BLUE} GC; Entered tests file:\n"

export DATA_FOLDER=$TRAVIS_BUILD_DIR/EXAMPLE
export EXAMPLE_FOLDER=$TRAVIS_BUILD_DIR/build/EXAMPLE
# export TEST_FOLDER=$TRAVIS_BUILD_DIR/build/TEST

# export PATH="$PATH:$TRAVIS_BUILD_DIR/installDir/openmpi-3.1.4/install/bin"
# export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$TRAVIS_BUILD_DIR/installDir/openmpi-3.1.4/install/lib/" 

case "${TEST_NUMBER}" in
1) cd $DATA_FOLDER/FULLMAT_DATA
   sh file_merge.sh Smatrix.mat
   cd $TRAVIS_BUILD_DIR/build
   mpirun --allow-run-as-root "-n" "8" "$EXAMPLE_FOLDER/smat" "-quant" "--data_dir" "$DATA_FOLDER/FULLMAT_DATA" "--explicitflag" "0" "-option" "--lrlevel" "0" ;;# test sampling-based construction of HODLR with LR for an EM scattering matrix stored in file
2) mpirun --allow-run-as-root "-n" "2" "$EXAMPLE_FOLDER/krr" "-quant" "--data_dir" "$DATA_FOLDER/KRR_DATA/susy_10Kn" "--dimn" "8" "--ntrain" "10000" "--ntest" "1000" "--sigma" "0.1" "--lambda" "1.0" "-option" "--xyzsort" "2" "--reclr_leaf" "2";; # test HODLR with LR for KRR
3) mpirun --allow-run-as-root "-n" "4" "$EXAMPLE_FOLDER/ie2d" "-quant" "--model2d" "10" "--nunk" "5000" "--wavelength" "0.08" "-option" "--lr_blk_num" "1" "--tol_comp" "1d-4" "--errfillfull" "0" "--reclr_leaf" "4" "--baca_batch" "16" "--lrlevel" "0" "--precon" "1" "--xyzsort" "2" "--nmin_leaf" "200" "--near_para" "0.01d0" "--pat_comp" "3" "--schulzlevel" "100" "--nbundle" "1" "--format" "1"
   mpirun --allow-run-as-root "-n" "4" "$EXAMPLE_FOLDER/ie2d" "-quant" "--model2d" "10" "--nunk" "5000" "--wavelength" "0.08" "-option" "--lr_blk_num" "1" "--tol_comp" "1d-4" "--errfillfull" "0" "--reclr_leaf" "4" "--baca_batch" "16" "--lrlevel" "100" "--precon" "1" "--xyzsort" "2" "--nmin_leaf" "200" "--near_para" "0.01d0" "--pat_comp" "3" "--schulzlevel" "100" "--nbundle" "1" "--format" "1";;   # test HODLR with BF and LR for 2d IE
4) mpirun --allow-run-as-root "-n" "4" "$EXAMPLE_FOLDER/ie2d" "-quant" "--model2d" "10" "--nunk" "5000" "--wavelength" "0.08" "-option" "--lr_blk_num" "1" "--tol_comp" "1d-4" "--errfillfull" "0" "--reclr_leaf" "4" "--baca_batch" "16" "--lrlevel" "0" "--precon" "1" "--xyzsort" "2" "--nmin_leaf" "200" "--near_para" "2.01d0" "--pat_comp" "3" "--schulzlevel" "100" "--nbundle" "1" "--format" "2"
   mpirun --allow-run-as-root "-n" "4" "$EXAMPLE_FOLDER/ie2d" "-quant" "--model2d" "10" "--nunk" "5000" "--wavelength" "0.08" "-option" "--lr_blk_num" "1" "--tol_comp" "1d-4" "--errfillfull" "0" "--reclr_leaf" "4" "--baca_batch" "16" "--lrlevel" "100" "--precon" "1" "--xyzsort" "2" "--nmin_leaf" "200" "--near_para" "2.01d0" "--pat_comp" "3" "--schulzlevel" "100" "--nbundle" "1" "--format" "2";;   # test HMAT with BF and LR for 2d IE
5) mpirun --allow-run-as-root "-n" "8" "$EXAMPLE_FOLDER/ie3d" "-quant" "--data_dir" "$DATA_FOLDER/EM3D_DATA/sphere_2300" "--wavelength" "2.0" "-option" "--lr_blk_num" "1" "--tol_comp" "1e-2" "--errfillfull" "0" "--reclr_leaf" "4" "--baca_batch" "1" "--lrlevel" "0" "--precon" "1" "--xyzsort" "2" "--nmin_leaf" "100" "--near_para" "0.01d0" "--pat_comp" "3" "--format" "1"
   mpirun --allow-run-as-root "-n" "8" "$EXAMPLE_FOLDER/ie3d" "-quant" "--data_dir" "$DATA_FOLDER/EM3D_DATA/sphere_2300" "--wavelength" "2.0" "-option" "--lr_blk_num" "1" "--tol_comp" "1e-2" "--errfillfull" "0" "--reclr_leaf" "4" "--baca_batch" "1" "--lrlevel" "100" "--precon" "1" "--xyzsort" "2" "--nmin_leaf" "100" "--near_para" "0.01d0" "--pat_comp" "3" "--format" "1";; # test HODLR with BF and LR for 3d IE
6) mpirun --allow-run-as-root "-n" "8" "$EXAMPLE_FOLDER/ie3d" "-quant" "--data_dir" "$DATA_FOLDER/EM3D_DATA/sphere_2300" "--wavelength" "2.0" "-option" "--lr_blk_num" "1" "--tol_comp" "1e-2" "--errfillfull" "0" "--reclr_leaf" "4" "--baca_batch" "1" "--lrlevel" "0" "--precon" "1" "--xyzsort" "2" "--nmin_leaf" "100" "--near_para" "2.01d0" "--pat_comp" "3" "--format" "2"
   mpirun --allow-run-as-root "-n" "8" "$EXAMPLE_FOLDER/ie3d" "-quant" "--data_dir" "$DATA_FOLDER/EM3D_DATA/sphere_2300" "--wavelength" "2.0" "-option" "--lr_blk_num" "1" "--tol_comp" "1e-2" "--errfillfull" "0" "--reclr_leaf" "4" "--baca_batch" "1" "--lrlevel" "100" "--precon" "1" "--xyzsort" "2" "--nmin_leaf" "100" "--near_para" "2.01d0" "--pat_comp" "3" "--format" "2";; # test HMAT with BF and LR for 3d IE   
7) mpirun --allow-run-as-root "-n" "2" "$EXAMPLE_FOLDER/ctest" "1" "$DATA_FOLDER/KRR_DATA/susy_10Kn" "8" "1" "0.1" "10.0" "200" "1e-2" "5" "0" "16" "10";; # test CPP interface with HODLR-LR for KRR with SUSY dataset
8) mpirun --allow-run-as-root "-n" "7" "$EXAMPLE_FOLDER/ctest" "2" "1000" "8" "1" "0.1" "10.0" "100" "1e-4" "5" "0" "16" "10";; # test CPP interface with HODLR-LR for KRR with random points
9) mpirun --allow-run-as-root "-n" "3" "$EXAMPLE_FOLDER/ctest" "3" "5000" "20" "100" "1e-4" "2" "0" "100";;  # test CPP interface with HODLR-LR for product of two random matrices 
10) mpirun --allow-run-as-root "-n" "5" "$EXAMPLE_FOLDER/full";; # test HODLR with LR for reading a KRR matrix from a file 
11) mpirun --allow-run-as-root "-n" "2" "$EXAMPLE_FOLDER/frontal" "-quant" "--data_dir" "$DATA_FOLDER/FULLMAT_DATA/Frontal_poisson_50.mat" "--nunk" "50" "--explicitflag" "0" "-option" "--lrlevel" "0" "--tol_comp" "1d-1" "--nmin_leaf" "100" # test sampling-based construction of HODLR with LR for an 3D poisson frontal matrix stored in file
   mpirun --allow-run-as-root "-n" "2" "$EXAMPLE_FOLDER/frontal" "-quant" "--data_dir" "$DATA_FOLDER/FULLMAT_DATA/Frontal_poisson_50.mat" "--nunk" "50" "--explicitflag" "0" "-option" "--lrlevel" "100" "--tol_comp" "1d-1" "--nmin_leaf" "100" # test sampling-based construction of HODLR with BF for an 3D poisson frontal matrix stored in file
   mpirun --allow-run-as-root "-n" "2" "$EXAMPLE_FOLDER/frontal" "-quant" "--data_dir" "$DATA_FOLDER/FULLMAT_DATA/Frontal_poisson_50.mat" "--nunk" "50" "--explicitflag" "1" "-option" "--lrlevel" "100" "--tol_comp" "1d-1" "--nmin_leaf" "100";; # test entry-evaluation-based construction of HODLR with BF for an 3D poisson frontal matrix stored in file   
12) mpirun --allow-run-as-root "-n" "4" "$EXAMPLE_FOLDER/ie2deigen" "-quant" "--model2d" "10" "--nunk" "5000" "--wavelength" "0.08" "--cmmode" "0" "--si" "0" "-option" "--lr_blk_num" "1" "--tol_comp" "1d-4" "--errfillfull" "0" "--reclr_leaf" "4" "--baca_batch" "16" "--lrlevel" "0" "--precon" "1" "--xyzsort" "2" "--nmin_leaf" "200" "--near_para" "2.01d0" "--pat_comp" "3" "--schulzlevel" "100" "--nbundle" "1" "--format" "1";; # test an eigen solver for HODLR with LR for 2d IE 
13) mpirun --allow-run-as-root "-n" "8" "$EXAMPLE_FOLDER/ie3deigen" "-quant" "--data_dir" "$DATA_FOLDER/EM3D_DATA/sphere_2300" "--wavelength" "2.0" "--si" "1" "--which" "LM" "--nev" "1" "--cmmode" "0" "-option" "--tol_comp" "1e-3" "--lrlevel" "0" "--precon" "1" "--xyzsort" "2" "--nmin_leaf" "100" "--near_para" "0.01d0" "--pat_comp" "3" "--format" "1";; # test an eigen solver for HODLR with LR for 3d IE
14) cd $DATA_FOLDER/FULLMAT_DATA
   sh file_merge.sh Frontal_elastic
   cd $TRAVIS_BUILD_DIR/build
    mpirun --allow-run-as-root "-n" "4" "$EXAMPLE_FOLDER/frontaldist" "-quant" "--nunk" "2500" "--data_dir" "$DATA_FOLDER/FULLMAT_DATA/Frontal_elastic/Frontal_elastic_2500" "--explicitflag" "0" "-option" "--nmin_leaf" "100" "--tol_comp" "1e-3" "--lrlevel" "100";; # test sampling-based construction of HODLR with BF for an 3D elastic Helmholtz frontal matrix stored parallel in file
*) printf "${RED} ###GC: Unknown test\n" ;;
esac

