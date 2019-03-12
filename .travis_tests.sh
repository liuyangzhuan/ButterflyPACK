#!/bin/sh
set -e

export RED="\033[31;1m"
export BLUE="\033[34;1m"
printf "${BLUE} GC; Entered tests file:\n"

export DATA_FOLDER=$TRAVIS_BUILD_DIR/EXAMPLE
export EXAMPLE_FOLDER=$TRAVIS_BUILD_DIR/build/EXAMPLE
# export TEST_FOLDER=$TRAVIS_BUILD_DIR/build/TEST

case "${TEST_NUMBER}" in
1) cd $DATA_FOLDER/FULLMAT_DATA
   sh file_merge.sh Smatrix.mat
   cd $TRAVIS_BUILD_DIR/build
   mpirun "-n" "8" "$EXAMPLE_FOLDER/smat" "$DATA_FOLDER/FULLMAT_DATA" "0" # test sampling-based construction of HODLR with LR for an EM scattering matrix stored in file
   mpirun "-n" "8" "$EXAMPLE_FOLDER/smat" "$DATA_FOLDER/FULLMAT_DATA" "100";; # test sampling-based construction of HODLR with BF for an EM scattering matrix stored in file
2) mpirun "-n" "2" "$EXAMPLE_FOLDER/krr" "$DATA_FOLDER/KRR_DATA/susy_10Kn" "8" "10000" "1000" "0.1" "1.0" "2" "2";; # test HODLR with LR for KRR
3) mpirun "-n" "4" "$EXAMPLE_FOLDER/ie2d" "1" "10" "5000" "0.08" "1d-4" "0" "4" "16" "0" "1" "2" "200" "0.01d0" "3" "100" "1" "1"
   mpirun "-n" "4" "$EXAMPLE_FOLDER/ie2d" "1" "10" "5000" "0.08" "1d-4" "0" "4" "16" "100" "1" "2" "200" "0.01d0" "3" "100" "1" "1";;   # test HODLR with BF and LR for 2d IE
4) mpirun "-n" "4" "$EXAMPLE_FOLDER/ie2d" "1" "10" "5000" "0.08" "1d-4" "0" "4" "16" "0" "1" "2" "200" "2.01d0" "3" "100" "1" "2"
   mpirun "-n" "4" "$EXAMPLE_FOLDER/ie2d" "1" "10" "5000" "0.08" "1d-4" "0" "4" "16" "100" "1" "2" "200" "2.01d0" "3" "100" "1" "2";;   # test HMAT with BF and LR for 2d IE
5) mpirun "-n" "8" "$EXAMPLE_FOLDER/ie3d" "$DATA_FOLDER/EM3D_DATA/sphere_2300" "2.0" "1" "1e-2" "0" "4" "1" "0" "1" "2" "100" "0.01d0" "3" "1"
   mpirun "-n" "8" "$EXAMPLE_FOLDER/ie3d" "$DATA_FOLDER/EM3D_DATA/sphere_2300" "2.0" "1" "1e-2" "0" "4" "1" "100" "1" "2" "100" "0.01d0" "3" "1";; # test HODLR with BF and LR for 3d IE
6) mpirun "-n" "8" "$EXAMPLE_FOLDER/ie3d" "$DATA_FOLDER/EM3D_DATA/sphere_2300" "2.0" "1" "1e-2" "0" "4" "1" "0" "1" "2" "100" "2.01d0" "3" "2"
   mpirun "-n" "8" "$EXAMPLE_FOLDER/ie3d" "$DATA_FOLDER/EM3D_DATA/sphere_2300" "2.0" "1" "1e-2" "0" "4" "1" "100" "1" "2" "100" "2.01d0" "3" "2";; # test HMAT with BF and LR for 3d IE   
7) mpirun "-n" "2" "$EXAMPLE_FOLDER/ctest" "1" "$DATA_FOLDER/KRR_DATA/susy_10Kn" "8" "1" "0.1" "10.0" "200" "1e-2" "3" "0" "100";; # test CPP interface with HODLR-LR for KRR with SUSY dataset
8) mpirun "-n" "7" "$EXAMPLE_FOLDER/ctest" "2" "1000" "8" "1" "0.1" "10.0" "100" "1e-4" "2" "0" "100";; # test CPP interface with HODLR-LR for KRR with random points
9) mpirun "-n" "3" "$EXAMPLE_FOLDER/ctest" "3" "5000" "20" "100" "1e-4" "2" "0" "100";;  # test CPP interface with HODLR-LR for product of two random matrices 
10) mpirun "-n" "5" "$EXAMPLE_FOLDER/full";; # test HODLR with LR for reading a KRR matrix from a file 
11) mpirun "-n" "2" "$EXAMPLE_FOLDER/frontal" "50" "$DATA_FOLDER/FULLMAT_DATA/Frontal_poisson_50.mat" "100" "1e-1" "0" # test sampling-based construction of HODLR with LR for an 3D poisson frontal matrix stored in file
   mpirun "-n" "2" "$EXAMPLE_FOLDER/frontal" "50" "$DATA_FOLDER/FULLMAT_DATA/Frontal_poisson_50.mat" "0" "1e-1" "0" # test sampling-based construction of HODLR with BF for an 3D poisson frontal matrix stored in file
   mpirun "-n" "2" "$EXAMPLE_FOLDER/frontal" "50" "$DATA_FOLDER/FULLMAT_DATA/Frontal_poisson_50.mat" "0" "1e-1" "1";; # test entry-evaluation-based construction of HODLR with BF for an 3D poisson frontal matrix stored in file   
12) mpirun "-n" "4" "$EXAMPLE_FOLDER/ie2deigen" "1" "10" "5000" "0.08" "1d-4" "0" "4" "16" "0" "1" "2" "200" "0.01d0" "3" "100" "1" "1" "0" "0"
*) printf "${RED} ###GC: Unknown test\n" ;;
esac

