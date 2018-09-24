#!/bin/sh
set -e

export RED="\033[31;1m"
export BLUE="\033[34;1m"
printf "${BLUE} GC; Entered tests file:\n"

export DATA_FOLDER=$TRAVIS_BUILD_DIR/EXAMPLE
export EXAMPLE_FOLDER=$TRAVIS_BUILD_DIR/build/EXAMPLE
# export TEST_FOLDER=$TRAVIS_BUILD_DIR/build/TEST

case "${TEST_NUMBER}" in
1) mpirun "-n" "2" "$EXAMPLE_FOLDER/fkerreg" "$DATA_FOLDER/SUSY/susy_10Kn" "8" "10000" "1000" "0.1" "1.0" "2" "2";;
2) mpirun "-n" "2" "$EXAMPLE_FOLDER/fem2d" ;;
*) printf "${RED} ###GC: Unknown test\n" ;;
esac
