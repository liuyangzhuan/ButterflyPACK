#!/bin/sh
set -e

export RED="\033[31;1m"
export BLUE="\033[34;1m"
printf "${BLUE} GC; Entered tests file:\n"

export DATA_FOLDER=$TRAVIS_BUILD_DIR/EXAMPLE
export EXAMPLE_FOLDER=$TRAVIS_BUILD_DIR/build/EXAMPLE
export TEST_FOLDER=$TRAVIS_BUILD_DIR/build/TEST

case "${TEST_NUMBER}" in
1)  mpirun "-n" "1" "$EXAMPLE_FOLDER/ctest" ;;
*) printf "${RED} ###GC: Unknown test\n" ;;
esac
