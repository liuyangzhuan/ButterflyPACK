#!/bin/bash

#modules:
module purge
module load gcc/9.1.0
module load openmpi/gcc-9.1.0/4.0.1
module load scalapack-netlib/gcc-9.1.0/2.2.0
module load cmake/3.19.2
module unload python
export PATH=/home/administrator/Desktop/Research/GPTune_master/env/bin/:$PATH
export BPACK_PYTHON_LIB_PATH=/home/administrator/Desktop/Research/ButterflyPACK_logdet/build/SRC_DOUBLE/
export PYTHONPATH=$BPACK_PYTHON_LIB_PATH:$PYTHONPATH
# export PYTHONPATH=$PWD/../EXAMPLE/:$PYTHONPATH
python --version




nmpi=1
export OMP_NUM_THREADS=1

# mpirun --allow-run-as-root -n $nmpi valgrind --tool=memcheck --leak-check=full --track-origins=yes python ../EXAMPLE/Test_python_allranks.py -option --xyzsort 1 --tol_comp 1e-6 --lrlevel 0 --reclr_leaf 5 --nmin_leaf 128
# mpirun --allow-run-as-root -n $nmpi python ../EXAMPLE/Test_python_allranks.py -option --xyzsort 1 --tol_comp 1e-6 --lrlevel 0 --reclr_leaf 5 --nmin_leaf 128 --errsol 1
# python ../EXAMPLE/Test_python_allranks.py -option --xyzsort 1 --tol_comp 1e-6 --lrlevel 0 --reclr_leaf 5 --nmin_leaf 128 --errsol 1




## The following sets the file names for the butterflypack file interface
################################################# 
export CONTROL_FILE="control.txt"  ## this file is used to pass flags and parameters between the master driver and butterflypack workers 
export DATA_FILE="data.bin" ## this file is used to pass covariance matrix and rhs from the master driver to butterflypack workers 
export RESULT_FILE="result.bin" ## this file is used to pass solution vector and logdet from butterflypack workers to the master driver 
export MAX_ID_FILE=10 ## this is the maximum number of BPACK instances 
#################################################

############## sequentially call the python driver Test_python_master.py, but parallelly launching the workers dPy_BPACK_worker.py 
for fid in $(seq 0 "$MAX_ID_FILE"); do
    rm -rf "$CONTROL_FILE.$fid" "$DATA_FILE.$fid" "$RESULT_FILE.$fid"
done
mpirun --allow-run-as-root -n $nmpi python -u ${BPACK_PYTHON_LIB_PATH}/dPy_BPACK_worker.py -option --xyzsort 1 --tol_comp 1e-10 --lrlevel 0 --reclr_leaf 5 --nmin_leaf 128 --errsol 1 | tee a.out_seperatelaunch_worker &
python -u ../EXAMPLE/Test_python_master.py | tee a.out_seperatelaunch_master 
# python -u ../EXAMPLE/Test_python_master_2mats.py | tee a.out_seperatelaunch_master 
python -c "from dPy_BPACK_wrapper import *; bpack_terminate()"




