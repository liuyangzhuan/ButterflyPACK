#!/bin/bash

#modules:
module load PrgEnv-gnu
module load cmake
module load python
module unload matlab

module unload cray-libsci
ulimit -s unlimited



export BPACK_PYTHON_LIB_PATH=$PWD/../build/SRC_DOUBLE/
export PYTHONPATH=$BPACK_PYTHON_LIB_PATH:$PYTHONPATH
# export PYTHONPATH=$PWD/../EXAMPLE/:$PYTHONPATH
python --version




nmpi=4
export OMP_NUM_THREADS=4
THREADS_PER_RANK=`expr $OMP_NUM_THREADS \* 2`	
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
# export OMP_PROC_BIND=true
# export OMP_PLACES=cores
export OMP_MAX_ACTIVE_LEVELS=1
# export OMP_NESTED=FALSE
export CRAY_LIBSCI_NUM_THREADS=1

# srun -n 32 -c 8 --cpu_bind=cores ../build/EXAMPLE/ie3d -quant --data_dir ../EXAMPLE/EM3D_DATA/preprocessor_3dmesh/sphere_9000 --wavelength 1.0 -option --lr_blk_num 1 --use_zfp 2 --elem_extract 2 --tol_comp 1d-4 --tol_rand 1d-4 --less_adapt 1 --errsol 0 --errfillfull 0 --reclr_leaf 5 --baca_batch 16 --lrlevel 0 --bp_cnt_lr 1 --precon 1 --xyzsort 1 --forwardN15flag 0 --nmin_leaf 1000 --near_para 0.01 --verbosity 0 --pat_comp 3 --nbundle 1 --format 1 --knn 10 --sample_para 2 --sample_para_outer 2 2>&1 | tee sphere_9000.out_precon_1_sort_1_mpi_32_format_1_less_adapt_1_tol_comp_1d-4_tol_rand_1d-4_lrlevel_0_nbundle_1_forwardN15flag_0_knn_10_omp4

# # mpirun --allow-run-as-root -n $nmpi valgrind --tool=memcheck --leak-check=full --track-origins=yes python ../EXAMPLE/Test_python_allranks.py -option --xyzsort 1 --tol_comp 1e-6 --lrlevel 0 --reclr_leaf 5 --nmin_leaf 128
srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores   python ../EXAMPLE/Test_python_allranks.py -option --xyzsort 1 --tol_comp 1e-6 --lrlevel 0 --reclr_leaf 5 --nmin_leaf 128 --errsol 1 | tee a.out_allranks




# ## The following sets the file names for the butterflypack file interface
# ################################################# 
# export CONTROL_FILE="control.txt"  ## this file is used to pass flags and parameters between the master driver and butterflypack workers 
# export DATA_FILE="data.bin" ## this file is used to pass covariance matrix and rhs from the master driver to butterflypack workers 
# export RESULT_FILE="result.bin" ## this file is used to pass solution vector and logdet from butterflypack workers to the master driver 
# export MAX_ID_FILE=10 ## this is the maximum number of BPACK instances 
# #################################################

# ############## sequentially call the python driver Test_python_master.py, but parallelly launching the workers dPy_BPACK_worker.py 
# for fid in $(seq 0 "$MAX_ID_FILE"); do
#     rm -rf "$CONTROL_FILE.$fid" "$DATA_FILE.$fid" "$RESULT_FILE.$fid"
# done
# srun -n ${nmpi} -c $THREADS_PER_RANK --cpu_bind=cores  python -u ${BPACK_PYTHON_LIB_PATH}/dPy_BPACK_worker.py -option --xyzsort 1 --tol_comp 1e-10 --lrlevel 0 --reclr_leaf 5 --nmin_leaf 128 --errsol 1 | tee a.out_seperatelaunch_worker &
# python -u ../EXAMPLE/Test_python_master.py | tee a.out_seperatelaunch_master 
# # python -u ../EXAMPLE/Test_python_master_2mats.py | tee a.out_seperatelaunch_master 
# python -c "from dPy_BPACK_wrapper import *; bpack_terminate()"




