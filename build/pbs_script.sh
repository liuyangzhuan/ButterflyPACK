#!/bin/sh
#PBS -S /bin/sh
#PBS -N HODLR_3d
#PBS -A vikramg_fluxoe
#PBS -q fluxoe
#PBS -l nodes=1:ppn=16:vikramg,pmem=3900mb,walltime=400:00:00,qos=flux
#PBS -M liuyangz@umich.edu
#PBS -m abe
#PBS -j oe
#PBS -V
# Let PBS handle your output
cd /home/liuyangz/my_research/DirectSolver/Han/DIRECT_SOLVER_3D_Bplus_Inverse_Scaling_sphere_04_22_17/160000_coarse_mkl2017_rerun_heaparray_CFIE
# Use mpirun to run with 7 nodes for 1 hour
export NUM_OMP_THREADS=8
mpirun -np 1 --map-by ppr:1:socket --bind-to none MLMDA_DIRECT_SOLVER_3D_EFIE_NEW | tee a.out
# export NUM_OMP_THREADS=1
# mpirun -np 1 --map-by ppr:1:node --bind-to none MLMDA_DIRECT_SOLVER_2D_EFIE_NEW | tee a.out
