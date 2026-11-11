#!/bin/bash -l
#SBATCH -A m2957
#SBATCH -q regular
#SBATCH -N 1
#SBATCH --constraint=cpu
#SBATCH -t 02:00:00
#SBATCH -J acoustic_cheb

# Sound-soft scattering from a sphere with the acoustic Chebyshev driver (ie3dacousticcheb) on one Perlmutter CPU node, solved
# by unpreconditioned GMRES with three compressed forward operators:
#   tbf  tensor H-BF: a forest of H-TENSOR blocks, one per pair of base patches on their parameter grids, butterflies at all
#        admissible levels (-quant --tensor 1, -option --lrlevel 100)
#   tlr  tensor LR: the same forest with low-rank (Tucker) admissible blocks (-quant --tensor 1, -option --lrlevel 0)
#   mbf  matrix H-BF: a forest of H-matrices with butterfly admissible blocks, one tree per base patch (-quant --forest 1,
#        -option --format 2)
# The sphere has 6*S^2 patches of N x N Chebyshev nodes (N=8, Nunk = 384*S^2); wavenum 5.79*S keeps about 6 points per
# wavelength. Run it from example_scripts, either with sbatch or inside an allocation:
#   salloc -A m2957 -C cpu -q interactive -N 1 -t 2:00:00
#   bash runit_acoustic_cheb_perlmutter.sh
# Settings (environment variables): FORMATS (default "tbf tlr mbf"), SPLITS ("8 16 32"), LAYOUTS ("8x16", MPI ranks x OpenMP
# threads per rank; tensor H-BF products are fastest with 4x32), BUILD_DIR (build), OUT_DIR ($PSCRATCH/acoustic_cheb).
# Each run writes OUT_DIR/<format>_S<S>_np<ranks>_t<threads>.log; a summary of all logs in OUT_DIR follows the runs
# ("bash runit_acoustic_cheb_perlmutter.sh summary" prints only the summary). Tensor LR at S=32 spends about 17 minutes in
# construction.

module load PrgEnv-gnu
module load cray-fftw
module unload cray-libsci
ulimit -s unlimited

FORMATS=${FORMATS:-"tbf tlr mbf"}
SPLITS=${SPLITS:-"8 16 32"}
LAYOUTS=${LAYOUTS:-"8x16"}
BUILD_DIR=${BUILD_DIR:-build}
OUT_DIR=${OUT_DIR:-$PSCRATCH/acoustic_cheb}
EXE=$(cd .. && pwd)/$BUILD_DIR/EXAMPLE/ie3dacousticcheb
if [ ! -x $EXE ]; then
	echo "$EXE not found: build it with run_cmake_build_gnu_perlmutter_acoustic_cheb.sh, and run this script from example_scripts"
	exit 1
fi
mkdir -p $OUT_DIR

# the library links the sequential OpenBLAS; the OpenMP OpenBLAS (about 4x faster here) and tcmalloc are preloaded into the
# application only, not into srun
TCMALLOC_ROOT=/global/cfs/cdirs/m2957/lib/lib/PrgEnv-gnu/gperftools-2.18.1
OMP_OPENBLAS=/global/cfs/cdirs/m2957/lib/lib/PrgEnv-gnu/OpenBLAS/build/install/lib/libopenblas.so.0
ZFP_LIB=/global/cfs/cdirs/m2957/liuyangz/my_research/zfp-1.0.0_gcc_perlmutter/install/lib64
export LD_LIBRARY_PATH="$TCMALLOC_ROOT/lib:$ZFP_LIB:${LD_LIBRARY_PATH:-}"
export TCMALLOC_RELEASE_RATE=1
export OMP_PLACES=cores OMP_PROC_BIND=spread OMP_MAX_ACTIVE_LEVELS=1

# run <format> <S> <ranks> <threads per rank>
run() {
	local fmt=$1 S=$2 nr=$3 nt=$4
	local k=$(awk "BEGIN {printf \"%.3f\", 5.79*$S}")
	local quant="-quant --geo sphere --split $S --N 8 --wavenum $k --onthefly 1 --gmres_restart 200 --eigtest 0"
	local opt="-option --tol_comp 1d-4 --tol_itersol 1d-6 --n_iter 600 --verbosity 1"
	case $fmt in
	tbf) quant="$quant --tensor 1"; opt="$opt --sample_para 0.5 --sample_para_outer 0.5 --htensor_mvp_level_batch 100 --lrlevel 100" ;;
	tlr) quant="$quant --tensor 1"; opt="$opt --sample_para 0.5 --sample_para_outer 0.5 --htensor_mvp_level_batch 100 --lrlevel 0" ;;
	mbf) quant="$quant --forest 1"; opt="$opt --format 2 --precon 2 --knn 30 --forwardN15flag 0 --less_adapt 1 --lrlevel 100" ;;
	*) echo "unknown format $fmt"; return ;;
	esac
	local log=$OUT_DIR/${fmt}_S${S}_np${nr}_t${nt}.log
	echo "$(date +%T) $fmt S=$S ranks=$nr threads=$nt -> $log"
	(cd $OUT_DIR && OMP_NUM_THREADS=$nt OPENBLAS_NUM_THREADS=$nt srun -n $nr -c $((2*nt)) --cpu-bind=cores \
		env LD_PRELOAD="$TCMALLOC_ROOT/lib/libtcmalloc.so:$OMP_OPENBLAS" $EXE $quant $opt > $log 2>&1)
}

if [ "$1" != "summary" ]; then
	if [ -z "$SLURM_JOB_ID" ]; then
		echo "no Slurm allocation: submit with sbatch or run inside salloc (see the header)"
		exit 1
	fi
	for S in $SPLITS; do
		for fmt in $FORMATS; do
			for lay in $LAYOUTS; do
				run $fmt $S ${lay%x*} ${lay#*x}
			done
		done
	done
fi

# summary: construction time, memory of the compressed (butterfly or low-rank) and dense blocks summed over ranks, GMRES
# iterations, time per product (the matrix forest has no product timer: GMRES time / (iterations + 1)) and far-field error
# against the Mie series
echo
printf "%-6s %4s %8s %7s %11s %15s %11s %6s %10s %11s\n" format S Nunk layout "constr (s)" "compressed (MB)" "dense (MB)" iters "s/product" "far error"
for log in $OUT_DIR/*_S*_np*_t*.log; do
	awk -v f=$(basename $log .log) '
	/ Nunk:/ { nunk = $2 }
	/tensor forest construction time:/ { con = $5 }
	/tensor forest memory:/ { mc = $4; md = $8 }
	/tensor forest solve time:/ { pp = $7 }
	/ Total construction time:/ { mcon = $4 }
	/MB costed for butterfly forward blocks/ { mmc = $1 }
	/MB costed for direct forward blocks/ { mmd = $1 }
	/ Solving:/ { sol = $2 }
	/ GMRES\(/ { it = $4 }
	/far-field error vs. Mie series/ { ff = $NF }
	END {
		split(f, a, "_"); fmt = a[1]; S = substr(a[2], 2); lay = substr(a[3], 3) "x" substr(a[4], 2)
		if (fmt == "mbf") { con = mcon; mc = mmc; md = mmd; if (it != "") pp = sol / (it + 1) }
		printf "%-6s %4s %8s %7s %11.1f %15.0f %11.0f %6s %10.4f %11s\n", fmt, S, nunk, lay, con, mc, md, it, pp, ff
	}' $log
done
