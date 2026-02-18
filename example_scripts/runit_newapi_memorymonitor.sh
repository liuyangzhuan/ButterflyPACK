#!/bin/bash -l

ulimit -s unlimited
ulimit -c unlimited


export OMP_NUM_THREADS=1
nmpi=2


python ../EXAMPLE/bpack_memory_monitor.py ./pid_list.txt &

logfile=a.out
mpirun --allow-run-as-root -n $nmpi /usr/bin/time -f "Time=%E MaxRSS=%M KB" ./EXAMPLE/ctest_simple_newapi 2>&1 | tee $logfile 

echo $(cut -d':' -f2 pid_list.txt | sort -u) > hostnames.txt
python ../EXAMPLE/bpack_memory_plotter.py

# Extract all numeric values after "MaxRSS="
values=$(grep -o 'MaxRSS=[0-9]*' "$logfile" | cut -d'=' -f2)
if [[ -z "$values" ]]; then
    echo "No MaxRSS entries found in $logfile" | tee -a "$logfile"
    exit 0
fi
# Sum all MaxRSS values (in KB)
total_kb=0
for v in $values; do
    total_kb=$((total_kb + v))
done
total_gb=$(echo "scale=6; $total_kb / 1048576" | bc -l)
printf "Total RSS: %.4f GB\n" "$total_gb" | tee -a "$logfile"
