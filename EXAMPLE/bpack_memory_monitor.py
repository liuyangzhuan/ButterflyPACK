import psutil
import time
import sys
import os
import numpy as np
import socket

def get_local_hostname():
    return socket.gethostname().split('.')[0]

def read_pids_from_file(filename, local_hostname, timeout=60, poll_interval=0.01):
    """
    Reads PIDs from file, returns only those on local hostname.
    Format: rank_0:hostname:545760
    """
    
    elapsed=0
    while not os.path.exists(filename):
        if elapsed >= timeout:
            print(f"Error: Timeout waiting for '{filename}'", file=sys.stderr)
            return False
        time.sleep(poll_interval)
        elapsed += poll_interval

    # print(f"Reading PIDs from: {filename}")
    processes = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line and ':' in line:
                    try:
                        parts = line.split(':')
                        if len(parts) == 3:
                            rank_str, hostname, pid_str = parts
                            hostname = hostname.split('.')[0]
                            if hostname == local_hostname:
                                processes.append({
                                    'rank': rank_str,
                                    'hostname': hostname,
                                    'pid': int(pid_str)
                                })
                        else:
                            # Old format - assume local
                            rank_str, pid_str = parts
                            processes.append({
                                'rank': rank_str,
                                'hostname': local_hostname,
                                'pid': int(pid_str)
                            })
                    except (ValueError, IndexError):
                        print(f"Warning: Could not parse line: {line}", file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: PID file not found at '{filename}'", file=sys.stderr)
        sys.exit(1)
    return processes

def write_time_stamps(processes, time_stamps, memory_stamps, nstamps, hostname):
    output_filename = f"bpack_monitor_{hostname}.txt"
    with open(output_filename, 'w') as outfile:
        # Header with rank names
        outfile.write("# time_sec : " + " ".join(p['rank'] for p in processes) + "\n")
        for istamp in range(nstamps):
            outfile.write("%10.1f : " % (time_stamps[istamp]))
            for ipid in range(len(processes)):
                outfile.write("%.6e " % (memory_stamps[istamp, ipid]))
            outfile.write("\n")
    # print(f"Written {nstamps} timestamps to {output_filename}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python bpack_monitor.py <pid_file_path>", file=sys.stderr)
        sys.exit(1)
    
    local_hostname = get_local_hostname()
    # print(f"Local hostname: {local_hostname}")

    time.sleep(5)

    pid_file = sys.argv[1]
    processes = read_pids_from_file(pid_file, local_hostname)
    
    if len(processes) == 0:
        print(f"No PIDs to monitor on this host ({local_hostname})")
        sys.exit(0)

    # print(f"Monitoring {len(processes)} local processes:")
    # for p in processes:
    #     print(f"  {p['rank']}:{p['pid']}")

    has_pss = hasattr(psutil.Process(os.getpid()).memory_full_info(), 'pss')
    if not has_pss:
        print("Warning: PSS not available, using RSS", file=sys.stderr)

    polling_interval = 1 # in seconds
    GB = 1024**3

    nstamps_max = int(48 * 3600 / polling_interval)
    # print(f'nstamps_max: {nstamps_max}')
    time_stamps = np.zeros(nstamps_max)
    memory_stamps = np.zeros((nstamps_max, len(processes)))

    time_beg_sec = time.time()
    istamp = 0

    while True:
        try:
            time_elapsed_sec = time.time() - time_beg_sec
            time_stamps[istamp] = time_elapsed_sec

            npid_active = 0
            for ipid, proc in enumerate(processes):
                pid = proc['pid']
                memory_res_gb = 0.0
                if psutil.pid_exists(pid):
                    try:
                        process = psutil.Process(pid)
                        mem_info = process.memory_full_info() if has_pss else process.memory_info()
                        memory_res_gb = mem_info.rss / GB
                        npid_active += 1
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
                memory_stamps[istamp, ipid] = memory_res_gb

            total_mem = np.sum(memory_stamps[istamp, :])
            # print(f'time_stamp: {istamp:8d} : {time_stamps[istamp]:10.1f}s, '
            #       f'total={total_mem:.3f} GB, active={npid_active}/{len(processes)}')

            if istamp >= nstamps_max - 1:
                print('Exceeding time stamps limit')
                write_time_stamps(processes, time_stamps, memory_stamps, istamp + 1, local_hostname)
                break

            if npid_active == 0 and istamp > 5:
                # print('No active PIDs remaining')
                write_time_stamps(processes, time_stamps, memory_stamps, istamp + 1, local_hostname)
                break

            time.sleep(polling_interval)
            istamp += 1

        except KeyboardInterrupt:
            print('\nKeyboard interrupt')
            write_time_stamps(processes, time_stamps, memory_stamps, istamp, local_hostname)
            sys.exit(0)

if __name__ == "__main__":
    main()

