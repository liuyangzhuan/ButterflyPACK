import numpy as np
import matplotlib.pyplot as plt
import os

BPACK_Algorithm = 'HODLR'

# Read hostnames from file
with open('hostnames.txt', 'r') as f:
    hostline = f.readline()
hosts = hostline.strip().split()
nHosts = len(hosts)
print(f'Found {nHosts} hosts: {", ".join(hosts)}')

# Read data from each host's monitor file
all_stamps = []
all_nRanks = []

for host in hosts:
    filename = f'bpack_monitor_{host}.txt'
    
    if not os.path.isfile(filename):
        print(f'Warning: {filename} not found, skipping')
        continue
    
    # Read first non-comment line to detect number of columns
    with open(filename, 'r') as f:
        firstline = f.readline()
        if firstline.startswith('#'):
            dataline = f.readline()
        else:
            dataline = firstline
    
    # Count columns by splitting on whitespace (ignore the ':')
    dataline = dataline.replace(':', ' ')
    cols = dataline.strip().split()
    nCols = len(cols)
    nRanksLocal = nCols - 1  # first column is time
    
    print(f'{filename}: {nRanksLocal} ranks')
    
    # Read data, skipping header if present
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    data_lines = []
    for line in lines:
        line = line.strip()
        if line and not line.startswith('#'):
            line = line.replace(':', ' ')
            values = [float(x) for x in line.split()]
            data_lines.append(values)
    
    stamps = np.array(data_lines)
    all_stamps.append(stamps)
    all_nRanks.append(nRanksLocal)

nHosts = len(all_stamps)  # update in case some files were missing
nRanksTotal = sum(all_nRanks)
print(f'Total ranks across all hosts: {nRanksTotal}')

# Combine data - assume all files have same timestamps
# Use first file's timestamps as reference
time_ref = all_stamps[0][:, 0]
nStamps = len(time_ref)

# Combine all rank data into one matrix
combined_data = np.zeros((nStamps, nRanksTotal))
col_offset = 0
for iHost in range(nHosts):
    nRanksLocal = all_nRanks[iHost]
    combined_data[:, col_offset:col_offset + nRanksLocal] = all_stamps[iHost][:, 1:]
    col_offset += nRanksLocal

# Plot
fig, ax = plt.subplots(figsize=(10, 6))

# Plot each rank
for iRank in range(nRanksTotal):
    ax.semilogy(time_ref, combined_data[:, iRank], alpha=0.5)

# Plot sum of all ranks
RSS = np.sum(combined_data, axis=1)
display_name = f'{BPACK_Algorithm} (sum {nRanksTotal} ranks)'
ax.semilogy(time_ref, RSS, 'r--', linewidth=2, label=display_name)

ax.grid(True, which='both', linestyle='--', alpha=0.5)
ax.set_xlabel('time (sec)')
ax.set_ylabel('memory (GB)')
ax.legend(loc='best')
strTitle = f'{BPACK_Algorithm} process memory profiles ({nHosts} nodes, {nRanksTotal} ranks)'
ax.set_title(strTitle)

plt.tight_layout()
plt.savefig('memory_profile.pdf', format='pdf', bbox_inches='tight')
print('Figure saved to memory_profile.pdf')
plt.show()

