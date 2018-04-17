import bisect
import numpy as np

bins = np.logspace(np.log10(1),np.log10(10000), 20)
bins =  [int(x) for x in bins] + [10000000000]
counts_alone = {}
for i in bins:
    counts_alone[i] = 0
f2 = "hg002_giab_distances.txt"
f3 = "HG002_distances.txt"
with open(f2, 'r') as fh:
    for line in fh:
        val = int(line.strip())
        if val == 0:
            counts_alone[bins[0]] += 1
            continue
        pos = bisect.bisect_left(bins, val)
        counts_alone[bins[pos]] += 1
counts_merged = {}
for i in bins:
    counts_merged[i] = 0

f2 = "giab_merge_distances.txt"
f3 = "rtg_distances.txt"
with open(f2, 'r') as fh:
    for line in fh:
        val = int(line.strip())
        if val == 0:
            counts_merged[bins[0]] += 1
            continue
        pos = bisect.bisect_left(bins, val)
        counts_merged[bins[pos]] += 1

for b in bins:
    print int(b), counts_alone[b], counts_merged[b]





