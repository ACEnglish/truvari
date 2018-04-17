import sys
import argparse
import vcf
from collections import defaultdict


MINSIZE = 50
SAMPLE = "HG002"
CHROM = "1"

def make_interval_tree(vcf_file, sizemin=10, passonly=False):
    """
    Return a dictonary of {chr:[start, end,..], ...}
    that can be queried with bisect to return a span to fetch variants against
    could stand to make a progress bar since this now takes so long
    """
    n_entries = 0
    cmp_entries = 0
    lookup = defaultdict(IntervalTree)
    for entry in vcf_file:
        if passonly and len(entry.FILTER):
            continue
        g = entry.genotype(SAMPLE)
        if not g.called:
            continue
        n_entries += 1
        start, end = get_vcf_boundaries(entry)
        if get_vcf_entry_size(entry) < sizemin:
            continue
        cmp_entries += 1
        lookup[entry.CHROM].addi(start, end, entry.start)
    return lookup, n_entries, cmp_entries

def get_vcf_entry_size(entry):
    """
    Calculate the size of the variant. Use SVLEN INFO tag if available. Otherwise inferr
    """
    if "SVLEN" in entry.INFO:
        if type(entry.INFO["SVLEN"]) is list:
            size = abs(entry.INFO["SVLEN"][0])
        else:
            size = abs(entry.INFO["SVLEN"])
    elif str(entry.ALT[0]).count("<"):
        start, end = get_vcf_boundaries(entry)
        size = end - start
    else:
        size = abs(len(entry.REF) - len(str(entry.ALT[0])))
    return size

v = vcf.Reader(filename=sys.argv[1])
used = {}
starts = defaultdict(list)
for entry in v:
    gt = entry.genotype(SAMPLE)
    #if not gt.called:
        #continue
    key = entry.CHROM + "_" + str(entry.start)
    if key in used:
        continue
    used[key] = True
    if entry.CHROM != "1":
        continue
    if entry.ALT[0].type == "BND":
        continue
    sz = get_vcf_entry_size(entry)
    if sz < MINSIZE:
        continue
    
    starts[entry.CHROM].append((entry.start, gt.called))

distances = []
for chrom in starts:
    starts[chrom].sort()
    for pos in range(1, len(starts[chrom]) - 1):
        if starts[chrom][pos][1]:
            mstart = starts[chrom][pos][0]
            dist = min(abs(mstart - starts[chrom][pos-1][0]), 
                       abs(mstart - starts[chrom][pos+1][0]),)
            distances.append(dist)
print "\n".join([str(x) for x in distances])

        
        
       
    
