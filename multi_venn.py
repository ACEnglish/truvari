import sys
import vcf
from collections import defaultdict
import itertools


"""
Over multiple vcfs, calculate their overlap.

Entries will match if a key of:
    CHROM:POS ID REF ALT matches. Doesn't do anything fancy with Alleles
"""
allVCFs = sys.argv[1:]

class hash_list(list):
    def __hash__(self):
        return hash(" ".join(self))
        
def entry_key(entry):
    key = "%s:%d %s %s %s" % (entry.CHROM, entry.POS, entry.ID, entry.REF, str(entry.ALT))
    return key

call_lookup = defaultdict(list)
file_abscnt = defaultdict(float)
for vcfn in allVCFs:
    v = vcf.Reader(filename=vcfn)
    for entry in v:
        call_lookup[entry_key(entry)].append(vcfn)
        file_abscnt[vcfn] += 1

count_lookup = {}
combo_gen = [x for l in range(1, len(allVCFs)+1) for x in itertools.combinations(allVCFs, l)]
for files_combo in combo_gen:
    files_combo = hash_list(files_combo)
    files_combo.sort()
    count_lookup[files_combo] = 0

for key in call_lookup:
    fkey = hash_list(sorted(call_lookup[key]))
    count_lookup[fkey] += 1

#1 I want to make a key "101010" so that they can be viz'd easier
#2 - I want to sort the count_lookup by their value so that we output them in order
print "#Files"
for i in allVCFs:
    print "#" + i

print "Group\tTotal\tFileName=pct_of_files_calls"
for combo, value in sorted(count_lookup.iteritems(), key=lambda (k,v): (v,k), reverse=True):
    #Key to id
    if value == 0:
        continue
    mid = ["0"]*len(allVCFs)
    for j in combo:
        mid[allVCFs.index(j)] = "1"
    print "%s\t%d" % ("".join(mid), value),
    for fkey in combo:
        if file_abscnt[fkey] > 0:
            print "%s=%.2f%%" % (fkey, count_lookup[combo] / file_abscnt[fkey] * 100),
        else:
            print "%s=0%%" % (fkey),
    print
    
