import sys
import pysam
import truvari
from itertools import combinations
from collections import defaultdict, Counter

# So this is an important parameter now.
# Essentially, this says how far a call can look for any match
# The thresholds will determine if its a TP or a FP
# CHUNKDIST must be >= REFDIST
CHUNKDIST = 500
SIZEMIN = 50 # Probably mean size_filt here...
# But this brings up an interesting question... I don't super like the sizefilt stuff.
# Actually, filt_call will just need to check if the BASE variant is >= sizeMin
# and COMP variant is >= sizefilt
class VCFChunker():
    def __init__(self, files, distance=500):
        self.files = {}
        for key, fn in files:
            self.files[key] = pysam.VariantFile(fn)

    def __iter__(self):
        last_pos = []
        stack = []
        cur = {}
        for key in self.files:
            cur[key] = next(self.files[key])
            last_pos.append((cur[key].chrom, cur[key].start))

        #while True:
            #The most upstream call in the last_pos will be added to stack
           
def file_zipper(*start_files):
    """
    Zip files to yield the entries in order.
    Files must be sorted
    Each parameter is a tuple of ('key', 'fn')
    where key is the identifier (so we know which file the yielded entry came from
    and fn is a VariantFile
    yields key, pysam.VariantRecord 
    """
    next_markers = []
    files = []
    names = []
    for name, i in start_files:
        try:
            next_markers.append(next(i))
            files.append(i)
            names.append(name)
        except StopIteration:
            # For when there are no variants in the file
            pass

    while next_markers:
        sidx = 0
        while not next_markers[sidx]:
            sidx += 1
        if sidx == len(files):
            return None
        for idx, i in enumerate(next_markers):
            if idx == sidx:
                continue
            if i.start < next_markers[sidx].start:
                sidx = idx
        ret = next_markers[sidx]
        key = names[sidx]
        try:
            next_markers[sidx] = next(files[sidx])
        except StopIteration:
            # This file is done
            files.pop(sidx)
            names.pop(sidx)
            next_markers.pop(sidx)
        yield key, ret

def chunker_test():
    """
    So the chunker is going to need to know how to filer the calls
    """
    f1 = pysam.VariantFile("repo_utils/test_files/input1.vcf.gz")
    f2 = pysam.VariantFile("repo_utils/test_files/input2.vcf.gz")
    #f2 = pysam.VariantFile("repo_utils/test_files/giab.vcf.gz")
    #f3 = pysam.VariantFile("repo_utils/test_files/input3.vcf.gz")


    #sys.stdout.write(str(f1.header))
    # This is the actual chunker. The above method is the zipper
    prev_pos = None
    """
    breakdown = defaultdict(list)
    for key, call in chunk:
        if key != 'file3': 
            breakdown[key].append(chunk)
    if len(breakdown) == 1:
        blank_chunks += 1
    """

    cur_chunk = defaultdict(list)
    for key, entry in file_zipper(('comp', f1), ('base', f2)):#, ('file3', f3)):
        if prev_pos and prev_pos + CHUNKDIST * 2 < entry.start:
            #cnt = Counter([_[0] for _ in cur_chunk])
            #print(len(cur_chunk), cnt)
            yield cur_chunk
            #for k, e in cur_chunk:
            #    sys.stdout.write(k + ' ' + str(e))
            #print('would chunk')
            cur_chunk = defaultdict(list)
            prev_pos = None # Gotta encounter a new SV to get a new prev_pos
        # truvari.bench.filter_call
        # which will also need to check that the entry is within includebed
        if truvari.entry_size(entry) >= SIZEMIN:
            prev_pos = entry.stop
            cur_chunk[key].append(entry)
    return cur_chunk
        
def runner_test():
    """
    """
    total_chunks = 0
    blank_chunks = 0
    total_pair_checks = 0
    for chunk in chunker_test():
        total_chunks += 1
        if len(chunk.keys()) == 1:
            blank_chunks += 1
        for i, j in combinations(chunk.keys(), 2):
            # here is where I would be making all of the pairs
            # combinations is for multiple files. Just pairs of files, I don't need to worry about it
            # truvari.bench.match_calls(a and b) 
            # MATCHRESULT needs to have an 'over threshold'
            sz1 = len(chunk[i]) 
            sz2 = len(chunk[j])
            woulddo = sz1 * sz2
            total_pair_checks += woulddo
            print('chunk %d - would be comparing %s (%d) -> %s (%d) = %d' % (total_chunks, i, sz1, j, sz2, woulddo))
        # And here's where we'd need to sort the match results
        # and do the logic to choose the match results
    print("%d chunks would be doing %s comparisons with %d blank chunks" % (total_chunks, total_pair_checks, blank_chunks))
            
if __name__ == '__main__':
    # parallelization
    # I think I got this covered.. like, the chunker is what's imapped from
    # and the runner will just need to be on a per-chunk basis instead of actually pulling
    # from the chunker like it currently does
    runner_test()
