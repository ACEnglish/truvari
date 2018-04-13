import sys
import vcf

import bisect
#Go through all the entries in a vcf and bin up the sizes
#Annotations, too? No. That's advanced features


v = vcf.Reader(filename=sys.argv[1])
#dynamic binsizes based on freqs?
bins = [0, 10, 20, 50, 100, 250, 500, 1000, 3000, 7000, 10000, 5e9, \
          -10,-20,-50,-100,-250,-500,-1000,-3000,-7000,-10000,-5e9]
bin_cnts = [0]*len(bins)
bins.sort()

def get_vcf_entry_size(entry):
    """
    Calculate the size of the variant. Use SVLEN INFO tag if available. Otherwise inferr
    """
    if "SVLEN" in entry.INFO:
        if type(entry.INFO["SVLEN"])is list:
            size = entry.INFO["SVLEN"][0]
        else:
            size = entry.INFO["SVLEN"]
    else:
        size = len(entry.REF) - len(str(entry.ALT[0]))
        if size == 0:
            size = len(entry.REF)
    return size

for entry in v:
    size = get_vcf_entry_size(entry)
    pos = bisect.bisect_left(bins, size)
    if size > 0:
        pos -= 1
    bin_cnts[pos] += 1
print "\n".join(["%d\t%d" % (x, y) for x,y in zip(bins, bin_cnts)])

