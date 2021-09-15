"""Pulls info fields from a VCF"""
import sys
import pysam

v = pysam.VariantFile(sys.argv[1])
infos = sys.argv[2].split(',')
for entry in v:
    out = []
    for i in infos:
        if i in entry.info:
            out.append(str(entry.info[i]))
        else:
            out.append('.')
    print("\t".join(out))
