import sys
import pysam
from truvari.phab import extract_haplotypes

ref_fn = "ref.fa"
vcf_fn = sys.argv[1]
v = pysam.VariantFile(vcf_fn)

for sample in v.header.samples:
    x = extract_haplotypes((vcf_fn, sample, False), ref_fn)
    print(x[0])
    print(x[1])
    

