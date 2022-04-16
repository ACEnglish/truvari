import pysam
import random
v = pysam.VariantFile("multi.vcf.gz")
n_header = v.header.copy()
n_header.add_line("##FORMAT=<ID=DP,Number=1,Type=Integer,Description='test'>")
o = pysam.VariantFile("/dev/stdout", 'w', header=n_header)
for entry in v:
    entry.translate(n_header)
    for samp in entry.samples:
        entry.samples[samp]["DP"] = random.randint(0, 30)
    o.write(entry)
