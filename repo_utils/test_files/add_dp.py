import pysam
import random
v = pysam.VariantFile("multi.vcf.gz")
n_header = v.header.copy()
n_header.add_line("##FORMAT=<ID=DP,Number=1,Type=Integer,Description='test'>")
n_header.add_line("##FORMAT=<ID=AD,Number=R,Type=Integer,Description='test'>")
o = pysam.VariantFile("/dev/stdout", 'w', header=n_header)
for entry in v:
    entry.translate(n_header)
    for samp in entry.samples:
        entry.samples[samp]["AD"] = (random.randint(5, 25), random.randint(5, 25))
        entry.samples[samp]["DP"] = sum(entry.samples[samp]["AD"])
    o.write(entry)
