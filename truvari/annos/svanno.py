import sys
import pysam
import truvari

v = pysam.VariantFile(sys.argv[1])
header = v.header
header.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SVTYPE">')
header.add_line('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SVLEN">')
o = pysam.VariantFile("/dev/stdout", 'w', header=header)

for entry in v:
    entry = truvari.copy_entry(entry, header)
    svtype = truvari.entry_variant_type(entry)
    entry.info["SVTYPE"] = svtype
    sz = truvari.entry_size(entry)
    entry.info["SVLEN"] = sz
    o.write(entry)

