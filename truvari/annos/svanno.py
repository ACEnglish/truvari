import sys
import pysam
import truvari

v = pysam.VariantFile(sys.argv[1])
header = v.header
header.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SVTYPE">')
o = pysam.VariantFile("/dev/stdout", 'w', header=header)

for entry in v:
    entry = truvari.copy_entry(entry, header)
    svtype = truvari.entry_variant_type(entry)
    entry.info["SVTYPE"] = svtype
    o.write(entry)

