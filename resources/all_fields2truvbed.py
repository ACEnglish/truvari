#print("""##INFO=<ID=SREP_periods,Number=.,Type=Float,Description="Simple Repeat period lengths">
##INFO=<ID=SREP_copies,Number=.,Type=Float,Description="Simple Repeat copy numbers">
##INFO=<ID=SREP_scores,Number=.,Type=Integer,Description="Simple Repeat scores">
##INFO=<ID=SREP_entropies,Number=.,Type=Float,Description="Simple Repeat entropies">
##INFO=<ID=SREP_repeats,Number=.,Type=String,Description="Simple Repeat sequences">""")

# Given an all_fields output from UCSC genome browser simple repeats, make a anno trf bed
import sys

for line in sys.stdin:
    if line.startswith("#"):
        index = {x:y for y,x in enumerate(line.strip()[1:].split('\t'))}
        continue
    data = line.strip().split('\t')
    print("\t".join([data[index['chrom']],
                     data[index['chromStart']],
                     data[index['chromEnd']],
                     data[index['period']],
                     data[index['copyNum']],
                     data[index['score']],
                     data[index['entropy']],
                     data[index['sequence']]]))
