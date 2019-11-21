
from truvari.annotation import BEDReader

trees = defaultdict(IntervalTree)
for entry in whatever:
    tree[entry.chrom].addi(entry.start, entry.stop, entry)

# Its all about the anno_data that I need to preserve, I reckon

# And then I need a get_header
