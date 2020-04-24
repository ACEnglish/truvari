


# 1 - get repeatmasker executable to run
# 2 - Parse the output
# 3 - Put the relevant information into the VCF
# 4 - Then I can work on finding which parameters should be played with
#     so that there can be better results - e.g. expsoing RMs params, 
#                                           what output lines are considered 'valid'
# 5 - Add parameters for 'simple' output or the full blown annotation. I'm going to default to simple

# Then I can reuse most of this for TRF.. hopefully. Where'
from truvari.annotation import BEDReader

trees = defaultdict(IntervalTree)
for entry in whatever:
    tree[entry.chrom].addi(entry.start, entry.stop, entry)

# Its all about the anno_data that I need to preserve, I reckon

# And then I need a get_header
