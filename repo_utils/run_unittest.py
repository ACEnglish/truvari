"""
One off unittests
"""
import os
import sys
import pysam
from collections import defaultdict
from intervaltree import IntervalTree

# Use the current truvari, not any installed libraries
sys.path.insert(0, os.getcwd())

import truvari

# Assume we're running in truvari root directory

"""
Boundary issues
"""
vcf_fn = "repo_utils/test_files/variants/boundary.vcf.gz"
bed_fn = "repo_utils/test_files/beds/boundary.bed"
region_start = 10
region_end = 20

vcf = pysam.VariantFile(vcf_fn)
for entry in vcf:
    state = entry.info['include'] == 'in'
    assert state == truvari.entry_within(entry, region_start, region_end), f"Bad Boundary {str(entry)}"

# removed
#regions = truvari.RegionVCFIterator(vcf, includebed=bed_fn)
#vcf.reset()
#truv_ans = defaultdict(lambda: False)
#for entry in regions.iterate(vcf):
    #truv_ans[truvari.entry_to_key(entry)] = True
#vcf.reset()
#for entry in vcf:
    #state = entry.info['include'] == 'in'
    #assert state == truv_ans[truvari.entry_to_key(entry)], f"Bad Boundary {str(entry)}"

"""
New Region Filtering
"""
vcf_fn = "repo_utils/test_files/variants/boundary_cpx.vcf.gz"
bed_fn = "repo_utils/test_files/beds/boundary_cpx.bed"

tree = defaultdict(IntervalTree)
with open(bed_fn, 'r') as fh:
    for line in fh:
        data = line.strip().split()
        tree[data[0]].addi(int(data[1]), int(data[2]) + 1)

vcf = pysam.VariantFile(vcf_fn)
for entry in truvari.region_filter(vcf, tree, True, False):
    assert entry.info['include'] == 'in', f"Bad in {str(entry)}"

vcf.reset()
for entry in truvari.region_filter(vcf, tree, False, False):
    assert entry.info['include'] == 'out', f"Bad out {str(entry)}"


vcf = pysam.VariantFile(vcf_fn)
for entry in truvari.region_filter_stream(vcf, tree, True, False):
    assert entry.info['include'] == 'in', f"Bad in {str(entry)}"

vcf.reset()
for entry in truvari.region_filter_stream(vcf, tree, False, False):
    assert entry.info['include'] == 'out', f"Bad out {str(entry)}"

vcf = pysam.VariantFile(vcf_fn)
for entry in truvari.region_filter_fetch(vcf, tree, False):
    assert entry.info['include'] == 'in', f"Bad in {str(entry)}"


"""
Filtering logic
"""
vcf = pysam.VariantFile("repo_utils/test_files/variants/filter.vcf")
matcher = truvari.Matcher()
matcher.params.sizemin = 0
matcher.params.sizefilt = 0
matcher.params.passonly = True
for entry in vcf:
    try:
        assert matcher.filter_call(entry), f"Didn't filter {str(entry)}"
    except ValueError as e:
        assert e.args[0].startswith("Cannot compare multi-allelic"), f"Unknown exception {str(entry)}"
