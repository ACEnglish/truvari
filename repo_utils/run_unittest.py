"""
One off unittests
"""
import os
import sys
from collections import defaultdict
from intervaltree import IntervalTree

# Use the current truvari, not any installed libraries
sys.path.insert(0, os.getcwd())

import truvari

from truvari.region_vcf_iter import region_filter_stream, region_filter_fetch
# Assume we're running in truvari root directory

"""
Boundary issues
"""
vcf_fn = "repo_utils/test_files/variants/boundary.vcf.gz"
bed_fn = "repo_utils/test_files/beds/boundary.bed"
region_start = 10
region_end = 20

vcf = truvari.VariantFile(vcf_fn)
for entry in vcf:
    state = entry.info['include'] == 'in'
    assert state == entry.within(region_start, region_end), f"Bad Boundary {str(entry)}"

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

vcf = truvari.VariantFile(vcf_fn)
for entry in vcf.regions_fetch(tree, True, False):
    assert entry.info['include'] == 'in', f"Bad in {str(entry)}"

vcf.reset()
for entry in vcf.regions_fetch(tree, False, False):
    assert entry.info['include'] == 'out', f"Bad out {str(entry)}"


vcf = truvari.VariantFile(vcf_fn)
for entry in region_filter_stream(vcf, tree, True, False):
    assert entry.info['include'] == 'in', f"Bad in {str(entry)}"

vcf.reset()
for entry in region_filter_stream(vcf, tree, False, False):
    assert entry.info['include'] == 'out', f"Bad out {str(entry)}"

vcf = truvari.VariantFile(vcf_fn)
for entry in region_filter_fetch(vcf, tree, False):
    assert entry.info['include'] == 'in', f"Bad in {str(entry)}"


"""
Filtering logic
"""
p = truvari.VariantParams(sizemin=0, sizefilt=0, passonly=True)
vcf = truvari.VariantFile("repo_utils/test_files/variants/filter.vcf", params=p)
for entry in vcf:
    try:
        assert entry.filter_call(), f"Didn't filter {str(entry)}"
    except ValueError as e:
        assert e.args[0].startswith("Cannot compare multi-allelic"), f"Unknown exception {str(entry)}"
