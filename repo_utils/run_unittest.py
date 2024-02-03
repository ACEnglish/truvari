"""
One off unittests
"""
import os
import sys
import pysam
from collections import defaultdict

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

v = pysam.VariantFile(vcf_fn)
regions = truvari.RegionVCFIterator(v, includebed=bed_fn)

truv_ans = defaultdict(lambda: False)
for entry in regions.iterate(v):
    truv_ans[truvari.entry_to_key(entry)] = True

v = pysam.VariantFile(vcf_fn)
for entry in v:
    state = entry.info['include'] == 'in'
    assert state == truv_ans[truvari.entry_to_key(entry)], f"Bad Boundary {str(entry)}"
