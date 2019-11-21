"""
The ``truvari` library provides methods for comparing vcf entries and calculating performance metrics for discovery
performance

User Manual:
    https://github.com/spiralgenetics/truvari
"""
from truvari.utils import (
    StatsBox,
    setup_progressbar,
    LogFileStderr,
    setup_logging
)

from truvari.comparisons import (
    entry_is_variant,
    vcf_to_key,
    var_sizesim,
    size_similarity,
    get_vcf_size_similarity,
    gt_comp,
    create_haplotype,
    var_pctsim_lev,
    overlaps,
    get_vcf_variant_type,
    same_variant_type,
    fetch_coords,
    get_vcf_boundaries,
    get_vcf_entry_size,
    get_rec_ovl,
    get_weighted_score,
    reciprocal_overlap,
    get_vcf_reciprocal_overlap,
    is_sv,
    type_match,
    filter_value,
    match_sorter
)

from truvari.genome_tree import *

