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
    sizesim,
    entry_size_similarity,
    entry_gt_comp,
    create_haplotype,
    entry_pctsim_lev,
    overlaps,
    get_vcf_variant_type,
    same_variant_type,
    fetch_coords,
    entry_boundaries,
    entry_size,
    reciprocal_overlap,
    weighted_score,
    reciprocal_overlap,
    entry_reciprocal_overlap,
    is_sv,
    #type_match,
    filter_value,
    match_sorter
)

from truvari.genome_tree import *
