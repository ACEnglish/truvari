"""
Truvari

User Manual:
    https://github.com/spiralgenetics/truvari
"""
__version__ = '3.0.0-dev'
from truvari.utils import (
    setup_progressbar,
    LogFileStderr,
    setup_logging,
    cmd_exe,
    HEADERMAT,
    MATCHRESULT
)

from truvari.comparisons import (
    entry_is_present,
    entry_to_key,
    sizesim,
    entry_size_similarity,
    entry_gt_comp,
    entry_create_haplotype,
    create_pos_haplotype,
    entry_pctsim,
    overlaps,
    entry_variant_type,
    entry_same_variant_type,
    fetch_coords,
    entry_boundaries,
    entry_size,
    weighted_score,
    reciprocal_overlap,
    entry_reciprocal_overlap,
    filter_value,
    match_sorter,
    copy_entry
)

from truvari.genome_tree import (
    GenomeTree,
    make_interval_tree
)

from truvari.annos.af_calc import (
    allele_freq_annos
)

from truvari.vcf2df import (
    vcf_to_df,
    SV,
    GT,
    SZBINTYPE,
    SVTYTYPE,
    SZBINS,
    SZBINMAX,
    QUALBINS,
    get_svtype,
    get_sizebin,
    get_gt,
    get_scalebin

)
