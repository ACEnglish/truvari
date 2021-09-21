"""
Full documentation at https://github.com/spiralgenetics/truvari/wiki

See help() of specific functions / objects for details

functions:
allele_freq_annos, cmd_exe, copy_entry, create_pos_haplotype, entry_boundaries,
entry_create_haplotype, entry_gt_comp, entry_is_present, entry_pctsim,
entry_reciprocal_overlap, entry_same_variant_type, entry_size,
entry_size_similarity, entry_to_key, entry_variant_type, fetch_coords,
filter_value, get_gt, get_scalebin, get_sizebin, get_svtype, make_interval_tree,
match_sorter, overlaps, reciprocal_overlap, restricted_float, setup_logging,
setup_progressbar, sizesim, vcf_to_df, weighted_score

dynamic objects:
GT, GenomeTree, LogFileStderr, SV

static objects:
HEADERMAT - regular expression of vcf header INFO/FORMAT fields with groups
MATCHRESULT - named tuple of match files
QUALBINS - 0-100 quality score bin strings (step size 10)
SVTYTYPE - pandas.CategoricalDtype of SV enum
SZBINMAX - integer list of maximum size for size bins
SZBINS - string list of size bins
SZBINTYPE - pandas.CategoricalDtype of SZBINS
"""

__version__ = '3.1.0-dev'

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

from truvari.annos.af_calc import (
    allele_freq_annos
)

from truvari.genome_tree import (
    GenomeTree,
    make_interval_tree
)


from truvari.utils import (
    setup_progressbar,
    LogFileStderr,
    setup_logging,
    cmd_exe,
    HEADERMAT,
    MATCHRESULT,
    restricted_float
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
