"""
Truvari - SV comparison and annotation toolkit

See `help()` of specific functions / objects for details

VariantRecord methods:
:meth:`entry_boundaries`
:meth:`entry_create_haplotype`
:meth:`entry_distance`
:meth:`entry_gt_comp`
:meth:`entry_is_filtered`
:meth:`entry_is_present`
:meth:`entry_pctsim`
:meth:`entry_reciprocal_overlap`
:meth:`entry_same_variant_type`
:meth:`entry_size`
:meth:`entry_size_similarity`
:meth:`entry_to_haplotype`
:meth:`entry_to_key`
:meth:`entry_variant_type`

Extra methods:
:meth:`allele_freq_annos`
:meth:`bed_ranges`
:meth:`build_anno_tree`
:meth:`calc_af`
:meth:`calc_hwe`
:meth:`compress_index_vcf`
:meth:`create_pos_haplotype`
:meth:`get_gt`
:meth:`get_scalebin`
:meth:`get_sizebin`
:meth:`get_svtype`
:meth:`msa2vcf`
:meth:`overlap_percent`
:meth:`overlaps`
:meth:`phab`
:meth:`reciprocal_overlap`
:meth:`ref_ranges`
:meth:`seqsim`
:meth:`sizesim`
:meth:`unroll_compare`
:meth:`vcf_ranges`
:meth:`weighted_score`

Dev methods:
:meth:`chunker`
:meth:`cmd_exe`
:meth:`file_zipper`
:meth:`help_unknown_cmd`
:meth:`make_temp_filename`
:meth:`opt_gz_open`
:meth:`optimize_df_memory`
:meth:`restricted_float`
:meth:`restricted_int`
:meth:`setup_logging`
:meth:`setup_progressbar`
:meth:`vcf_to_df`

Objects:
:class:`GT`
:class:`RegionVCFIterator`
:class:`LogFileStderr`
:class:`MatchResult`
:class:`Matcher`
:class:`SV`

Data:
:data:`truvari.HEADERMAT`
:data:`truvari.QUALBINS`
:data:`truvari.SVTYTYPE`
:data:`truvari.SZBINMAX`
:data:`truvari.SZBINS`
:data:`truvari.SZBINTYPE`
"""

__version__ = '3.5.0'


from truvari.annotations.af_calc import (
    allele_freq_annos,
    calc_af,
    calc_hwe
)

from truvari.comparisons import (
    create_pos_haplotype,
    entry_boundaries,
    entry_create_haplotype,
    entry_distance,
    entry_gt_comp,
    entry_is_filtered,
    entry_is_present,
    entry_pctsim,
    entry_reciprocal_overlap,
    entry_same_variant_type,
    entry_size,
    entry_size_similarity,
    entry_to_haplotype,
    entry_to_key,
    entry_variant_type,
    overlap_percent,
    overlaps,
    reciprocal_overlap,
    seqsim,
    sizesim,
    unroll_compare,
    weighted_score,
)

from truvari.matching import (
    MatchResult,
    Matcher,
    chunker,
    file_zipper
)

from truvari.msa2vcf import (
    msa2vcf
)

from truvari.phab import (
    phab
)

from truvari.region_vcf_iter import (
    RegionVCFIterator,
    build_anno_tree
)

from truvari.utils import (
    HEADERMAT,
    LogFileStderr,
    bed_ranges,
    cmd_exe,
    compress_index_vcf,
    help_unknown_cmd,
    make_temp_filename,
    opt_gz_open,
    ref_ranges,
    restricted_float,
    restricted_int,
    setup_logging,
    setup_progressbar,
    vcf_ranges,
)

from truvari.vcf2df import (
    GT,
    QUALBINS,
    SV,
    SVTYTYPE,
    SZBINMAX,
    SZBINS,
    SZBINTYPE,
    get_gt,
    get_scalebin,
    get_sizebin,
    get_svtype,
    optimize_df_memory,
    vcf_to_df,
)
