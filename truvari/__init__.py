"""

Truvari - SV comparison and annotation toolkit

See `help()` of specific functions / objects for details

VariantRecord methods:

:meth:`entry_boundaries`
:meth:`entry_distance`
:meth:`entry_gt_comp`
:meth:`entry_is_filtered`
:meth:`entry_is_present`
:meth:`entry_reciprocal_overlap`
:meth:`entry_same_variant_type`
:meth:`entry_shared_ref_context`
:meth:`entry_seq_similarity`
:meth:`entry_size`
:meth:`entry_size_similarity`
:meth:`entry_to_hash`
:meth:`entry_to_key`
:meth:`entry_within`
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

Dev methods:

:meth:`benchdir_count_entries`
:meth:`chunker`
:meth:`cmd_exe`
:meth:`consolidate_phab_vcfs`
:meth:`coords_within`
:meth:`count_entries`
:meth:`file_zipper`
:meth:`help_unknown_cmd`
:meth:`make_temp_filename`
:meth:`opt_gz_open`
:meth:`optimize_df_memory`
:meth:`performance_metrics`
:meth:`region_filter`
:meth:`restricted_float`
:meth:`restricted_int`
:meth:`setup_logging`
:meth:`vcf_to_df`

Objects:

:class:`Bench`
:class:`BenchOutput`
:class:`GT`
:class:`RegionVCFIterator`
:class:`LogFileStderr`
:class:`MatchResult`
:class:`Matcher`
:class:`StatsBox`
:class:`SV`

Data:

:data:`truvari.HEADERMAT`
:data:`truvari.QUALBINS`
:data:`truvari.SVTYTYPE`
:data:`truvari.SZBINMAX`
:data:`truvari.SZBINS`
:data:`truvari.SZBINTYPE`
"""

__version__ = '4.2.1'


from truvari.annotations.af_calc import (
    allele_freq_annos,
    calc_af,
    calc_hwe
)

from truvari.bench import (
    Bench,
    BenchOutput,
    StatsBox,
)

from truvari.comparisons import (
    coords_within,
    create_pos_haplotype,
    entry_boundaries,
    entry_distance,
    entry_gt_comp,
    entry_is_filtered,
    entry_is_present,
    entry_reciprocal_overlap,
    entry_same_variant_type,
    entry_shared_ref_context,
    entry_seq_similarity,
    entry_size,
    entry_size_similarity,
    entry_to_hash,
    entry_to_key,
    entry_variant_type,
    entry_within,
    overlap_percent,
    overlaps,
    reciprocal_overlap,
    seqsim,
    sizesim,
    unroll_compare,
)

from truvari.matching import (
    MatchResult,
    Matcher,
    chunker,
    file_zipper
)

from truvari.msatovcf import (
    msa2vcf
)

from truvari.phab import (
    phab,
)

from truvari.region_vcf_iter import (
    RegionVCFIterator,
    build_anno_tree,
    region_filter,
)

from truvari.stratify import (
    count_entries,
    benchdir_count_entries,
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
    performance_metrics,
    ref_ranges,
    restricted_float,
    restricted_int,
    setup_logging,
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
