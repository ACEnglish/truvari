"""
Truvari - SV comparison and annotation toolkit

See `help()` of specific functions / objects for details

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
:class:`VariantRecord`

Extra methods:

:meth:`allele_freq_annos`
:meth:`bed_ranges`
:meth:`best_seqsim`
:meth:`build_anno_tree`
:meth:`calc_af`
:meth:`calc_hwe`
:meth:`compress_index_vcf`
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
:meth:`roll_compare`
:meth:`seqsim`
:meth:`sizesim`
:meth:`unroll_seqsim`
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

Data:

:data:`truvari.HEADERMAT`
:data:`truvari.QUALBINS`
:data:`truvari.SVTYTYPE`
:data:`truvari.SZBINMAX`
:data:`truvari.SZBINS`
:data:`truvari.SZBINTYPE`
"""

__version__ = '5.0.0'


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
    best_seqsim,
    coords_within,
    overlap_percent,
    overlaps,
    reciprocal_overlap,
    roll_seqsim,
    seqsim,
    sizesim,
    unroll_seqsim,
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
    build_region_tree,
    build_anno_tree,
    merge_region_tree_overlaps,
    extend_region_tree,
    region_filter,
    region_filter_fetch,
    region_filter_stream,
)

from truvari.stratify import (
    count_entries,
    benchdir_count_entries,
)

from truvari.utils import (
    HEADERMAT,
    LogFileStderr,
    bed_ranges,
    check_vcf_index,
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

from truvari.variants import (
    VariantFile,
    VariantRecord,
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
