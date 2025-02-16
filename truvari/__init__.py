"""
Truvari - SV comparison and annotation toolkit

See `help()` of specific functions / objects for details

Objects:

:class:`Bench`
:class:`BenchOutput`
:class:`GT`
:class:`LogFileStderr`
:class:`MatchResult`
:class:`StatsBox`
:class:`SV`
:class:`VariantFile`
:class:`VariantRecord`
:class:`VariantParams`

Extra methods:

:meth:`bed_ranges`
:meth:`benchdir_count_entries`
:meth:`best_seqsim`
:meth:`build_region_tree`
:meth:`check_vcf_index`
:meth:`chunker`
:meth:`cmd_exe`
:meth:`compress_index_vcf`
:meth:`coords_within`
:meth:`count_entries`
:meth:`extend_region_tree`
:meth:`file_zipper`
:meth:`help_unknown_cmd`
:meth:`get_gt`
:meth:`get_scalebin`
:meth:`get_sizebin`
:meth:`get_svtype`
:meth:`make_temp_filename`
:meth:`merge_region_tree_overlaps`
:meth:`msa2vcf`
:meth:`opt_gz_open`
:meth:`optimize_df_memory`
:meth:`overlap_percent`
:meth:`overlaps`
:meth:`performance_metrics`
:meth:`phab`
:meth:`read_bed_tree`
:meth:`reciprocal_overlap`
:meth:`restricted_float`
:meth:`restricted_int`
:meth:`ref_ranges`
:meth:`roll_seqsim`
:meth:`seqsim`
:meth:`setup_logging`
:meth:`sizesim`
:meth:`unroll_seqsim`
:meth:`vcf_ranges`
:meth:`vcf_to_df`

Data:

:data:`truvari.HEADERMAT`
:data:`truvari.QUALBINS`
:data:`truvari.SVTYTYPE`
:data:`truvari.SZBINMAX`
:data:`truvari.SZBINS`
:data:`truvari.SZBINTYPE`
"""

__version__ = '5.2.0'


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
    read_bed_tree,
    merge_region_tree_overlaps,
    extend_region_tree,
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

from truvari.variant_file import VariantFile
from truvari.variant_params import VariantParams
from truvari.variant_record import VariantRecord

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
