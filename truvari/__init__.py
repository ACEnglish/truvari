"""
Truvari - SV comparison and annotation toolkit

See `help()` of specific functions / objects for details.

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
# pylint: disable=undefined-all-variable

__version__ = '5.4.0'

__all__ = [
    # Classes
    "Bench", "BenchOutput", "StatsBox", "MatchResult", "GT", "LogFileStderr",
    "SV", "VariantFile", "VariantParams", "VariantRecord",

    # Methods
    "bed_ranges", "benchdir_count_entries", "best_seqsim", "build_region_tree", "check_vcf_index",
    "chunker", "cmd_exe", "compress_index_vcf", "coords_within", "count_entries", "extend_region_tree",
    "file_zipper", "help_unknown_cmd", "get_gt", "get_scalebin", "get_sizebin", "get_svtype",
    "make_temp_filename", "merge_region_tree_overlaps", "msa2vcf", "opt_gz_open", "optimize_df_memory",
    "overlap_percent", "overlaps", "performance_metrics", "read_bed_tree", "reciprocal_overlap",
    "restricted_float", "restricted_int", "ref_ranges", "roll_seqsim", "seqsim", "setup_logging",
    "sizesim", "unroll_seqsim", "vcf_ranges", "vcf_to_df",

    # Data
    "HEADERMAT", "QUALBINS", "SVTYTYPE", "SZBINMAX", "SZBINS", "SZBINTYPE"
]

import importlib

_lazy_map = {
    # Classes
    "Bench": ("truvari.bench", "Bench"),
    "BenchOutput": ("truvari.bench", "BenchOutput"),
    "StatsBox": ("truvari.bench", "StatsBox"),
    "MatchResult": ("truvari.matching", "MatchResult"),
    "GT": ("truvari.vcf2df", "GT"),
    "LogFileStderr": ("truvari.utils", "LogFileStderr"),
    "SV": ("truvari.vcf2df", "SV"),
    "VariantFile": ("truvari.variant_file", "VariantFile"),
    "VariantParams": ("truvari.variant_params", "VariantParams"),
    "VariantRecord": ("truvari.variant_record", "VariantRecord"),

    # Functions (grouped by modules)
    # comparisons
    "best_seqsim": ("truvari.comparisons", "best_seqsim"),
    "coords_within": ("truvari.comparisons", "coords_within"),
    "overlap_percent": ("truvari.comparisons", "overlap_percent"),
    "overlaps": ("truvari.comparisons", "overlaps"),
    "reciprocal_overlap": ("truvari.comparisons", "reciprocal_overlap"),
    "roll_seqsim": ("truvari.comparisons", "roll_seqsim"),
    "seqsim": ("truvari.comparisons", "seqsim"),
    "sizesim": ("truvari.comparisons", "sizesim"),
    "unroll_seqsim": ("truvari.comparisons", "unroll_seqsim"),

    # matching
    "chunker": ("truvari.matching", "chunker"),
    "file_zipper": ("truvari.matching", "file_zipper"),

    # msatovcf
    "msa2vcf": ("truvari.msatovcf", "msa2vcf"),

    # region_vcf_iter
    "build_region_tree": ("truvari.region_vcf_iter", "build_region_tree"),
    "read_bed_tree": ("truvari.region_vcf_iter", "read_bed_tree"),
    "merge_region_tree_overlaps": ("truvari.region_vcf_iter", "merge_region_tree_overlaps"),
    "extend_region_tree": ("truvari.region_vcf_iter", "extend_region_tree"),

    # stratify
    "count_entries": ("truvari.stratify", "count_entries"),
    "benchdir_count_entries": ("truvari.stratify", "benchdir_count_entries"),

    # utils
    "HEADERMAT": ("truvari.utils", "HEADERMAT"),
    "bed_ranges": ("truvari.utils", "bed_ranges"),
    "check_vcf_index": ("truvari.utils", "check_vcf_index"),
    "cmd_exe": ("truvari.utils", "cmd_exe"),
    "compress_index_vcf": ("truvari.utils", "compress_index_vcf"),
    "help_unknown_cmd": ("truvari.utils", "help_unknown_cmd"),
    "make_temp_filename": ("truvari.utils", "make_temp_filename"),
    "opt_gz_open": ("truvari.utils", "opt_gz_open"),
    "performance_metrics": ("truvari.utils", "performance_metrics"),
    "ref_ranges": ("truvari.utils", "ref_ranges"),
    "restricted_float": ("truvari.utils", "restricted_float"),
    "restricted_int": ("truvari.utils", "restricted_int"),
    "setup_logging": ("truvari.utils", "setup_logging"),
    "vcf_ranges": ("truvari.utils", "vcf_ranges"),

    # vcf2df
    "QUALBINS": ("truvari.vcf2df", "QUALBINS"),
    "SVTYTYPE": ("truvari.vcf2df", "SVTYTYPE"),
    "SZBINMAX": ("truvari.vcf2df", "SZBINMAX"),
    "SZBINS": ("truvari.vcf2df", "SZBINS"),
    "SZBINTYPE": ("truvari.vcf2df", "SZBINTYPE"),
    "get_gt": ("truvari.vcf2df", "get_gt"),
    "get_scalebin": ("truvari.vcf2df", "get_scalebin"),
    "get_sizebin": ("truvari.vcf2df", "get_sizebin"),
    "get_svtype": ("truvari.vcf2df", "get_svtype"),
    "optimize_df_memory": ("truvari.vcf2df", "optimize_df_memory"),
    "vcf_to_df": ("truvari.vcf2df", "vcf_to_df"),
}

def __getattr__(name):
    if name in _lazy_map:
        module_name, attr_name = _lazy_map[name]
        module = importlib.import_module(module_name)
        value = getattr(module, attr_name)
        globals()[name] = value  # Cache to avoid repeat lookups
        return value
    raise AttributeError(f"module 'truvari' has no attribute '{name}'")
