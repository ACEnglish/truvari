"""
Full documentation at https://github.com/spiralgenetics/truvari/wiki

See `help()` of specific functions / objects for details

methods:
:meth:`allele_freq_annos`
:meth:`bed_ranges`
:meth:`cmd_exe`
:meth:`copy_entry`
:meth:`create_pos_haplotype`
:meth:`entry_boundaries`
:meth:`entry_create_haplotype`
:meth:`entry_distance`
:meth:`entry_gt_comp`
:meth:`entry_is_present`
:meth:`entry_pctsim`
:meth:`entry_reciprocal_overlap`
:meth:`entry_same_variant_type`
:meth:`entry_size`
:meth:`entry_size_similarity`
:meth:`entry_to_key`
:meth:`entry_variant_type`
:meth:`fetch_coords`
:meth:`filter_value`
:meth:`get_gt`
:meth:`get_scalebin`
:meth:`get_sizebin`
:meth:`get_svtype`
:meth:`help_unknown_cmd`
:meth:`make_interval_tree`
:meth:`match_sorter`
:meth:`optimize_df_memory`
:meth:`overlaps`
:meth:`reciprocal_overlap`
:meth:`ref_ranges`
:meth:`restricted_float`
:meth:`setup_logging`
:meth:`setup_progressbar`
:meth:`sizesim`
:meth:`vcf_to_df`
:meth:`weighted_score`

objects:
:class:`GT`
:class:`GenomeTree`
:class:`LogFileStderr`
:class:`SV`

data:
:data:`truvari.HEADERMAT`
:data:`truvari.MATCHRESULT`
:data:`truvari.QUALBINS`
:data:`truvari.SVTYTYPE`
:data:`truvari.SZBINMAX`
:data:`truvari.SZBINS`
:data:`truvari.SZBINTYPE`
"""

__version__ = '3.1.0-dev'

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

from truvari.annos.af_calc import (
    allele_freq_annos
)

from truvari.genome_tree import (
    GenomeTree,
    make_interval_tree
)


from truvari.utils import (
    HEADERMAT,
    LogFileStderr,
    MATCHRESULT,
    bed_ranges,
    cmd_exe,
    copy_entry,
    help_unknown_cmd,
    ref_ranges,
    restricted_float,
    setup_logging,
    setup_progressbar,
)

from truvari.comparisons import (
    create_pos_haplotype,
    entry_boundaries,
    entry_create_haplotype,
    entry_distance,
    entry_gt_comp,
    entry_is_present,
    entry_pctsim,
    entry_reciprocal_overlap,
    entry_same_variant_type,
    entry_size,
    entry_size_similarity,
    entry_to_key,
    entry_variant_type,
    fetch_coords,
    filter_value,
    match_sorter,
    overlaps,
    reciprocal_overlap,
    sizesim,
    weighted_score,
)
