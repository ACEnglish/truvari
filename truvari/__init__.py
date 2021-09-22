"""
Full documentation at https://github.com/spiralgenetics/truvari/wiki

See `help()` of specific functions / objects for details

methods:
:meth:`allele_freq_annos`
:meth:`cmd_exe`
:meth:`copy_entry`
:meth:`create_pos_haplotype`
:meth:`entry_boundaries`
:meth:`entry_create_haplotype`
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
:meth:`make_interval_tree`
:meth:`match_sorter`
:meth:`overlaps`
:meth:`reciprocal_overlap`
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
    restricted_float,
    copy_entry
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
    match_sorter
)
