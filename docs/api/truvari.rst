truvari package
===============

Overview
--------
.. automodule:: truvari
   :members:
   :undoc-members:
   :show-inheritance:

VariantRecord Methods
---------------------
entry_boundaries
^^^^^^^^^^^^^^^^
.. autofunction:: entry_boundaries

entry_distance
^^^^^^^^^^^^^^
.. autofunction:: entry_distance

entry_gt_comp
^^^^^^^^^^^^^
.. autofunction:: entry_gt_comp

entry_is_filtered
^^^^^^^^^^^^^^^^^
.. autofunction:: entry_is_filtered

entry_is_present
^^^^^^^^^^^^^^^^
.. autofunction:: entry_is_present

entry_reciprocal_overlap
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: entry_reciprocal_overlap

entry_same_variant_type
^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: entry_same_variant_type

entry_seq_similarity
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: entry_seq_similarity

entry_size
^^^^^^^^^^
.. autofunction:: entry_size

entry_size_similarity
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: entry_size_similarity

entry_to_hash
^^^^^^^^^^^^^
.. autofunction:: entry_to_hash

entry_to_key
^^^^^^^^^^^^
.. autofunction:: entry_to_key

entry_variant_type
^^^^^^^^^^^^^^^^^^
.. autofunction:: entry_variant_type

shared_reference_context
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: shared_reference_context

Extra Methods
-------------
allele_freq_annos
^^^^^^^^^^^^^^^^^
.. autofunction:: allele_freq_annos  

bed_ranges
^^^^^^^^^^
.. autofunction:: bed_ranges

build_anno_tree
^^^^^^^^^^^^^^^
.. autofunction:: build_anno_tree

calc_af
^^^^^^^
.. autofunction:: calc_af

calc_hwe
^^^^^^^^
.. autofunction:: calc_hwe

compress_index_vcf
^^^^^^^^^^^^^^^^^^
.. autofunction:: compress_index_vcf

create_pos_haplotype
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: create_pos_haplotype

get_gt
^^^^^^
.. autofunction:: get_gt

get_scalebin
^^^^^^^^^^^^
.. autofunction:: get_scalebin

get_sizebin
^^^^^^^^^^^
.. autofunction:: get_sizebin

get_svtype
^^^^^^^^^^
.. autofunction:: get_svtype

msa2vcf
^^^^^^^
.. autofunction:: msa2vcf

overlap_percent
^^^^^^^^^^^^^^^
.. autofunction:: overlap_percent

overlaps
^^^^^^^^
.. autofunction:: overlaps

phab
^^^^
.. autofunction:: phab

phab_multi
^^^^^^^^^^
.. autofunction:: phab_multi

reciprocal_overlap
^^^^^^^^^^^^^^^^^^
.. autofunction:: reciprocal_overlap

ref_ranges
^^^^^^^^^^
.. autofunction:: ref_ranges

seqsim
^^^^^^
.. autofunction:: seqsim

sizesim
^^^^^^^
.. autofunction:: sizesim

unroll_compare
^^^^^^^^^^^^^^
.. autofunction:: unroll_compare

vcf_ranges
^^^^^^^^^^
.. autofunction:: vcf_ranges

Dev methods
-----------
benchdir_count_entries
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: benchdir_count_entries

chunker
^^^^^^^
.. autofunction:: chunker

cmd_exe
^^^^^^^
.. autofunction:: cmd_exe

consolidate_phab_vcfs
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: consolidate_phab_vcfs

count_entries
^^^^^^^^^^^^^
.. autofunction:: count_entries

fchain
^^^^^^
.. autofunction:: fchain

file_zipper
^^^^^^^^^^^
.. autofunction:: file_zipper

help_unknown_cmd
^^^^^^^^^^^^^^^^
.. autofunction:: help_unknown_cmd

make_temp_filename
^^^^^^^^^^^^^^^^^^
.. autofunction:: make_temp_filename

opt_gz_open
^^^^^^^^^^^
.. autofunction:: opt_gz_open

optimize_df_memory
^^^^^^^^^^^^^^^^^^
.. autofunction:: optimize_df_memory

performance_metrics
^^^^^^^^^^^^^^^^^^^
.. autofunction:: performance_metrics

restricted_float
^^^^^^^^^^^^^^^^
.. autofunction:: restricted_float

restricted_int
^^^^^^^^^^^^^^^^
.. autofunction:: restricted_int

setup_logging
^^^^^^^^^^^^^
.. autofunction:: setup_logging

setup_progressbar
^^^^^^^^^^^^^^^^^
.. autofunction:: setup_progressbar

vcf_to_df
^^^^^^^^^
.. autofunction:: vcf_to_df

Objects
-------

Bench
^^^^^
.. autoclass:: Bench
   :members:

BenchOutput
^^^^^^^^^^^
.. autoclass:: Bench
   :members: BenchOutput

GT
^^
.. autoclass:: GT
   :members:

RegionVCFIterator
^^^^^^^^^^^^^^^^^
.. autoclass:: RegionVCFIterator
   :members:

LogFileStderr
^^^^^^^^^^^^^
.. autoclass:: LogFileStderr
   :members:

MatchResult
^^^^^^^^^^^
.. autoclass:: MatchResult
   :members:

Matcher
^^^^^^^
.. autoclass:: Matcher
   :members:

SV
^^
.. autoclass:: SV
   :members:

Data
----
HEADERMAT
^^^^^^^^^
regular expression of vcf header INFO/FORMAT fields with groups

.. autodata:: HEADERMAT


QUALBINS
^^^^^^^^
0-100 quality score bin strings (step size 10)

.. autodata:: QUALBINS

SVTYTYPE
^^^^^^^^
:class:`pandas.CategoricalDtype` of :class:`truvari.SV`

SZBINMAX
^^^^^^^^
integer list of maximum size for size bins

.. autodata:: SZBINMAX

SZBINS
^^^^^^
string list of size bins

.. autodata:: SZBINS

SZBINTYPE
^^^^^^^^^
:class:`pandas.CategoricalDtype` of :data:`truvari.SZBINS`
