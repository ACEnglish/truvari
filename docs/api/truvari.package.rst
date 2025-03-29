Truvari package
===============

Overview
--------
.. automodule:: truvari
   :members:
   :undoc-members:
   :show-inheritance:

.. _variant_handling:

Variant Handling
----------------

.. autoclass:: VariantFile
   :members:

.. autoclass:: VariantRecord
   :members:

.. _variant_params:

.. autoclass:: VariantParams
   :members:

Objects
-------

.. _match_result:

.. autoclass:: MatchResult
   :members:

.. autoclass:: GT
   :members:

.. autoclass:: SV
   :members:

.. autoclass:: Bench
   :members:

.. autoclass:: BenchOutput
   :members:

.. autoclass:: StatsBox
   :members:

.. autoclass:: LogFileStderr
   :members:

Extra Methods
-------------
.. autofunction:: bed_ranges

.. autofunction:: benchdir_count_entries

.. autofunction:: best_seqsim

.. autofunction:: build_region_tree

.. autofunction:: read_bed_tree

.. autofunction:: check_vcf_index

.. autofunction:: chunker

.. autofunction:: cmd_exe

.. autofunction:: compress_index_vcf

.. autofunction:: coords_within

.. autofunction:: count_entries

.. autofunction:: extend_region_tree

.. autofunction:: file_zipper

.. autofunction:: help_unknown_cmd

.. autofunction:: get_gt

.. autofunction:: get_scalebin

.. autofunction:: get_sizebin

.. autofunction:: get_svtype

.. autofunction:: make_temp_filename

.. autofunction:: merge_region_tree_overlaps

.. autofunction:: msa2vcf

.. autofunction:: opt_gz_open

.. autofunction:: optimize_df_memory

.. autofunction:: overlap_percent

.. autofunction:: overlaps

.. autofunction:: performance_metrics

.. autofunction:: reciprocal_overlap

.. autofunction:: restricted_float

.. autofunction:: restricted_int

.. autofunction:: ref_ranges

.. autofunction:: roll_seqsim

.. autofunction:: seqsim

.. autofunction:: setup_logging

.. autofunction:: sizesim

.. autofunction:: unroll_seqsim

.. autofunction:: vcf_ranges

.. autofunction:: vcf_to_df

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
