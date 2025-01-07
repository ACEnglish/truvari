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

Extra Methods
-------------
.. autofunction:: bed_ranges

.. autofunction:: best_seqsim

.. autofunction:: read_bed_tree

.. autofunction:: check_vcf_index

.. autofunction:: compress_index_vcf

.. autofunction:: get_gt

.. autofunction:: get_scalebin

.. autofunction:: get_sizebin

.. autofunction:: get_svtype

.. autofunction:: msa2vcf

.. autofunction:: overlap_percent

.. autofunction:: overlaps

.. autofunction:: phab

.. autofunction:: reciprocal_overlap

.. autofunction:: ref_ranges

.. autofunction:: roll_seqsim

.. autofunction:: seqsim

.. autofunction:: sizesim

.. autofunction:: unroll_seqsim

.. autofunction:: vcf_ranges

Dev methods
-----------
.. autofunction:: benchdir_count_entries

.. autofunction:: chunker

.. autofunction:: cmd_exe

.. autofunction:: coords_within

.. autofunction:: count_entries

.. autofunction:: file_zipper

.. autofunction:: help_unknown_cmd

.. autofunction:: make_temp_filename

.. autofunction:: opt_gz_open

.. autofunction:: optimize_df_memory

.. autofunction:: performance_metrics

.. autofunction:: restricted_float

.. autofunction:: restricted_int

.. autofunction:: setup_logging

.. autofunction:: vcf_to_df

Objects
-------

.. autoclass:: Bench
   :members:

.. autoclass:: BenchOutput
   :members:

.. autoclass:: StatsBox
   :members:

.. autoclass:: GT
   :members:

.. autoclass:: LogFileStderr
   :members:

.. autoclass:: MatchResult
   :members:

.. autoclass:: Matcher
   :members:

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
