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

VariantFile
^^^^^^^^^^^
.. autoclass:: VariantFile
   :members:

VariantRecord
^^^^^^^^^^^^^
.. autoclass:: VariantRecord
   :members:

Extra Methods
-------------
bed_ranges
^^^^^^^^^^
.. autofunction:: bed_ranges

best_seqsim
^^^^^^^^^^^
.. autofunction:: best_seqsim

read_bed_tree
^^^^^^^^^^^^^^^
.. autofunction:: read_bed_tree

check_vcf_index
^^^^^^^^^^^^^^^
.. autofunction:: check_vcf_index

compress_index_vcf
^^^^^^^^^^^^^^^^^^
.. autofunction:: compress_index_vcf

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

reciprocal_overlap
^^^^^^^^^^^^^^^^^^
.. autofunction:: reciprocal_overlap

ref_ranges
^^^^^^^^^^
.. autofunction:: ref_ranges

roll_seqsim
^^^^^^^^^^^
.. autofunction:: roll_seqsim

seqsim
^^^^^^
.. autofunction:: seqsim

sizesim
^^^^^^^
.. autofunction:: sizesim

unroll_seqsim
^^^^^^^^^^^^^
.. autofunction:: unroll_seqsim

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

coords_within
^^^^^^^^^^^^^
.. autofunction:: coords_within

count_entries
^^^^^^^^^^^^^
.. autofunction:: count_entries

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
.. autoclass:: BenchOutput
   :members:

StatsBox
^^^^^^^^
.. autoclass:: StatsBox
   :members:

GT
^^
.. autoclass:: GT
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
