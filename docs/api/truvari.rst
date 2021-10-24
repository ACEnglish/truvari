truvari package
===============

Overview
--------
.. automodule:: truvari
   :members:
   :undoc-members:
   :show-inheritance:

Objects
-------
GT
^^
.. autoclass:: GT
   :members:

GenomeTree
^^^^^^^^^^
.. autoclass:: GenomeTree
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

VariantRecord Methods
---------------------
copy_entry
^^^^^^^^^^
.. autofunction:: copy_entry

entry_boundaries
^^^^^^^^^^^^^^^^
.. autofunction:: entry_boundaries

entry_create_haplotype
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: entry_create_haplotype

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

entry_pctsim
^^^^^^^^^^^^
.. autofunction:: entry_pctsim

entry_reciprocal_overlap
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: entry_reciprocal_overlap

entry_same_variant_type
^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: entry_same_variant_type

entry_size
^^^^^^^^^^
.. autofunction:: entry_size

entry_size_similarity
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: entry_size_similarity

entry_to_key
^^^^^^^^^^^^
.. autofunction:: entry_to_key

entry_variant_type
^^^^^^^^^^^^^^^^^^
.. autofunction:: entry_variant_type

Extra Methods
---------------------
allele_freq_annos
^^^^^^^^^^^^^^^^^
.. autofunction:: allele_freq_annos  

bed_ranges
^^^^^^^^^^
.. autofunction:: bed_ranges

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

overlaps
^^^^^^^^
.. autofunction:: overlaps

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

weighted_score
^^^^^^^^^^^^^^
.. autofunction:: weighted_score

Dev methods
-----------
cmd_exe
^^^^^^^
.. autofunction:: cmd_exe

help_unknown_cmd
^^^^^^^^^^^^^^^^
.. autofunction:: help_unknown_cmd

optimize_df_memory
^^^^^^^^^^^^^^^^^^
.. autofunction:: optimize_df_memory

restricted_float
^^^^^^^^^^^^^^^^
.. autofunction:: restricted_float

setup_logging
^^^^^^^^^^^^^
.. autofunction:: setup_logging

setup_progressbar
^^^^^^^^^^^^^^^^^
.. autofunction:: setup_progressbar

vcf_to_df
^^^^^^^^^
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
