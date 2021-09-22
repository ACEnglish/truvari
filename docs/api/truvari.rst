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

SV
^^
.. autoclass:: SV
   :members:

Methods
-------
allele_freq_annos
^^^^^^^^^^^^^^^^^
.. autofunction:: allele_freq_annos  

bed_ranges
^^^^^^^^^^
.. autofunction:: bed_ranges

cmd_exe
^^^^^^^
.. autofunction:: cmd_exe

copy_entry
^^^^^^^^^^
.. autofunction:: copy_entry

create_pos_haplotype
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: create_pos_haplotype

entry_boundaries
^^^^^^^^^^^^^^^^
.. autofunction:: entry_boundaries

entry_create_haplotype
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: entry_create_haplotype

entry_gt_comp
^^^^^^^^^^^^^
.. autofunction:: entry_gt_comp

entry_is_present
^^^^^^^^^^^^^^^^
.. autofunction:: entry_is_present

entry_pctsim
^^^^^^^^^^^^
.. autofunction:: entry_pctsim

entry_reciprocal_overlap
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: entry_reciprocal_overlap

ref_ranges
^^^^^^^^^^
.. autofunction:: ref_ranges

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

fetch_coords
^^^^^^^^^^^^
.. autofunction:: fetch_coords

filter_value
^^^^^^^^^^^^
.. autofunction:: filter_value

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

make_interval_tree
^^^^^^^^^^^^^^^^^^
.. autofunction:: make_interval_tree

match_sorter
^^^^^^^^^^^^
.. autofunction:: match_sorter

overlaps
^^^^^^^^
.. autofunction:: overlaps

reciprocal_overlap
^^^^^^^^^^^^^^^^^^
.. autofunction:: reciprocal_overlap

restricted_float
^^^^^^^^^^^^^^^^
.. autofunction:: restricted_float

setup_logging
^^^^^^^^^^^^^
.. autofunction:: setup_logging

setup_progressbar
^^^^^^^^^^^^^^^^^
.. autofunction:: setup_progressbar

sizesim
^^^^^^^
.. autofunction:: sizesim

vcf_to_df
^^^^^^^^^
.. autofunction:: vcf_to_df

weighted_score
^^^^^^^^^^^^^^
.. autofunction:: weighted_score

Data
----
HEADERMAT
^^^^^^^^^
regular expression of vcf header INFO/FORMAT fields with groups

.. autodata:: HEADERMAT

MATCHRESULT
^^^^^^^^^^^
named tuple of match files

.. autodata:: MATCHRESULT

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
