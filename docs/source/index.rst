.. Truvari documentation master file, created by
   sphinx-quickstart on Thu May 14 20:15:26 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Truvari
=======

Blurb - What is this thing about

See Updates for information on new versions

Installation
------------
Truvari uses python 3.7+ and can be installed with pip:

:: 

  pip install truvari

To build from the repository, simply run:

::
  
  python setup.py sdist bdist_wheel
  pip install dist/Truvari-<version>.tar.gz

Where ``<version>`` is which version you just built.

Quick Start
-----------

Each sub-command contains help documentation. Start with ``truvari -h`` to see available commands.

The current most common use case is for structural variation benchmarking:

::

  truvari bench -b base.vcf.gz -c comp.vcf.gz -r reference.fasta -o output_dir/

More Information
----------------
See the documentation for more detail

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. automodule:: truvari.bench
   :members:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
