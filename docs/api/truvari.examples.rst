Truvari API QuickStart
======================

Truvari provides functionality to facilitate comparison of SVs in VCF variant records.
These functions are available to developers by simply replacing calls to `pysam.VariantFile` with `truvari.VariantFile`.
All of the `pysam` functionality is still available in this object with the exception of context managers (i.e. `with`
statements).

.. code-block:: python

    import truvari
    vcf = truvari.VariantFile("input.vcf.gz")
    for entry in vcf:
        print(entry.info['SVTYPE']) # pysam access to the variant's INFO fields
        print(entry.allele_freq_annos()) # truvari calculation of a variant's allele frequency


Details of all available functions are in :ref:`package documentation <variant_handling>`
