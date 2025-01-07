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

        # pysam GT access
        if 'GT' in entry.samples['SAMPLE']:
            print(entry.samples[0]['GT'])
        # truvari GT access
        print(entry.gt('SAMPLE'))


Details of all available functions are in :ref:`package documentation <variant_handling>`

Besides some helpful accession functions, the `truvari.VariantRecord` also makes comparing two VCF entries easy.

.. code-block:: python

    # Given two `truvari.VariantRecords`, entry1 and entry2
    match = entry1.match(entry2)
    print("Entries' Sequence Similarity:", match.seqsim)
    print("Entries' Size Similarity:", match.sizesim)
    print("Is the match above thresholds:", match.state)

This returns a `truvari.MatchResult`. We can customize the thresholds for matching by giving the `truvari.VariantFile` a
`truvari.Matcher`.

.. code-block:: python

    # Turn off sequence and size similarity. Turn on reciprocal overlap
    matcher = truvari.Matcher(seqsim=0, sizesim=0, recovl=0.5, ...)
    vcf = truvari.VariantFile("input.vcf.gz", matcher=matcher)
    entry1 = next(vcf)
    entry2 = next(vcf)
    match = entry1.match(entry2)

Another useful function is a quick way filter variants. The `truvari.Matcher` has parameters for e.g. minimum or maximum
size of SVs one wants to analyze which can be leveraged via:

.. code-block:: python

    matcher = truvari.Matcher(sizemin=200, sizemax=500)
    vcf = truvari.VariantFile("input.vcf.gz", matcher=matcher)
    # Grab all of the variant records between sizemin and sizemax
    results = [entry for entry in vcf if not entry.size_filter()]

Additional filtering for things such as monomorphic reference sites or single-end BNDs are available by calling `entry.filter_call()`
