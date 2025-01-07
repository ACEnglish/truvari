Truvari API QuickStart
======================

Truvari provides functionality to facilitate the comparison of structural variants (SVs) in VCF variant records. Developers can easily leverage this functionality by replacing calls to `pysam.VariantFile` with `truvari.VariantFile`. The `truvari.VariantFile` retains all `pysam` functionality except for context managers (i.e., `with` statements).

.. code-block:: python

    import truvari
    vcf = truvari.VariantFile("input.vcf.gz")
    for entry in vcf:
        print(entry.info['SVTYPE'])  # Access variant's INFO fields using pysam
        print(entry.allele_freq_annos())  # Calculate variant's allele frequency with truvari

        # Access genotype (GT) using pysam
        if 'GT' in entry.samples['SAMPLE']:
            print(entry.samples['SAMPLE']['GT'])
        # Access genotype (GT) using truvari
        print(entry.gt('SAMPLE'))

Details of all available functions are in :ref:`package documentation <variant_handling>`.

Comparing Variants
------------------

The `truvari.VariantRecord` simplifies comparing two VCF entries.

.. code-block:: python

    # Given two `truvari.VariantRecords`, entry1 and entry2
    match = entry1.match(entry2)
    print("Entries' Sequence Similarity:", match.seqsim)
    print("Entries' Size Similarity:", match.sizesim)
    print("Is the match above thresholds:", match.state)

This returns a `truvari.MatchResult`. You can customize matching thresholds by providing a `truvari.VariantParams` to the `truvari.VariantFile`.

.. code-block:: python

    # Disable sequence and size similarity; enable reciprocal overlap
    matcher = truvari.VariantParams(seqsim=0, sizesim=0, recovl=0.5)
    vcf = truvari.VariantFile("input.vcf.gz", matcher=matcher)
    entry1 = next(vcf)
    entry2 = next(vcf)
    match = entry1.match(entry2)

Filtering Variants
------------------

The `truvari.VariantParams` provides parameters for filtering variants, such as minimum or maximum SV sizes.

.. code-block:: python

    matcher = truvari.VariantParams(sizemin=200, sizemax=500)
    vcf = truvari.VariantFile("input.vcf.gz", matcher=matcher)
    # Retrieve all variant records within sizemin and sizemax
    results = [entry for entry in vcf if not entry.size_filter()]

Additional filters, such as excluding monomorphic reference sites or single-end BNDs, can be applied using `entry.filter_call()`.

Subsetting to Regions
---------------------

To subset a VCF to regions specified in a BED file, use:

.. code-block:: python

    for entry in vcf.bed_fetch("regions.bed"):
        print("Entry's variant type:", entry.var_type())
        print("Entry's variant size:", entry.var_size())

If your regions of interest are stored in an in-memory object instead of a BED file, use the `.regions_fetch` method:

.. code-block:: python

    from collections import defaultdict
    from pyintervaltree import IntervalTree
    tree = defaultdict(IntervalTree)
    tree['chr1'].addi(10, 100)
    tree['chr2'].addi(2000, 2200)
    count = 0
    for entry in vcf.regions_fetch(tree):
        count += 1
    print(f"Total of {count} variants")

To iterate over variants that are not within the regions, use `vcf.regions_fetch(tree, within=False)`.

Parsing BND Information
-----------------------

Truvari simplifies parsing BND information from VCF entries:

.. code-block:: python

    # Example entry:
    # chr1  23272628  SV_1  G  G]chr5:52747359]  .  PASS  SVTYPE=BND;EVENTTYPE=TRA:UNBALANCED;SUBCLONAL=n;COMPLEX=n;MATEID=SV_171  GT:PSL:PSO  0/1:.:.
    print(entry.bnd_position())
    # ('chr5', 52747359)
    print(entry.bnd_direction_strand())
    # ('right', 'direct')

