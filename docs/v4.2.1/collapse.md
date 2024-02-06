`collapse` is Truvari's approach to SV merging. After leveraging `bcftools` to merge VCFs, `truvari collapse` can then iterate over the calls and create clusters of SVs that match over the [provided thresholds](https://github.com/spiralgenetics/truvari/wiki/bench#matching-parameters). This is also useful when removing potentially redundant calls within a single sample. 

Example
=======
To start, we merge multiple VCFs (each with their own sample) and ensure there are no multi-allelic entries via: 
```bash
bcftools merge -m none one.vcf.gz two.vcf.gz | bgzip > merge.vcf.gz
```

This will `paste` SAMPLE information between vcfs when calls have the exact same chrom, pos, ref, and alt.
For example, consider two vcfs:

    >> one.vcf:
    chr1 1 ... GT 0/1
    chr1 5 ... GT 1/1
    >> two.vcf:
    chr1 1 ... GT 1/1
    chr1 7 ... GT 0/1

`bcftools merge` creates:

    >> merge.vcf:
    chr1 1 ... GT 0/1 1/1
    chr1 5 ... GT 1/1 ./.
    chr1 7 ... GT ./. 0/1    

This VCF can then be collapsed to allow 'fuzzier' matching than the exact merge just performed.

```bash
truvari collapse -i merge.vcf.gz -o truvari_merge.vcf -c truvari_collapsed.vcf -f /path/to/reference.fa
```

For example, if we collapsed our example merge.vcf by matching any calls within 3bp, we'd create:

    >> truvari_merge.vcf
    chr1 1 ... GT 0/1 1/1
    chr1 5 ... GT 1/1 0/1
    >> truvari_collapsed.vcf
    chr1 7 ... GT ./. 0/1    

--choose behavior
=================
When collapsing, the default `--choose` behavior is to take the first variant from a cluster to
be written to the output while the others will be placed in the collapsed output. 
Other choosing options are `maxqual` (the call with the highest quality score) or `common` (the call with the highest minor allele count).

Samples with no genotype information in the kept variant will be filled by the first
collapsed variant containing genotype information.                                                                                    

--gt
====
For some results, one may not want to collapse variants with conflicting genotypes from a single sample. With the `--gt all` parameter, variants which are present (non `0/0` or `./.`) in the same sample are not collapsed. With the `-gt het` parameter, only variants which are both heterozygous in a sample (e.g. `0/1` and `0/1`) are prevented from collapsing. The `--gt het` is useful for some SV callers which will redundantly call variants and typically genotype them all as `1/1`.

--intra
=======
When a single sample is run through multiple SV callers, one may wish to consolidate those results. After the `bcftools merge` of the VCFs, there will be one SAMPLE column per-input. With `--intra`, collapse will consolidate the sample information so that only a single sample column is present in the output. Since the multiple callers may have different genotypes or other FORMAT fields with conflicting information, `--intra` takes the first column from the VCF, then second, etc. For example, if we have an entry with:
```
FORMAT    RESULT1     RESULT2
GT:GQ:AD  ./.:.:3,0  1/1:20:0,30
```
The `--intra` output would be:
```
FORMAT    RESULT1
GT:GQ:AD  1/1:20:3,0
```
As you can see in this example, 1) The first sample name is the only one preserved. 2) conflicting FORMAT fields can be consolidated in a non-useful way (here the AD of `3,0` isn't informative to a `1/1` genotype). We're working to provide an API to help users write custom intra-sample consolidation scripts.

--hap mode
==========
When using `--hap`, we assume phased variants from a single individual. Only the
single best matching call from the other haplotype will be collapsed,
and the consolidated genotype will become 1/1

For example, if we collapse anything at the same position:

    chr1 1 .. GT 0|1
    chr1 1 .. GT 1|0
    chr1 2 .. GT 1|0

will become:

    chr1 1 .. GT 1/1
    chr1 2 .. GT 1|0

--chain mode
============
Normally, every variant in a set of variants that are collapsed together matches every other variant in the set. However, when using `--chain` mode, we allow 'transitive matching'. This means that all variants match to only at least one other variant in the set. In situations where a 'middle' variant has two matches that don't match each other, without `--chain` the locus will produce two variants whereas using `--chain` will produce one.
For example, if we have

    chr1 5 ..
    chr1 7 ..
    chr1 9 ..

When we collapse anything within 2bp of each other, without `--chain`, we output:

    chr1 5 ..
    chr1 9 ..

With `--chain`, we would collapse `chr1 9` as well, producing

    chr1 5 ..

Annotations
===========
`collapse` produces two files. The output file has kept variants along with unanalyzed (< sizemin) variants. The collapsed file contains the variants that were collapsed into the kept variants. 

The output file has only two annotations added to the `INFO`. 
- `CollapseId` - Identifier of the variant when comparing to the collapse outputs. 
- `NumCollapsed` - Number of variants collapsed into this variant
- `NumConsolidated` - Number of samples' genotypes consolidated into this call's genotypes

The collapsed file has all of the annotations added by [[bench|bench#definition-of-annotations-added-to-tp-vcfs]]. Note that `MatchId` is tied to the output file's `CollapseId`. See [MatchIds](https://github.com/spiralgenetics/truvari/wiki/MatchIds) for details.

```
usage: collapse [-h] -i INPUT [-o OUTPUT] [-c COLLAPSED_OUTPUT] [-f REFERENCE] [-k {first,maxqual,common}] [--debug]
                [-r REFDIST] [-p PCTSIM] [-B MINHAPLEN] [-P PCTSIZE] [-O PCTOVL] [-t] [--use-lev] [--hap] [--chain]
                [--no-consolidate] [--null-consolidate NULL_CONSOLIDATE] [-s SIZEMIN] [-S SIZEMAX] [--passonly]

Structural variant collapser

Will collapse all variants within sizemin/max that match over thresholds

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Comparison set of calls
  -o OUTPUT, --output OUTPUT
                        Output vcf (stdout)
  -c COLLAPSED_OUTPUT, --collapsed-output COLLAPSED_OUTPUT
                        Where collapsed variants are written (collapsed.vcf)
  -f REFERENCE, --reference REFERENCE
                        Indexed fasta used to call variants
  -k {first,maxqual,common}, --keep {first,maxqual,common}
                        When collapsing calls, which one to keep (first)
  --debug               Verbose logging
  --hap                 Collapsing a single individual's haplotype resolved calls (False)
  --chain               Chain comparisons to extend possible collapsing (False)
  --no-consolidate      Skip consolidation of sample genotype fields (True)
  --null-consolidate NULL_CONSOLIDATE
                        Comma separated list of FORMAT fields to consolidate into the kept entry by taking the first non-null
                        from all neighbors (None)

Comparison Threshold Arguments:
  -r REFDIST, --refdist REFDIST
                        Max reference location distance (500)
  -p PCTSIM, --pctsim PCTSIM
                        Min percent allele sequence similarity. Set to 0 to ignore. (0.95)
  -B MINHAPLEN, --minhaplen MINHAPLEN
                        Minimum haplotype sequence length to create (50)
  -P PCTSIZE, --pctsize PCTSIZE
                        Min pct allele size similarity (minvarsize/maxvarsize) (0.95)
  -O PCTOVL, --pctovl PCTOVL
                        Min pct reciprocal overlap (0.0) for DEL events
  -t, --typeignore      Variant types don't need to match to compare (False)
  --use-lev             Use the Levenshtein distance ratio instead of edlib editDistance ratio (False)

Filtering Arguments:
  -s SIZEMIN, --sizemin SIZEMIN
                        Minimum variant size to consider for comparison (50)
  -S SIZEMAX, --sizemax SIZEMAX
                        Maximum variant size to consider for comparison (50000)
  --passonly            Only consider calls with FILTER == PASS
```