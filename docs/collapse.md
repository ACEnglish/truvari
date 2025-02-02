`collapse` is Truvari's approach to SV merging. After leveraging `bcftools` to merge VCFs, `truvari collapse` can then iterate over the calls and create clusters of SVs that match over the [provided thresholds](https://github.com/ACEnglish/truvari/wiki/bench#matching-parameters). This is also useful when removing potentially redundant calls within a single sample. 

Since `collapse` uses the same variant comparison engine as `bench`, there is a lot of overlap in the behavior. Be sure to also be familiar with the [bench documentation](https://github.com/ACEnglish/truvari/wiki/bench).

Example
=======
To start, we merge multiple VCFs (each with their own sample) and ensure there are no multi-allelic entries via: 
```bash
bcftools merge -m none one.vcf.gz two.vcf.gz | bgzip > merge.vcf.gz
```
WARNING! If you have symbolic variants, see [the below section](https://github.com/ACEnglish/truvari/wiki/collapse#symbolic-variants) on using bcftools.

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
truvari collapse -i merge.vcf.gz -o truvari_merge.vcf -c truvari_collapsed.vcf
```

For example, if we collapsed our example merge.vcf by matching any calls within 3bp, we'd create:

    >> truvari_merge.vcf
    chr1 1 ... GT 0/1 1/1
    chr1 5 ... GT 1/1 0/1
    >> truvari_collapsed.vcf
    chr1 7 ... GT ./. 0/1    

Symbolic Variants
=================
bcftools may not handle symbolic variants correctly since it doesn't consider their END position. To correct for this, ensure that every input variant has a unique ID and use `bcftools merge -m id`. For example:
```
# A.vcf
chr1	147022730	SV1	N	<DEL>	.	PASS	SVLEN=-570334;END=147593064
# B.vcf
chr1	147022730	SV2	N	<DEL>	.	PASS	SVLEN=-990414;END=148013144

# bcftools merge -m none A.vcf B.vcf
# Premature collapse
chr1	147022730	SV1;SV2	N	<DEL>	.	PASS	SVLEN=-570334;END=147593064

# bcftools merge -m id A.vcf B.vcf
chr1	147022730	SV1	N	<DEL>	.	PASS	SVLEN=-570334;END=147593064
chr1	147022730	SV2	N	<DEL>	.	PASS	SVLEN=-990414;END=148013144
```

This bug has been replicated with bcftools 1.18 and 1.21.

--choose behavior
=================
When collapsing, the default `--choose` behavior is to take the `first` variant by position from a cluster to
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

    chr1 1 ..
    chr1 4 ..
    chr1 7 ..
    chr1 10 ..

We take the `chr1 1` variant and find all its matches. When we collapse anything within 5bp of each other, without `--chain`, we output:

    chr1 1 ..
    chr1 7 ..

With `--chain`, we would allow one level of transitive matching. This means that after finding the `chr1 1 -> chr1 4` match, we check `chr1 4` against all the remaining variants and would output

    chr1 1 ..
    chr1 10 ..

Note that this leaves `chr1 10` because we don't do multiple levels of transitive matching, meaning we never compare `chr1 7` to `chr1 10`. This is preferred because otherwise variants which have a continuous range of similarity could all be collapsed into a single variant. e.g., if the position in this example were sizes and, we wouldn't want the 1bp variant being a kept representation for all the variants.

Annotations
===========
`collapse` produces two files. The output file has kept variants along with unanalyzed (< sizemin) variants. The collapsed file contains the variants that were collapsed into the kept variants. 

The output file has only two annotations added to the `INFO`. 
- `CollapseId` - Identifier of the variant when comparing to the collapse outputs. 
- `NumCollapsed` - Number of variants collapsed into this variant
- `NumConsolidated` - Number of samples' genotypes consolidated into this call's genotypes

The collapsed file has all of the annotations added by [[bench|bench#definition-of-annotations-added-to-tp-vcfs]]. Note that `MatchId` is tied to the output file's `CollapseId`. See [MatchIds](https://github.com/spiralgenetics/truvari/wiki/MatchIds) for details.