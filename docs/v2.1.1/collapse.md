# truvari collapse

Will collapse all variants within sizemin/max that match over thresholds. All variants outside size boundaries will be placed into the output.

Collapse will iterate over a single VCF and 'collapse' variants that match one another above the provided thresholds. This is useful when removing potentially redundant calls or merging calls between samples/software. 

Merging SVs is a tricky business. Most approaches will simply use distances or reciprocal overlap of calls to merge them. There may be a more informed way to merge calls, especially when they're sequence-resolved breakpoint exact calls. This command helps do a more sophisticated merge.

To start, we can merge multiple VCFs (each with their own sample) using `bcftools merge -m none one.vcf two.vcf | bgzip > merge.vcf`.
This will `paste` SAMPLE information between vcfs when calls have the exact same chrom, pos, ref, and alt.
For example:

    >> one.vcf:
    chr1 1 ... GT 0/1
    chr1 5 ... GT 1/1
    >> two.vcf:
    chr1 1 ... GT 1/1
    chr1 7 ... GT 0/1

Creates:

    >> merge.vcf:
    chr1 1 ... GT 0/1 1/1
    chr1 5 ... GT 1/1 ./.
    chr1 7 ... GT ./. 0/1    

This single VCF can then be collapsed to allow 'fuzzier' matching than the exact merge just performed.
For example, if we collapsed the above example merge.vcf by matching any calls within 3bp, we'd create:

    >> merge.collapse.vcf
    chr1 1 ... GT 0/1 1/1
    chr1 5 ... GT 1/1 0/1

When collapsing, the first variant from a matching set of variants will
be written to the output while the others will be placed in the collapsed output.

Samples with no genotype information in the first variant will be filled by the first
collapsed variant containing genotype information.                                                                                    

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

When using `--chain` mode, instead of collapsing all variants matching the first variant
together, we'll collapse all variants in a matching set together.
For example, if we have

    chr1 5 ..
    chr1 7 ..
    chr1 9 ..

When we collapse anything within 2bp of each other, without --chain, we output:

    chr1 5 ..
    chr1 9 ..

With `--chain`, we would collapse `chr1 9` as well, producing

    chr1 5 ..


```
usage: collapse [-h] -i INPUT [-o OUTPUT] [-c COLLAPSED_OUTPUT] [-f REFERENCE] [--debug] [-r REFDIST] [-p PCTSIM] [-B BUFFER] [-P PCTSIZE]
                [-O PCTOVL] [-t] [--hap] [--chain] [-s SIZEMIN] [-S SIZEMAX] [--passonly]

Structural variant collapser

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Comparison set of calls
  -o OUTPUT, --output OUTPUT
                        Output vcf (stdout)
  -c COLLAPSED_OUTPUT, --collapsed-output COLLAPSED_OUTPUT
                        Where collapsed variants are written (collapsed.vcf)
  -f REFERENCE, --reference REFERENCE
                        Indexed fasta used to call variants
  --debug               Verbose logging
  --hap                 Collapsing a single individual's haplotype resolved calls (False)
  --chain               Chain comparisons to extend possible collapsing (False)

Comparison Threshold Arguments:
  -r REFDIST, --refdist REFDIST
                        Max reference location distance (500)
  -p PCTSIM, --pctsim PCTSIM
                        Min percent allele sequence similarity. Set to 0 to ignore. (0.95)
  -B BUFFER, --buffer BUFFER
                        Percent of the reference span to buffer the haplotype sequence created
  -P PCTSIZE, --pctsize PCTSIZE
                        Min pct allele size similarity (minvarsize/maxvarsize) (0.95)
  -O PCTOVL, --pctovl PCTOVL
                        Minimum pct reciprocal overlap (0.0) for DEL events
  -t, --typeignore      Variant types don't need to match to compare (False)

Filtering Arguments:
  -s SIZEMIN, --sizemin SIZEMIN
                        Minimum variant size to consider for comparison (50)
  -S SIZEMAX, --sizemax SIZEMAX
                        Maximum variant size to consider for comparison (50000)
  --passonly            Only consider calls with FILTER == PASS
```