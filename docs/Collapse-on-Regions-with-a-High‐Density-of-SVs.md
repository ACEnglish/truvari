Some regions of the reference genome can give rise to a huge number of SVs. These regions are typically in telomeres or centromeres, especially on chm13. Furthermore, some multi-sample VCFs might contain samples with a high genetic distance from the reference which can also create an unreasonable number of SVs. These 'high-density' regions can cause difficulties for `truvari collapse` when there are too many comparisons that need to be made.

If you find `truvari collapse` 'freezing' during processing where it is no longer writing variants to the output VCFs, you should investigate if there are regions which should be excluded from analysis. To do this, first run:

```
truvari anno chunks input.vcf.gz > counts.bed
```

The `counts.bed` will have the chromosome, start, and end position of each chunk in the VCF as well as three additional columns. The 4th column has the number of SVs inside the chunk while the 5th and 6th have a comma-separated counts of the number of variants in 'sub-chunks'. These sub-chunk counts correspond to advanced separation of the variants beyond just their window. The 5th is after separating by svtype and svlen, the 6th is re-chunking the svtype/svlen separation by distance.

If you find spans in the `counts.bed` with a huge number of SVs, these are prime candidates for exclusion to speed up `truvari collapse`. To exclude them, you first subset to the regions of interest and then use `bedtools` to create a bed file that will skip them.

```
# exclude regions with more than 30k SVs. You should test this threshold.
awk '$4 >= 30000' counts.bed > to_exclude.bed
bedtools complement -g genome.bed -i to_exclude.bed > to_analyize.bed

truvari collapse --bed to_analyze.bed -i input.vcf.gz
```

When considering what threshold you would like to use, just looking at the 4th column may not be sufficient as the 'sub-chunks' may be smaller and will therefore run faster. There also is no guarantee that a region with high SV density will be slow. For example, if all SVs in a chunk could collapse, it would only take O(N - 1) comparisons to complete the collapse. The problems arise when the SVs have few redundant SVs and therefore requires O(N**2) comparisons.

A conservative workflow for figuring out which regions to exclude would run `truvari collapse` without a `--bed`, wait for regions to 'freeze', and then check if the last variant written out by `truvari collapse` is just before a high-density chunk in the `counts.bed`. That region could then be excluded and collapsing repeated until success. 