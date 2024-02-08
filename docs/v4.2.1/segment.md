Segmentation: Normalization of SVs into disjointed genomic regions

For SVs with a span in the genome (currently only DELs), split the overlaps into disjointed regions. This is an experimental tool that explores the possibility of assisting SV association analysis.

This tool adds an INFO field `SEGCNT` which holds the number of original SVs that overlap the newly reported region. It also adds a FORMAT field `SEG`, which is the 'allele coverage' per-sample. For example, if a sample has two overlapping heterozygous deletions, the shared region will have `SEG=2`. If the two deletions were homozygous then `SEG=4`.

In the below example, we have three deletions found across three samples.

![](https://github.com/spiralgenetics/truvari/blob/develop/imgs/segment_example.png)

The `segment` added annotations for the regions would then be:
| Region | INFO/SEGCNT | S1/SEG | S2/SEG | S3/SEG |
|--------|-------------|--------|--------|--------|  
| A | 1 | 2 | 0 | 0 |
| B | 2 | 2 | 1 | 0 |
| C | 3 | 2 | 2 | 2 |
| D | 2 | 2 | 1 | 0 |
| E | 1 | 0 | 1 | 0 |