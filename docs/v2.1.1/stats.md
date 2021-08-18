
# Running

```
usage: stats [-h] [-d DATAFRAME] [-o OUT] [--qmin QMIN] [--qmax QMAX]
             VCF [VCF ...]

VCF Stats for SVs with sizebins and svtypes

positional arguments:
  VCF                   VCFs to annotate (stdin)

optional arguments:
  -h, --help            show this help message and exit
  -d DATAFRAME, --dataframe DATAFRAME
                        Write dataframe joblib to file (None)
  -o OUT, --out OUT     Stats output file (/dev/stdout
  --qmin QMIN           Minimum QUAL score found in VCF (0)
  --qmax QMAX           Maximum QUAL score found in VCF (100)
  ```
  
# Example
`$ truvari stats tp-call.vcf.gz`

## How it works
Counts per SVTYPE, SZBIN, QUALBIN, Genotype the number of events found.

The text output makes multiple tables for 'total' and per-sample. Per-sample is the subset of total where the sample has a HET or HOM genotype.

Optionally outputs a joblib dictionary (`--dataframe`) with keys for the total, and per-sample. The total array doesn't have the GT dimension. See this [notebook](https://github.com/spiralgenetics/truvari/blob/develop/TruvariStatsPlots.ipynb) for examples of how to use this data. 
