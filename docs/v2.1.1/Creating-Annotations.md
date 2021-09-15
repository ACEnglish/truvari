## Bed-like reference annotations 
Making a reference range annotation

1) Download the track for your reference. UCSC Genome Browser is pretty good
2) A truvari annotation track has the following columns:
- CHROM : chromosome position of the reference range
- START : start position of the reference range
- END : end position of the reference range
- ANNOS* : One or more annotation columns on each line. These are the values that will be added to the VCF INFO fields

These should be zero-based half-open intervals

3) Add the header information
For each ANNOS*, create a header line. These headerlines are the exact same as VCF INFO header lines. See
https://samtools.github.io/hts-specs/VCFv4.2.pdf for details.
These header lines should be in the same order as the columns in each row in the annotation file

Example
```
##INFO=<ID=SREP_periods,Number=.,Type=Float,Description="Simple Repeat period lengths">
##INFO=<ID=SREP_copies,Number=.,Type=Float,Description="Simple Repeat copy numbers">
1       10000   10468   6       77.2
1       10627   10800   29      6
```

So, any entry that intersects with 1:1000-10468 will have the following added to their INFO:

```
SREP_periods=6;SREP_copies=77.2
```

4) If some columns are not populated on certain lines, use '.' for Null