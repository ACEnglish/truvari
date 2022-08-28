
In addition to looking at performance of a single set of variation against a baseline, one may wish to measure the consistency between multiple sets of variation. The tool `truvari consistency` can automatically create that result. 

Running
=======

```
usage: consistency [-h] [-j] VCFs [VCFs ...]

Over multiple vcfs, calculate their intersection/consistency.

Calls will match between VCFs if they have a matching key of:
    CHROM:POS ID REF ALT

positional arguments:
  VCFs        VCFs to intersect

optional arguments:
  -h, --help  show this help message and exit
  -j, --json  Output report in json format
```
Example
=======

```bash
truvari consistency fileA.vcf fileB.vcf fileC.vcf
```

Matching Entries
================

VCF entries will be considered matching if and only if they have an exact same key of `CHROM:POS ID REF ALT`. Because of this stringency, it is recommend that you compare the tp-base.vcf or fn.vcf results from each individual VCF's Truvari output.

Output Report
=============

Below is an example report:

```text
#
# Total 5534 calls across 3 VCFs
#
#File	NumCalls
fileA.vcf	4706
fileB.vcf	4827
fileC.vcf	4882
#
# Summary of consistency
#
#VCFs	Calls	Pct
3	3973	71.79%
2	935	16.90%
1	626	11.31%
#
# Breakdown of VCFs' consistency
#
#Group	Total	TotalPct	PctOfFileCalls
111	3973	71.79%	84.42% 82.31% 81.38%
011	351	6.34%	7.27% 7.19%
101	308	5.57%	6.54% 6.31%
110	276	4.99%	5.86% 5.72%
001	250	4.52%	5.12%
010	227	4.10%	4.70%
100	149	2.69%	3.17%
```

At the top we see that we compared 5,534 unique variants between the 3 files, with fileC.vcf having the most calls at 4,882.

The "Summary of consistency" shows us that 3,973 (71.79%) of all the calls are shared between the 3 VCFs, while 626 (11.31%) are only found in one of the VCFs.

Reading the "Breakdown of VCFs' consistency", a `Group` is a unique key for presence (1) or absence (0) of a call within each of the listed `#Files`. For example: `Group 111` is calls present in all VCFs; `Group 011` is calls present in only the 2nd and 3rd VCFs (i.e. fileB.vcf and fileC.vcf).

We see that `Group 101` has calls belonging to the 1st and 3rd `#Files` (i.e. fileA.vcf and fileC.vcf). This group has a total of 308 calls that intersect, or 5.57%  of all calls in all VCFs. This 308 represents 6.54% of calls in fileA.vcf and 6.31% of calls in fileC.vcf. 

Finally, we see that fileA.vcf has the least amount of calls unique to it on the `Group 100` line.