
In addition to looking at performance of a single set of variation against a baseline, one may wish to measure the consistency between multiple sets of variation. The tool `truvari consistency` can automatically create that result. 

Running
=======

```
usage: consistency [-h] [-d] [-j] [-o OUTPUT] VCFs [VCFs ...]

Over multiple vcfs, calculate their intersection/consistency.

Calls will match between VCFs if they have a matching key of:
    CHROM POS ID REF ALT

positional arguments:
  VCFs                  VCFs to intersect

optional arguments:
  -h, --help            show this help message and exit
  -d, --no-dups         Disallow duplicate SVs
  -j, --json            Output report in json format
  -o OUTPUT, --output OUTPUT
                        Write tsv of variant keys and their flag
```

Example
=======

```bash
truvari consistency fileA.vcf fileB.vcf fileC.vcf
```

Matching Entries
================

VCF entries will be considered matching if and only if they have an exact same key of `CHROM POS ID REF ALT`. Because of this stringency, it is recommend that you compare the tp-base.vcf or fn.vcf results from each individual VCF's Truvari output. The key and flags can be written with the `--output` option.

Duplicates
==========

If there are VCFs with duplicate keys, they are handled by appending a e.g. `.1` to the key. If you'd like to ignore duplicates, just add `--no-dups`

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

Json
====
Below is a consistency report in json format.
```json
{
    "vcfs": [
        "repo_utils/test_files/variants/input1.vcf.gz",
        "repo_utils/test_files/variants/input2.vcf.gz",
        "repo_utils/test_files/variants/input3.vcf.gz"
    ],
    "total_calls": 3513,
    "num_vcfs": 3,
    "vcf_counts": {
        "repo_utils/test_files/variants/input1.vcf.gz": 2151,
        "repo_utils/test_files/variants/input2.vcf.gz": 1783,
        "repo_utils/test_files/variants/input3.vcf.gz": 2065
    },
    "shared": [
        {
            "vcf_count": 3,
            "num_calls": 701,
            "call_pct": 0.1995445488186735
        },
        {
            "vcf_count": 2,
            "num_calls": 1084,
            "call_pct": 0.3085681753487048
        },
        {
            "vcf_count": 1,
            "num_calls": 1728,
            "call_pct": 0.4918872758326217
        }
    ],
    "detailed": [
        {
            "group": "111",
            "total": 701,
            "total_pct": 0.1995445488186735,
            "repo_utils/test_files/variants/input1.vcf.gz": 0.32589493258949326,
            "repo_utils/test_files/variants/input2.vcf.gz": 0.393157599551318,
            "repo_utils/test_files/variants/input3.vcf.gz": 0.3394673123486683
        },
        {
            "group": "001",
            "total": 645,
            "total_pct": 0.18360375747224594,
            "repo_utils/test_files/variants/input1.vcf.gz": 0,
            "repo_utils/test_files/variants/input2.vcf.gz": 0,
            "repo_utils/test_files/variants/input3.vcf.gz": 0.31234866828087166
        },
        {
            "group": "100",
            "total": 598,
            "total_pct": 0.17022487902077996,
            "repo_utils/test_files/variants/input1.vcf.gz": 0.2780102278010228,
            "repo_utils/test_files/variants/input2.vcf.gz": 0,
            "repo_utils/test_files/variants/input3.vcf.gz": 0
        },
        {
            "group": "101",
            "total": 487,
            "total_pct": 0.1386279533162539,
            "repo_utils/test_files/variants/input1.vcf.gz": 0.22640632264063226,
            "repo_utils/test_files/variants/input2.vcf.gz": 0,
            "repo_utils/test_files/variants/input3.vcf.gz": 0.2358353510895884
        },
        {
            "group": "010",
            "total": 485,
            "total_pct": 0.13805863933959578,
            "repo_utils/test_files/variants/input1.vcf.gz": 0,
            "repo_utils/test_files/variants/input2.vcf.gz": 0.27201346045989905,
            "repo_utils/test_files/variants/input3.vcf.gz": 0
        },
        {
            "group": "110",
            "total": 365,
            "total_pct": 0.10389980074010817,
            "repo_utils/test_files/variants/input1.vcf.gz": 0.1696885169688517,
            "repo_utils/test_files/variants/input2.vcf.gz": 0.2047111609646663,
            "repo_utils/test_files/variants/input3.vcf.gz": 0
        },
        {
            "group": "011",
            "total": 232,
            "total_pct": 0.06604042129234272,
            "repo_utils/test_files/variants/input1.vcf.gz": 0,
            "repo_utils/test_files/variants/input2.vcf.gz": 0.13011777902411667,
            "repo_utils/test_files/variants/input3.vcf.gz": 0.11234866828087167
        }
    ]
}
```

Output TSV
==========

When given a `--output`, a TSV is written which holds the variant key used for comparison (first 5 columns of the VCF entry), the FLAG corresponding to the 'group' in the main output (e.g. FLAG 3 == group 0011), and a COUNT of how many VCFs the key was found in.