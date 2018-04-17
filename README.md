<center>

```
████████╗██████╗ ██╗   ██╗██╗   ██╗ █████╗ ██████╗ ██╗
╚══██╔══╝██╔══██╗██║   ██║██║   ██║██╔══██╗██╔══██╗██║
   ██║   ██████╔╝██║   ██║██║   ██║███████║██████╔╝██║
   ██║   ██╔══██╗██║   ██║╚██╗ ██╔╝██╔══██║██╔══██╗██║
   ██║   ██║  ██║╚██████╔╝ ╚████╔╝ ██║  ██║██║  ██║██║
   ╚═╝   ╚═╝  ╚═╝ ╚═════╝   ╚═══╝  ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝
```

</center>

Structural variant caller comparison tool for VCF

Given a benchmark and callset, calculate the sensitivity/specificity/f-measure.

Spiral Genetics, 2018


Installation
============

Truvari requires the following modules:

  $ pip install pyvcf Levenshtein swalign intervaltree progressbar2


Quick start
===========

  $ ./truvari.py -b base_calls.vcf -c compare_calls.vcf -o output_dir/


Outputs
=======

  * tp-call.vcf -- annotated true positive calls from the CALLS
  * tp-base.vcf -- anotated true positive calls form the BASE
  * fn.vcf -- false negative calls from BASE
  * fp.vcf -- false positive calls from CALLS
  * base-filter.vcf -- size filtered calls from BASE
  * call-filter.vcf -- size filtered calls from CALLS
  * summary.txt -- json output of performance stats
  * log.txt -- run log


Methodology
===========

```
Input:
    BaseCall - Benchmark TruthSet of SVs
    CompCalls - Comparison SVs from another program
For each BaseCall, fetch CompCalls overlapping within ±500 buffer. 
If LevDistRatio >= 0.7 and SizeRatio >= 0.7:
    Match BaseCall with best CompCall as TP
Only use a CompCall once. 
```


Comparing VCFs without sequence resolved calls
----------------------------------------------

If the base or call vcfs do not have sequence resolved calls, simply set `--pctsim=0` to turn off
sequence comparison.

Definition of annotations added to TP vcfs
--------------------------------------------
<table>
<tr><th>Anno                   </th><th>Definition</th></tr>
<tr><td>PctSimilarity          </td><td>Pct sequence similarity between this variant and its closest match</td></tr>
<tr><td>StartDistance          </td><td>Distance of this call's start from comparison call's start</td></tr>
<tr><td>EndDistance            </td><td>Distance of this call's end from comparison call's end</td></tr>
<tr><td>PctRecOverlap          </td><td>Percent reciprocal overlap of the two calls' coordinates</td></tr>
<tr><td>SizeDiff               </td><td>Difference in size(basecall) and size(evalcall)</td></tr>
<tr><td>NumNeighbors           </td><td>Number of calls in B that were in the neighborhood (REFDIST) of this call</td></tr>
<tr><td>NumThresholdNeighbors  </td><td>Number of calls in B that are within threshold distances of this call</td></tr>
</table>

Difference between --sizemin and --sizefilt
-------------------------------------------

`--sizemin` is the minimum size of a base to be considered a TP. Those removed are counted as base
size filtered. It is also the minimum size a call to be considered a FP. Those removed are counted
as call size filtered.

`--sizefilt` is the minimum size of a call to be added into the IntervalTree for searching. It should
be less than sizemin for edge case variants around sizemin.

For example: sizemin is 50 and sizefilt is 30. A 50bp base is 98% similar to a 49bp call at the
same position.

These two calls should be considered identical. If we instead removed calls less than sizemin, we'd
incorrectly classify the base as a false negative.

This does have the side effect of artificially inflating specificity. If that same 49bp call in the
above were below the similarity threshold, it would not be classified as a FP due to the sizemin
threshold. So we're giving the call a better chance to be useful and less chance to be detrimental
to final statistics. SVs are fuzzy like this.


Using the GIAB Report
---------------------

When running against the 
[GIAB SV v0.5 benchmark](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_UnionSVs_12122017/), 
you can create a detailed report of calls summarized by the GIAB VCF's SVTYPE, SVLEN, Technology, 
and Repeat annotations.

To do this, run Truvari with the flag `--giabreport`.

In your output directory, you will find a file named `giab_report.txt`.

Next, make a copy of the 
[Truvari Report Template Google Sheet](https://docs.google.com/spreadsheets/d/1T3EdpyLO1Kq-bJ8SDatqJ5nP_wwFKCrH0qhxorvTVd4/edit?usp=sharing).

Finally, paste ALL of the information inside `giab_report.txt` into the "RawData" tab. Be careful not 
to alter the report text in any way. If successul, the "Formatted" tab you will have a fully formated report.

This currently only works with GIAB SV v0.5. Work will need to be done to ensure Truvari can parse future 
GIAB SV releases.
