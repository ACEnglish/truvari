
  _______      __      __        _ 
 |__   __|     \ \    / /       (_)
    | |_ __ _   \ \  / /_ _ _ __ _ 
    | | '__| | | \ \/ / _` | '__| |
    | | |  | |_| |\  / (_| | |  | |
    |_|_|   \__,_| \/ \__,_|_|  |_|
                                   
                                   

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

Input:
    BaseCall - Benchmark TruthSet of SVs
    CompCalls - Comparison SVs from another program
For each BaseCall, fetch CompCalls overlapping within ±500 buffer. 
If LevDistRatio >= 0.7 and SizeRatio >= 0.7:
    Match BaseCall with best CompCall as TP
Only use a CompCall once. 


Comparing VCFs without sequence resolved calls
----------------------------------------------

If the base or call vcfs do not have sequence resolved calls, simply set `--pctsim=0` to turn of
sequence comparison.

Definition of annotations** added to TP vcfs
--------------------------------------------

PctSimilarity          "Pct sequence similarity between this variant and its closest match"
StartDistance          "Distance of this call's start from comparison call's start"
EndDistance            "Distance of this call's end from comparison call's end"
PctRecOverlap          "Percent reciprocal overlap of the two calls' coordinates"
SizeDiff               "Difference in size(basecall) and size(evalcall)"
NumNeighbors           "Number of calls in B that were in the neighborhood (REFDIST) of this call"
NumThresholdNeighbors  "Number of calls in B that are within threshold distances of this call"


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
