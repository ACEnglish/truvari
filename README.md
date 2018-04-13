    __                              _
   / /________  ___   ______ ______(_)
  / __/ ___/ / / / | / / __ `/ ___/ /
 / /_/ /  / /_/ /| |/ / /_/ / /  / /
 \__/_/   \__,_/ |___/\__,_/_/  /_/


Structural variant caller comparison tool for VCF

Given a benchmark and callset, calculate the sensitivity/specificity/f-measure.

Spiral Genetics, 2018


Installation
============

Truvari requires the following modules:

  $ pip install edlib pyvcf progressbar2 intervaltree


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

Everything in all caps is a parameter.

For varA in BASE greater than SIZEMIN:
    For varB in CALLS within varA.start - REFDIST to varA.end + REFDIST greater than SIZEMIN
        if already matched varB:
            continue
        if similarity* >= PCTSIM
            add to candidates
        if more than 0 candidates:
            matchVar is highest similarity score of candidates
            add to TP-base and TP-call
            add annotations** to varA and matchVar
            mark matchVar as already matched
        else:
            add to FN

For every varB in CALLS greater than SIZEMIN:
    if varB not already match:
        add to FP


How similarity* is calculated
-----------------------------

TYPEMATCH is a flag to force same type (default False)

if TYPEMATCH and (len(entryA.REF) < len(entryB.ALT)) != (len(entryA.REF) < len(entryB.ALT)):
    return False
if ref and alt sequences are identical:
    return True

ref_dist is Levenshtein distance of the two BASE/CALL reference allele sequences
alt_dist is Levenshtein distance of the two BASE/CALL alternate allele sequences
max_ref is max length of two ref allele sequences
max_alt is max length of two alt allele sequences

similarity = 1 - ((ref_dist + alt_dist) / (max_ref + max_alt)


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

--sizemin is the minimum size of a base to be considered a TP. Those removed are counted as base
size filtered. It is also the minimum size a call to be considered a FP. Those removed are counted
as call size filtered.

--sizefilt is the minimum size of a call to be added into the IntervalTree for searching. It should
be less than sizemin for edge case variants around sizemin.

For example: sizemin is 50 and sizefilt is 30. A 50bp base is 98% similar to a 49bp call at the
same position.

These two calls should be considered identical. If we instead removed calls less than sizemin, we'd
incorrectly classify the base as a false negative.

This does have the side effect of artificially inflating specificity. If that same 49bp call in the
above were below the similarity threshold, it would not be classified as a FP due to the sizemin
threshold. So we're giving the call a better chance to be useful and less chance to be detrimental
to final statistics. SVs are fuzzy like this.
