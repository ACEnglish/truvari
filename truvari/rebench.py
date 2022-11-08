"""
A simple Truvari bench limitation comes about with its use of `--includebed`. Currently, we expect callers to produce
calls over every locus covered by the 'tier' bed file. However, this isn't always the case, particularally when
considering tandem repeat callers with a pre-defined set of regions over which it attempts to discover variants. If
there are regions inside the `--includebed` that aren't in the caller's pre-defined set, the caller will be penalized
with FNs.

A user could solve this problem by intersecting their `--includebed` with their caller's pre-defined set. But this may
create ambiguities when users compare their results to other reported results. The user would have to run once with the
full tier list to see what percent of TRs they're getting overall and run again to get the precision/recall of sites
they're actually trying to find. While this isn't an intractable problem, there is an opportunity for a more elegant
solution.

Truvari bench has a major limitation that it expects variants to have consistent representations.
However, depending on the alignment parameters, there can be effects such as split variants. For example, a truthset may
report an insertion as a single event whereas a caller can place the event over multiple VCF entries.

There are a couple of ways around this. One possibility is to perform variant harmonization which is already provided by
`truvari phab`. By re/co-aligning the variants from the base and comparison VCFs together, they get more consistent
representations. This idea of realignment with identical parameters has been demonstrated with similar techniques.
Another possibility is something like what's performed by rtg (I think) and hap-eval by senteion. These approaches
effectively achieve the same goal by reconstructing haplotypes created by groups of variants and then directly comparing
the base VCF's haplotypes with those produced by the comparison VCFs variants.

In the interest of usability and consistency and wanting to standardize things, we'll address these two limitations by
building a new tool `truvari rebench`.

This tool takes as input a `truvari bench` output directory and attempts to reclassify true/false calls with a more
comprehensive comparison.

First, it will load the first result's `--includebed` and intersect it with a tool's pre-defined set (name TBD..
something like `--subset`) and exclude any FNs covered by `--subset`.

Second, for each region from the union of `--subset` and `--includebed` it counts the number of TP/FN/FP
per-region. If there are any discrepancies which would indicate the FNs/FPs would benefit from variant harmonization,
the region is reprocessed.

This identification of discrepancies is currently magic. As a start, rebench will simply check if there are FN > 0 and
FP > 0 in a region. I think this can be refined improved, but this first conditional should be a strong enough start.

As for performing the variant harmonization, we have two options. The first will just be to feed the region into truvari
phab. The second will be to recreate the haplotype construction/comparison of hap-eval. 

Each of these methods has limitations. First, phab is slow and it assumes that variants are phased. This isn't always
true.

Second, hap-eval *may* have two different problems. 
The first is that it may over simplify regions in that all variants in a region become True/False together. Imagine
there are two distinct SVs in a region produced by a caller and one SV inside the truth-set. And let's say one of the
caller's SV is an exact match to the truth-set SV and that this variant is large. The second SV from the caller is small
and a FP. Because the first variant is large and 100% identical, the penalty form the mismatched bases of the second
call could be 'washed out'. To say that a different way, a 1000bp deletion in the base/comp VCF when considered with a
50bp FP gives a similarity of 95%, which is over the (default) threshold and the FP is misclassified as a TP.

Furthermore,
hap-eval's approach of aligning the haplotypes created by variants over an entire region (if that's what it's doing) may
be succeptable to incorrectly matching variants. For example, imagine an insertion of 50bp at position X and an
insertion of 50bp at position X+maxdist. The haplotypes each creates over the region automatically shared maxdist
base-pairs. If we assume maxdist == 1000 (the default) the WORST that their sequence similarity can be is `1 - (50 /
1000)` or 95%.

These two problems need to be addressed. But let's assume we come up with some solution to these problems and we have
harmonized variants that we can then compare. Truvari rebench will the update the FN/FP numbers appropriately.

We finish up by simply creating new files inside of the output directory. For now, I'll assume that we're making all the
same files (e.g. `summary.json` and `tp-base.vcf` etc), but we'll just prefix the filenames with `rebench.` so they're
identifiable. Therefore, we'll have two of each file inside a truvari benchmark result (e.g. `summary.json` and
`rebench.summary.json`)

After we get this stuff done. There are other features we'll be able to add into the tool. We can allow a
`--stratifications` bed file which has information such as TR repeat region score (an idea we're pursuing for the TR
region annotation file). As we're going through all the variants and regions already, we can build output files that
give informative reports such as precision/recall as a function of region score. I don't want to over-engineer the
reports here (yet). Instead, let's just think about what the output file would look like so that people can take it and
do whatever they need to downstream. Maybe something as simple as - stratifications is a bed file with a 4th column
having some annotation (TR motif length, let's say). We then output `rebench.strats.bed` with three more columns added of
TP/FN/FP counts. That simple tsv can easily be processed by other tools/scripts/notebooks.

Need to builds
--------------
- truvari bench --no-compress to optionally skip compression/indexing of the output VCFs. Will help with fetching in
  rebench
- save a config.json in truvari bench output (helps usability/reproducibility)
- rename summary.txt to summary.json

User workflow:
--------------

1. Run truvari bench on your calls and the TR specific benchmark. Use the TR regions as your includebed
2. Run truvari rebench on the output directory as well as a few more parameters.

tr-bench tries to harmonize variants in TR regions where - for example - there's a high number of FN and FP. These
regions may have more matching than regular truvari bench is reporting and we can therefore reduce the number of FN/FP

tr-bench also allows you to refine the summaries based on regions. In step 1 of the workflow, we ask the user to run
using the TR region bed which contains all the tandem repeats. However, many tr callers are built on pre-defined regions
which are a subset of the TR regions. Therefore, the callers shouldn't be penalized for FNs on regions over which they
didn't even try to discover variants.

Remaing questions:
------------------
- Percent of variant bases over the regions post-phab over threshold - If that's 
- Do we want an 'unresolved' mode - Like, we should end up with copy numbers so we can just compare the copy numbers
  directly? But that's pretty much the same as the lengths and stuff... Ugh. This would just be a different tool
  entirely

Inputs:
    - Truvari bench output directory
        Note: We're going to assume that truvari bench was run with the TR catalog. That will have annotations we'll
        want to reuse in tr-bench
    - phab parameters
    - tr-bench parameters:
        - Bed file of regions
            Is this all TR regions or just regions where the caller tries to discover?
            Probably the latter because we want to allow people to correct their recall so it isn't penalized for
            regions they don't even attempt to call
        - Number of regions to sample
            If there are a bunch of candidate phab regions, probably don't want to process all of them, so
            we'll allow a random selection of sites. You could then 'estimate' based on those percents how many might
            become updated
        - Comparison thesholds
            How many variants or what percent of the sequence needs to be matching for the sites to be considered
            equal
        - ???

Outputs:
    We'll write results to the truvari bench output directory
    - phab.* mainly phab.tp-base phab.tp-call phab.summary.txt
    - Also, let's refactor the summary.txt to summary.json, that's always annoyed me
    - For each of the bed files of regions, we can report how many base/comp variants there are per-region before/after phab

** Question 1: are there any parameters I should try to pull from the bench? IDK.
    But if there are I'll need to make a `confg` that's saved into truvari bench directories

Operations:
    1. Given the bed file - Do counts/evaluations on a per-region basis.
    2. Figure out a way to select which regions should be phab'd
    3. Re-count the variants in those regions to see if they should have their state change
        - one way change of false becoming true.
    4. Update the summary.txt with the variants that became true `phab.summary.txt`
    5. Need to figure out how to report regions..
    6. Eventually would think about ROC curves, maybe based on the score

Random ideas:
    - So, MSA is slow. But I think you can get more similar results if you just align haplotypes with consistent
      alignment parameters... So maybe I should build a phab mode that just runs edlib on the haplotypes and makes new
      VCFs. 
"""
