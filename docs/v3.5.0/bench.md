
Quick start
===========
Run this command where base is your 'truth set' SVs and comp is the comparison set of SVs.
```bash
truvari bench -b base_calls.vcf -c comp_calls.vcf -f reference.fa -o output_dir/
```

Outputs
=======
Truvari bench creates an output directory with the following files.
<table><tr><th>File</th><th>Description</th>
<tr><td>tp-call.vcf</td><td>Annotated true positive calls from the comparison VCF</td></tr>
<tr><td>tp-base.vcf</td><td>Annotated true positive calls form the base VCF</td></tr>
<tr><td>fp.vcf</td><td>False positive calls from comparison</td></tr>
<tr><td>fn.vcf</td><td>False negative calls from base</td></tr>
<tr><td>summary.txt</td><td>Json output of performance stats</td></tr>
<tr><td>log.txt</td><td>Run's logging output</td></tr>
<tr><td>giab_report.txt</td><td>(optional) Summary of GIAB benchmark calls. See "Using the GIAB Report" below.</td></tr>
<tr><td>giab.jl</td><td>(optional) DataFrame of benchmark results.</td></tr>
</table>

summary.txt
===========
The following stats are generated for benchmarking your call set and written to summary.txt as a json dump.

<table><tr><th>Metric</th>  <th>Definition</th>
<tr><td>TP-base</td>        <td>Number of matching calls from the base vcf</td></tr>
<tr><td>TP-call</td>        <td>Number of matching calls from the comp vcf</td></tr>
<tr><td>FP</td>             <td>Number of non-matching calls from the comp vcf</td></tr>
<tr><td>FN</td>             <td>Number of non-matching calls from the base vcf</td></tr>
<tr><td>precision</td>      <td>TP-call / (TP-call + FP)</td></tr>
<tr><td>recall</td>         <td>TP-base / (TP-base + FN)</td></tr>
<tr><td>f1</td>             <td>2 * ((recall * precision) / (recall + precision))</td></tr>
<tr><td>base cnt</td>       <td>Number of calls in the base vcf</td></tr>
<tr><td>call cnt</td>       <td>Number of calls in the comp vcf</td></tr>
<tr><td>TP-call_TP-gt</td>  <td>TP-call with genotype match</td></tr>
<tr><td>TP-call_FP-gt</td>  <td>TP-call without genotype match</td></tr>
<tr><td>TP-base_TP-gt</td>  <td>TP-base with genotype match</td></tr>
<tr><td>TP-base_FP-gt</td>  <td>TP-base without genotype match</td></tr>
<tr><td>gt_concordance</td> <td>TP-call_TP-gt / (TP-call_TP-gt + TP-call_FP-gt)</td></tr>
<tr><td>gt_matrix</td>      <td>Base GT x Comp GT Matrix of all Base calls' best, TP match</td></tr>
</table>

Matching Parameters
===================
Picking matching parameters can be more of an art than a science. It really depends on the precision of your callers and the tolerance you wish to allow them such that it is a fair comparison.

For example, depth of coverage callers (such as CNVnator) will have very 'fuzzy' boundaries, and don't report the exact deleted sequence but only varying regions. So thresholds of `pctsim=0`, `pctsize=.5`, `pctovl=.5`, `refdist=1000` may seem fair.

[BioGraph](https://github.com/spiralgenetics/biograph) and many long-read callers report precise breakpoints and full alternate allele sequences. When benchmarking those results, we want to ensure our accuracy by using the stricter default thresholds.

If you're still having trouble picking thresholds, it may be beneficial to do a few runs of Truvari bench over different values. Start with the strict defaults and gradually increase the leniency. From there, you can look at the performance metrics and manually inspect differences between the runs to find out what level you find acceptable. Truvari is meant to be flexible for comparison. More importantly, Truvari helps one clearly report the thresholds used for reproducibility (see the json at the top of your log.txt).

Here is a rundown of each matching parameter.
<table><tr><th>Parameter</th><th>Default</th><th>Definition</th>
<tr><td>refdist</td>         <td>500</td>    <td>Maximum distance comparison calls must be within from base call's start/end</td></tr>
<tr><td>pctsim</td>          <td>0.7</td>    <td>Edit distance ratio between the REF/ALT haplotype sequences of base and comparison call. See "Comparing Sequences of Variants" below.</td></tr>
<tr><td>unroll</td>          <td>False</td>  <td>Use the unroll technique for sequence comparison. See "Comparing Sequences of Variants" below.</td></tr>
<tr><td>pctsize</td>         <td>0.7</td>    <td>Ratio of min(base_size, comp_size)/max(base_size, comp_size)</td></tr>
<tr><td>pctovl</td>          <td>0.0</td>    <td>Ratio of two calls' (overlapping bases)/(longest span)</td></tr>
<tr><td>typeignore</td>      <td>False</td>  <td>Types don't need to match to compare calls.</td></tr>
</table>

Matching Parameter Diagrams
===========================
Below are matching parameter diagrams to illustrate (approximately) how they work.

```
 █ = Deletion ^ = Insertion

--refdist REFDIST (500)
  Max reference location distance

    ACTGATCATGAT
     |--████--|    
          █████      
  
  Calls are within reference distance of 2

--pctsize PCTSIZE (0.7)
  Min pct allele size similarity

    ACTGATCATGA    sizes
      █████     -> 5bp
        ████    -> 4bp

  variants have 0.8 size similarity


--pctovl PCTOVL (0.0)
  Min pct reciprocal overlap

    ACTGATCATGA  ranges
      █████      [2,7)
        ████     [4,8)

  variants have 0.6 reciprocial overlap


--pctsim PCTSIM (0.7)
  Min percent allele sequence similarity

    A-CTG-ACTG
     ^   ^       haplotypes
     |   └ACTG -> CTGACTGA
     └CTGA     -> CTGACTGA

  haplotypes have 100% sequence similarity
```

Comparing VCFs without sequence resolved calls
==============================================

If the base or comp vcfs do not have sequence resolved calls, simply set `--pctsim=0` to turn off
sequence comparison. The `--reference` does not need to be set when not using sequence comparison.

Difference between --sizemin and --sizefilt
===========================================

`--sizemin` is the minimum size of a base call to be considered.  

`--sizefilt` is the minimum size of a comparison call that will be matched to base calls. It should
be less than `sizemin` for edge case variants.

For example: `sizemin` is 50 and `sizefilt` is 30. A 50bp base call is 98% similar to a 49bp comparison
call at the same position.

These two calls should be considered matching. If we instead removed comparison calls less than `sizemin`, 
we'd incorrectly classify the 50bp base call as a false negative.

This does have the side effect of artificially inflating specificity. If that same 49bp call in the
above were below the similarity threshold, it would not be classified as a FP due to the `sizemin`
threshold. So we're giving the call a better chance to be useful and less chance to be detrimental
to final statistics.

Definition of annotations added to vcfs
==========================================
<table>
<tr><th>Anno                   </th><th>Definition</th></tr>
<tr><td>TruScore	       </td><td>Truvari score for similarity of match. `((pctsim + pctsize + pctovl) / 3 * 100)`</td></tr>
<tr><td>PctSeqSimilarity       </td><td>Pct sequence similarity between this variant and its closest match</td></tr>
<tr><td>PctSizeSimilarity      </td><td>Pct size similarity between this variant and its closest match</td></tr>
<tr><td>PctRecOverlap          </td><td>Percent reciprocal overlap percent of the two calls</td></tr>
<tr><td>StartDistance          </td><td>Distance of the base call's start from comparison call's start</td></tr>
<tr><td>EndDistance            </td><td>Distance of the base call's end from comparison call's end</td></tr>
<tr><td>SizeDiff               </td><td>Difference in size of base and comp calls</td></tr>
<tr><td>GTMatch                </td><td>Base/comp calls' Genotypes match</td></tr>
<tr><td>MatchId                </td><td>Id to help tie base/comp calls together {chunkid}.{baseid}.{compid}</td></tr>
<tr><td>Multi                  </td><td>Call is false due to non-multimatching</td></tr>
</table>

See [[MatchIds wiki|MatchIds]] for details.

Using the GIAB Report
=====================

When running against the GIAB SV benchmark (link below), you can automatically add a text `giab_report.txt` to the output directory containing calls summarized by the GIAB VCF's SVTYPE, SVLEN, Technology, and Repeat annotations. This also adds a DataFrame `giab.jl` (similar to running `truvari vcf2df --bench-dir`).

See the [Recipes Discussions](https://github.com/spiralgenetics/truvari/discussions/categories/vcf2df-recipes) on how make interesting summaries and plots.

While Truvari can use other benchmark sets, this formatted report currently only works with GIAB SV v0.5 and v0.6. Work
will need to be done to ensure Truvari can parse future GIAB SV releases.

[GIAB v0.6 Download Link](https://bit.ly/2SlWYbB)

ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/

Include Bed & VCF Header Contigs 
================================

If an `--includebed` is provided, only base and comp calls contained within the defined regions are used 
for comparison. This is similar to pre-filtering your base/comp calls using:

```bash
(zgrep "#" my_calls.vcf.gz && bedtools intersect -u -a my_calls.vcf.gz -b include.bed) | bgzip > filtered.vcf.gz
```

with the exception that Truvari requires the start and the end to be contained in the same includebed region 
whereas `bedtools intersect` does not.

If an `--includebed` is not provided, the comparison is restricted to only the contigs present in the base VCF
header. Therefore, any comparison calls on contigs not in the base calls will not be counted toward summary 
statistics and will not be present in any output vcfs.

Extending an Include Bed
========================

The option `--extend` extends the regions of interest (set in `--includebed` argument) by the given number of bases on each side, allowing truth set variants to match call set variants that are just outside of the original region. If a call set variant is in the extended regions it can potentially match a truth set variant that is in the original regions turning it to TP, but call set variants in the extended regions that don't have a match are not counted as FP. This strategy is similar to the one implemented for size matching - only the truth set variants longer than sizemin (equal to 50 by default) are considered, but they are allowed to match shorter comparison set variants of size sizefilt (30 by default) or longer.

See this [discussion](https://github.com/ACEnglish/truvari/discussions/99)for details.

Comparing Sequences of Variants
=========================================

While many SV overlapping approaches match variants with metrics such as reciprocal overlap and
start/end distances, few account for the possibility to represent the same sequence change in multiple ways.

Truvari has implemented two approaches to compare variant sequences. The original method is to use the reference context of a pair of variants . The new method is called 'unroll'.

## Reference context
For the reference context method, consider a hypothetical tandem repeat expansion of the reference sequence 'AB'. Here, we'll represent the 'insertion' sequence as lower-case 'ab', though it should be considered equivalent to 'AB'. Three equally valid descriptions of this 
variant would be:

```text
#POS INS  Haplotype
  0  ab   abAB
  1  ba   AbaB
  2  ab   ABab
```

Therefore, to compare the sequence similarity, Truvari builds the haplotypes over the range of a pair of calls'
`min(starts):max(ends)` before making the the sequence change introduced by the variants. In python, this line 
looks like:

``` python
hap1_seq = ref.get_seq(a1_chrom, start + 1, a1_start).seq + a1_seq + ref.get_seq(a1_chrom, a1_end + 1, end).seq
```

Where `a1_seq1` is the longer of the REF or ALT allele.

Future versions of Truvari may take this concept further to reconstruct haplotypes based on phasing information
between variants in order to allow for more flexibility in complex multi-allelic representations.


## Unroll
The method of giving a pair of calls the same reference context can be achieved using an 'unroll' method. For a formal description, see [this gist](https://gist.github.com/ACEnglish/1e7421c46ee10c71bee4c03982e5df6c). 

The main idea is that in order to move variants upstream/downstream, the reference sequence flanking the variant will need to be moved downstream/upstream respectively. Or, to say this another way, we can think of the alternate sequences as being circular instead of linear. This means that in order to move the variant e.g. 1bp downstream for an INS, we could remove the first base from the ALT and append it to the end. So in our 'ab' example above, we only need to unroll the insertion at a position by the distance between it and another variant e.g. the INS `ab` at POS 2 becomes identical to the INS `ba` at POS 1 by rolling `2-1 = 1` bases from the start to the end.

This unroll method has a number of benefits and a couple of considerations.
* Unroll does not need a `--reference` for comparison, which saves I/O time.
* In the limited testing performed so far, unroll increase the number of matching SVs
* Unroll seems to decrease the number of 'suspect' matches in smaller size regimes
* The pattern between PctSizeSimilarity and PctSeqSimilarity is simpler
* The PctSeqSimilarity of variant pairs reported by the reference context method is higher than the unroll method. This is mainly due to reference context 'inflating' the similarity (i.e. adding identical sequence to both variants in the pair naturally increases their identity).  However, a majority of pairs seem to have identical similarity reported by either method. 


Methodology
===========
Here is a high-level pseudocode description of the steps Truvari bench conducts to compare the two VCFs.
```
* zip the Base and Comp calls together in sorted order
* create chunks of all calls overlapping within ±`--chunksize` basepairs
* make a |BaseCall| x |CompCall| match matrix for each chunk
* build a Match for each call pair in the chunk - annotate as TP if >= all thresholds 
* if the chunk has no Base or Comp calls
** return them all as FNs/FPs 
* if multimatch:
** return the max of each row (base) and of each column (call) from the match matrix
* else:
** ravel and sort the match matrix
** return the first use of each call in its current state (True or False)
** return the second use of each call with state=False and set Multi=True
```
![](https://github.com/spiralgenetics/truvari/blob/develop/imgs/TruvariBenchMethod.png)

--dup-to-ins
============

Most SV benchmarks only report DEL and INS SVTYPEs. The flag `--dup-to-ins` will interpret SVs with SVTYPE == DUP to SVTYPE == INS. Note that DUPs generally aren't sequence resolved (i.e. the ALT isn't a sequence) like INS. Therefore, `--dup-to-ins` typically should be used without sequence comparison via `--pctsim 0`