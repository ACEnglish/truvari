```
████████╗██████╗ ██╗   ██╗██╗   ██╗ █████╗ ██████╗ ██╗
╚══██╔══╝██╔══██╗██║   ██║██║   ██║██╔══██╗██╔══██╗██║
   ██║   ██████╔╝██║   ██║██║   ██║███████║██████╔╝██║
   ██║   ██╔══██╗██║   ██║╚██╗ ██╔╝██╔══██║██╔══██╗██║
   ██║   ██║  ██║╚██████╔╝ ╚████╔╝ ██║  ██║██║  ██║██║
   ╚═╝   ╚═╝  ╚═╝ ╚═════╝   ╚═══╝  ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝
```

Structural variant comparison tool for VCFs

Given benchmark and comparsion sets of SVs, calculate the recall, precision, and f-measure.

[Spiral Genetics](https:www.spiralgenetics.com)

[Motivation](https://docs.google.com/presentation/d/17mvC1XOpOm7khAbZwF3SgtG2Rl4M9Mro37yF2nN7GhE/edit)

UPDATES
=======

Truvari has some big changes. In order to keep up with the retirement of Python 2.7 https://pythonclock.org/
We're now only supporting Python 3.

Additionally, we now package Truvari so it and its dependencies can be installed directly. See Installation 
below. This will enable us to refactor the code for easier maintenance and reusability.

Finally, we now automatically report genotype comparisons in the summary stats.

Installation
============

Truvari uses Python 3.7 and can be installed with pip:

  $ pip install Truvari 


Quick start
===========

  $ truvari -b base_calls.vcf -c compare_calls.vcf -o output_dir/

Outputs
=======

  * tp-call.vcf -- annotated true positive calls from the COMP
  * tp-base.vcf -- anotated true positive calls form the BASE
  * fn.vcf -- false negative calls from BASE
  * fp.vcf -- false positive calls from COMP
  * base-filter.vcf -- size filtered calls from BASE
  * call-filter.vcf -- size filtered calls from COMP
  * summary.txt -- json output of performance stats
  * log.txt -- run log
  * giab_report.txt -- (optional) Summary of GIAB benchmark calls. See "Using the GIAB Report" below.

summary.txt
===========

The following stats are generated for benchmarking your call set.
<table><tr><th>Metric</th><th>Definition</th>
<tr><td>TP-base</td><td>Number of matching calls from the base vcf</td></tr>
<tr><td>TP-call</td><td>Number of matching calls from the comp vcf</td></tr>
<tr><td>FP</td><td>Number of non-matching calls from the comp vcf</td></tr>
<tr><td>FN</td><td>Number of non-matching calls from the base vcf</td></tr>
<tr><td>precision</td><td>TP-call / (TP-call +  FP)</td></tr>
<tr><td>recall</td><td>TP-base / (TP-base + FN)</td></tr>
<tr><td>f1</td><td>(recall * precision) / (recall + precision)</td></tr>
<tr><td>base cnt</td><td>Number of calls in the base vcf</td></tr>
<tr><td>call cnt</td><td>Number of calls in the comp vcf</td></tr>
<tr><td>base size filtered</td><td>Number of base vcf calls outside of (sizemin, sizemax)</td></tr>
<tr><td>call size filtered</td><td>Number of comp vcf calls outside of (sizemin, sizemax)</td></tr>
<tr><td>base gt filtered</td><td>Number of base calls not passing the no-ref parameter filter</td></tr>
<tr><td>call gt filtered</td><td>Number of comp calls not passing the no-ref parameter filter</td></tr>
<tr><td>TP-call_TP-gt</td><td>TP-call's with genotype match</td></tr>
<tr><td>TP-call_FP-gt</td><td>TP-call's without genotype match</td></tr>
<tr><td>TP-base_TP-gt</td><td>TP-base's with genotype match</td></tr>
<tr><td>TP-base_FP-gt</td><td>TP-base's without genotype match</td></tr>
<tr><td>gt_precision</td><td>TP-call_TP-gt / (TP-call_TP-gt + FP + TP-call_FP-gt)</td></tr>
<tr><td>gt_recall</td><td>TP-base_TP-gt / (TP-base_TP-gt / FN)</td></tr>
<tr><td>gt_f1</td><td>(gt_recall * gt_precision) / (gt_recall + gt_precision)</td></tr>
</table>

Methodology
===========

```
Input:
    BaseCall - Benchmark TruthSet of SVs
    CompCalls - Comparison SVs from another program
Build IntervalTree of CompCalls
For each BaseCall:
  Fetch CompCalls overlapping within *refdist*. 
    If typematch and LevDistRatio >= *pctsim* \
    and SizeRatio >= *pctsize* and PctRecOvl >= *pctovl*: 
      Add CompCall to list of Neighbors
  Sort list of Neighbors by TruScore ((2*sim + 1*size + 1*ovl) / 3.0)
  Take CompCall with highest TruScore and BaseCall as TPs
  Only use a CompCall once if not --multimatch
  If no neighbors: BaseCall is FN
For each CompCall:
  If not used: mark as FP
```

Matching Parameters
--------------------
<table><tr><th>Parameter</th><th>Default</th><th>Definition</th>
<tr><td>refdist</td><td>500</td>
<td>Maximum distance comparison calls must be within from base call's start/end</td></tr>
<tr><td>pctsim</td><td>0.7</td>
<td>Levenshtein distance ratio between the REF/ALT haplotype sequences of base and comparison call.
See "Comparing Haplotype Sequences of Variants" below.
<tr><td>pctsize</td><td>0.7</td>
<td>Ratio of min(base_size, comp_size)/max(base_size, comp_size)</td></tr>
<tr><td>pctovl</td><td>0.0</td>
<td>Ratio of two calls' (overlapping bases)/(longest span)</td></tr>
<tr><td>typeignore</td><td>False</td>
<td>Types don't need to match to compare calls.</td>
</table>

Comparing VCFs without sequence resolved calls
----------------------------------------------

If the base or comp vcfs do not have sequence resolved calls, simply set `--pctsim=0` to turn off
sequence comparison.

Difference between --sizemin and --sizefilt
-------------------------------------------

`--sizemin` is the minimum size of a base call or comparison call to be considered.  

`--sizefilt` is the minimum size of a call to be added into the IntervalTree for searching. It should
be less than `sizemin` for edge case variants.

For example: `sizemin` is 50 and `sizefilt` is 30. A 50bp base call is 98% similar to a 49bp call at 
the same position.

These two calls should be considered matching. If we instead removed calls less than `sizemin`, we'd
incorrectly classify the 50bp base call as a false negative.

This does have the side effect of artificially inflating specificity. If that same 49bp call in the
above were below the similarity threshold, it would not be classified as a FP due to the `sizemin`
threshold. So we're giving the call a better chance to be useful and less chance to be detrimental
to final statistics.

Definition of annotations added to TP vcfs
--------------------------------------------
<table>
<tr><th>Anno                   </th><th>Definition</th></tr>
<tr><td>TruScore	       </td><td>Truvari score for similarity of match. `((2*sim + 1*size + 1*ovl) / 3.0)`</td></tr>
<tr><td>PctSeqSimilarity       </td><td>Pct sequence similarity between this variant and its closest match</td></tr>
<tr><td>PctSizeSimilarity      </td><td>Pct size similarity between this variant and it's closest match</td></tr>
<tr><td>PctRecOverlap          </td><td>Percent reciprocal overlap of the two calls' coordinates</td></tr>
<tr><td>StartDistance          </td><td>Distance of this call's start from matching  call's start</td></tr>
<tr><td>EndDistance            </td><td>Distance of this call's end from matching  call's end</td></tr>
<tr><td>SizeDiff               </td><td>Difference in size(basecall) and size(compcall)</td></tr>
<tr><td>NumNeighbors           </td><td>Number of comparison calls that were in the neighborhood (REFDIST) of the base call</td></tr>
<tr><td>NumThresholdNeighbors  </td><td>Number of comparison calls that passed threshold matching of the base call</td></tr>
</table>

NumNeighbors and NumThresholdNeighbors are also added to the FN vcf.

Using the GIAB Report
---------------------

When running against the GIAB SV benchmark (link below), you can create a detailed report of 
calls summarized by the GIAB VCF's SVTYPE, SVLEN, Technology, and Repeat annotations.

To create this report.

1. Run truvari with the flag `--giabreport`.
2. In your output directory, you will find a file named `giab_report.txt`.
3. Next, make a copy of the 
[Truvari Report Template Google Sheet](https://docs.google.com/spreadsheets/d/1T3EdpyLO1Kq-bJ8SDatqJ5nP_wwFKCrH0qhxorvTVd4/edit?usp=sharing).
4. Finally, paste ALL of the information inside `giab_report.txt` into the "RawData" tab. Be careful not 
to alter the report text in any way. If successul, the "Formatted" tab you will have a fully formated report.

While Truvari can use other benchmark sets, this formatted report currently only works with GIAB SV v0.5 and v0.6. Work
will need to be done to ensure Truvari can parse future GIAB SV releases.

<a href="ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/">GIAB v0.6 Download Link</a>

Include Bed & VCF Header Contigs 
--------------------------------

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


Comparing Haplotype Sequences of Variants
---------------------------------------

While many SV overlapping approaches match changes within regions with metrics such as reciprocal overlap and
start/end distances, few account for the possibility to represent the same sequence change in multiple ways.

For example, consider a tandem repeat expansion of the reference sequence 'AB'. Here, we'll represent the 'insertion'
sequence as lower-case 'ab', though it should be considered equivalent to 'AB'. Three equally valid descriptions of this 
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
between variants in order to allow for more flexibilty in complex multi-allelic representations.

More Information
----------------

Find more details and discussions about Truvari on the [WIKI page](https://github.com/spiralgenetics/truvari/wiki).



<a href="https://www.spiralgenetics.com" rel="Spiral Genetics" style="width:400px;">![Spiral Genetics](http://static1.squarespace.com/static/5a81ef7629f187c795c973c3/t/5a986ab453450a17fc3003e8/1533115866026/?format=1500w)</a>
