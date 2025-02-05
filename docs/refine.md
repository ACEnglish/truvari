As described in the [[phab|phab]] documentation, a constraint on Truvari `bench` finding matches is that there needs to be some consistency in how the variants are represented. To help automate the process of running Truvari `phab` on a benchmarking result and recomputing benchmarking performance on harmonized variants, we present the tool `refine`.

Quick Start
===========

Basic
-----
After making a `bench` result:
```bash
truvari bench -b base.vcf.gz -c comp.vcf.gz -f ref.fa -o result/
```
Use `refine` on the `result/`
```bash
truvari refine result/
```

Optionally, you may automatically hook into refine via `truvari bench --refine ...`.

Output
======
* `refine.base.vcf.gz` - ga4gh refinement result with keys such as FORMAT/BD annotating TP/FN/FP from the baseline variants.
* `refine.comp.vcf.gz - ga4gh refinement result with keys such as FORMAT/BD annotating TP/FN/FP from the comparison variants.
* `refine.variant_summary.json` - result with nearly the same structure as [[summary.json|bench#summaryjson]]
* `refine.regions.txt` - Tab-delimited file with variant counts of refined regions
* `refine.region_summary.json` - Per-region performance metrics
* `phab_bench/` - Bench results on the subset of variants harmonized

To see an example output, look at [test data](https://github.com/ACEnglish/truvari/tree/develop/answer_key/refine/refine_output_three)

Recounting Variants
===================

`refine` will try to turn FN and FP variants into TP variants. It does this by running `phab` to harmonize variant representations and then rerunning `bench`. In order for an originally FN/FP variant to become TP, all variants in the `phab` harmonized region must become TP. Otherwise, variants are left with their original `bench` state.

The option `--write-phab` will ignore the original variant representations and will count/write the `phab` harmonized variant representations. Note that this may cause some variants to split/combine in such a manner that the variant counts will shift. For example, two technical replicate comparison variants benchmarked against the same baseline with `refine --write-phab` my end up with different "base cnt" SVs.

`--align`
=========
By default, Truvari will make the haplotypes and use an external call `mafft` to perform a multiple sequence alignment between them and the reference to harmonize the variants. While this is the most accurate alignment technique, it isn't fast. If you're willing to sacrifice some accuracy for a huge speed increase, you can use `--align wfa`, which also doesn't require an external tool. Another option is `--align poa` which performs a partial order alignment which is faster than mafft but less accurate and slower than wfa but more accurate. However, `poa` appears to be non-deterministic which is not ideal for some benchmarking purposes.

`--use-original-vcfs`
=====================

By default, `refine` will use the base/comparison variants from the `bench` results `tp-base.vcf.gz`, `fn.vcf.gz`, `tp-comp.vcf.gz`, and `fp.vcf.gz` as input for `phab`. However, this contains a filtered subset of variants originally provided to `bench` since it removes variants e.g. below `--sizemin` or not `--passonly`. 

With the `--use-original-vcfs` parameter, all of the original calls from the input vcfs are fetched. This parameter is useful in recovering matches in situations when variants in one call set are split into two variants which are smaller than the minimum size analyzed by `bench`. For example, imagine a base VCF with a 20bp DEL, a comp VCF with two 10bp DEL, and `bench --sizemin 20` was used. `--use-original-vcfs` will consider the two 10bp comp variants during phab harmonization with the 20bp base DEL.

`--regions`
===========

This parameter specifies which regions to re-evaluate. If this is not provided, the original `bench` result's `candidate.refine.bed` is used. 

Note that the larger these regions are the slower MAFFT (used by `phab`) will run. 

`--coords`
==========

The default of, `refine --coords R` will harmonize the bed file coordinates in `--regions` ± `--buffer` base-pairs. If the original `bench --includebed` has coordinates which span the regions that should be harmonized (as is the case in the GIABTR benchmark), then `--coords O` should be used.

`--subset`
==========

By default, refine will still count all variants from the original `bench` results. One may count/report only the variants that are part of the `refine --regions` (default candidate.refine.bed) by specifying `--subset`. This is useful for situations such as the TR benchmarking where the set of regions analyzed by a tool isn't equal to the set of variants inside the benchmark.  In this case, we only want to benchmark the regions that are both in the benchmark and analyzed by the tool, therefore we use `--subset`.

Using `refine.regions.txt`
==========================
| Column            | Description                             |
| ----------------- | --------------------------------------- |
| chrom             | Region's chromosome                     |
| start             | Region's start                          |
| end               | Region's end                            |
| in_tpbase         | Input's True Positive base count        |
| in_tp             | Input's True Positive comparison count  |
| in_fp             | Input's false positive count            |
| in_fn             | Input's false negative count            |
| refined           | Boolean for if region was refined       |
| out_tpbase        | Output's true positive base count       |
| out_tp            | Output's true positive comparison count |
| out_fn            | Outputs false positive count            |
| out_fp            | Output's false negative count           |
| state             | True/False state of the region          |

Note that the `out_*` holds variants counts either from the `phab` harmonized variants when `refined == True` or will be equal to their `in_*` counterpart for regions which were not refined.

Performance by Regions
======================

Because `truvari phab` can alter variant counts during harmonization, one may wish to assess the performance on a per-region basis rather than the per-variant basis.  In the `refine.regions.txt`, a column `state` will have a TP/FN/FP value as defined by the following rules:

```python
false_pos = (data['out_fp'] != 0)
false_neg = (data['out_fn'] != 0)
any_false = false_pos | false_neg

true_positives = (data['out_tp'] != 0) & (data['out_tpbase'] != 0) & ~any_false

true_negatives = (data[['out_tpbase', 'out_tp', 'out_fn', 'out_fp']] == 0).all(axis=1)

baseP = (data['out_tpbase'] != 0) | (data['out_fn'] != 0)
compP = (data['out_tp'] != 0) | (data['out_fp'] != 0)
```

This logic has two edge cases to consider. 1) a region with at least one false-positive and one false-negative will be counted as both a false-positive and a false-negative. 2) Regions within `--refdist` may experience 'variant bleed' where they e.g. have an out_tp, but no other variants because a neighboring region actually contains the the corresponding `out_tpbase`. For the first case, we simply count the region twice and set its state in `refine.regions.txt` to "FP,FN". For the second case, we set the state to 'UNK' and ignore it when calculating the region summary. Future versions may figure out exactly how to handle (prevent?) 'UNK' regions.

These by-region state counts are summarized and written to `refine.region_summary.json`. The definition of metrics inside this json are:
| Key    | Definition                                   | Formula                         |
|--------|----------------------------------------------|---------------------------------|
| TP     | True Positive region count                   |                                 |
| TN     | True Negative region count                   |                                 |
| FP     | False Positive region count                  |                                 |
| FN     | False Negative region count                  |                                 |
| base P | Regions with base variant(s)                 |                                 |
| base N | Regions without base variant(s)              |                                 |
| comp P | Regions with comparison variant(s)           |                                 |
| comp N | Regions without comparison variant(s)        |                                 |
| PPV    | Positive Predictive Value (a.k.a. precision) | TP / comp P                     |
| TPR    | True Positive Rate (a.k.a. recall)           | TP / base P                     |
| TNR    | True Negative Rate (a.k.a. specificity)      | TN / base N                     |
| NPV    | Negative Predictive Value                    | TN / comp N                     |
| ACC    | Accuracy                                     | (TP + TN) / (base P + base N)   |
| BA     | Balanced Accuracy                            | (TPR + TNR) / 2                 |
| F1     | f1 score                                     | 2 * ((PPV * TPR) / (PPV + TPR)) |
| UND    | Regions without an undetermined state        |                                 |

Even though PPV is synonymous with precision, we use these abbreviated names when dealing with per-region performance in order to help users differentiate from the by-variant performance reports.