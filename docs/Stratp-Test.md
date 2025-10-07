The StratP test identifies variant stratifications that significantly impact a tool’s performance. The test evaluates whether a caller’s accuracy (e.g., recall) is dependent on a stratification by computing a directional enrichment score using expected counts under independence. P-values of the significance of this enrichment is then calculated via permutation-based resampling. This is necessary because variant features can co-occur non-randomly. For example, tandem repeats may more frequently contain insertions than deletions, creating imbalanced subgroups that could bias naive stratification analyses. By evaluating enrichment across groups defined by combinations of feature values and assessing significance through permutation, StratP helps mitigate such confounding effects without requiring explicit modeling of feature dependencies. The resulting p-values provide a robust, interpretable measure for prioritizing stratifications that are impactful to a tool.

The simplest way to start using the StratP test is to compare to SV benchmarks with pre-set stratifications. Currently available presets are from GIAB v1.1 SVs and the SMaHT MIMS SV benchmarks.

```bash
truvari stratp --preset mims truvari_output/
```

Output
------

The columns in the stratp report are:

| Column   | Definition                                                         |
| -------- | ------------------------------------------------------------------ |
| feature  | Stratification feature (e.g. `svtype`)                             |
| value    | Stratification value (e.g. `DEL`)                                  |
| acc      | Accuracy of calls in this stratification (e.g. `TP / (TP+FN)`)     |
| obs      | Number of variants observed in this stratification                 |
| delta    | Observed mean ranking minus permutation tests' mean ranking        |
| q1       | 1% Quantile of permutation tests                                   |
| q50      | 50% Quantile of permutation tests                                  |
| q99      | 99% Quantile of permutation tests                                  |
| pval     | Raw permutation test p-value                                       |
| adj_pval | Benjamini/Hochberg multiple test correction adjusted pvalue        |
| reject   | True or False based on rejection of null hypothesis at given alpha |

#### How to interpret stratp_test output
As an example, consider the below, truncated StratP test report

| feature   | value           | acc    | obs   | delta     | q1       | q50     | q99     | pval   | adj_pval | reject |
| --------- | --------------- | ------ | ----- | --------- | -------- | ------- | ------- | ------ | -------- | ------ |
| VAF_bin   | SomaticLow(<5%) | 0.140  | 13184 | -111.7586 | -22.653  | -0.0418 | 22.6734 | 0.0    | 0.0      | True   |
| isolated  | False           | 0.391  | 12401 | -49.437   | -20.6976 | 0.0522  | 21.1852 | 0.0    | 0.0      | True   |
| svtype    | INS             | 0.4985 | 20545 | -18.8122  | -21.0539 | -0.0179 | 21.2513 | 0.0182 | 0.1213   | False  |
| szbin     | [200,300)       | 0.5299 | 2546  | 2.3111    | -31.8288 | -0.0405 | 31.6662 | 0.5651 | 1.0      | False  |
| svtype    | DEL             | 0.6264 | 13595 | 18.8122   | -21.1796 | -0.0359 | 21.0722 | 0.9806 | 1.0      | False  |
| VAF_bin   | GermlineHom     | 0.9446 | 5308  | 72.259    | -24.9471 | -0.0244 | 24.533  | 1.0    | 1.0      | False  |

This report tells us that from top-to-bottom (lowest to highest p-value), we have ranked features on which the caller's
accuracy is most dependent. For example, the top hit of `VAF_bin SomaticLow(<5%)` indicates that the
average accuracy of variants with this feature is `0.140`, which is significantly lower than the average performance of
the caller across other values in this stratification. Furthermore, because the multiple-test correction adjusted
p-value (`adj_pval`) is still below our alpha (`0.01`), we can reject the null hypothesis that the caller's performance is
independent of this stratification.

While few stratifications will reach the level of significance, we can still use the sorted p-values to roughly assess a
callers performance on all stratifications. For example, the third lowest p-value (`svtype INS, p=0.1213`) has
an `acc=0.4985`, which is lower than other stratifications with higher p-values (e.g. `svtype DEL, p=0.9806, acc=0.62`).
This ordering tells is that this example caller performs better on DEL than INS, and performs really well on GermlineHom
SVs.


**NOTE!** A stratification being found as significant to a result's accuracy is not strictly-speaking an indication of the tool's ability to perform on said stratification. For example, one SV caller (toolA) may find that e.g. DEL are not significantly reducing its recall while for toolB DELs do significantly reduce its recall. However, toolA may have a recall of 0.50 on DELs while toolB has 0.80 recall on DELs. This is possible if toolA has an average recall of 0.50 across all stratifications while toolB averages 0.95 recall. Stratp performs its comparison of stratification feature:values per-result and compares each stratification to each tool's average recall. Therefore, the reasonable interpretation of a stratp significant stratification is that the particular stratification is an outlier relative to a tool's overall accuracy. These outliers are then places where tool developers can focus their attention.

# Custom `--features`

StratP starts by converting a truvari output directory into a pandas DataFrame via `truvari vcf2df -i -f`. Any of the columns that are inside of that DataFrame can be turned into a stratification for the StratP test to run on. However, the columns MUST be categorical. The simplest stratifications are custom fields made by the `vcf2df` of `svtype,szbin`. If a column is specified that has a continuous variable, each unique value will be treated as its own categorical value. StratP makes it easy to bin continuous values into categorical values. For example, say you want to make stratifications of the field INFO/TruScore. There are three ways to bin that continuous value.

1. `--feature TruScore:w3`

This will look for the minimum/maximum values of TruScore and then make 3 equal-width bins across the domain with `pd.cut(df['TruScore`], bins=3)`. The value 3 can be any integer that corresponds to the number of bins made.

2. `--feature TruScore:f5` 

This will look at the distribution of TruScore values and use `pd.qcut(df['TruScore'], q=5)` to make bins that have equal frequency of values in each. The value 5 can be any integer that corresponds to the number of bins made.

3. `--feature TruScore:c0-70-80-95-100`

This makes custom bins around the `-` delimited integers with `pd.cut(df['TruScore`], bins=[0, 70, 80, 95, 100])` The values for the bins can be any set of integers or floats.

Multiple features are specified by a comma-separated list e.g. `--feature svtype,szbin,TruScore:f4`. 

There are a couple of considerations about picking your features. First, some feature/values (stratifications) may be highly redundant. For example, if we run `truvari anno trf` and `truvari anno remap`, the INFO/TRF flag and INFO/REMAP=tandem, are likely to almost always the same. You can run with `--check-co 0.50` to do a chi-squared test and then a phi measurement (essentially a correlation coefficient) to get warnings anytime phi is ≥0.50. Note, that the `--check-co` requires `scipy` to be installed, which isn't done automatically by truvari.

The second consideration is how many stratifications you build. If there are too few stratifications, it will be difficult to perform the permutation test. For example, if you only provide `--feature svtype,szbin`, there's only DEL,INS and about 8 size bins. When testing all the groups that have DEL compared to all the groups that aren't DEL, you can only have the 8 INS size bins to compare against the 8 DEL size bins. This is less than the `--min-obs`, so stratp won't perform the permutation test and the output p-values will be blank. However, if you provide too many feature/value combinations, groups will become too small to perform the observed/expected initial scoring/ranking of groups. StratP will warn and drop all groups with fewer than `--min-obs` SVs.

## --tail

Left tail corresponds to testing for stratifications that reduce accuracy whereas right tail tests for stratifications that increase accuracy relative to the overall average accuracy.

## --states

By default, strap will check the states for baseline SVs to identify the true/false variants per-group. This means it will only analyze the baseline SVs the the reported `acc` reflects recall. For `--states comp`, it will only analyze the comparison SVs and the reported `acc` reflects precision. 

