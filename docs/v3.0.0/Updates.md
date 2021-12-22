# Truvari 3.0
*Coming soon*

As Truvari's adoption and functionality grows, we decided to spend time working on sustainability and performance of the tool.  Multiple [Actions](https://github.com/spiralgenetics/truvari/actions) for CI/CD have been added. Many components have been refactored for speed, and other 'cruft' code has been removed.  Some of these changes (particularly the switch to using edlib for sequence similarity) affects the results. Therefore, we've bumped to a new major release version.

* Working on speed improvements
* Added edlib as the default when calculating pctseq_sim, keeping Levenstein as an option (`--use-lev`).
* `truvari bench` summary's gt_precision/gt_recall are replaced by gt_concordance, which is just the percent of TP-comp calls with a concordant genotype. `--no-ref` has better functionality. `--giabreport` is different.
* Added `—keep common` to `truvari collapse`, which allows one to choose to keep the allele with the highest MAC.
* `truvari collapse --hap` wasn't working correctly. The assumptions about the calls being phased wasn't being 
properly used (e.g. don't collapse 1|1) and the NumCollapsed was being populated before the single-best 
match was chosen. The latter is a reporting problem, but the former had an effect on the results with 
~3% of collapsed calls being mis-collapsed.
* `truvari anno trf` is now faster and simpler in its approach and whats reported.. and hopefully more useful.
* `truvari anno grm` has min_size and regions arguments added.
* truv2df has become `truvari vcf2df` where the default is vcf conversion with options to run on a `truvari bench` output directory. It also allows a specific sample to be parsed with `--format` and better Number=A handling.
* NeighId added to `truvari anno numneigh`, which works like bedtools cluster.
* The method af_calc now makes MAC/AC.
* Added 'partial' to `truvari anno remap`.
* Added `truvari anno svinfo`.
* Removed `truvari stats` as `truvari vcf2df` is better and began building [community-driven summaries](https://github.com/spiralgenetics/truvari/discussions/categories/vcf2df-recipes).
* Ubiquitous single version.
* Added a Dockerfile and instructions for making a Truvari docker container.
* Code and repository cleaning.
* Github actions for automated pylint, testing, and releases to pypi.
* Preserving per-version documentation from the wiki in `docs/`.


# Truvari 2.1
*Jan 27, 2021*

We've expanded and improved Truvari's [annotations](https://github.com/spiralgenetics/truvari/wiki/anno). We've added an [SV "collapsing" tool](https://github.com/spiralgenetics/truvari/wiki/collapse). And we've added a way to [turn VCFs into pandas DataFrames](https://github.com/spiralgenetics/truvari/wiki/truv2df) easily for downstream analysis/QC. 

# Truvari 2.0
*May 14, 2020*

After performing a drastic code refactor, we were able to create several helper methods from Truvari's core functionality around SV comparisons and VCF manipulations. This reusable code gave us an opportunity to create tools relevant for SV analysis. 

Truvari now contains multiple subcommands. In addition to the original benchmarking functionality (`truvari bench`), Truvari can generate SV relevant summary statistics, compute consistency of calls within VCFs, and we've begun to develop annotations for SVs. Details on these tools are on the [WIKI](https://github.com/spiralgenetics/truvari/wiki).

We are committed to continually improving Truvari with the hopes of advancing the study and analysis of structural variation.

# Truvari 1.3
*September 25th, 2019*

Truvari has some big changes. In order to keep up with the o deement of Python 2.7 https://pythonclock.org/
We're now only supporting Python 3.

Additionally, we now package Truvari so it and its dependencies can be installed directly. See Installation 
below. This will enable us to refactor the code for easier maintenance and reusability.

Finally, we now automatically report genotype comparisons in the summary stats.