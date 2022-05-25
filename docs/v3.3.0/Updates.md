# Truvari 3.3
*In Progress*

* New utilities `vcf_ranges` and `make_temp_filename`
* New annotations `dpcnt` and `lcr`
* Fixed a bug in `truvari collapse --keep` that prevented the `maxqual` or `common` options from working
* Increased determinism for `truvari collapse` so that in cases of tied variant position the longer allele is returned. If the alleles also have the same length, they are sorted alphabetically by the REF
* New `truvari bench --extend` functionality. See [discussion](https://github.com/ACEnglish/truvari/discussions/99) for details

# Truvari 3.2
*Apr 1, 2022*

* Removed `truvari.copy_entry` for `pysam.VariantRecord.translate` a 10x faster operation
* Faster `truvari collapse`  ([@c8b319b](https://github.com/ACEnglish/truvari/commit/c8b319b0e717a9e342f52e4a5e927f154eeb0e4a))
* When building `MatchResult` between variants with shared start/end positions, we save processing work by skipping haplotype creation and just compare REFs/ALTs directly.
* Updated documentation to reference the paper https://doi.org/10.1101/2022.02.21.481353
* New `truvari anno density` for identifying regions with 'sparse' and 'dense' overlapping SVs ([details](https://github.com/spiralgenetics/truvari/wiki/anno#truvari-anno-density))
* Better `bench` genotype reporting with `summary.txt` having a `gt_matrix` of Base GT x Comp GT for all Base calls' best, TP match.
* New `truvari anno bpovl` for intersecting against tab-delimited files ([details](https://github.com/spiralgenetics/truvari/wiki/anno#truvari-anno-bpovl))
* New `truvari divide` command to split VCFs into independent parts ([details](https://github.com/ACEnglish/truvari/wiki/divide))
* Replaced `--buffer` parameter with `--minhaplen` for slightly better matching specificity
* Bugfix - `truvari anno trf` no longer duplicates entries spanning multple parallelization regions
* Bugfix - `collapse` MatchId/CollapseId annotation wasn't working
* Bugfixes - from [wwliao](https://github.com/wwliao) ([@4dd9968](https://github.com/ACEnglish/truvari/commit/4dd99683912236f433166889bb0b5667e9fa936d) [@ef2cfb3](https://github.com/ACEnglish/truvari/commit/ef2cfb366b60a5af4671d65d3ed987b08da72227))
* Bugfixes - Issues [#107](https://github.com/ACEnglish/truvari/issues/107), [#108](https://github.com/ACEnglish/truvari/issues/108)

# Truvari 3.1
*Dec 22, 2021*

* `bench` now annotates FPs by working a little differently. See [[bench|bench#methodology]] for details.
* Recalibrated TruScore and new reciprocal overlap measurement for sequence resolved `INS` ([details](https://github.com/spiralgenetics/truvari/discussions/92))
* Match objects are now usable via the SDK. See [#94](https://github.com/spiralgenetics/truvari/discussions/94) for an example of using Truvari programmatically
* `file_zipper` VCF iteration strategy (`GenomeTree` -> `RegionVCFIterator`) that improves speed, particularly when using `--includebed`
* `collapse` refactored to use Match object and for prettier code, cleaner output.
* `anno remap` now optionally adds `INFO` field of the location of the top N hits.
* An experimental tool `truvari segment` added to help SV association analysis.
* `vcf2df` now supports pulling `FORMAT` fields from multiple samples.
* `vcf2df` now adds `('_ref', '_alt')`, or `('_ref', '_het', '_hom')` for `INFO,Number=[R|G]` fields, respectively.
* Improved documentation, including http://truvari.readthedocs.io/ for developers.
* Increasing/diversifying test coverage exposed minor bugs which were fixed.
* `bench --no-ref --cSample` bug fixes.
* Minor usability feature implemented in `help_unknown_cmd`.

# Truvari 3.0
*Sep 15, 2021*

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