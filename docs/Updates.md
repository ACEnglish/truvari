# Truvari 5.0
*in progress*

* Reference context sequence comparison is now deprecated and sequence similarity calculation improved by also checking lexicographically minimum rotation's similarity. [details](https://github.com/ACEnglish/truvari/wiki/bench#comparing-sequences-of-variants)
* Symbolic variants (`<DEL>`, `<INV>`, `<DUP>`) can now be resolved for sequence comparison when a `--reference` is provided. The function for resolving the sequences is largely similar to [this discussion](https://github.com/ACEnglish/truvari/discussions/216)
* Symbolic variants can now match to resolved variants, even with `--pctseq 0`, with or without the new sequence resolving procedure.
* Symbolic variant sub-types are ignored e.g. `<DUP:TANDEM> == <DUP>`
* `--sizemax` now default to `-1`, meaning all variant ≥ `--sizemin / --sizefilt` are compared
* Redundant variants which are collapsed into kept (a.k.a. removed) variants now more clearly labeled (`--removed-output` instead of `--collapsed-output`)
* Fixed 'Unknown error' caused by unset TMPDIR ([#229](https://github.com/ACEnglish/truvari/issues/229) and [#245](https://github.com/ACEnglish/truvari/issues/245))
* Fixes to minor record keeping bugs in refine/ga4gh better ensure all variants are counted/preserved
* BND variants are now compared by bench ([details](https://github.com/ACEnglish/truvari/wiki/bench#bnd-comparison))
* Cleaner outputs by not writing matching annotations (e.g. `PctSeqSimilarity`) that are `None`
* Major refactor of Truvari package API for easy reuse of SV comparison functions ([details](https://truvari.readthedocs.io/en/latest/))

# Truvari 4.3.1
*September 9, 2024*
* `bench`
  * Correctly filtering `ALT=*` alleles ([details](https://github.com/ACEnglish/truvari/discussions/219)) and monomorphic reference
    * including test coverage this time
* `stratify`
  * Default behavior is to count variants within ([#221](https://github.com/ACEnglish/truvari/issues/221))
* `collapse`
  * Faster sub-chunking operations by dropping use of pyintervaltree
* `anno chunks`
  * New command for identifying windows with a high number of SVs ([details](https://github.com/ACEnglish/truvari/wiki/Collapse-on-Regions-with-a-High%E2%80%90Density-of-SVs))


# Truvari 4.3.0
*July 31, 2024*
* `refine` & `stratify`
  * Fixed variant and bed boundary overlapping issue 
* general
  * Definition of variants within a region now includes replacement style from TR callers having an anchor base 1bp upstream of catalog/includebed regions
  * Propagating MAFFT errors ([#204](https://github.com/ACEnglish/truvari/issues/204))
  * FIPS compliance ([#205](https://github.com/ACEnglish/truvari/issues/205))
  * Allow csi-indexed vcf files ([#209](https://github.com/ACEnglish/truvari/issues/209))
  * bcftools sort tempfile ([#213](https://github.com/ACEnglish/truvari/issues/213))

# Truvari 4.2.2
*March 28, 2024*

* `collapse`
  * Fewer comparisons needed per-chunk on average
  * Fixed `--chain` functionality ([details](https://github.com/ACEnglish/truvari/issues/196))
  * Fixed `--gt` consolidation of format fields
* `bench`
  * Faster result at the cost of less complete annotations with `--short` flag
* `refine`
  * Assures variants are sequence resolved before incorporating into consensus
  * `bench --passonly --sizemax` parameters are used when building consensus for a region. Useful for `refine --use-original-vcfs`
  * When a refined region has more than 5k variants, it is skipped and a warning of the region is written to the log
  * Flag `--use-region-coords` now expands `--region` coordinates by 100bp (`phab --buffer` default) to allow variants to harmonize out of regions.
* general
  * Dynamic bed/vcf parsing tries to choose faster of streaming/fetching variants
 
# Truvari 4.2.1
*February 6, 2024*
* `collapse`
  * Faster handling of genotype data for `--gt` and `--keep common`
* general
  * Fix to bed end position bug for including variants ([details](https://github.com/ACEnglish/truvari/issues/193))
  * Fix to Dockerfile
* `refine`
  * Changes to `--recount` that accompany the fix to bed end positions.
* New command `ga4gh` to convert Truvari results into GA4GH truth/query VCFs with intermediates tags

# Truvari 4.2
*January 12, 2024*
* `collapse`
  * New parameter `--gt` disallows intra-sample events to collapse ([details](https://github.com/ACEnglish/truvari/wiki/collapse#--gt))
  * New parameter `--intra` for consolidating SAMPLE information during intra-sample collapsing ([details](https://github.com/ACEnglish/truvari/wiki/collapse#--intra))
  * Preserve phasing information when available
  * Faster O(n-1) algorithm instead of O(n^2)
  * Faster sub-chunking strategy makes smaller chunks of variants needing fewer comparisons
  * Fixed rare non-determinism error in cases where multiple variants are at the same position and equal qual/ac could be ordered differently.
* `phab`
  * Correct sample handling with `--bSamples` `--cSamples` parameters
  * Faster generation of consensus sequence
  * Resolved 'overlapping' variant issue causing variants to be dropped
  * New `poa` approach to harmonization. Faster than mafft but less accurate. Slower than wfa but more accurate.
* `bench`
  * New, easier `MatchId` field to track which baseline/comparison variants match up [details](https://github.com/ACEnglish/truvari/wiki/MatchIds)
  * `entry_is_present` method now considers partial missing variants (e.g. `./1`) as present
  * Removed the 'weighted' metrics from `summary.json`
* `consistency`
  * Fixed issue with counting duplicate records
  * Added flag to optionally ignore duplicate records
* `anno svinfo` now overwrites existing SVLEN/SVTYPE info fields
* general
  * Reduced fn matches for unroll sequence similarity by reporting maximum of multiple manipulations of variant sequence (roll up/down/none). Comes at a small, but reasonable, expense of some more fp matches.
  * Bump pysam version
  * Fixed bug in `unroll` sequence similarity that sometimes rolled from the wrong end 
  * Fixed bug for handling of None in ALT field
  * `truvari.compress_index_vcf` forces overwriting of tabix index to prevent annoying crashes


# Truvari 4.1
*August 7, 2023*

* `bench`
  * Creates `candidate.refine.bed` which hooks into `refine` on whole-genome VCFs [details](https://github.com/ACEnglish/truvari/wiki/bench#refining-bench-output)
  * `--recount` for correctly assessing whole-genome refinement results
  * experimental 'weighted' summary metrics [details](https://github.com/ACEnglish/truvari/wiki/bench#weighted-performance)
  * Unresolved SVs (e.g. `ALT == <DEL>`) are filtered when `--pctseq != 0`
* `phab`
  * ~2x faster via reduced IO from operating in stages instead of per-region
  * Removed most external calls (e.g. samtools doesn't need to be in the environment anymore)
  * new `--align wfa` allows much faster (but slightly less accurate) variant harmonization
  * increased determinism of results [detals](https://github.com/ACEnglish/truvari/commit/81a9ab85b91b0c530f9faeedfa4e7e0d68a5e8c2)
* `refine`
  * Faster bed file intersection of `--includebed` and `--regions`
  * Refine pre-flight check
  * Correct refine.regions.txt end position from IntervalTree correction
  * Better refine region selection with `--use-original`
  * `--use-includebed` switched to `--use-region-coords` so that default behavior is to prefer the includebed's coordinates
  * `--use-original-vcfs` to use the original pre-bench VCFs
  * `refine.variant_summary.json` is cleaned of uninformative metrics
* `stratify`
  * parallel parsing of truvari directory to make processing ~4x faster
* `msa2vcf` Fixed REPL decomposition bug to now preserve haplotypes
* `anno grpaf` - expanded annotation info fields
* `anno density` - new parameter `--stepsize` for sliding windows
* `collapse`
  * New optional `--median-info` fields [#146](https://github.com/ACEnglish/truvari/issues/146)
* Minor updates
  * Fix some `anno` threading on macOS [#154](https://github.com/ACEnglish/truvari/issues/154)
  * Monomorphic/multiallelic check fix in `bench`
  * `PHAB_WRITE_MAFFT` environment variable to facilitate updating functional test answer key
  * Slightly slimmer docker container

# Truvari 4.0
*March 13, 2023*

As part of the GIAB TR effort, we have made many changes to Truvari's tooling to enable comparison of variants in TR regions down to 5bp. Additionally, in order to keep Truvari user friendly we have made changes to the UI. Namely, we've updated some default parameters, some command-line arguments, and some outputs. There are also a few new tools and how a couple of tools work has changed. Therefore, we decided to bump to a new major release. If you're using Truvari in any kind of production capacity, be sure to test your pipeline before moving to v4.0.

* New `refine` command for refining benchmarking results. [Details](refine)
* `bench`
  * [Unroll](bench#unroll) is now the default sequence comparison approach. 
  * New `--pick` parameter to control the number of matches a variant can participate in [details](bench#controlling-the-number-of-matches)
  * The `summary.txt` is now named `summary.json`
  * Outputs parameters to `params.json`
  * Output VCFs are sorted, compressed, and indexed
  * Ambiguous use of 'call' in outputs corrected to 'comp' (e.g. `tp-call.vcf.gz` is now `tp-comp.vcf.gz`)
  * Renamed `--pctsim` parameter to `--pctseq`
  * Fixed bug where FP/FN weren't getting the correct, highest scoring match reported
  * Fixed bug where `INFO/Multi` wasn't being properly applied
  * Fixed bug where variants spanning exactly one `--includebed` region were erroneously being counted.
  * Removed parameters: `--giabreport`, `--gtcomp`,`--multimatch`, `--use-lev`, `--prog`, `--unroll`
* `collapse`
  * Renamed `--pctsim` parameter to `--pctseq`
  * Runtime reduction by ~40% with short-circuiting during `Matcher.build_match`
  * Better output sorting which may allow pipelines to be a little faster.
* `vcf2df`
  * More granular sizebins for `[0,50)` including better handling of SNPs
  * `--multisample` is removed. Now automatically add all samples with `--format`
  * key index column removed and replaced by chrom, start, end. Makes rows easier to read and easier to work with e.g. pyranges
* `anno`
  * Simplified ui. Commands that work on a single VCF and can stream (stdin/stdout) no longer use `--input` but a positional argument.
  * Added `addid`
* `consistency`
  * Slight speed improvement
  * Better json output format
* `segment`
  * Added `--passonly` flag
  * Changed UI, including writing to stdout by default
  * Fixed END and 1bp DEL bugs, now adds N to segmented variants' REF, and info fields SVTYPE/SVLEN
* API
  * Began a focused effort on improving re-usability of Truvari code.
  * Entry point to run benchmarking programmatically with [Bench object](https://truvari.readthedocs.io/en/latest/truvari.html#bench).
  * Better development version tracking. [details](https://github.com/ACEnglish/truvari/commit/4bbf8d9a5be3b6a3f935afbd3a9b323811b676a0)
  * Improved developer documentation. See [readthedocs](https://truvari.readthedocs.io/)
* general
  * msa2vcf now left-trims and decomposes variants into indels 
  * Functional tests reorganization
  * Fix for off-by-one errors when using pyintervaltree. See [ticket](https://github.com/ACEnglish/truvari/issues/137)
  * Removed progressbar and Levenshtein dependencies as they are no longer used.

# Truvari 3.5
*August 27, 2022*

* `bench`
  * `--dup-to-ins` flag automatically treats SVTYPE==DUP as INS, which helps compare some programs/benchmarks
  * New `--unroll` sequence comparison method for `bench` and `collapse` ([details](bench#unroll))
* Major `anno trf` refactor (TODO write docs) including:
  * annotation of DEL is fixed (was reporting the ALT copy numbers, not the sample's copy numbers after incorporating the ALT
  * allow 'denovo' annotation by applying any TRF annotations found, not just those with corresponding annotations
* New `anno grpaf` annotates vcf with allele frequency info for groups of samples
* New `phab` for variant harmonization ([details](../phab))
* backend
  * `truvari.entry_size` returns the length of the event in the cases where len(REF) == len(ALT) (e.g. SNPs entry_size is 1)
  * New key utility for `truvari.build_anno_trees`
* general
  * Float metrics written to the VCF (e.g. PctSizeSimilarity) are rounded to precision of 4
  * Nice colors in some `--help` with [rich](https://github.com/Textualize/rich/)
* `divide` 
  * output shards are now more easily sorted (i.e. `ls divide_result/*.vcf.gz` will return the shards in the order they were made)
  * compression/indexing of sub-VCFs in separate threads, reducing runtime
* user issues
  * Monomorphic reference ALT alleles no longer throw an error in `bench` ([#131](https://github.com/ACEnglish/truvari/issues/131))
  * `SVLEN Number=A` fix ([#132](https://github.com/ACEnglish/truvari/issues/132))

# Truvari 3.4
*July 7, 2022*

* Improved performance of `consistency` (see [#127](https://github.com/ACEnglish/truvari/pull/127))
* Added optional json output of `consistency` report
* Allow GT to be missing, which is allowed by VCF format specification
* TRF now uses `truvari.entry_variant_type` instead of trying to use `pysam.VariantRecord.info["SVLEN"]`
directly which allows greater flexibility.
* vcf2df now parses fields with `Number=\d` (e.g. 2+), which is a valid description
* `truvari.seqsim` is now case insensitive (see [#128](https://github.com/ACEnglish/truvari/issues/128))
* Collapse option to skip consolidation of genotype information so kept alleles are unaltered
* `truvari anno dpcnt --present` will only count the depths of non ./. variants
* New collapse annotation `NumConsolidate` records how many FORMATs were consolidated
* Official [conda](https://anaconda.org/bioconda/truvari) support

# Truvari 3.3
*May 25, 2022*

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