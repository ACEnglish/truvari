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