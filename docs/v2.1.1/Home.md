(Work in progress.. some links may not work.. Some details may be missing)

# Before you start
VCFs aren't always created with a strong adherence to the format's specification. To make matters worse, many tools don't require input VCFs to strictly adhere to the format. 

Truvari expects input VCFs to be valid so that it will only output valid VCFs. 

We've developed a separate tool that runs multiple validation programs and standard VCF parsing libraries in order to validate a VCF. 

Run [this program](https://github.com/acenglish/usable_vcf) over any VCFs that are giving Truvari trouble. 

Furthermore, Truvari expects 'resolved' SVs (e.g. DEL/INS) and will not interpret BND signals across SVTYPEs (e.g. combining two BND lines to match a DEL call). A brief description Truvari bench methodology is linked below.

# Index:

- [Updates](https://github.com/spiralgenetics/truvari/wiki/Updates)
- Truvari Commands:
  - [bench](https://github.com/spiralgenetics/truvari/wiki/bench)
    - [The multimatch parameter](https://github.com/spiralgenetics/truvari/wiki/The--multimatch-parameter)
    - [Edit Distance Ratio vs Sequence Similarity](https://github.com/spiralgenetics/truvari/wiki/Edit-Distance-Ratio-vs-Sequence-Similarity)
    - [Multi-allelic VCFs](https://github.com/spiralgenetics/truvari/wiki/Multi-allelic-VCFs)
    - [stratifications.sh](https://github.com/spiralgenetics/truvari/wiki/stratifications.sh)
    - [Comparing two SV programs](https://github.com/spiralgenetics/truvari/wiki/Comparing-two-SV-programs)
  - [consistency](https://github.com/spiralgenetics/truvari/wiki/consistency)
  - [anno](https://github.com/spiralgenetics/truvari/wiki/anno)
    - [Creating an annotation source](https://github.com/spiralgenetics/truvari/wiki/making-annotations)
  - [collapse](https://github.com/spiralgenetics/truvari/wiki/collapse)
  - [vcf2df](https://github.com/spiralgenetics/truvari/wiki/vcf2df)