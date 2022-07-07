The wiki holds documentation most relevant for develop. For information on a specific version of Truvari, see [`docs/`](https://github.com/spiralgenetics/truvari/tree/develop/docs)

Truvari: Refined Structural Variant Comparison Preserves Allelic Diversity  
doi: https://doi.org/10.1101/2022.02.21.481353

# Before you start
VCFs aren't always created with a strong adherence to the format's specification. 

Truvari expects input VCFs to be valid so that it will only output valid VCFs. 

We've developed a separate tool that runs multiple validation programs and standard VCF parsing libraries in order to validate a VCF. 

Run [this program](https://github.com/acenglish/usable_vcf) over any VCFs that are giving Truvari trouble. 

Furthermore, Truvari expects 'resolved' SVs (e.g. DEL/INS) and will not interpret BND signals across SVTYPEs (e.g. combining two BND lines to match a DEL call). A brief description of Truvari bench methodology is linked below.

# Index

- [[Updates|Updates]]
- [[Installation|Installation]]
- Truvari Commands:
  - [[bench|bench]]
    - [[The multimatch parameter|The--multimatch-parameter]]
    - [[Edit Distance Ratio vs Sequence Similarity|Edit-Distance-Ratio-vs-Sequence-Similarity]]
    - [[Multi-allelic VCFs|Multi-allelic-VCFs]]
    - [[stratifications.sh|GIAB-stratifications.sh]]
    - [[Comparing two SV programs|Comparing-two-SV-programs]]
  - [[consistency|consistency]]
  - [[anno|anno]]
  - [[collapse|collapse]]
  - [[vcf2df|vcf2df]]
  - [[segment|segment]]
- [[Development|Development]]
- [[Citations|Citations]]
- [[Resources|Resources]]