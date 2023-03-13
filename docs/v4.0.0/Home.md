The wiki holds documentation most relevant for develop. For information on a specific version of Truvari, see [`docs/`](https://github.com/spiralgenetics/truvari/tree/develop/docs)

Citation:  
English, A.C., Menon, V.K., Gibbs, R.A. et al. Truvari: refined structural variant comparison preserves allelic diversity. Genome Biol 23, 271 (2022). https://doi.org/10.1186/s13059-022-02840-6

# Before you start
VCFs aren't always created with a strong adherence to the format's specification. 

Truvari expects input VCFs to be valid so that it will only output valid VCFs. 

We've developed a separate tool that runs multiple validation programs and standard VCF parsing libraries in order to validate a VCF. 

Run [this program](https://github.com/acenglish/usable_vcf) over any VCFs that are giving Truvari trouble. 

Furthermore, Truvari expects 'resolved' SVs (e.g. DEL/INS) and will not interpret BND signals across SVTYPEs (e.g. combining two BND lines to match a DEL call). A brief description of Truvari bench methodology is linked below.

Finally, Truvari does not handle multi-allelic VCF entries and as of v4.0 will throw an error if multi-allelics are encountered. Please use `bcftools norm` to split multi-allelic entries. 

# Index

- [[Updates|Updates]]
- [[Installation|Installation]]
- Truvari Commands:
  - [[anno|anno]]
  - [[bench|bench]]
  - [[collapse|collapse]]
  - [[consistency|consistency]]
  - [[divide|divide]]
  - [[phab|phab]]
  - [[refine|refine]]
  - [[segment|segment]]
  - [[stratify|stratify]]
  - [[vcf2df|vcf2df]]
- [[Development|Development]]
- [[Citations|Citations]]
- [[Resources|Resources]]