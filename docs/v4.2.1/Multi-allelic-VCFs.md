Truvari only compares the first alternate allele in VCFs. If a VCF contains multi-allelic sites such as:

```
chr2     1948201     .     T     TACAACACGTACGATCAGTAGAC,TCAACACACAACACGTACGATCAGTAGAC     ....
```

Then pre-process the VCFs with bcftools:

```bash
bcftools norm -m-any base_calls.vcf.gz | bgzip > base_calls_split.vcf.gz
```