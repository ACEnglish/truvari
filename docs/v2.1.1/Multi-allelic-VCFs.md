Truvari only runs on comparing the first alternate allele in the base and comparison calls. If VCFs containing multi-allelic sites:

```
chr2     1948201     .     T     TACAACACGTACGATCAGTAGAC,TCAACACACAACACGTACGATCAGTAGAC     ....
```

Pre-process the VCFs with bcftools to get a more complete result:

```bash
bcftools norm -m-any base_calls.vcf.gz | bgzip > base_calls_split.vcf.gz
```