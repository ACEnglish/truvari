truvari anno hompct -i multi.vcf.gz \
        | truvari anno remap -H 10 -r reference.fa \
        | truvari anno gcpct -r reference.fa \
        | truvari anno gtcnt \
        | truvari anno numneigh \
        | truvari anno svinfo \
        | bgzip > hold.vcf.gz
