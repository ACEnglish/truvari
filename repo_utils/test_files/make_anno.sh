truvari anno hompct -i multi.vcf.gz \
        | truvari anno remap -r reference.fa \
        | truvari anno gcpct -r reference.fa \
        | truvari anno gtcnt \
        | truvari anno numneigh \
        | bgzip > hold.vcf.gz
