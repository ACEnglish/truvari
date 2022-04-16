truvari anno hompct -i multi_dp.vcf.gz \
        | truvari anno remap -H 10 -r reference.fa \
        | truvari anno gcpct -r reference.fa \
        | truvari anno gtcnt \
        | truvari anno numneigh \
        | truvari anno svinfo \
        | truvari anno dpcnt \
        | truvari anno lcr \
        | vcf-sort \
        | bgzip > hold.vcf.gz
