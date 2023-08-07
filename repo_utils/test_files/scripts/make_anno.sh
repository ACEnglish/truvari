in_file=$1
ref=$2
truvari anno hompct -i $in_file \
        | truvari anno remap -H 10 -r $ref \
        | truvari anno gcpct -r $ref \
        | truvari anno gtcnt \
        | truvari anno numneigh \
        | truvari anno svinfo \
        | truvari anno dpcnt \
        | truvari anno lcr \
        | bcftools sort -O z -o hold.vcf.gz
