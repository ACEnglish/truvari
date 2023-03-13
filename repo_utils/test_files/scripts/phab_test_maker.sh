# record keeping for making phab tests
vcf_region=chr1:26399065-26401053
fa_region=chr1:26399066-26401053
offset=26399065

bcftools view -r $vcf_region ~/scratch/code/adotto/variants/data/GRCh38.variants.sqoff.vcf.gz \
    | python translate.py $offset | bgzip > phab_base.vcf.gz
tabix phab_base.vcf.gz

bcftools view -r $vcf_region ~/scratch/code/adotto/variants/benchmarking/truthsets/HPRC-cur.20211005-align2-GRCh38.dip.vcf.gz \
    | python translate.py $offset | bgzip > phab_comp.vcf.gz
tabix phab_comp.vcf.gz
samtools faidx ~/scratch/insertion_ref/msru/data/reference/grch38/GRCh38_1kg_mainchrs.fa $fa_region \
    | sed 's/>chr1.*/>chr1/' > phab_ref.fa

truvari phab -o m_test \
    -b phab_base.vcf.gz \
    -c phab_comp.vcf.gz \
    -f phab_ref.fa \
    -r chr1:700-900


