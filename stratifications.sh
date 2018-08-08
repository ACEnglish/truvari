## Run the GIAB recommended tandem-repeat stratifications analysis
## Arg1 = the Truvari result directory to operate over
## Arg2 = The repeats annotation file

repeats=`readlink -f $2`
cd $1
# FN/INS/non-TR
(grep "#" fn.vcf && (grep 'SVTYPE=INS' fn.vcf | grep 'TRgt100=FALSE')) | vcf-sort | bgzip > fn_ins_nonTR.vcf.gz
tabix -p vcf fn_ins_nonTR.vcf.gz
# FN/INS/TR
(grep "#" fn.vcf && (grep 'SVTYPE=INS' fn.vcf | grep 'TRgt100=TRUE')) | vcf-sort | bgzip > fn_ins_TR.vcf.gz
tabix -p vcf fn_ins_TR.vcf.gz

# FP/INS/non-TR
# grep -v "#" fp.vcf | awk '{if (length($4) - length($5) < 0) print $0}' | bedtools intersect -a stdin \
(grep "#" fp.vcf && (grep "#" fp.vcf && grep 'SVTYPE=INS' fp.vcf) | bedtools intersect -a stdin \
            -b $repeats -f 0.2 -v) | vcf-sort | bgzip > fp_ins_nonTR.vcf.gz
tabix -p vcf fp_ins_nonTR.vcf.gz

# FP/INS/TR
# grep -v "#" fp.vcf | awk '{if (length($4) - length($5) < 0) print $0}' | bedtools intersect -a stdin \
(grep "#" fp.vcf && (grep "#" fp.vcf && grep 'SVTYPE=INS' fp.vcf) | bedtools intersect -a stdin \
            -b $repeats -f 0.2 -u) | vcf-sort | bgzip > fp_ins_TR.vcf.gz
tabix -p vcf fp_ins_TR.vcf.gz

# FN/DEL/non-TR
(grep "#" fn.vcf && (grep 'SVTYPE=DEL' fn.vcf | grep 'TRgt100=FALSE')) | vcf-sort | bgzip > fn_del_nonTR.vcf.gz
tabix -p vcf fn_del_nonTR.vcf.gz
# FN/DEL/TR
(grep "#" fn.vcf && (grep 'SVTYPE=DEL' fn.vcf | grep 'TRgt100=TRUE')) | vcf-sort | bgzip > fn_del_TR.vcf.gz
tabix -p vcf fn_del_TR.vcf.gz
# FP/DEL/non-TR
# grep -v "#" fp.vcf | awk '{if (length($4) - length($5) > 0) print $0}' | bedtools intersect -a stdin \
(grep "#" fp.vcf && (grep "#" fp.vcf && grep 'SVTYPE=DEL' fp.vcf) | bedtools intersect -a stdin \
            -b $repeats -f 0.2 -v) | vcf-sort | bgzip > fp_del_nonTR.vcf.gz
tabix -p vcf fp_del_nonTR.vcf.gz
# FP/DEL/TR
# grep -v "#" fp.vcf | awk '{if (length($4) - length($5) > 0) print $0}' | bedtools intersect -a stdin \
(grep "#" fp.vcf && (grep "#" fp.vcf && grep 'SVTYPE=DEL' fp.vcf) | bedtools intersect -a stdin \
            -b $repeats -f 0.2 -u) | vcf-sort | bgzip > fp_del_TR.vcf.gz
tabix -p vcf fp_del_TR.vcf.gz
cd - 
