## Run the GIAB recommended tandem-repeat stratifications analysis
## Arg1 = the Truvari result directory to operate over
## Arg2 = The repeats annotation file

cd $1
repeats=$2 #../repeats_trfSimplerepLowcomplex_merged50_slop50_gt100.bed.gz
# FN/INS/non-TR
grep 'SVTYPE=INS' fn.vcf | grep 'TRgt100=FALSE' > fn_ins_nonTR.vcf
# FN/INS/TR
grep 'SVTYPE=INS' fn.vcf | grep 'TRgt100=TRUE' > fn_ins_TR.vcf

# FP/INS/non-TR
# grep -v "#" fp.vcf | awk '{if (length($4) - length($5) < 0) print $0}' | bedtools intersect -a stdin \
(grep "#" fp.vcf && grep 'SVTYPE=INS' fp.vcf) | bedtools intersect -a stdin \
            -b $repeats -f 0.2 -v > fp_ins_nonTR.vcf

# FP/INS/TR
# grep -v "#" fp.vcf | awk '{if (length($4) - length($5) < 0) print $0}' | bedtools intersect -a stdin \
(grep "#" fp.vcf && grep 'SVTYPE=INS' fp.vcf) | bedtools intersect -a stdin \
            -b $repeats -f 0.2 -u > fp_ins_TR.vcf

# FN/DEL/non-TR
grep 'SVTYPE=DEL' fn.vcf | grep 'TRgt100=FALSE' > fn_del_nonTR.vcf
# FN/DEL/TR
grep 'SVTYPE=DEL' fn.vcf | grep 'TRgt100=TRUE' > fn_del_TR.vcf
# FP/DEL/non-TR
# grep -v "#" fp.vcf | awk '{if (length($4) - length($5) > 0) print $0}' | bedtools intersect -a stdin \
(grep "#" fp.vcf && grep 'SVTYPE=DEL' fp.vcf) | bedtools intersect -a stdin \
            -b $repeats -f 0.2 -v > fp_del_nonTR.vcf
# FP/DEL/TR
# grep -v "#" fp.vcf | awk '{if (length($4) - length($5) > 0) print $0}' | bedtools intersect -a stdin \
(grep "#" fp.vcf && grep 'SVTYPE=DEL' fp.vcf) | bedtools intersect -a stdin \
            -b $repeats -f 0.2 -u > fp_del_TR.vcf
cd - 
