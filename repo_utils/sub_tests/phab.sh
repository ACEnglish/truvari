# ------------------------------------------------------------
#                                 phab - disabled for now (tools)
# ------------------------------------------------------------
#run test_phab $truv phab -o $OD/phab_result \
#        -b $INDIR/phab_base.vcf.gz \
#        -c $INDIR/phab_comp.vcf.gz \
#        -f $INDIR/phab_ref.fa \
#        -r chr1:700-900
#assert_exit_code 0
#
#run test_phab_result
#assert_equal $(fn_md5 $ANSDIR/phab_result.vcf.gz) $(fn_md5 $OD/phab_result/output.vcf.gz)
