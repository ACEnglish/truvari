# ------------------------------------------------------------
#                                 phab
# ------------------------------------------------------------
run test_phab $truv phab -o $OD/phab_result.vcf.gz \
                         -b $INDIR/variants/phab_base.vcf.gz \
                         -c $INDIR/variants/phab_comp.vcf.gz \
                         -f $INDIR/references/phab_ref.fa \
                         -r chr1:700-900 --align mafft

if [ $test_phab ]; then
    dump_logs
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/phab/phab_result.vcf.gz) $(fn_md5 $OD/phab_result.vcf.gz)
fi

run test_phab_wfa $truv phab -o $OD/phab_result_wfa.vcf.gz \
                         -b $INDIR/variants/phab_base.vcf.gz \
                         -c $INDIR/variants/phab_comp.vcf.gz \
                         -f $INDIR/references/phab_ref.fa \
                         -r chr1:700-900 --align wfa
if [ $test_phab_wfa ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/phab/phab_result_wfa.vcf.gz) $(fn_md5 $OD/phab_result_wfa.vcf.gz)
fi

run test_phab_badparams $truv phab -o $OD/phab_result\
    -b doesntexist.vcf \
    -c alsobad.vcf \
    -f noref \
    -r chrP:4802-99292
if [ $test_phab_badparams ]; then
    assert_exit_code 100
fi
