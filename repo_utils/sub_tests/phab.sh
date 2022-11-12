# ------------------------------------------------------------
#                                 phab
# ------------------------------------------------------------
export PATH=$INDIR/external/fake_mafft/:$PATH 
run test_phab $truv phab -o $OD/phab_result \
        -b $INDIR/variants/phab_base.vcf.gz \
        -c $INDIR/variants/phab_comp.vcf.gz \
        -f $INDIR/references/phab_ref.fa \
        -r chr1:700-900
assert_exit_code 0

if [ $test_rebench ]; then
    assert_exit code 0
    run test_phab_result
    assert_equal $(fn_md5 $ANSDIR/phab/phab_result/output.vcf.gz) $(fn_md5 $OD/phab_result/output.vcf.gz)
fi
