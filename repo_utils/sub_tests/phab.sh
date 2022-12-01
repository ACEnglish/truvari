# ------------------------------------------------------------
#                                 phab
# ------------------------------------------------------------
export PATH=$INDIR/external/fake_mafft/:$PATH 

run test_phab $truv phab -k $OD/phab_result \
                         -o $OD/phab_result.vcf.gz \
                         -b $INDIR/variants/phab_base.vcf.gz \
                         -c $INDIR/variants/phab_comp.vcf.gz \
                         -f $INDIR/references/phab_ref.fa \
                         -r chr1:700-900

if [ $test_phab ]; then
    assert_exit_code 0
    run test_phab_result
    assert_equal $(fn_md5 $ANSDIR/phab/phab_result/output.vcf.gz) $(fn_md5 $OD/phab_result/chr1:700-900/output.vcf.gz)
fi

run test_phab_badparams $truv phab -o $OD/phab_result\
    -b doesntexist.vcf \
    -c alsobad.vcf \
    -f noref \
    -r chrP:4802-99292
if [ $test_phab_badparams ]; then
    assert_exit_code 100
fi
