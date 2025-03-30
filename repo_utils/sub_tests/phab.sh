# ------------------------------------------------------------
#                                 phab
# ------------------------------------------------------------
run phab $truv phab -o $OD/phab_result.vcf.gz \
                    -r chr1:700-900 --align mafft \
                    --no-dedup \
                    -f $INDIR/references/phab_ref.fa \
                    $INDIR/variants/phab_base.vcf.gz \
                    $INDIR/variants/phab_comp.vcf.gz

if [ $phab ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/phab/phab_result.vcf.gz) $(fn_md5 $OD/phab_result.vcf.gz)
fi

run phab_wfa $truv phab -o $OD/phab_result_wfa.vcf.gz \
                        -r chr1:700-900 --align wfa \
                        -f $INDIR/references/phab_ref.fa \
                        $INDIR/variants/phab_base.vcf.gz \
                        $INDIR/variants/phab_comp.vcf.gz
if [ $phab_wfa ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/phab/phab_result_wfa.vcf.gz) $(fn_md5 $OD/phab_result_wfa.vcf.gz)
fi

run phab_poa $truv phab -o $OD/phab_result_poa.vcf.gz \
                        -r chr1:700-900 --align poa \
                        -f $INDIR/references/phab_ref.fa \
                        --dedup $INDIR/variants/phab_base.vcf.gz \
                        $INDIR/variants/phab_comp.vcf.gz
if [ $phab_poa ]; then
    assert_exit_code 0
    # poa isn't deterministic
    #assert_equal $(fn_md5 $ANSDIR/phab/phab_result_poa.vcf.gz) $(fn_md5 $OD/phab_result_poa.vcf.gz)
fi


run phab_badparams $truv phab -o $INDIR/variants/filter.vcf -f noref -r chrP:4802-99292 \
                              doesntexist.vcf $INDIR/variants/filter.vcf -s A,A
if [ $phab_badparams ]; then
    assert_exit_code 100
fi
