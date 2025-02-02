# ------------------------------------------------------------
#                                 ga4gh
# ------------------------------------------------------------
run ga4gh $truv ga4gh -i $ANSDIR/refine/refine_output_three/ -o $OD/ga4gh_norefine
if [ $ga4gh ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/ga4gh/ga4gh_norefine_truth.vcf.gz) $(fn_md5 $OD/ga4gh_norefine_truth.vcf.gz)
    assert_equal $(fn_md5 $ANSDIR/ga4gh/ga4gh_norefine_query.vcf.gz) $(fn_md5 $OD/ga4gh_norefine_query.vcf.gz)
fi

run ga4gh_refine $truv ga4gh -b 0 -w -i $ANSDIR/refine/refine_output_three/ -o $OD/ga4gh_withrefine
if [ $ga4gh_refine ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/ga4gh/ga4gh_withrefine_truth.vcf.gz) $(fn_md5 $OD/ga4gh_withrefine_truth.vcf.gz)
    assert_equal $(fn_md5 $ANSDIR/ga4gh/ga4gh_withrefine_query.vcf.gz) $(fn_md5 $OD/ga4gh_withrefine_query.vcf.gz)
fi

run ga4gh_badparam $truv ga4gh -w -i notreal.file1234 -o $ANSDIR/ga4gh/ga4gh_withrefine
if [ $ga4gh_badparam ]; then
    assert_exit_code 1
fi
