# ------------------------------------------------------------
#                                 ga4gh
# ------------------------------------------------------------
run ga4gh $truv ga4gh -r -i $ANSDIR/refine/refine_output_three/ -o $OD/ga4gh_norefine
if [ $ga4gh ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/ga4gh/ga4gh_norefine.base.vcf.gz) $(fn_md5 $OD/ga4gh_norefine.base.vcf.gz)
    assert_equal $(fn_md5 $ANSDIR/ga4gh/ga4gh_norefine.comp.vcf.gz) $(fn_md5 $OD/ga4gh_norefine.comp.vcf.gz)
fi

run ga4gh_refine $truv ga4gh -b 0 -w -i $ANSDIR/refine/refine_output_three/ -o $OD/ga4gh_withrefine
if [ $ga4gh_refine ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/ga4gh/ga4gh_withrefine.base.vcf.gz) $(fn_md5 $OD/ga4gh_withrefine.base.vcf.gz)
    assert_equal $(fn_md5 $ANSDIR/ga4gh/ga4gh_withrefine.comp.vcf.gz) $(fn_md5 $OD/ga4gh_withrefine.comp.vcf.gz)
fi

run ga4gh_refine_orig $truv ga4gh -b 0 -i $ANSDIR/refine/refine_output_three/ -o $OD/ga4gh_withrefine_orig
if [ $ga4gh_refine_orig ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/ga4gh/ga4gh_withrefine_orig.base.vcf.gz) $(fn_md5 $OD/ga4gh_withrefine_orig.base.vcf.gz)
    assert_equal $(fn_md5 $ANSDIR/ga4gh/ga4gh_withrefine_orig.comp.vcf.gz) $(fn_md5 $OD/ga4gh_withrefine_orig.comp.vcf.gz)
fi

run ga4gh_badparam $truv ga4gh -w -i notreal.file1234 -o $ANSDIR/ga4gh/ga4gh_withrefine
if [ $ga4gh_badparam ]; then
    assert_exit_code 1
fi
