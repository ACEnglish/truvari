# ------------------------------------------------------------
#                                 segment
# ------------------------------------------------------------
run test_segment $truv segment --no-info $ANSDIR/anno_answers.vcf.gz -o $OD/segment.vcf
if [ $test_segment ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/segment/segment.vcf) $(fn_md5 $OD/segment.vcf)
fi

run test_segment2 $truv segment --passonly $INDIR/variants/real_small_comp.vcf.gz -o $OD/segment2.vcf
if [ $test_segment2 ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/segment/segment2.vcf) $(fn_md5 $OD/segment2.vcf)
fi
