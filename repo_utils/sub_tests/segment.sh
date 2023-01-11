# ------------------------------------------------------------
#                                 segment
# ------------------------------------------------------------
run test_segment $truv segment $ANSDIR/anno_answers.vcf.gz -o $OD/segment.vcf
if [ $test_segment ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/segment/segment.vcf) $(fn_md5 $OD/segment.vcf)
fi
