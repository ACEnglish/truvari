# ------------------------------------------------------------
#                                 segment
# ------------------------------------------------------------
run segment $truv segment --passonly $INDIR/variants/multi.vcf.gz -o $OD/segment.vcf
if [ $segment ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/segment/segment.vcf) $(fn_md5 $OD/segment.vcf)
fi
