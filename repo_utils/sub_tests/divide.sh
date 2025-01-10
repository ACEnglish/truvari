# ------------------------------------------------------------
#                                 divide
# ------------------------------------------------------------
run divide $truv divide -m 100 $INDIR/variants/multi.vcf.gz $OD/divided/
if [ $divide ]; then
    assert_exit_code 0

    run divide_result
    for i in $ANSDIR/divide/divided/*.vcf.gz
    do
        assert_equal $(fn_md5 $i) $(fn_md5 $OD/divided/$(basename $i))
    done
fi
