# ------------------------------------------------------------
#                                 divide
# ------------------------------------------------------------
run test_divide $truv divide $INDIR/multi.vcf.gz $OD/divided/ --no-compress
if [ $test_divide ]; then
    assert_exit_code 0

    run test_divide_result
    for i in $ANSDIR/divided/*.vcf
    do
        assert_equal $(fn_md5 $i) $(fn_md5 $OD/divided/$(basename $i))
    done
fi
