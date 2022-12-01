# ------------------------------------------------------------
#                                 stratify
# ------------------------------------------------------------
run test_stratify $truv stratify -w \
        $INDIR/beds/include.bed \
        $INDIR/variants/input1.vcf.gz \
        -o $OD/stratify.txt
if [ $test_stratify ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/stratify/stratify.txt) $(fn_md5 $OD/stratify.txt)
fi
