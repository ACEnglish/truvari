# ------------------------------------------------------------
#                                 consistency
# ------------------------------------------------------------

run test_consistency $truv consistency $INDIR/variants/input*.vcf.gz
if [ $test_consistency ]; then
    assert_equal $(fn_md5 $ANSDIR/consistency/dup_consistency.txt) $(fn_md5 $STDOUT_FILE)
fi

run test_consistency_json $truv consistency -d -o $OD/consistency.txt --json $INDIR/variants/input*.vcf.gz
if [ $test_consistency_json ]; then
    assert_equal $(fn_md5 $ANSDIR/consistency/consistency.json) $(fn_md5 $STDOUT_FILE)
fi

