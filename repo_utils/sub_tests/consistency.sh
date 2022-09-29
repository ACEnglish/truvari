# ------------------------------------------------------------
#                                 consistency
# ------------------------------------------------------------
run test_consistency $truv consistency $INDIR/input*.vcf.gz
if [ $test_consistency ]; then
    assert_exit_code 0
fi

run test_consistency_results $truv consistency $INDIR/input*.vcf.gz
if [ $test_consistency_results ]; then
    assert_equal $(fn_md5 $ANSDIR/consistency.txt) $(fn_md5 $STDOUT_FILE)
fi

run test_consistency_results_json $truv consistency --json $INDIR/input*.vcf.gz
if [ $test_consistency_results_json ]; then
    assert_equal $(fn_md5 $ANSDIR/consistency.json) $(fn_md5 $STDOUT_FILE)
fi
