# ------------------------------------------------------------
#                                 consistency
# ------------------------------------------------------------

run consistency $truv consistency $INDIR/variants/input*.vcf.gz
if [ $consistency ]; then
    assert_equal $(fn_md5 $ANSDIR/consistency/dup_consistency.txt) $(fn_md5 $STDOUT_FILE)
fi

run consistency_json $truv consistency -d -o $OD/consistency.txt --json $INDIR/variants/input*.vcf.gz
if [ $consistency_json ]; then
    assert_equal $(fn_md5 $ANSDIR/consistency/consistency.json) $(fn_md5 $STDOUT_FILE)
    assert_equal $(fn_md5 $ANSDIR/consistency/consistency.txt) $(fn_md5 $OD/consistency.txt)
fi

