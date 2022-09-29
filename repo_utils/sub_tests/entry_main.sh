# ------------------------------------------------------------
#                                 entry help
# ------------------------------------------------------------
run test_help $truv
if [ $test_help ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $STDERR_FILE) $(fn_md5 $ANSDIR/help.txt)
fi

run test_help_error $truv banch
if [ $test_help ]; then
    assert_exit_code 2
fi
