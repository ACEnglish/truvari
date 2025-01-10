# ------------------------------------------------------------
#                                 entry help
# ------------------------------------------------------------
run help $truv
if [ $help ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $STDERR_FILE) $(fn_md5 $ANSDIR/help.txt)
fi

run help_error $truv banch
if [ $help ]; then
    assert_exit_code 2
fi
