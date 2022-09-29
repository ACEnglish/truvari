# ------------------------------------------------------------
#                                 version
# ------------------------------------------------------------
run test_version $truv version
if [ $test_version ]; then
    assert_in_stdout "$(python3 -c 'import truvari; print(f"Truvari v{truvari.__version__}")')"
fi
