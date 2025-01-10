# ------------------------------------------------------------
#                                 version
# ------------------------------------------------------------
run version $truv version
if [ $version ]; then
    assert_in_stdout "$(python3 -c 'import truvari; print(f"Truvari v{truvari.__version__}")')"
fi
