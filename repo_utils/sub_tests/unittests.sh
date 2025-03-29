# ------------------------------------------------------------
#                                 unittests
# ------------------------------------------------------------
run unittests coverage run --concurrency=multiprocessing repo_utils/run_unittest.py
if [ $unittests ]; then
    assert_exit_code 0 
fi
