# ------------------------------------------------------------
#                                 unittests
# ------------------------------------------------------------
run unittests coverage run --concurrency=multiprocessing -p repo_utils/run_unittest.py
if [ $unittest ]; then
    assert_exit_code 0 
fi
