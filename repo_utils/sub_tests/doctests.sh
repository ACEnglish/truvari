# ------------------------------------------------------------
#                                 doctests
# ------------------------------------------------------------
run test_doctests coverage run --concurrency=multiprocessing -p repo_utils/run_doctests.py
if [ $test_doctests ]; then
    assert_exit_code 0 
fi
