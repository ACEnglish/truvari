# ------------------------------------------------------------
#                                 doctests
# ------------------------------------------------------------
run doctests coverage run --concurrency=multiprocessing -p repo_utils/run_doctests.py
if [ $doctests ]; then
    assert_exit_code 0 
fi
