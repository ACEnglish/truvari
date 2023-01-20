INDIR=repo_utils/test_files
ANSDIR=repo_utils/answer_key
OD=test_results
COVERAGE_RCFILE=.coveragerc

# truvari in code coverage
truv="coverage run --concurrency=multiprocessing,thread -p -m truvari.__main__"

# ------------------------------------------------------------
#                                 test helpers
# ------------------------------------------------------------
sort_vcf() {
    fn=$1
    python3 -c "import sys; from pysam import bcftools; sys.stdout.write(bcftools.sort(sys.argv[1]))" $fn > tmp
    mv tmp $fn
}

fn_md5() {
    fn=$1
    # simple md5sum checking
    md5sum $fn | cut -f1 -d\  
}

info_tests() {
    # pull and query info fields from vcfs
    name=$1
    base_vcf=$ANSDIR/anno_answers.vcf.gz
    comp_vcf=$2
    infos=$3
    python3 repo_utils/info_puller.py $base_vcf ${infos} | sort > $OD/answer.txt
    python3 repo_utils/info_puller.py $comp_vcf ${infos} | sort > $OD/result.txt
    assert_equal $(fn_md5 $OD/answer.txt) $(fn_md5 $OD/result.txt)
}

df_check() {
    # check if joblib saved pandas dataframes are equivalent
    test_name=$1
    base_df=$2
    comp_df=$3
    run $test_name python3 -c """
import joblib;
a = joblib.load(\"$base_df\")
b = joblib.load(\"$comp_df\")
assert a.equals(b), \"$base_df != $comp_df\";
"""
    assert_exit_code 0
}

dump_logs() {
    # For debugging, put this after any failing command
    echo '---- STDERR ----'
    cat $STDERR_FILE
    echo '---- STDOUT ----'
    cat $STDOUT_FILE
}
