test -e ssshtest || curl -O https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

source ssshtest

# Work inside of the repo folder
cd "$( dirname "${BASH_SOURCE[0]}" )"/../
INDIR=repo_utils/test_files
ANSDIR=$INDIR/answer_key
OD=test_results
COVERAGE_RCFILE=.coveragerc

# Reset test results
rm -rf $OD
mkdir -p $OD

truv="coverage run --concurrency=multiprocessing -p -m truvari.__main__"
# ------------------------------------------------------------
#                                 test helpers
# ------------------------------------------------------------
fn_md5() {
    fn=$1
    # simple md5sum checking
    md5sum $fn | cut -f1 -d\  
}

bench() {
    # run and test truvari bench
    f1=$1
    f2=$2
    k=${f1}${f2}
    rm -rf $OD/bench${k}
    run test_bench_${k} $truv bench -b $INDIR/input${f1}.vcf.gz \
                                      -c $INDIR/input${f2}.vcf.gz \
                                      -f $INDIR/reference.fa \
                                      -o $OD/bench${k}/
    assert_exit_code 0

    for i in $ANSDIR/bench${k}/*.vcf
    do
        bname=$(basename $i | sed 's/[\.|\-]/_/g')
        run test_bench${k}_${bname}
        assert_equal $(fn_md5 $i) $(fn_md5 $OD/bench${k}/$(basename $i))
    done
}

collapse_hap() {
    # run and test truvari collapse
    run test_collapse_$1 $truv collapse -f $INDIR/reference.fa \
                     -i $INDIR/input${1}.vcf.gz \
                     -o $OD/input${1}_collapsed.vcf \
                     -c $OD/input${1}_removed.vcf \
                     --hap
    assert_exit_code 0

    run test_collapse_${1}_collapsed
    assert_equal $(fn_md5 $ANSDIR/input${1}_collapsed.vcf) $(fn_md5 $OD/input${1}_collapsed.vcf)

    run test_collapse_${1}_removed
    assert_equal $(fn_md5 $ANSDIR/input${1}_removed.vcf) $(fn_md5 $OD/input${1}_removed.vcf)
}

collapse_multi() {
    # tests multi sample collapse with provided keep method
    keep=$1
    run test_collapse_multi_$keep $truv collapse -f $INDIR/reference.fa \
                                             -i $INDIR/multi.vcf.gz \
                                             -o $OD/multi_collapsed_${keep}.vcf \
                                             -c $OD/multi_removed_${keep}.vcf \
                                             --keep $keep
    assert_exit_code 0

    run test_collapse_multi_${keep}_collapsed
    assert_equal $(fn_md5 $ANSDIR/multi_collapsed_${keep}.vcf) $(fn_md5 $OD/multi_collapsed_${keep}.vcf)

    run test_collapse_multi_${keep}_removed
    assert_equal $(fn_md5 $ANSDIR/multi_removed_${keep}.vcf) $(fn_md5 $OD/multi_removed_${keep}.vcf)
}

info_tests() {
    # pull and query info fields from vcfs
    name=$1
    base_vcf=$ANSDIR/anno_answers.vcf.gz
    comp_vcf=$2
    infos=$3
    python repo_utils/info_puller.py $base_vcf ${infos} | sort > $OD/answer.txt
    python repo_utils/info_puller.py $comp_vcf ${infos} | sort > $OD/result.txt
    run test_infos_${name}
    assert_equal $(fn_md5 $OD/answer.txt) $(fn_md5 $OD/result.txt)
}

df_check() {
    # check if joblib saved pandas dataframes are equivalent
    test_name=$1
    base_df=$2
    comp_df=$3
    run $test_name python -c """
import joblib;
a = joblib.load(\"$base_df\")
b = joblib.load(\"$comp_df\")
assert a.equals(b), \"$base_df != $comp_df\";
"""
    assert_exit_code 0
}

dump_logs() {
    echo '---- STDERR ----'
    cat $STDERR_FILE
    echo '---- STDOUT ----'
    cat $STDOUT_FILE
}
# ------------------------------------------------------------
#                                 entry help
# ------------------------------------------------------------
run test_help $truv
assert_exit_code 0

run test_help2 $truv
assert_equal $(fn_md5 $STDERR_FILE) $(fn_md5 $ANSDIR/help.txt)

# ------------------------------------------------------------
#                                 version
# ------------------------------------------------------------
run test_version $truv version
assert_in_stdout "$(python -c 'import truvari; print(f"Truvari v{truvari.__version__}")')"

# ------------------------------------------------------------
#                                 bench
# ------------------------------------------------------------
bench 1 2
bench 1 3
bench 2 3

rm -rf $OD/bench_giab
run test_bench_giab $truv bench -b $INDIR/giab.vcf.gz \
                                -c $INDIR/input1.vcf.gz \
                                -f $INDIR/reference.fa \
                                -o $OD/bench_giab/ \
                                --includebed $INDIR/giab.bed \
                                --giabreport
assert_exit_code 0

run test_bench_giab_report
assert_equal $(fn_md5 $ANSDIR/bench_giab_report.txt) $(fn_md5 $OD/bench_giab/giab_report.txt)

# ------------------------------------------------------------
#                                 collapse
# ------------------------------------------------------------
collapse_hap 1
collapse_hap 2
collapse_hap 3

collapse_multi first
collapse_multi common
collapse_multi maxqual

# ------------------------------------------------------------
#                                 consistency
# ------------------------------------------------------------
run test_consistency $truv consistency $INDIR/input*.vcf.gz
assert_exit_code 0


run test_consistency_results $truv consistency $INDIR/input*.vcf.gz
assert_equal $(fn_md5 $ANSDIR/consistency.txt) $(fn_md5 $STDOUT_FILE)

# ------------------------------------------------------------
#                                 anno
# ------------------------------------------------------------
VCF=$INDIR/multi.vcf.gz
REF=$INDIR/reference.fa

#                                 hompct
run test_anno_hompct $truv anno hompct -i $VCF -o $OD/anno_hompct.vcf
assert_exit_code 0
info_tests hompct $OD/anno_hompct.vcf HOMPCT

#                                 remap
run test_anno_remap $truv anno remap -i $VCF -r $REF -o $OD/anno_remap.vcf
assert_exit_code 0
info_tests remap $OD/anno_remap.vcf REMAP

#                                 gcpct
run test_anno_gcpct $truv anno gcpct -i $VCF -r $REF -o $OD/anno_gcpct.vcf
assert_exit_code 0
info_tests gcpct $OD/anno_gcpct.vcf GCPCT

#                                 gtcnt
run test_anno_gtcnt $truv anno gtcnt -i $VCF -o $OD/anno_gtcnt.vcf
assert_exit_code 0
info_tests gtcnt $OD/anno_gtcnt.vcf GTCNT

#                                 numneigh
run test_anno_numneigh $truv anno numneigh -i $VCF -o $OD/anno_numneigh.vcf
assert_exit_code 0
info_tests numneigh $OD/anno_numneigh.vcf NumNeighbors,NeighId

#                                 grm
run test_anno_grm $truv anno grm -i $INDIR/input2.vcf.gz -r $REF -o $OD/grm.jl
assert_exit_code 0

df_check test_grm_result $ANSDIR/grm.jl $OD/grm.jl

#                               af
# waiting on unit tests. might make an anno tool out of this eventually, though

# requires external executables which are too large to put into the repository repo_utils/test_files/
#                                 trf
run test_anno_trf $truv anno trf -i $INDIR/input1.vcf.gz \
                                 -s $INDIR/simplerepeat.txt.gz \
                                 -f $INDIR/reference.fa \
                                 -e $INDIR/external/trf  \
                                 -o $OD/trf.vcf
assert_exit_code 0

cat $OD/trf.vcf
run test_anno_trf_result
assert_equal $(fn_md5 $ANSDIR/trf.vcf) $(fn_md5 $OD/trf.vcf)

#                                 repmask
#run test_anno_repmask $truv anno repmask -i $INDIR/input1.vcf.gz -e $INDIR/external/RepeatMasker -T 1 -o $OD/repmask.vcf
#assert_exit_code 0

#run test_anno_repmask_result
#assert_equal $(fn_md5 $ANSDIR/repmask.vcf) $(fn_md5 $OD/repmask.vcf)


# ------------------------------------------------------------
#                                 vcf2df
# ------------------------------------------------------------
run test_vcf2df $truv vcf2df -f -i $INDIR/input1.vcf.gz $OD/vcf2df.jl
assert_exit_code 0

df_check test_vcf2df_result $ANSDIR/vcf2df.jl $OD/vcf2df.jl

run test_vcf2df_dir $truv vcf2df -f -i -b $OD/bench23/ $OD/truv2df.jl
assert_exit_code 0
df_check test_vcf2df_dir_result $ANSDIR/truv2df.jl $OD/truv2df.jl

# ------------------------------------------------------------
#                                 coverage.py
# ------------------------------------------------------------
coverage combine
coverage report --include=truvari/*
coverage html --include=truvari/* -d $OD/htmlcov/
coverage json --include=truvari/* -o $OD/coverage.json
python repo_utils/coverage_maker.py $OD/coverage.json

