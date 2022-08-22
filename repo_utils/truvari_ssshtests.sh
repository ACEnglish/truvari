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

bench() {
    # run and test truvari bench
    f1=$1
    f2=$2
    k=$3
    rm -rf $OD/bench${k}
    run test_bench_${k} $truv bench -T 2 -b $INDIR/input${f1}.vcf.gz \
                                      -c $INDIR/input${f2}.vcf.gz \
                                      -f $INDIR/reference.fa \
                                      -o $OD/bench${k}/ ${4}
    assert_exit_code 0

    for i in $ANSDIR/bench${k}/*.vcf
    do
        bname=$(basename $i | sed 's/[\.|\-]/_/g')
        result=$OD/bench${k}/$(basename $i)
        sort_vcf $result
        run test_bench${k}_${bname}
        assert_equal $(fn_md5 $i) $(fn_md5 $result)
    done
}

collapse() {
    # run and test truvari collapse
    run test_collapse_$1 $truv collapse -f $INDIR/reference.fa \
                     -i $INDIR/input${1}.vcf.gz \
                     -o $OD/input${1}_collapsed.vcf \
                     -c $OD/input${1}_removed.vcf \
                     ${2}

    sort_vcf $OD/input${1}_collapsed.vcf 
    sort_vcf $OD/input${1}_removed.vcf 

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
    python3 repo_utils/info_puller.py $base_vcf ${infos} | sort > $OD/answer.txt
    python3 repo_utils/info_puller.py $comp_vcf ${infos} | sort > $OD/result.txt
    run test_infos_${name}
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
assert_equal $(fn_md5 $STDERR_FILE) $(fn_md5 $ANSDIR/help.txt)

run test_help_error $truv banch
assert_exit_code 2

# ------------------------------------------------------------
#                                 version
# ------------------------------------------------------------
run test_version $truv version
assert_in_stdout "$(python3 -c 'import truvari; print(f"Truvari v{truvari.__version__}")')"

# ------------------------------------------------------------
#                                 bench
# ------------------------------------------------------------
bench 1 2 12
bench 1 3 13
bench 2 3 23 --multimatch

# Testing --includebed
bench 1 3 13_includebed "--includebed $INDIR/include.bed"

# Testing --extend
bench 1 3 13_extend "--includebed $INDIR/include.bed --extend 500"

rm -rf $OD/bench_giab
run test_bench_giab $truv bench -b $INDIR/giab.vcf.gz \
                                -c $INDIR/input1.vcf.gz \
                                -f $INDIR/reference.fa \
                                -o $OD/bench_giab/ \
                                -T 1 \
                                --includebed $INDIR/giab.bed \
                                --multimatch \
                                --giabreport \
                                --prog
assert_exit_code 0

run test_bench_giab_report
assert_equal $(fn_md5 $ANSDIR/bench_giab_report.txt) $(fn_md5 $OD/bench_giab/giab_report.txt)

run test_bench_badparams $truv bench -b nofile.vcf -c nofile.aga -f notref.fa -o $OD
assert_exit_code 100

# ------------------------------------------------------------
#                                 collapse
# ------------------------------------------------------------
collapse 1 --null-consolidate=PL,DP
collapse 2 --hap
collapse 3 --chain

collapse_multi first
collapse_multi common
collapse_multi maxqual

run test_collapse_badparams $truv collapse -i nofile.vcf -c nofile.aga -f notref.fa -o $OD --hap --chain --keep common
assert_exit_code 100

# ------------------------------------------------------------
#                                 consistency
# ------------------------------------------------------------
run test_consistency $truv consistency $INDIR/input*.vcf.gz
assert_exit_code 0

run test_consistency_results $truv consistency $INDIR/input*.vcf.gz
assert_equal $(fn_md5 $ANSDIR/consistency.txt) $(fn_md5 $STDOUT_FILE)

run test_consistency_results $truv consistency --json $INDIR/input*.vcf.gz
assert_equal $(fn_md5 $ANSDIR/consistency.json) $(fn_md5 $STDOUT_FILE)

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
run test_anno_remap $truv anno remap -H 10 -i $VCF -r $REF -o $OD/anno_remap.vcf
assert_exit_code 0
info_tests remap $OD/anno_remap.vcf REMAP,REMAPHits

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

#                                 svinfo
run test_anno_svinfo $truv anno svinfo -i $VCF -o $OD/anno_svinfo.vcf
assert_exit_code 0
info_tests svinfo $OD/anno_svinfo.vcf SVTYPE,SVLEN

#                                 grm
run test_anno_grm $truv anno grm -i $INDIR/input2.vcf.gz -r $REF -o $OD/grm.jl
assert_exit_code 0

df_check test_grm_result $ANSDIR/grm.jl $OD/grm.jl

#                                 trf
run test_anno_trf $truv anno trf -i $INDIR/multi.vcf.gz \
                                 -r $INDIR/repeats.adotto.bed.gz \
                                 -f $INDIR/reference.fa \
                                 -e $INDIR/external/trf  \
                                 -m 5 \
                                 -o $OD/trf.vcf
assert_exit_code 0

run test_anno_trf_reg $truv anno trf -i $INDIR/multi.vcf.gz \
                                 -r $INDIR/repeats.adotto.bed.gz \
                                 -f $INDIR/reference.fa \
                                 -e $INDIR/external/trf  \
                                 -m 5 -R \
                                 -o $OD/trf.reg.vcf

assert_exit_code 0

run test_anno_badparam $truv anno trf -i $INDIR/input_null.vcf \
                                 -r $INDIR/simplerepeat_null.bed \
                                 -f $INDIR/reference.fa \
                                 -e $INDIR/external/trf_dne  \
                                 -o $OD/trf_null.vcf
assert_exit_code 1


# TRF isn't deterministic for some reason, so it gives a different answer in  action
# run test_anno_trf_result
# assert_equal $(fn_md5 $ANSDIR/trf.vcf) $(fn_md5 $OD/trf.vcf)

#                                 repmask
run test_anno_repmask $truv anno repmask -i $INDIR//multi.vcf.gz -o $OD/repmask.vcf -e $INDIR/external/fakeRM.py
assert_exit_code 0
assert_equal $(fn_md5 $ANSDIR/repmask.vcf) $(fn_md5 $OD/repmask.vcf)

run test_anno_repmask_err $truv anno repmask -i $INDIR/input1.vcf.gz -o $OD/repmask.vcf -e $INDIR/external/fakeRM.py
assert_exit_code 1

#                                 bpovl
run test_anno_bpovl $truv anno bpovl -i $INDIR/input1.vcf.gz \
                    -o $OD/anno_bpovl.jl \
                    -a $INDIR/anno.gtf.gz -p gff --sizemin 2
assert_exit_code 0
df_check test_anno_bpovl_result $ANSDIR/anno_bpovl.jl $OD/anno_bpovl.jl

#                                 density
run truvari_anno_density $truv anno density -i $INDIR/input3.vcf.gz \
                    -o $OD/anno_density.jl \
                    -g $INDIR/genome.bed -m $INDIR/mask.bed
assert_exit_code 0
df_check test_anno_density_result $ANSDIR/anno_density.jl $OD/anno_density.jl

#                                 dpcnt
run test_anno_dpcnt $truv anno dpcnt -i $VCF -o $OD/anno_dpcnt.vcf
assert_exit_code 0
info_tests dpcnt $OD/anno_dpcnt.vcf DPCNT,ADCNT

#                                 lcr
run test_anno_lcr $truv anno lcr -i $VCF -o $OD/anno_lcr.vcf
assert_exit_code 0
info_tests lcr $OD/anno_lcr.vcf LCR

#                                 grpaf
run test_anno_grpaf $truv anno grpaf -i $INDIR/grpaf.vcf.gz -l $INDIR/grpaf.labels.txt -o $OD/anno_grpaf.vcf
assert_exit_code 0
assert_equal $(fn_md5 $ANSDIR/anno_grpaf.vcf) $(fn_md5 $OD/anno_grpaf.vcf)

run test_anno_grpaf_strict $truv anno grpaf --strict -i $INDIR/grpaf.vcf.gz -l $INDIR/grpaf.labels.txt -o $OD/anno_grpaf.vcf
assert_exit_code 1

run test_anno_grpaf_subset $truv anno grpaf --tags AF,HWE,ExcHet -i $INDIR/grpaf.vcf.gz -l $INDIR/grpaf.labels.txt -o $OD/anno_grpaf.subtags.vcf
assert_exit_code 0
assert_equal $(fn_md5 $ANSDIR/anno_grpaf.subtags.vcf) $(fn_md5 $OD/anno_grpaf.subtags.vcf)

# ------------------------------------------------------------
#                                 vcf2df
# ------------------------------------------------------------
run test_vcf2df_bare $truv vcf2df $INDIR/input1.vcf.gz $OD/vcf2df_bare.jl
assert_exit_code 0
df_check test_vcf2df_bare_result $ANSDIR/vcf2df_bare.jl $OD/vcf2df_bare.jl

run test_vcf2df $truv vcf2df -f -i $INDIR/input1.vcf.gz $OD/vcf2df.jl
assert_exit_code 0
df_check test_vcf2df_result $ANSDIR/vcf2df.jl $OD/vcf2df.jl

run test_vcf2df_dir $truv vcf2df -f -i -b $OD/bench23/ $OD/truv2df.jl
assert_exit_code 0
df_check test_vcf2df_dir_result $ANSDIR/truv2df.jl $OD/truv2df.jl

run test_vcf2df_single $truv vcf2df -f -i -s NA12878 $INDIR/input2.vcf.gz $OD/single_vcf2df.jl
assert_exit_code 0
df_check test_vcf2df_single $ANSDIR/single_vcf2df.jl $OD/single_vcf2df.jl

run test_vcf2df_subset $truv vcf2df -f -i -s HG00733,NA12878 -m $INDIR/multi.vcf.gz $OD/subset_vcf2df.jl
assert_exit_code 0
df_check test_vcf2df_subset $ANSDIR/subset_vcf2df.jl $OD/subset_vcf2df.jl

run test_vcf2df_multi $truv vcf2df -f -i -m $INDIR/multi.vcf.gz $OD/multi_vcf2df.jl
assert_exit_code 0
df_check test_vcf2df_multi $ANSDIR/multi_vcf2df.jl $OD/multi_vcf2df.jl


# ------------------------------------------------------------
#                                 segment
# ------------------------------------------------------------
run test_segment $truv segment $ANSDIR/anno_answers.vcf.gz $OD/segment.vcf
assert_exit_code 0

run test_segment_result
assert_equal $(fn_md5 $ANSDIR/segment.vcf) $(fn_md5 $OD/segment.vcf)

# ------------------------------------------------------------
#                                 divide
# ------------------------------------------------------------
run test_divide $truv divide $INDIR/multi.vcf.gz $OD/divided/ --no-compress
assert_exit_code 0

run test_divide_result
for i in $ANSDIR/divided/*.vcf
do
    assert_equal $(fn_md5 $i) $(fn_md5 $OD/divided/$(basename $i))
done

# ------------------------------------------------------------
#                                 doctests
# ------------------------------------------------------------

run test_doctests coverage run --concurrency=multiprocessing -p repo_utils/run_doctests.py
assert_exit_code 0 

# ------------------------------------------------------------
#                                 coverage.py
# ------------------------------------------------------------

printf "\n${BOLD}generating test coverage reports${NC}\n"
coverage combine
coverage report --include=truvari/*
coverage html --include=truvari/* -d $OD/htmlcov/
coverage json --include=truvari/* -o $OD/coverage.json
python3 repo_utils/coverage_maker.py $OD/coverage.json

