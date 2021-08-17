
test -e ssshtest || curl -O https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

source ssshtest

BASDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
INDIR=$BASDIR/test_files
ANSDIR=$INDIR/answer_key/

fn_md5() {
    md5sum $1 | cut -f1 -d\  
}

# ------------------------------------------------------------
#                                 version
# ------------------------------------------------------------
run test_version truvari version
assert_in_stdout "$(python -c 'import truvari; print(f"Truvari v{truvari.__version__}")')"

# ------------------------------------------------------------
#                                 bench
# ------------------------------------------------------------
bench() {
    f1=$1
    f2=$2
    k=${f1}${f2}
    rm -rf bench${k}
    run test_bench_${k} truvari bench -b $INDIR/input${f1}.vcf.gz \
                                      -c $INDIR/input${f2}.vcf.gz \
                                      -f $INDIR/reference.fa \
                                      -o bench${k}/
    assert_exit_code $? 0

    for i in $ANSDIR/bench${k}/*.vcf
    do
        bname=$(basename $i | sed 's/[\.|\-]/_/g')
        run test_bench${k}_${bname}
        assert_equal $(fn_md5 $i) $(fn_md5 bench${k}/$(basename $i))
    done
}

bench 1 2
bench 1 3
bench 2 3

# ------------------------------------------------------------
#                                 collapse
# ------------------------------------------------------------

collapse_hap() {
    run test_collapse_$1 truvari collapse -f $INDIR/reference.fa \
                     -i $INDIR/input${1}.vcf.gz \
                     -o input${1}_collapsed.vcf \
                     -c input${1}_removed.vcf \
                     --hap
    assert_exit_code $? 0
    
    
    run test_collapse_${1}_collapsed
    assert_equal $(fn_md5 $ANSDIR/input${1}_collapsed.vcf) $(fn_md5 input${1}_collapsed.vcf)

    run test_collapse_${1}_removed
    assert_equal $(fn_md5 $ANSDIR/input${1}_removed.vcf) $(fn_md5 input${1}_removed.vcf)
}

collapse_hap 1
collapse_hap 2
collapse_hap 3

run test_collapse_multi truvari collapse -f $INDIR/reference.fa \
                                         -i $INDIR/multi.vcf.gz \
                                         -o multi_collapsed.vcf \
                                         -c multi_removed.vcf
assert_exit_code $? 0

run test_collapse_multi_collapsed
assert_equal $(fn_md5 $ANSDIR/multi_collapsed.vcf) $(fn_md5 multi_collapsed.vcf)

run test_collapse_multi_removed
assert_equal $(fn_md5 $ANSDIR/multi_removed.vcf) $(fn_md5 multi_removed.vcf)

# ------------------------------------------------------------
#                                 consistency
# ------------------------------------------------------------

run test_consistency truvari consistency $INDIR/input1.vcf.gz \
                                         $INDIR/input2.vcf.gz \
                                         $INDIR/input3.vcf.gz
assert_exit_code $? 0

run test_consistency_results
assert_equal $(fn_md5 $ANSDIR/consistency.txt) $(fn_md5 consistency.txt)

# ------------------------------------------------------------
#                                 anno
# ------------------------------------------------------------

#So for these, I can just do a bcftools query to compare. 
#I can run them all into a single base VCF with the answer,
#Then when I run this test, `bcftools query` compare the two,
#If there are any differences then I fail the test. 
#These tests would just be validating that
#1 - the commands work without throwing errors
#2 - the command outputs stay the same

#If a value needs to be changed, it will be changed in the answer key


#                                 gcpct
#vcf_anno_compare answer_key.vcf.gz quiz.vcf.gz INFOA INFOB INFOC etc

#                                 gtcnt
#vcf_anno_compare INFO/A INFO/B INFO/C
#                                 remap
#                                 numneigh
#                                 hompct

#This will need a custom checker
#                                 grm


#These can also be tested with the rest of the anno framework. 
#However I need extra tools to make it happen.

#                                 repmask
#                                 trf



# ------------------------------------------------------------
#                                 truv2df
# ------------------------------------------------------------
#I can probably use the same custom checker from grm for this



# ------------------------------------------------------------
#                                 stats
# ------------------------------------------------------------
#This needs to be deprecated


