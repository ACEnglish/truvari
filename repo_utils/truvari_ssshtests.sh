
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
test $(which bcftools) || sudo apt-get install bcftools 
info_tests() {
    name=$1
    base_vcf=$ANSDIR/anno_answers.vcf.gz
    comp_vcf=$2
    infos=$3
    bcftools query -f "%CHROM\t%POS\t${INFOS}\n" $base_vcf | sort > answer.txt
    bcftools query -f "%CHROM\t%POS\t${INFOS}\n" $comp_vcf | sort > result.txt
    run test_infos_${name}
    assert_equal $(fn_md5 answer.txt) $(fn_md5 result.txt)
}

VCF=$INDIR/multi.vcf.gz
REF=$INDIR/reference.fa
#                                 hompct
run test_anno_hompct truvari anno hompct -i $VCF -o anno_hompct.vcf
assert_exit_code $? 0
info_tests hompct anno_hompct.vcf INFO/HOMPCT

#                                 remap
run test_anno_remap truvari anno remap -i $VCF -r $REF -o anno_remap.vcf
assert_exit_code $? 0
info_tests remap anno_remap.vcf INFO/REMAP

#                                 gcpct
run test_anno_gcpct truvari anno gcpct -i $VCF -r $REF -o anno_gcpct.vcf
assert_exit_code $? 0
info_tests gcpct anno_gcpct.vcf INFO/GCPCT

#                                 gtcnt
run test_anno_gtcnt truvari anno gtcnt -i $VCF -o anno_gtcnt.vcf
assert_exit_code $? 0
info_tests gtcnt anno_gtcnt.vcf INFO/GTCNT

#                                 numneigh
run test_anno_numneigh truvari anno numneigh -i $VCF -o anno_numneigh.vcf
assert_exit_code $? 0
info_tests numneigh anno_numneigh.vcf INFO/NumNeighbors,INFO/NeighId

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


