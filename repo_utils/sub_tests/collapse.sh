# ------------------------------------------------------------
#                                 collapse
# ------------------------------------------------------------
collapse() {
    $truv collapse -f $INDIR/references/reference.fa \
                   -i $INDIR/variants/input${1}.vcf.gz \
                   -o $OD/input${1}_collapsed.vcf \
                   -c $OD/input${1}_removed.vcf \
                   ${2}
}

collapse_assert() {
    sort_vcf $OD/input${1}_collapsed.vcf 
    sort_vcf $OD/input${1}_removed.vcf 

    assert_exit_code 0

    assert_equal $(fn_md5 $ANSDIR/collapse/input${1}_collapsed.vcf) $(fn_md5 $OD/input${1}_collapsed.vcf)
    assert_equal $(fn_md5 $ANSDIR/collapse/input${1}_removed.vcf) $(fn_md5 $OD/input${1}_removed.vcf)
}

collapse_multi() {
    # tests multi sample collapse with provided keep method
    keep=$1
    run collapse_multi_$keep $truv collapse -f $INDIR/references/reference.fa \
                                             -i $INDIR/variants/multi.vcf.gz \
                                             -o $OD/multi_collapsed_${keep}.vcf \
                                             -c $OD/multi_removed_${keep}.vcf \
                                             --keep $keep
}

collapse_multi_assert() {
    assert_exit_code 0
    keep=$1
    sort_vcf $OD/multi_collapsed_${keep}.vcf
    sort_vcf $OD/multi_removed_${keep}.vcf

    assert_equal $(fn_md5 $ANSDIR/collapse/multi_collapsed_${keep}.vcf) $(fn_md5 $OD/multi_collapsed_${keep}.vcf)
    assert_equal $(fn_md5 $ANSDIR/collapse/multi_removed_${keep}.vcf) $(fn_md5 $OD/multi_removed_${keep}.vcf)
}

run collapse_1 collapse 1 "--null-consolidate=PL,DP --write-resolved"
if [ $collapse_1 ]; then
    collapse_assert 1
fi

run collapse_2 collapse 2 --hap
if [ $collapse_2 ]; then
    collapse_assert 2
fi

run collapse_3 collapse 3 --chain
if [ $collapse_3 ]; then
    collapse_assert 3
fi

run collapse_multi_first collapse_multi first
if [ $collapse_multi_first ]; then
    collapse_multi_assert first
fi

run collapse_multi_common collapse_multi common
if [ $collapse_multi_common ]; then
    collapse_multi_assert common
fi

run collapse_multi_maxqual collapse_multi maxqual
if [ $collapse_multi_maxqual ]; then
    collapse_multi_assert maxqual
fi

run collapse_badparams $truv collapse -i nofile.vcf -c nofile.aga -f notref.fa -o $OD --hap --chain --keep common
if [ $collapse_badparams ]; then
    assert_exit_code 100
fi

# Lower collapse sub-chunk threshold
export COLLAP_SUB=1
run collapse_median $truv collapse -f $INDIR/references/reference.fa \
                   -i $INDIR/variants/input1.vcf.gz \
                   -o $OD/input1_median_collapsed.vcf \
                   -c $OD/input1_median_removed.vcf \
                   --median-info
if [ $collapse_median ]; then
    collapse_assert 1_median
fi
unset COLLAP_SUB

run collapse_intragt $truv collapse -i $INDIR/variants/bcftools_merged.vcf.gz \
                        -o $OD/inputintragt_collapsed.vcf \
                        -c $OD/inputintragt_removed.vcf \
                        --intra --gt all --keep maxqual
if [ $collapse_intragt ]; then
    collapse_assert intragt
fi

run collapse_intragth $truv collapse -i $INDIR/variants/bcftools_merged.vcf.gz \
                        -o $OD/inputintragth_collapsed.vcf \
                        -c $OD/inputintragth_removed.vcf \
                        --intra --gt het --keep maxqual
if [ $collapse_intragth ]; then
    collapse_assert intragth
fi


run collapse_chain $truv collapse -i $INDIR/variants/issue196_chain.vcf.gz \
                        -o $OD/inputissue196_collapsed.vcf \
                        -c $OD/inputissue196_removed.vcf \
                        --chain --pctseq 0 --pctsize 0.35
if [ $collapse_chain ]; then
    collapse_assert issue196
fi

