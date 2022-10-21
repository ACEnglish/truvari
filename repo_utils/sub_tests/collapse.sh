# ------------------------------------------------------------
#                                 collapse
# ------------------------------------------------------------
collapse() {
    $truv collapse -f $INDIR/reference.fa \
                   -i $INDIR/input${1}.vcf.gz \
                   -o $OD/input${1}_collapsed.vcf \
                   -c $OD/input${1}_removed.vcf \
                   ${2}
}

collapse_assert() {
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
}

collapse_multi_assert() {
    assert_exit_code 0
    keep=$1
    sort_vcf $OD/multi_collapsed_${keep}.vcf
    sort_vcf $OD/multi_removed_${keep}.vcf

    run test_collapse_multi_${keep}_collapsed
    assert_equal $(fn_md5 $ANSDIR/multi_collapsed_${keep}.vcf) $(fn_md5 $OD/multi_collapsed_${keep}.vcf)

    run test_collapse_multi_${keep}_removed
    assert_equal $(fn_md5 $ANSDIR/multi_removed_${keep}.vcf) $(fn_md5 $OD/multi_removed_${keep}.vcf)
}

run test_collapse_1 collapse 1 --null-consolidate=PL,DP
if [ $test_collapse_1 ]; then
    collapse_assert 1
fi

run test_collapse_2 collapse 2 --hap
if [ $test_collapse_2 ]; then
    collapse_assert 2
fi

run test_collapse_3 collapse 3 --chain
if [ $test_collapse_3 ]; then
    collapse_assert 3
fi

run test_collapse_multi_first collapse_multi first
if [ $test_collapse_multi_first ]; then
    collapse_multi_assert first
fi

run test_collapse_multi_common collapse_multi common
if [ $test_collapse_multi_common ]; then
    collapse_multi_assert common
fi

run test_collapse_multi_maxqual collapse_multi maxqual
if [ $test_collapse_multi_maxqual ]; then
    collapse_multi_assert maxqual
fi

run test_collapse_badparams $truv collapse -i nofile.vcf -c nofile.aga -f notref.fa -o $OD --hap --chain --keep common
if [ $test_collapse_badparams ]; then
    assert_exit_code 100
fi
