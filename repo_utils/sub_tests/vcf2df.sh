# ------------------------------------------------------------
#                                 vcf2df
# ------------------------------------------------------------
run test_vcf2df_bare $truv vcf2df -n $INDIR/input1.vcf.gz $OD/vcf2df_bare.jl
if [ $test_vcf2df_bare ]; then
    assert_exit_code 0
    df_check test_vcf2df_bare_result $ANSDIR/vcf2df_bare.jl $OD/vcf2df_bare.jl
fi

run test_vcf2df $truv vcf2df -f -i -n $INDIR/input1.vcf.gz $OD/vcf2df.jl
if [ $test_vcf2df ]; then
    assert_exit_code 0
    df_check test_vcf2df_result $ANSDIR/vcf2df.jl $OD/vcf2df.jl
fi

run test_vcf2df_dir $truv vcf2df -f -i -b -n $OD/bench23/ $OD/truv2df.jl
if [ $test_vcf2df_dir ]; then
    assert_exit_code 0
    df_check test_vcf2df_dir_result $ANSDIR/truv2df.jl $OD/truv2df.jl
fi

run test_vcf2df_single $truv vcf2df -f -i -s NA12878 $INDIR/input2.vcf.gz $OD/single_vcf2df.jl
if [ $test_vcf2df_single ]; then
    assert_exit_code 0
    df_check test_vcf2df_single $ANSDIR/single_vcf2df.jl $OD/single_vcf2df.jl
fi

run test_vcf2df_subset $truv vcf2df -f -i -s HG00733,NA12878 $INDIR/multi.vcf.gz $OD/subset_vcf2df.jl
if [ $test_vcf2df_subset ]; then
    assert_exit_code 0
    df_check test_vcf2df_subset $ANSDIR/subset_vcf2df.jl $OD/subset_vcf2df.jl
fi

run test_vcf2df_multi $truv vcf2df -f -i $INDIR/multi.vcf.gz $OD/multi_vcf2df.jl
if [ $test_vcf2df_multi ]; then
    assert_exit_code 0
    df_check test_vcf2df_multi $ANSDIR/multi_vcf2df.jl $OD/multi_vcf2df.jl
fi
