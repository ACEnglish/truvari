# ------------------------------------------------------------
#                                 vcf2df
# ------------------------------------------------------------
run test_vcf2df_bare $truv vcf2df -n $INDIR/variants/input1.vcf.gz $OD/vcf2df_bare.jl
if [ $test_vcf2df_bare ]; then
    assert_exit_code 0
    df_check test_vcf2df_bare $ANSDIR/vcf2df/vcf2df_bare.jl $OD/vcf2df_bare.jl
fi

run test_vcf2df $truv vcf2df -f -i -n $INDIR/variants/input1.vcf.gz $OD/vcf2df.jl
if [ $test_vcf2df ]; then
    assert_exit_code 0
    df_check test_vcf2df $ANSDIR/vcf2df/vcf2df.jl $OD/vcf2df.jl
fi

run test_vcf2df_dir $truv vcf2df -f -i -b -n $ANSDIR/bench/bench23/ $OD/truv2df.jl
if [ $test_vcf2df_dir ]; then
    assert_exit_code 0
    df_check test_vcf2df_dir $ANSDIR/vcf2df/truv2df.jl $OD/truv2df.jl
fi

run test_vcf2df_single $truv vcf2df -f -i -s NA12878 $INDIR/variants/input2.vcf.gz $OD/single_vcf2df.jl
if [ $test_vcf2df_single ]; then
    assert_exit_code 0
    df_check test_vcf2df_single $ANSDIR/vcf2df/single_vcf2df.jl $OD/single_vcf2df.jl
fi

run test_vcf2df_subset $truv vcf2df -f -i -s HG00733,NA12878 $INDIR/variants/multi.vcf.gz $OD/subset_vcf2df.jl
if [ $test_vcf2df_subset ]; then
    assert_exit_code 0
    df_check test_vcf2df_subset $ANSDIR/vcf2df/subset_vcf2df.jl $OD/subset_vcf2df.jl
fi

run test_vcf2df_multi $truv vcf2df -f -i $INDIR/variants/multi.vcf.gz $OD/multi_vcf2df.jl
if [ $test_vcf2df_multi ]; then
    assert_exit_code 0
    df_check test_vcf2df_multi $ANSDIR/vcf2df/multi_vcf2df.jl $OD/multi_vcf2df.jl
fi
