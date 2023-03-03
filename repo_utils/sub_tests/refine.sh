# ------------------------------------------------------------
#                                 refine
# ------------------------------------------------------------
# export PATH=$INDIR/external/fake_mafft/:$PATH 

run test_refine_one $truv bench -b $INDIR/refine_data/hg002_base.vcf.gz \
                                -c $INDIR/refine_data/hg002_comp.vcf.gz \
                                --includebed $INDIR/refine_data/h1_hc_tr_hg002.bed \
                                -s 5 -o $OD/refine_output_one

run test_refine_one $truv refine -k -u -f $INDIR/refine_data/chr20.fa.gz $OD/refine_output_one

if [ $test_refine_one ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/refine/refine_output_one/refine.regions.txt) \
                 $(fn_md5 $OD/refine_output_one/refine.regions.txt)
    assert_equal $(fn_md5 $ANSDIR/refine/refine_output_one/refine.variant_summary.json) \
                 $(fn_md5 $OD/refine_output_one/refine.variant_summary.json)
    assert_equal $(fn_md5 $ANSDIR/refine/refine_output_one/refine.region_summary.json) \
                 $(fn_md5 $OD/refine_output_one/refine.region_summary.json)
fi

run test_refine_two $truv bench -b $INDIR/refine_data/hg002_base.vcf.gz \
              -c $INDIR/refine_data/hg002_comp.vcf.gz \
              --includebed $INDIR/refine_data/h1_hc_tr_hg002.bed \
              -s 5 -o $OD/refine_output_two

run test_refine_two $truv refine -k -I -r $INDIR/refine_data/h2_hc_tr_hg002.bed \
               -f $INDIR/refine_data/chr20.fa.gz \
               $OD/refine_output_two

if [ $test_refine_two ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/refine/refine_output_two/refine.regions.txt) \
                 $(fn_md5 $OD/refine_output_two/refine.regions.txt)
    assert_equal $(fn_md5 $ANSDIR/refine/refine_output_two/refine.variant_summary.json) \
                 $(fn_md5 $OD/refine_output_two/refine.variant_summary.json)
    assert_equal $(fn_md5 $ANSDIR/refine/refine_output_two/refine.region_summary.json) \
                 $(fn_md5 $OD/refine_output_two/refine.region_summary.json)
fi
