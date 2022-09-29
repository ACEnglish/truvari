# ------------------------------------------------------------
#                                 anno
# ------------------------------------------------------------
VCF=$INDIR/multi.vcf.gz
REF=$INDIR/reference.fa

#                                 hompct
run test_anno_hompct \
    $truv anno hompct -i $VCF -o $OD/anno_hompct.vcf
if [ $test_anno_hompct ]; then
    assert_exit_code 0
    info_tests hompct $OD/anno_hompct.vcf HOMPCT
fi

#                                 remap
run test_anno_remap \
    $truv anno remap -H 10 -i $VCF -r $REF -o $OD/anno_remap.vcf
if [ $test_anno_remap ]; then
    assert_exit_code 0
    info_tests remap $OD/anno_remap.vcf REMAP,REMAPHits
fi

#                                 gcpct
run test_anno_gcpct \
    $truv anno gcpct -i $VCF -r $REF -o $OD/anno_gcpct.vcf
if [ $test_anno_gcpct ]; then
    assert_exit_code 0
    info_tests gcpct $OD/anno_gcpct.vcf GCPCT
fi

#                                 gtcnt
run test_anno_gtcnt \
    $truv anno gtcnt -i $VCF -o $OD/anno_gtcnt.vcf
if [ $test_anno_gtcnt ]; then
    assert_exit_code 0
    info_tests gtcnt $OD/anno_gtcnt.vcf GTCNT
fi

#                                 numneigh
run test_anno_numneigh \
    $truv anno numneigh -i $VCF -o $OD/anno_numneigh.vcf
if [ $test_anno_numneigh ]; then
    assert_exit_code 0
    info_tests numneigh $OD/anno_numneigh.vcf NumNeighbors,NeighId
fi

#                                 svinfo
run test_anno_svinfo \
    $truv anno svinfo -i $VCF -o $OD/anno_svinfo.vcf
if [ $test_anno_svinfo ]; then
    assert_exit_code 0
    info_tests svinfo $OD/anno_svinfo.vcf SVTYPE,SVLEN
fi

#                                 grm
run test_anno_grm \
    $truv anno grm -i $INDIR/input2.vcf.gz -r $REF -o $OD/grm.jl
if [ $test_anno_grm ]; then
    assert_exit_code 0
    df_check test_grm_result $ANSDIR/grm.jl $OD/grm.jl
fi

#                                 trf
run test_anno_trf \
    $truv anno trf -i $INDIR/multi.vcf.gz \
                   -r $INDIR/repeats.adotto.bed.gz \
                   -f $INDIR/reference.fa \
                   -e $INDIR/external/trf  \
                   -m 5 \
                   -o $OD/trf.vcf
if [ $test_anno_trf ]; then
    assert_exit_code 0
fi

run test_anno_trf_reg \
    $truv anno trf -i $INDIR/multi.vcf.gz \
                   -r $INDIR/repeats.adotto.bed.gz \
                   -f $INDIR/reference.fa \
                   -e $INDIR/external/trf  \
                   -m 5 -R \
                   -o $OD/trf.reg.vcf
if [ $test_anno_trf_reg ]; then
    assert_exit_code 0
fi

run test_anno_badparam \
    $truv anno trf -i $INDIR/input_null.vcf \
                   -r $INDIR/simplerepeat_null.bed \
                   -f $INDIR/reference.fa \
                   -e $INDIR/external/trf_dne  \
                   -o $OD/trf_null.vcf
if [ $test_anno_badparam ]; then
    assert_exit_code 1
fi


# TRF isn't deterministic for some reason, so it gives a different answer in  action
# run test_anno_trf_result
# assert_equal $(fn_md5 $ANSDIR/trf.vcf) $(fn_md5 $OD/trf.vcf)

#                                 repmask
run test_anno_repmask \
    $truv anno repmask -i $INDIR//multi.vcf.gz -o $OD/repmask.vcf -e $INDIR/external/fakeRM.py
if [ $test_anno_repmask ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/repmask.vcf) $(fn_md5 $OD/repmask.vcf)
fi

run test_anno_repmask_err \
    $truv anno repmask -i $INDIR/input1.vcf.gz -o $OD/repmask.vcf -e $INDIR/external/fakeRM.py
if [ $test_anno_repmask_err ]; then
    assert_exit_code 1
fi

#                                 bpovl
run test_anno_bpovl \
    $truv anno bpovl -i $INDIR/input1.vcf.gz \
                     -o $OD/anno_bpovl.jl \
                     -a $INDIR/anno.gtf.gz -p gff --sizemin 2
if [ $test_anno_bpovl ]; then
    assert_exit_code 0
    df_check test_anno_bpovl_result $ANSDIR/anno_bpovl.jl $OD/anno_bpovl.jl
fi

#                                 density
run truvari_anno_density \
    $truv anno density -i $INDIR/input3.vcf.gz \
                       -o $OD/anno_density.jl \
                       -g $INDIR/genome.bed -m $INDIR/mask.bed
if [ $truvari_anno_density ]; then
    assert_exit_code 0
    df_check test_anno_density_result $ANSDIR/anno_density.jl $OD/anno_density.jl
fi

#                                 dpcnt
run test_anno_dpcnt \
    $truv anno dpcnt -i $VCF -o $OD/anno_dpcnt.vcf
if [ $test_anno_dpcnt ]; then
    assert_exit_code 0
    info_tests dpcnt $OD/anno_dpcnt.vcf DPCNT,ADCNT
fi

#                                 lcr
run test_anno_lcr \
    $truv anno lcr -i $VCF -o $OD/anno_lcr.vcf
if [ $test_anno_lcr ]; then
    assert_exit_code 0
    info_tests lcr $OD/anno_lcr.vcf LCR
fi

#                                 grpaf
run test_anno_grpaf \
    $truv anno grpaf -i $INDIR/grpaf.vcf.gz -l $INDIR/grpaf.labels.txt -o $OD/anno_grpaf.vcf
if [ $test_anno_grpaf ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/anno_grpaf.vcf) $(fn_md5 $OD/anno_grpaf.vcf)
fi

run test_anno_grpaf_strict \
    $truv anno grpaf --strict -i $INDIR/grpaf.vcf.gz -l $INDIR/grpaf.labels.txt -o $OD/anno_grpaf.vcf
if [ $test_anno_grpaf_strict ]; then
    assert_exit_code 1
fi

run test_anno_grpaf_subset \
    $truv anno grpaf --tags AF,HWE,ExcHet \
                     -i $INDIR/grpaf.vcf.gz \
                     -l $INDIR/grpaf.labels.txt \
                     -o $OD/anno_grpaf.subtags.vcf
if [ $test_anno_grpaf_subset ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/anno_grpaf.subtags.vcf) $(fn_md5 $OD/anno_grpaf.subtags.vcf)
fi
