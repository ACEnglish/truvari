# ------------------------------------------------------------
#                                 anno
# ------------------------------------------------------------
VCF=$INDIR/variants/multi.vcf.gz
REF=$INDIR/references/reference.fa

#                                 hompct
run anno_hompct \
    $truv anno hompct -i $VCF -o $OD/anno_hompct.vcf
if [ $anno_hompct ]; then
    assert_exit_code 0
    info_tests hompct $OD/anno_hompct.vcf HOMPCT
fi

#                                 remap
run anno_remap \
    $truv anno remap -H 10 $VCF -r $REF -o $OD/anno_remap.vcf
if [ $anno_remap ]; then
    assert_exit_code 0
    info_tests remap $OD/anno_remap.vcf REMAP,REMAPHits
fi

#                                 gcpct
run anno_gcpct \
    $truv anno gcpct $VCF -r $REF -o $OD/anno_gcpct.vcf
if [ $anno_gcpct ]; then
    assert_exit_code 0
    info_tests gcpct $OD/anno_gcpct.vcf GCPCT
fi

#                                 gtcnt
run anno_gtcnt \
    $truv anno gtcnt $VCF -o $OD/anno_gtcnt.vcf
if [ $anno_gtcnt ]; then
    assert_exit_code 0
    info_tests gtcnt $OD/anno_gtcnt.vcf GTCNT
fi

#                                 numneigh
run anno_numneigh \
    $truv anno numneigh $VCF -o $OD/anno_numneigh.vcf
if [ $anno_numneigh ]; then
    assert_exit_code 0
    info_tests numneigh $OD/anno_numneigh.vcf NumNeighbors,NeighId
fi

#                                 svinfo
run anno_svinfo \
    $truv anno svinfo $VCF -o $OD/anno_svinfo.vcf
if [ $anno_svinfo ]; then
    assert_exit_code 0
    info_tests svinfo $OD/anno_svinfo.vcf SVTYPE,SVLEN
fi

#                                 grm
run anno_grm \
    $truv anno grm -i $INDIR/variants/input2.vcf.gz -r $REF -o $OD/grm.jl -t 2
if [ $anno_grm ]; then
    assert_exit_code 0
    df_check anno_grm $ANSDIR/anno/grm.jl $OD/grm.jl
fi

#                                 trf
run anno_trf \
    $truv anno trf -i $INDIR/variants/multi.vcf.gz \
                   -r $INDIR/beds/repeats.adotto.bed.gz \
                   -f $INDIR/references/reference.fa \
                   -e $INDIR/external/trf  \
                   -m 5 \
                   -o $OD/trf.vcf
if [ $anno_trf ]; then
    assert_exit_code 0
fi

run anno_trf_reg \
    $truv anno trf -i $INDIR/variants/multi.vcf.gz \
                   -r $INDIR/beds/repeats.adotto.bed.gz \
                   -f $INDIR/references/reference.fa \
                   -e $INDIR/external/trf  \
                   -m 5 -R \
                   -o $OD/trf.reg.vcf
if [ $anno_trf_reg ]; then
    assert_exit_code 0
fi

run anno_badparam \
    $truv anno trf -i $INDIR/input_null.vcf \
                   -r $INDIR/simplerepeat_null.bed \
                   -f $INDIR/reference.fa \
                   -e $INDIR/external/trf_dne  \
                   -o $OD/trf_null.vcf
if [ $anno_badparam ]; then
    assert_exit_code 1
fi

# TRF isn't deterministic for some reason, so it can give a different answer in action?
# run anno_trf_result
# assert_equal $(fn_md5 $ANSDIR/trf.vcf) $(fn_md5 $OD/trf.vcf)

#                                 repmask
run anno_repmask \
    $truv anno repmask -i $INDIR/variants/multi.vcf.gz -o $OD/repmask.vcf -e $INDIR/external/fakeRM.py
if [ $anno_repmask ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/anno/repmask.vcf) $(fn_md5 $OD/repmask.vcf)
fi

run anno_repmask_err \
    $truv anno repmask -i $INDIR/variants/input1.vcf.gz -o $OD/repmask.vcf -e $INDIR/external/fakeRM.py
if [ $anno_repmask_err ]; then
    assert_exit_code 1
fi

#                                 bpovl
run anno_bpovl \
    $truv anno bpovl $INDIR/variants/input1.vcf.gz \
                     -o $OD/anno_bpovl.jl \
                     -a $INDIR/misc/anno.gtf.gz -p gff --sizemin 2
if [ $anno_bpovl ]; then
    assert_exit_code 0
    df_check anno_bpovl $ANSDIR/anno/anno_bpovl.jl $OD/anno_bpovl.jl
fi

#                                 density
run anno_density \
    $truv anno density $INDIR/variants/input3.vcf.gz \
                       -o $OD/anno_density.jl \
                       -g $INDIR/beds/genome.bed -m $INDIR/beds/mask.bed
if [ $anno_density ]; then
    assert_exit_code 0
    df_check anno_density $ANSDIR/anno/anno_density.jl $OD/anno_density.jl
fi

#                                 dpcnt
run anno_dpcnt \
    $truv anno dpcnt $VCF -o $OD/anno_dpcnt.vcf
if [ $anno_dpcnt ]; then
    assert_exit_code 0
    info_tests dpcnt $OD/anno_dpcnt.vcf DPCNT,ADCNT
fi

#                                 lcr
run anno_lcr \
    $truv anno lcr $VCF -o $OD/anno_lcr.vcf
if [ $anno_lcr ]; then
    assert_exit_code 0
    info_tests lcr $OD/anno_lcr.vcf LCR
fi

#                                 grpaf
run anno_grpaf \
    $truv anno grpaf $INDIR/variants/grpaf.vcf.gz -l $INDIR/misc/grpaf.labels.txt -o $OD/anno_grpaf.vcf
if [ $anno_grpaf ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/anno/anno_grpaf.vcf) $(fn_md5 $OD/anno_grpaf.vcf)
fi

run anno_grpaf_strict \
    $truv anno grpaf --strict $INDIR/variants/grpaf.vcf.gz -l $INDIR/misc/grpaf.labels.txt -o $OD/anno_grpaf.vcf
if [ $anno_grpaf_strict ]; then
    assert_exit_code 1
fi

run anno_grpaf_subset \
    $truv anno grpaf --tags AF,HWE,ExcHet \
                     $INDIR/variants/grpaf.vcf.gz \
                     -l $INDIR/misc/grpaf.labels.txt \
                     -o $OD/anno_grpaf.subtags.vcf
if [ $anno_grpaf_subset ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/anno/anno_grpaf.subtags.vcf) $(fn_md5 $OD/anno_grpaf.subtags.vcf)
fi

#                           add id 
run anno_addid \
    $truv anno addid $INDIR/variants/multi.vcf.gz -o $OD/addid.vcf
if [ $anno_addid ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/anno/addid.vcf) $(fn_md5 $OD/addid.vcf)
fi


#                           chunks
run anno_chunks $truv anno chunks $INDIR/variants/multi.vcf.gz -o $OD/chunks.bed
if [ $anno_chunks ]; then
    assert_exit_code 0
    assert_equal $(fn_md5 $ANSDIR/anno/chunks.bed) $(fn_md5 $OD/chunks.bed)
fi
