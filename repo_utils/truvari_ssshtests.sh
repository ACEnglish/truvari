test -e ssshtest || curl -O https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

source ssshtest

First, need to make a reference of the first 1Mbp of chr20 from grch38

# ------------------------------------------------------------
#                                 version
# ------------------------------------------------------------


# ------------------------------------------------------------
#                                 bench
# ------------------------------------------------------------

On second thought, Im not going to actually test that all the
parameters work. I am just going to test that the base command works
truvari bench --base input.vcf.gz
              --comp input.vcf.gz
              --output test_output/bench_out
              --reference reference.fa
check that files in bench_out are identical to answer_key_bench
There should be unit tests if I really parameters to work.
And as there become tickets/errors/issues, more ssshtests can be
added

# ------------------------------------------------------------
#                                 collapse
# ------------------------------------------------------------

truvari collapse --input
                 --output
                 --collapsed-output
                 --reference

# ------------------------------------------------------------
#                                 consistency
# ------------------------------------------------------------
Just need multiple VCFs.. this one will be pretty easy

# ------------------------------------------------------------
#                                 anno
# ------------------------------------------------------------

So for these, I can just do a bcftools query to compare. 
I can run them all into a single base VCF with the answer,
Then when I run this test, `bcftools query` compare the two,
If there are any differences then I fail the test. 
These tests would just be validating that
1 - the commands work without throwing errors
2 - the command outputs stay the same

If a value needs to be changed, it will be changed in the answer key

Downfalls 
- you can just change the answer_key to force passing.

#                                 gcpct
vcf_anno_compare answer_key.vcf.gz quiz.vcf.gz INFOA INFOB INFOC etc

#                                 gtcnt
vcf_anno_compare INFO/A INFO/B INFO/C
#                                 remap
#                                 numneigh
#                                 hompct

These can also be tested with the rest of the anno framework. 
However I need extra tools to make it happen.

#                                 repmask
#                                 trf

This will need a custom checker
#                                 grm


# ------------------------------------------------------------
#                                 truv2df
# ------------------------------------------------------------
I can probably use the same custom checker from grm for this



# ------------------------------------------------------------
#                                 stats
# ------------------------------------------------------------
This needs to be deprecated


