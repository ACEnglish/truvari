Fake MAFFT

We need to run MAFFT from inside of truvari when doing functional tests. However, MAFFT is a larger executable and isn't
that easy to install, especially on the instances that the github actions run on. Therefore, we have a fake mafft.

This works by saving real mafft results and using this fake mafft to simply lookup what the correct result should be.

Build the results by first soft-linking the directory where the phab results sit with `ln -s`

Next, run `python make_fake_mafft_results.py soft-link-dir/haps.fa soft-link-dir/*/haps.fa <etc> > mafft_results.json`

This script will warn you when there are conflicting hashes. I assume that's okay because identical inputs make
identical outputs? 

Also, be sure in the ssshtest to update the $PATH to point to this directory so that `mafft` is found in the
environment. See `sub_tests/phab.sh` for an example
