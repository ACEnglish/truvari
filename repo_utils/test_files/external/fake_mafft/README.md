Fake MAFFT

We need to run MAFFT from inside of truvari when doing functional tests. However, MAFFT is a larger executable and isn't
that easy to install, especially on the instances that the github actions run on. Therefore, we have a fake mafft.

This works by saving real mafft results and using this fake mafft to simply lookup what the correct result should be.

Build the results by running the ssshtest by setting the environment variable `PHAB_WRITE_MAFFT=1`. 
This will make files named `fm_<md5sum>.msa` into `repo_utils/test_files/external/fake_mafft/lookup/`.
