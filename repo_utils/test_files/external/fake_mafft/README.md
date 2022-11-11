Fake MAFFT

We need to run MAFFT from inside of truvari when doing functional tests. However, MAFFT is a larger executable and isn't
that easy to install, especially on the instances that the github actions run on. Therefore, we have a fake mafft.

This works by saving real mafft results and using this fake mafft to simply lookup what the correct result should be.

For example, say we build a test for `truvari phab`. We start by actually running the command.
Next, we point `make_fake_mafft_results.py` to the `haps.fa`. `haps.fa` is the input given to MAFFT by truvari.
`make_fake_mafft_results.py` will hash the `haps.fa` so that when that exact same input is seen later during testing, it
knows which `msa.fa` to output. 

The lookup between md5sums and out files is stored in `mafft_results.json`. Be sure that the values in the json are
pointing relative to the fake `mafft` script via `ln -s`. This allows us to run fake mafft from any directory.

Also, be sure in the ssshtest to update the $PATH to point to this directory so that `mafft` is found in the
environment.
