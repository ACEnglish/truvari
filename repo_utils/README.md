Truvari Testing Infrastructure
===============================

Truvari uses [ssshtests](https://github.com/ryanlayer/ssshtest/) to manage functional tests, [pylint](https://pylint.pycqa.org/en/latest/) to keep the code pretty. Test coverage is automatically generated using [coverage.py](https://coverage.readthedocs.io/en/6.5.0/)

CI/CD is performed with github actions, one for functional tests and one for pylint.

Functional Tests
=================
There are two main scripts. The first is `truvari_ssshtest.sh` which is the main entry point for running the tests. This file uses `source` to call the other pieces. The second, `setup_test.sh`, can be thought of a config file and also holds shared functions to be used by the individual tests. The main entry point script runs from the base directory of the code so that it can automatically call the truvari source code directly without needing to install the package.

The actual tests of the sub-commands are organized under the directory `sub_tests/`.

Input files for testing (e.g. VCFs, references) are stored inside of `test_files/`. The output files for comparison are inside of `test_files/answer_key/`. The `setup_test.sh` sets two variables to point to these directories, `$INDIR` and `$ANSDIR`. Each run of the functional tests should write its outputs to the `$OD` variable, which points to `test_results/`. Note that the `$OD` is removed and remade at the beginning of each run of `truvari_ssshtest.sh`

Doctests
--------
Doctests are run inside of functional tests, but are kind of a special case. See `repo_utils/run_doctests.py`, which is run by `repo_utils/sub_tests/doctests.sh`.

Coverage
--------
The `.coveragerc` file is the config. You can see how the `coverage` and `truvari` commands are combined to create the `$truv` variable which serves as the command run by nearly all the functional tests. At the bottom of the `truvari_ssshtest.sh` are the commands to generate the coverage reports and automatically make the coverage badge for `README.md`

Pylint
======

The `.pylintrc` file is the config. The script `repo_utils/pylint_maker.py` will do all the pylint work and automatically make the pylint score badge for `README.md` 

Tricks
======
The hardest part of maintaining functional tests is creating/maintaining the inputs and outputs. There are a few trick scripts I've made to help.

- `test_files/make_anno.sh` makes a VCF with all INFO fields of the multiple `truvari anno` fields.
- `test_files/add_dp.py` was used to add fake depth information in VCFs.

When functional tests fail because of files not matching up (e.g. md5sum aren't equal), I usually manually inspect with something like `vimdiff` on the `test_results/` and the `repo_utils/test_files/answer_key/`. If I'm CERTAIN that the `test_results/` is correct, then I copy that into `repo_utils/test_files/answer_key/`.
