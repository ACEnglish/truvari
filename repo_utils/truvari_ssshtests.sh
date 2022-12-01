# Functional tests
test -e ssshtest || curl -O https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
source ssshtest

# Run from truvari's base directory
cd "$( dirname "${BASH_SOURCE[0]}" )"/../

TESTSRC=repo_utils/
source $TESTSRC/setup_test.sh

# Reset test results
rm -rf $OD
mkdir -p $OD

source $TESTSRC/sub_tests/anno.sh
source $TESTSRC/sub_tests/bench.sh
source $TESTSRC/sub_tests/collapse.sh
source $TESTSRC/sub_tests/consistency.sh
source $TESTSRC/sub_tests/divide.sh
source $TESTSRC/sub_tests/doctests.sh
source $TESTSRC/sub_tests/entry_main.sh
source $TESTSRC/sub_tests/phab.sh
source $TESTSRC/sub_tests/refine.sh
source $TESTSRC/sub_tests/segment.sh
source $TESTSRC/sub_tests/stratify.sh
source $TESTSRC/sub_tests/vcf2df.sh
source $TESTSRC/sub_tests/version.sh


# Don't generate coverage when doing subset of tests
if [ -z "$1" ]; then
    printf "\n${BOLD}generating test coverage reports${NC}\n"
    coverage combine
    coverage report --include=truvari/*
    coverage html --include=truvari/* -d $OD/htmlcov/
    coverage json --include=truvari/* -o $OD/coverage.json
    python3 repo_utils/coverage_maker.py $OD/coverage.json
fi
