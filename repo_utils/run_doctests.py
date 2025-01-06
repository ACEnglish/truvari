"""
Runs doctests. Only works when run from the repository's root directory
"""
import os
import sys
import doctest

# Use the current truvari, not any installed libraries
sys.path.insert(0, os.getcwd())

from truvari import (
    bench,
    comparisons,
    msatovcf,
    utils,
    variants,
    vcf2df,
)
    

def tester(module):
    """
    Runs test on module and will set the exit status if there are any failures
    """
    ret = doctest.testmod(module)
    print(f"{ret.attempted} tests for {module.__name__}. {ret.failed} failed")
    return ret.failed

fails = 0
fails += tester(variants)
fails += tester(comparisons)
fails += tester(utils)
fails += tester(vcf2df)
fails += tester(msatovcf)

os.remove("log.txt")
sys.exit(fails)
