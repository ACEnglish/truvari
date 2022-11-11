"""
Makes a lookup from a haps.fa to the msa.fa output

This allows us to 'save' the mafft results and fake its functionality
"""
import os
import sys
import glob
import json
import hashlib
# Provide all the haps files
all_haps = sys.argv[1:]

lookup = {}
for haps in all_haps:
    # For every haps.fa, read it in, get a hash key, save hash key to file's msa.fa
    with open(haps, 'r') as fh:
        data = fh.read()
        m_key = hashlib.md5(data.encode('utf-8')).hexdigest()
        if m_key in lookup:
            print('error', haps)
        lookup[m_key] = os.path.join(os.path.dirname(haps), 'msa.fa')

print(json.dumps(lookup, indent=4))
    
