

"""
Use a truvari benchmark results that's been turned to a TinyDB
to create a pretty report
"""
import sys
from tinydb import *
from collections import Counter
db = TinyDB(sys.argv[1])
q = Query()

def count_by(key, docs):
    """
    get a count for all the keys
    """
    cnt = Counter()
    for x in docs:
        cnt[x[key]] += 1
    for key in cnt:
        print key, cnt[key]

def count_by_multi(keys, docs):
    """
    Go through all the docs and count by multiple keys
    """
    cnt = Counter()
    for x in docs:
        mkey = []
        for k in keys:
            mkey.append(x[k])
        cnt[":".join(mkey)] += 1

    for key in cnt:
        print key, cnt[key]

tp_base = db.search(q.BENCH == 'tp-base')
fn = db.search(q.BENCH == 'fn')

#OverallNumbers 
print "TP", len(tp_base)
print "FN", len(fn)

print "TP Illumina",  len(db.search((q.BENCH == 'tp-base') & (q.Illcalls != 0)))
print "FN Illumina",  len(db.search((q.BENCH == 'fn') & (q.Illcalls != 0)))

print "TP by size"
count_by("sizecat", tp_base)
print "FN by size"
count_by("sizecat", fn)

print "TP by type"
count_by("SVTYPE", tp_base)
print "FN by type", 
count_by("SVTYPE", fn)

print "TP by Type + Size"
count_by_multi(["sizecat", "SVTYPE"], tp_base)
print "FN by Type + Size"
count_by_multi(["sizecat", "SVTYPE"], fn)

exit(0)

#Breakdowns by Tech
#Breakdown by Size 
#Breakdown by Type
#Breakdown by Tech + Size
#Breakdown by Tech + Type
#Breakdown by Tech + Size + Type

base_total = len(db.search((q.BENCH == 'tp-base') | (q.BENCH == 'fn')))
call_total = len(db.search((q.BENCH == 'tp-call') | (q.BENCH == 'fn')))
