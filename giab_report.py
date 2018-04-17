

"""
Use a truvari benchmark results that's been turned to a TinyDB
to create a pretty report

You know what, If I just pop all of this stuff into truvari
I can pretty quickly add the work to make this into a report.
Just turn the two vcfs into lists of dicts
add an option for -g/--giabreport
Then without tinydb, I can make these same table.
"""
import sys
from tinydb import *
from collections import Counter, defaultdict
db = TinyDB(sys.argv[1])
q = Query()

def count_by(key, docs):
    """
    get a count for all the keys
    key is a dictionary of count and their order
    """
    main_key = key.keys()[0]
    cnt = Counter()
    for x in docs:
        cnt[x[main_key]] += 1
    for k in key[main_key]:
        print "%s\t%s" % (k, cnt[k])

def twoxtable(key1, key2, docs):
    """
    Parse any set of docs and create a 2x2 table
    """
    main_key1 = key1.keys()[0]
    main_key2 = key2.keys()[0]
    cnt = defaultdict(Counter)
    for x in docs:
        cnt[x[main_key1]][x[main_key2]] += 1
    
    print ".\t"+"\t".join([str(i) for i in key1[main_key1]])
    for y in key2[main_key2]:
        o = [str(y)]
        for x in key1[main_key1]:
            o.append(str(cnt[x][y]))
        print "\t".join(o)

def print_table(cnt, default=0):
    """
    print a 2d table
    """
    x = cnt.keys()
    y = set()
    for i in x:
        y.update(cnt[i].keys())
    y = list(y)
    for j in y:
        o = [str(j)]
        for i in x:
            o.append(str(cnt[i][j]))
        print "\t".join(o)
        
def collapse_techs(docs):
    """
    Make a new annotation about the presence of techs inplace
    called techs
    Illcalls
    PBcalls
    CGcalls
    TenXcalls
    """
    calls = ["Illcalls", "PBcalls", "CGcalls", "TenXcalls"]
    for d in docs:
        new_anno = []
        for i in calls:
            if d[i] > 0:
                new_anno.append(i.rstrip("calls"))
        d["techs"] = "+".join(new_anno)
        
tp_base = db.search(q.BENCH == 'tp-base')
fn = db.search(q.BENCH == 'fn')

size_keys = {"sizecat": ["50to99", "100to299", "300to999", "gt1000"]}
svtype_keys = {"SVTYPE": ["DEL", "INS", "COMPLEX"]}
tech_keys = {"techs": ["I+PB+CG+TenX", "I+PB+CG", "I+PB+TenX", "PB+CG+TenX",
                       "I+PB", "I+CG", "I+TenX", "PB+CG", "PB+TenX", "CG+TenX",
                       "I", "PB", "CG","TenX"]}
rep_keys = {"REPTYPE": [ "SIMPLEDEL", "SIMPLEINS", "DUP", "SUBSDEL", "SUBSINS", "CONTRAC"]}

#OverallNumbers 
print "TP\t%s" % (len(tp_base))
print "FN\t%s" % (len(fn))
print
print "TP_size"
count_by(size_keys, tp_base)
print "FN_size"
count_by(size_keys, fn)
print
print "TP_type"
count_by(svtype_keys, tp_base)
print "FN_type"
count_by(svtype_keys, fn)
print
print "TP_Type+Size"
twoxtable(svtype_keys, size_keys, tp_base)
print "FN_Type+Size"
twoxtable(svtype_keys, size_keys, fn)
print
print "TP_REPTYPE"
count_by(rep_keys, tp_base)
print "FN_REPTYPE"
count_by(rep_keys, fn)
print
print "TP_size+REPTYPE"
twoxtable(size_keys, rep_keys, tp_base)
print "FN_size+REPTYPE"
twoxtable(size_keys, rep_keys, fn)
print
#Need to add an annotation of collapsed techs
collapse_techs(tp_base)
collapse_techs(fn)
print "TP_Tech"
count_by(tech_keys, tp_base)
print "FN_Tech"
count_by(tech_keys, fn)
print
print "TP_Size+Tech"
twoxtable(size_keys, tech_keys, tp_base)
print "FN_Size+Tech"
twoxtable(size_keys, tech_keys, fn)
print
print "TP_Type+Tech"
twoxtable(svtype_keys, tech_keys, tp_base)
print "FN_Type+Tech"
twoxtable(svtype_keys, tech_keys, fn)



#the harder one - size by type by tech

exit(0)

#Breakdowns by Tech
#Breakdown by Size 
#Breakdown by Type
#Breakdown by Tech + Size
#Breakdown by Tech + Type
#Breakdown by Tech + Size + Type

base_total = len(db.search((q.BENCH == 'tp-base') | (q.BENCH == 'fn')))
call_total = len(db.search((q.BENCH == 'tp-call') | (q.BENCH == 'fn')))
