import sys
import argparse
from tinydb import TinyDB, Query

USAGE = """Query a TinyDB
Use 'q' for the query. Be sure to wrap the query in \"quotes\"

Example:
((q.SVTYPE == 'DEL') & (abs(q.SVLEN) >= 100))
"""
parser = argparse.ArgumentParser(prog="query_tinydb", description=USAGE,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-d", "--database", type=str, required=True,
                    help="Name of the database")
parser.add_argument("query", metavar="QUERY", type=str,
                    help="Query to execute")
args = parser.parse_args()
#default returns counts - could add a parameter to report properties XYZ of each
db = TinyDB(args.database)
q = Query()
print args.query
try:
    expand = "v = len(db.search(%s))" % args.query
    exec(expand)
except Exception:
    print 'invalid query'
finally:
    print "Count", v

