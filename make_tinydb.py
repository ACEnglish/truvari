import os
import sys
import argparse

import vcf
from tinydb import TinyDB, Query


#By default, don't include the filtered calls - just to speed 
#things along

#check valid directory first
def make_entries(vcf_reader):
    """
    Turn all vcf entries into documents for the database
    """
    for entry in vcf_reader:
        data = {}
        for key in entry.INFO:
            data[key] = entry.INFO[key]
        for sample in entry.samples:
            for fmt in sample.data._fields:
                data[sample.sample + "_" + fmt] = sample[fmt]
        data["CHROM"] = entry.CHROM
        data["POS"] = entry.POS
        data["ID"] = entry.ID
        data["FILTER"] = ";".join(entry.FILTER)
        data["QUAL"] = entry.QUAL
        yield data
USAGE = "Turn vcfs in a truvari result into a TinyDB for easier queries"
parser = argparse.ArgumentParser(prog="make_tinydb", description=USAGE,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-d", "--database", type=str, required=True,
                    help="Name of the database")
parser.add_argument("-i", "--input", type=str, required=True,
                    help="Truvari results directory to parse")
args = parser.parse_args()
if os.path.exists(args.database):
    print "%s already exists!" % args.database
    exit(1)

all_files = ["fn.vcf", "fp.vcf", "tp-base.vcf", "tp-call.vcf"]
for i in all_files:
    if not os.path.exists(os.path.join(args.input, i)):
        print "Directory %s missing %s" % (args.input, i)
        exit(1)

db = TinyDB(args.database)

#I could do some stuff around the lookup id to return - for this call look up it's comparison call..
# But I don't want to, yet

all_keys = set()
all_data = []
for fn in all_files:
    v = vcf.Reader(filename=os.path.join(args.input, fn))
    for data in make_entries(v):
        data["BENCH"] = fn.split('.')[0]
        all_keys.update(data.keys())
        all_data.append(data)
db.insert_multiple(all_data)
db.close()
