"""
For GIAB calls, parse the formats and tags to make summaries about:

PerSizeBin Info/sizecat
    Filter flag
    INFO:
        SVType  :   count
        ClusterIDs (Indiv_Tech_Method_id) : Count per
        PBcalls : count : count
        Illcalls : count : count
        TenXcalls : count : count
        CGcalls: count count
        HG2count: count : count
        HG3count: count : count
        HG4count: count : count
        NumTechs: count : count
        MultiTech: count
        MendelianError: count
        REPTYPE: type : count
        HG2_GT : type : count
        HG3_GT : type : count
        HG4_GT : type : count
        HG2_GTcons1 : type : count
        HG3_GTcons1 : type : count
        HG4_GTcons1 : type : count

in this size regime, the number of calls with 5 techs is 10
50-100bp : NumTechs : 5 : 10

Use Truvari's --giabreport for a better looking output
"""
import sys
import vcf
from collections import defaultdict, Counter

v = vcf.Reader(open(sys.argv[1]))
summary = defaultdict()
for entry in v:
    sizecat = entry.INFO["sizecat"]
    if sizecat not in summary:
        summary[sizecat] = {}
        summary[sizecat]["Total"] = 0
        summary[sizecat]["Program"] = Counter()
        summary[sizecat]["Filter"] = Counter()
        summary[sizecat]["SVType"] = Counter()
        summary[sizecat]["PBcalls"] = Counter()
        summary[sizecat]["Illcalls"] = Counter()
        summary[sizecat]["TenXcalls"] = Counter()
        summary[sizecat]["CGcalls"] = Counter()
        summary[sizecat]["HG2count"] = Counter()
        summary[sizecat]["HG3count"] = Counter()
        summary[sizecat]["HG4count"] = Counter()
        summary[sizecat]["NumTechs"] = Counter()
        summary[sizecat]["MultiTech"] = 0
        summary[sizecat]["MendelianError"] = 0
        summary[sizecat]["REPTYPE"] = Counter()
        summary[sizecat]["HG2_GT"] = Counter()
        summary[sizecat]["HG3_GT"] = Counter()
        summary[sizecat]["HG4_GT"] = Counter()
        summary[sizecat]["HG2_GTcons1"] = Counter()
        summary[sizecat]["HG3_GTcons1"] = Counter()
        summary[sizecat]["HG4_GTcons1"] = Counter()

    summary[sizecat]["Total"] += 1
    if len(entry.FILTER):
        filt = ",".join(sorted(entry.FILTER))
    else:
        filt = "PASS"
    # This assumes only allows a program to contribute to a call once
    # Could stand to do some cleaning of this, also
    for prog in set([x.split('_')[2] for x in entry.INFO["ClusterIDs"].split(':')]):
        summary[sizecat]["Program"][prog] += 1
    summary[sizecat]["Filter"][filt] += 1
    summary[sizecat]["SVType"][entry.INFO["SVTYPE"]] += 1
    summary[sizecat]["PBcalls"][entry.INFO["PBcalls"]] += 1
    summary[sizecat]["Illcalls"][entry.INFO["Illcalls"]] += 1
    summary[sizecat]["TenXcalls"][entry.INFO["TenXcalls"]] += 1
    summary[sizecat]["CGcalls"][entry.INFO["CGcalls"]] += 1
    summary[sizecat]["HG2count"][entry.INFO["HG2count"]] += 1
    summary[sizecat]["HG3count"][entry.INFO["HG3count"]] += 1
    summary[sizecat]["HG4count"][entry.INFO["HG4count"]] += 1
    summary[sizecat]["NumTechs"][entry.INFO["NumTechs"]] += 1

    if entry.INFO["MultiTech"] == "TRUE":
        summary[sizecat]["MultiTech"] += 1

    if entry.INFO["MendelianError"] == "TRUE":
        summary[sizecat]["MendelianError"] += 1
    summary[sizecat]["REPTYPE"][entry.INFO["REPTYPE"]] += 1

    summary[sizecat]["HG2_GT"][entry.genotype("HG002")["GT"]] += 1
    summary[sizecat]["HG3_GT"][entry.genotype("HG003")["GT"]] += 1
    summary[sizecat]["HG4_GT"][entry.genotype("HG004")["GT"]] += 1
    summary[sizecat]["HG2_GTcons1"][entry.genotype("HG002")["GTcons1"]] += 1
    summary[sizecat]["HG3_GTcons1"][entry.genotype("HG003")["GTcons1"]] += 1
    summary[sizecat]["HG4_GTcons1"][entry.genotype("HG004")["GTcons1"]] += 1

for szc in ["20to49", "50to99", "100to299", "300to999", "gt1000"]:
    if szc not in summary:
        continue
    print szc, "\tTotal\t", summary[szc]["Total"], "\tMultiTech\t", summary[szc]["MultiTech"], "\tMendelianError\t", summary[szc]["MendelianError"]
    for key in ["Filter", "SVType", "PBcalls", "Illcalls", "TenXcalls", "CGcalls", "HG2count", "HG3count", "HG4count",
                "NumTechs", "REPTYPE", "HG2_GT", "HG3_GT", "HG4_GT", "HG2_GTcons1", "HG3_GTcons1", "HG4_GTcons1",
                "Program"]:
        print key
        nkeys = sorted(summary[szc][key].keys())
        print "\t", "\t".join([str(x) for x in nkeys])
        o = []
        for i in nkeys:
            o.append(str(summary[szc][key][i]))
        print "\t", "\t".join(o)
