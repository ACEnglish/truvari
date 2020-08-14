"""
Over multiple vcfs, calculate their intersection/consistency.

Calls will match between VCFs if they have a matching key of:
    CHROM:POS ID REF ALT
"""
import io
import gzip
import bisect
import argparse
import itertools

from collections import defaultdict, namedtuple, Counter


def parse_vcf(fn):
    """
    Simple vcf reader
    """
    VCFLine = namedtuple("VCFline", "CHROM POS ID REF ALT QUAL FILT INFO FORMAT SAMPLES")
    if fn.endswith(".gz"):
        fh = io.TextIOWrapper(gzip.open(fn))
    else:
        fh = open(fn, 'r')
    for line in fh:
        if line.startswith("#"):
            continue
        data = line.strip().split('\t')
        yield VCFLine(*data[:9], SAMPLES=data[9:])


class hash_list(list):

    """
    A list that's hashable
    """

    def __hash__(self):
        """
        Only method needed
        """
        return hash(" ".join(self))


def entry_key(entry):
    """
    Turn a vcf entry into a key
    """
    key = "%s:%s %s %s %s" % (entry.CHROM, entry.POS, entry.ID, entry.REF, str(entry.ALT))
    return key


def read_files(allVCFs):
    """
    Load all vcfs and count their number of entries
    """
    # call exists in which files
    call_lookup = defaultdict(list)
    # total number of calls in a file
    file_abscnt = defaultdict(float)
    for vcfn in allVCFs:
        v = parse_vcf(vcfn)
        # disallow intra vcf duplicates
        seen = {}
        for entry in v:
            key = entry_key(entry)
            if key in seen:
                continue
            seen[key] = True
            bisect.insort(call_lookup[key], vcfn)
            file_abscnt[vcfn] += 1

    return call_lookup, file_abscnt


def create_file_intersections(allVCFs):
    """
    Generate all possible intersections of vcfs
    """
    count_lookup = {}
    combo_gen = [x for l in range(1, len(allVCFs) + 1) for x in itertools.combinations(allVCFs, l)]
    for files_combo in combo_gen:
        files_combo = hash_list(files_combo)
        files_combo.sort()
        count_lookup[files_combo] = 0
    return count_lookup

def parse_args(args):
    """ parse args """
    parser = argparse.ArgumentParser(prog="consistency_report", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("allVCFs", metavar='VCFs', nargs='+',
                        help="VCFs to intersect")
    args = parser.parse_args(args)
    return args

def make_consistency_overlap(count_lookup, file_abscnt, allVCFs):
    """
    # 1 I want to make a key "101010" so that they can be viz'd easier
    # 2 - I want to sort the count_lookup by their value so that we output them in order
    # The group
    """
    all_consistency = Counter()
    all_overlap = []

    for combo, value in sorted(count_lookup.items(), key=lambda i: (i[1], i[0]), reverse=True):
        # There are no calls here, so we just ignore it... But I think I want to keep that information
        cur_data = []
        if value == 0:
            continue

        my_group = ["0"] * len(allVCFs)
        m_cnt = 0
        for j in combo:
            my_group[allVCFs.index(j)] = "1"
            m_cnt += 1
        cur_data.append("".join(my_group))
        cur_data.append(value)
        cur_data.append([])

        all_consistency[m_cnt] += value
        for fkey in combo:
            if file_abscnt[fkey] > 0:
                cur_data[-1].append("%.2f%%" % (count_lookup[combo] / file_abscnt[fkey] * 100))
        all_overlap.append(cur_data)
    return all_consistency, all_overlap

def write_report(total_unique_calls, allVCFs, file_abscnt, all_consistency, all_overlap):
    """
    Write the report
    """
    print("#\n# Total %d calls across %d VCFs\n#" % (total_unique_calls, len(allVCFs)))
    print("#File\tNumCalls")
    for fn in allVCFs:
        print("%s\t%d" % (fn, file_abscnt[fn]))

    print("#\n# Summary of consistency\n#")
    print("#VCFs\tCalls\tPct")

    for i in sorted(all_consistency.keys(), reverse=True):
        print("%d\t%d\t%.2f%%" % (i, all_consistency[i], all_consistency[i] / total_unique_calls * 100))

    print("#\n# Breakdown of VCFs' consistency\n#")
    print("#Group\tTotal\tTotalPct\tPctOfFileCalls")
    for my_group, value, combo in all_overlap:
        c_text = ""
        pos = 0
        for i in my_group:
            if i == '1':
                c_text += combo[pos] + " "
                pos += 1
            else:
                c_text += "0% "
        print("%s\t%d\t%.2f%%\t%s" % (my_group, value, value / total_unique_calls * 100, c_text))

def consistency_main(args):
    """
    Run the program
    """
    args = parse_args(args)

    call_lookup, file_abscnt = read_files(args.allVCFs)

    count_lookup = create_file_intersections(args.allVCFs)

    for key in call_lookup:
        count_lookup[hash_list(call_lookup[key])] += 1

    all_consistency, all_overlap = make_consistency_overlap(count_lookup, file_abscnt, args.allVCFs)

    total_unique_calls = sum(all_consistency.values())

    write_report(total_unique_calls, args.allVCFs, file_abscnt, all_consistency, all_overlap)
