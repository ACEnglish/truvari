""" Wrapper around TRF to annotate a VCF """
import sys
import logging
import argparse
import tempfile
from collections import defaultdict

import pysam
import pyfaidx
import truvari
from acebinf import cmd_exe, setup_logging

# Start with just insertions and that sequence
# Eventually you can intersect with known tandem repeat regions as well
# and use those regions (with the edit) to better predict a 'difference' in number of copies
# Or you can just annotate the absolute number of the individual... probably bot

# Also, running single threaded isn't great, but also I don't care at this point...
# It'll be a heavy step, but I'm not judging the software on speed for a bit...


class TRFAnno():
    """ Class for trf annotation """
    TRNAME = tempfile.NamedTemporaryFile().name # Need to remove this
    FANAME = tempfile.NamedTemporaryFile().name # Need to remove this
    TRFCOLS = [("TRF_starts", int),
               ("TRF_ends", int),
               ("TRF_periods", int),
               ("TRF_copies", float),
               ("TRF_consizes", int),
               ("TRF_pctmats", int),
               ("TRF_pctindels", int),
               ("TRF_scores", int),
               ("TRF_As", int),
               ("TRF_Cs", int),
               ("TRF_Gs", int),
               ("TRF_Ts",  int),
               ("TRF_entropies", float),
               ("TRF_repeats", str),
               ("unk1", None),
               ("unk2", None),
               ("unk3", None)]

    def __init__(self, in_vcf, out_vcf="/dev/stdout", executable="trf409.linux64",
                 min_length=50, threshold=0.8, full=False,
                 trf_params="2 7 7 80 10 50 500 -m -f -h -d -ngs",
                 refanno=None, ref=None, pctovl=0.2):
        """ The setup """
        self.in_vcf = in_vcf
        self.out_vcf = out_vcf
        self.executable = executable
        self.min_length = min_length
        self.threshold = threshold
        self.full = full

        if "-ngs" not in trf_params:
            trf_params = trf_params + " -ngs "
        self.trf_params = trf_params

        if refanno and not ref:
            logging.error("-R requires -r")
            exit(1)

        self.refanno = truvari.make_bedanno_tree(refanno) if refanno else None
        self.ref = pyfaidx.Fasta(ref) if ref else None
        self.pctovl = pctovl

        self.n_header = self.edit_header()
        self.cmd = f"{self.executable} {TRFAnno.FANAME} {self.trf_params} > {TRFAnno.TRNAME}"

    def edit_header(self):
        """
        New VCF INFO fields
        """
        header = None
        with pysam.VariantFile(self.in_vcf, 'r') as fh:
            header = fh.header.copy()
        # if intersect_only: Do I want to sometimes turn this off?
        # Probably, actully
        header.add_line(('##INFO=<ID=TRF_repeats,Number=.,Type=String,'
                         'Description="TRF repeat sequences">'))
        header.add_line(('##INFO=<ID=TRF_periods,Number=.,Type=Integer,'
                         'Description="TRF periods">'))
        header.add_line(('##INFO=<ID=TRF_copies,Number=.,Type=Float,'
                         'Description="TRF repeat copies">'))
        header.add_line(('##INFO=<ID=TRF_scores,Number=.,Type=Integer,'
                         'Description="TRF repeat scores">'))
      
        # Take this out
        if self.refanno and self.ref:
            header.add_line(('##INFO=<ID=TRF_diffs,Number=.,Type=Float,'
                             'Description="Alternate allele copy number difference from '
                             'its reference SREP_repeat copy number">'))
        if self.refanno:
            for headinfo in self.refanno[1]:
                header.add_line(headinfo)

        if self.full:
            header.add_line(('##INFO=<ID=TRF_starts,Number=.,Type=Integer,'
                             'Description="TRF starts">'))
            header.add_line(('##INFO=<ID=TRF_ends,Number=.,Type=Integer,'
                             'Description="TRF ends">'))
            header.add_line(('##INFO=<ID=TRF_periods,Number=.,Type=Integer,'
                             'Description="TRF periods">'))
            header.add_line(('##INFO=<ID=TRF_consizes,Number=.,Type=Integer,'
                             'Description="TRF consizes">'))
            header.add_line(('##INFO=<ID=TRF_pctmats,Number=.,Type=Integer,'
                             'Description="TRF pctmats">'))
            header.add_line(('##INFO=<ID=TRF_pctindels,Number=.,Type=Integer,'
                             'Description="TRF pctindels">'))
            header.add_line(('##INFO=<ID=TRF_As,Number=.,Type=Integer,'
                             'Description="TRF As">'))
            header.add_line(('##INFO=<ID=TRF_Cs,Number=.,Type=Integer,'
                             'Description="TRF Cs">'))
            header.add_line(('##INFO=<ID=TRF_Gs,Number=.,Type=Integer,'
                             'Description="TRF Gs">'))
            header.add_line(('##INFO=<ID=TRF_Ts,Number=.,Type=Integer,'
                             'Description="TRF Ts">'))
            header.add_line(('##INFO=<ID=TRF_entropies,Number=.,Type=Float,'
                             'Description="TRF entropies">'))
            # TRF HIT - if refanno, just report if we hit a reference TRF region
            # Might want some pct/ovl -- this is handled... nowhere
            pass
        # TODO: Need to put a source line that says this thing was run with whatever parameters
        return header

    def run_trf(self, seq):
        """
        Runs trf on the ref/alt sequences
        returns list of hits
        """
        def parse_output():
            """
            Parse the outputs from trf, turn to a dictionary
            """
            hits = []
            with open(TRFAnno.TRNAME, 'r') as fh:
                name = fh.readline()
                if name == "": # no hits
                    return hits
                name = name.strip()[1:]
                while True:
                    # If there are multiple, need to parameters for 'take best' or take top N or something
                    # Will need name now that there's ref/alt seq
                    data = fh.readline()
                    if data == "":
                        break
                    if data.startswith("@"):
                        name = data.strip()[1:]
                        continue
                    data = data.strip().split(' ')
                    data = {x[0]: x[1](y) for x, y in zip(TRFAnno.TRFCOLS, data) if not x[0].startswith("unk")}
                    # don't really need until parallel
                    data["TRF_scores"] = int(data["TRF_scores"])
                    hits.append(data)
            return hits

        with open(TRFAnno.FANAME, 'w') as fout:
            #for seq in seqs:
            fout.write(">a\n%s\n" % (seq))
                
        ret = cmd_exe(self.cmd)
        if ret.ret_code != 0:
            logging.error("Couldn't run trf")
            logging.error(str(ret))
            exit(ret.ret_code)
        return parse_output()

    def __old_srep(self):
        """
        legacy code
        """
        srep_hits = None
        # Make srep_hits if we have the bed
        if self.refanno:
            srep_hits = self.refanno[0][entry.chrom].overlap(entry.start, entry.stop)


    def annotate(self, entry, altseq):
        """
        Annotates an Insertion, returns the annotation information
        """
        trf_hits = self.run_trf(altseq)
        best_hit = None
        entry_size = truvari.entry_size(entry)
        for hit in trf_hits:
            size_trf = abs(hit["TRF_ends"] - hit["TRF_starts"]) + 1
            pct = size_trf / entry_size  # The TR that covers the most of the insertion
            if pct < self.threshold:
                continue
            if best_hit is None:
                best_hit = hit
                continue
            if best_hit is None \
               or (hit["TRF_periods"] < best_hit["TRF_periods"] \
               and hit["TRF_copies"] > best_hit["TRF_copies"]):
                best_hit = hit
        return best_hit

    def annotate_entry(self, entry, altseq):
        """
        Annotates an insertion. Insertions are assumed to have a 
        refspan of a single anchor base and full ALT sequence
        Returns the edited entry
        """
        return self.edit_entry(entry, [self.annotate(entry, altseq)])

    def edit_entry(self, entry, trf_annos):
        """
        puts the annos in vcf entry
        return the edited entry
        """
        if not trf_annos:
            return entry
       
        # Let's assume we're only doing the alt allele
        srep_hits = None
        try:
            entry = truvari.copy_entry(entry, self.n_header)
        except TypeError:
            return entry
        n_dat = defaultdict(list)
        for i in trf_annos:
            if i is None:
                continue
            if self.full:
                for key, cnvt in TRFAnno.TRFCOLS:
                    if cnvt:
                        n_dat[key].append(cnvt(i[key]))
            else:
                n_dat["TRF_repeats"].append(i["TRF_repeats"])
                n_dat["TRF_periods"].append(int(i["TRF_periods"]))
                n_dat["TRF_copies"].append(float(i["TRF_copies"]))
                n_dat["TRF_scores"].append(int(i["TRF_scores"]))
        
        # Adding srep hits to n_dat also
        # this is getting kinda choppy with the converters
        # may have over engineered it early
        if srep_hits is not None:
            for i in srep_hits:
                for k,v in i.data.items():
                    n_dat[k].extend(v)
        
        # Can calculate the diffs
        if self.refanno and self.ref and srep_hits and alt_annos: 
            lookup = {}
            diffs = []
            # make a lookup of the trf_annos from the reference
            for i in srep_hits:
                i = i.data
                lookup[i["SREP_repeats"][0]] = i["SREP_copies"][0]
            
            has_diff = False
            # Subtract each alt_repeat from its corresponding ref_repeat
            for alt_repeat, alt_copy in zip(n_dat["TRF_repeats"], n_dat["TRF_copies"]):
                if alt_repeat in lookup:
                    has_diff = True
                    diffs.append(alt_copy - lookup[alt_repeat])
                else:
                    diffs.append(None)
            
            if has_diff:
                n_dat["TRF_diffs"] = diffs

        for key in n_dat:
            entry.info[key] = n_dat[key]
        
                    
        return entry

    def run(self):
        """
        The work
        """
        logging.info("Annotating VCF")
        # should probably 'batch' it instead of running individually
        with pysam.VariantFile(self.in_vcf) as vcf, \
                pysam.VariantFile(self.out_vcf, 'w', header=self.n_header) as output:
            for entry in vcf:
                sz = truvari.entry_size(entry)
                if sz >= self.min_length \
                   and truvari.entry_variant_type(entry) == "INS":
                    entry = self.annotate_entry(entry, entry.alts[0])
                output.write(entry)

def parse_args(args):
    """
    Pull the command line parameters
    """
    def restricted_float(x):
        x = float(x)
        if x < 0.0 or x > 1.0:
            raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
        return x
    parser = argparse.ArgumentParser(prog="trf", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, default="/dev/stdin",
                        help="VCF to annotate (stdin)")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output filename (stdout)")
    parser.add_argument("-e", "--executable", type=str, default="trf409.linux64",
                        help="Path to tandem repeat finder (%(default)s)")
    parser.add_argument("-f", "--full", action="store_true",
                        help="Write full trf output to entries")
    parser.add_argument("-m", "--min-length", type=int, default=50,
                        help="Minimum size of entry to annotate (%(default)s)")
    parser.add_argument("-t", "--threshold", type=restricted_float, default=.8,
                        help="Threshold for pct of allele covered (%(default)s)")
    parser.add_argument("-T", "--trf-params", type=str, default="2 7 7 80 10 50 500 -m -f -h -d -ngs",
                        help="Default parameters to send to trf (%(default)s)")
    parser.add_argument("-R", "--ref-bed", type=str, default=None,
                        help="Reference bed of simple repeats")
    parser.add_argument("-r", "--ref", type=str, default=None,
                        help="Reference fasta file (use with -R)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    setup_logging(args.debug)
    return args


def trf_main(cmdargs):
    """ Main """
    args = parse_args(cmdargs)
    anno = TRFAnno(in_vcf=args.input,
                   out_vcf=args.output,
                   executable=args.executable,
                   min_length=args.min_length,
                   threshold=args.threshold,
                   full=args.full,
                   trf_params=args.trf_params,
                   refanno=args.ref_bed,
                   ref=args.ref)
    anno.run()
    logging.info("Finished trf")


if __name__ == '__main__':
    trf_main(sys.argv[1:])


"""
1) I can't guarantee that TRF alt seq hits are going to happend
    But I'm returning nulls - not good. need to remove I think
    | 

So 1- you can give up on the reference, totally un-needed unti you get to 'denovo mode'
Which at this point you should just abandon until it beocmes a feature request
"""
