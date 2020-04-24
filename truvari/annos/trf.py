""" Wrapper around TRF to annotate a VCF """
import sys
import logging
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
    MFNAME = "output.trf"
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
                 min_length=50, bestn=5, layered=False, full=False, 
                 trf_params = "2 7 7 80 10 50 500 -m -f -h -d -ngs",
                 refanno=None, ref=None, pctovl=0.2):
        """ The setup """
        self.in_vcf = in_vcf
        self.out_vcf = out_vcf
        self.executable = executable
        self.min_length = min_length
        self.bestn = bestn
        self.layered = layered
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
        self.cmd = self.executable + ' input.fasta ' + self.trf_params + " > " + TRFAnno.MFNAME
    
    def run_trf(self, altseq, refseq=None):
        """
        Runs trf on the ref/alt sequences
        returns {'a': [althitsdict,..], 'r':[refhitsdict]}
        """
        def parse_output():
            """
            Parse the outputs from trf, turn to a dictionary
            """
            hits = defaultdict(list)
            with open(TRFAnno.MFNAME, 'r') as fh:
                while True:
                    # If there are multiple, need to parameters for 'take best' or take top N or something
                    # Will need name now that there's ref/alt seq
                    name = fh.readline()
                    data = fh.readline()
                    if data == "":
                        break
                    name = name.strip()[1:]
                    data = data.strip().split(' ')
                    data = {x[0]:y for x,y in zip(TRFAnno.TRFCOLS, data) if not x[0].startswith("unk")}
                    # don't really need until parallel
                    data["TRF_scores"] = int(data["TRF_scores"])
                    hits[name].append(data)
            return hits
        
        # UNSAFE - need to pick a NamedTemporary file and stick with it
        # same for TRFAnno.MFNAME
        with open("input.fasta", 'w') as fout:
            fout.write(">a\n%s\n" % (altseq))
            # and write refseq...
        ret = cmd_exe(self.cmd)
        if ret.ret_code != 0:
            logging.error("Couldn't run trf")
            logging.error(str(ret))
            exit(ret.ret_code)
        return parse_outputs()

    def insertion_anno(self, entry):
        """
        Annotates an insertion. Insertions are assumed to have a 
        refspan of a single anchor base and full ALT sequence
        """
        trf_hits = None
        srep_hits = None
        # Make srep_hits if we have the bed
        if self.refanno:
            srep_hits = self.refanno[entry.chrom].overlap(entry.start)
        
        # If we also have the reference, we want to do the fancy diff comparison
        if self.refanno and self.ref:
            #This has a whole other thing for making ref/alt
            trf_hits = run_trf(entry.alts[0], and_refseq_whatever)
        else:
            trf_hits = run_trf(entry.alts[0])
        entry = self.annotate_entry(entry, trf_hits, srep_hits)

    def run(self):
        """
        The work
        """
        # should probably 'batch' it instead of running individually
        with pysam.VariantFile(self.in_vcf) as vcf, \
             pysam.VariantFile("/dev/stdout", 'w', header=self.n_header) as output:
            for entry in vcf:
                sz = truvari.entry_size(entry)
                if sz < self.min_size:
                    output.write(entry)
                    continue
                svtype = truvari.entry_variant_type(entry)
                if svtype == "INS":
                    entry = self.insertion_anno(entry)
                elif svtype == "DEL"
                    entry = self.deletion_anno(entry)
                output.write(entry)

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
        header.add_line(('##INFO=<ID=TRF_copies,Number=.,Type=Float,'
                         'Description="TRF repeat copies">'))
        header.add_line(('##INFO=<ID=TRF_scores,Number=.,Type=Integer,'
                         'Description="TRF repeat scores">'))
        if self.refanno and self.ref:
            header.add_line(('##INFO=<ID=TRF_diffs,Number=.,Type=Integer,'
                             'Description="When an entry hits a reference STR, '
                             'calculate the length difference">'))
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
    
    def annotate_entry(self, entry, trf_annos=None, srep_hits=None):
        """
        puts the annos in vcf entry
        """
                
        if not trf_annos and not srep_hits:
            return entry
        
        if len(annos) > self.bestn:
            annos = sorted(annos, key=lambda x: x["TRF_scores"], reverse=True)[:self.bestn]
        
        entry = truvari.copy_entry(entry, self.n_header)
        n_dat = defaultdict(list)
        for i in annos:
            if self.full:
                for key, cnvt in TRFAnno.TRFCOLS:
                    if cnvt:
                        n_dat[key].append(cnvt(i[key]))
            else:
                n_dat["TRF_repeats"].append(i["TRF_repeats"])
                n_dat["TRF_copies"].append(float(i["TRF_copies"]))
                n_dat["TRF_scores"].append(int(i["TRF_scores"]))
        
        for key in n_dat:
            entry.info[key] = n_dat[key]

        return entry

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
    parser.add_argument("-b", "--bestn", type=int, default=5,
                        help="Top trf hits to report (%(default)s)")
    parser.add_argument("-t", "--trf-params", type=str, default="2 7 7 80 10 50 500 -m -f -h -d -ngs",
                        help="Default parameters to send to tr (%(default)s)")
    parser.add_argument("-R", "--ref-bed", type=str, default=None,
                        help="Reference bed of tandem repeat regions")
    parser.add_argument("-r", "--ref", type=str, default=None,
                        help="Reference fasta file (use with -R)")
    parser.add_argument("-p", "--pctovl", type=restricted_float, default=0.2,
                        help="Minimum percent reciprocal overlap (use with -R) (%(default)s)")
    parser.add_argument("-l", "--layered", action="store_true",
                        help="(non-functional)")
    parser.add_argument("--debug", action="store_true", 
                        help="Verbose logging")
    args = parser.parse_args()
    setup_logging(args.debug)
    return args
 
def main(cmdargs):
    """ Main """
    args = parse_args(cmdargs)
    setup_logging(True)
    anno = TRFAnno(in_vcf = args.input,
                    out_vcf = args.output,
                    executable = args.executable,
                    min_length = args.min_length
                    bestn = args.bestn,
                    layered = args.layered,
                    full = args.full,
                    trf_params = args.trf_params,
                    refanno = args.ref_bed,
                    ref = args.ref,
                    pctovl = args.pctovl)
    anno.run()

if __name__ == '__main__':
    main()
