"""
Annotator that can consume the 1kg phase3 integrated SVs
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz
"""
import re
import sys
import logging

from collections import defaultdict
from intervaltree import IntervalTree

import pysam
import truvari

sv_alt_match = re.compile("\<(?P<SVTYPE>.*)\>")


class OneKg():

    def __init__(self, anno_file, chrom, start, end):
        self.anno_file = anno_file
        self.tree = IntervalTree()
        try:
            for entry in self.anno_file.fetch(chrom, start, end):
                self.tree.addi(entry.start, entry.stop, entry)
        except Exception:
            pass
        self.tree_bts = list(self.tree.boundary_table.keys())
        # Everything past the end would need to be pumped through
        self.tree_bts.append(sys.maxsize)
        self.n_header = None

    def load_header(self, in_vcf):
        """
        Returns the header of the information we'll be adding to the annotated vcfs
        """

        ret = in_vcf.header.copy()
        ret.add_line(('##INFO=<ID=OKG_MSTART,Number=1,Type=Integer,Description="Mitochondrial '
                      'start coordinate of inserted sequence">'))
        ret.add_line(('##INFO=<ID=OKG_MLEN,Number=1,Type=Integer,Description="Estimated length '
                      'of mitochondrial insert">'))
        ret.add_line(('##INFO=<ID=OKG_MEND,Number=1,Type=Integer,Description="Mitochondrial end'
                      ' coordinate of inserted sequence">'))
        ret.add_line(('##INFO=<ID=OKG_MEINFO,Number=4,Type=String,Description="Mobile element '
                      'info of the form NAME,START,END<POLARITY; If there is only 5\' OR 3\' '
                      'support for this call, will be NULL NULL for START and END">'))
        ret.add_line(('##INFO=<ID=OKG_AF,Number=.,Type=Float,Description="Estimated allele '
                      'frequency in the range (0,1)">'))
        ret.add_line(('##INFO=<ID=OKG_EAS_AF,Number=.,Type=Float,Description="Allele frequency '
                      'in the EAS populations calculated from AC and AN, in the range (0,1)">'))
        ret.add_line(('##INFO=<ID=OKG_EUR_AF,Number=.,Type=Float,Description="Allele frequency '
                      'in the EUR populations calculated from AC and AN, in the range (0,1)">'))
        ret.add_line(('##INFO=<ID=OKG_AFR_AF,Number=.,Type=Float,Description="Allele frequency '
                      'in the AFR populations calculated from AC and AN, in the range (0,1)">'))
        ret.add_line(('##INFO=<ID=OKG_AMR_AF,Number=.,Type=Float,Description="Allele frequency '
                      'in the AMR populations calculated from AC and AN, in the range (0,1)">'))
        ret.add_line(('##INFO=<ID=OKG_SAS_AF,Number=.,Type=Float,Description="Allele frequency '
                      'in the SAS populations calculated from AC and AN, in the range (0,1)">'))
        ret.add_line(('##INFO=<ID=OKG_SVTYPE,Number=1,Type=String,Description="OneThousandGenome'
                      'ALT Type">'))
        self.n_header = ret

    def annotate(self, entry, refdist=500, size_min=50, size_max=50000):
        """
        Given an pyvcf Variant entry do the matching
        """
        # Biggest shortcut, only annotate SVs
        if "SVLEN" not in entry.info:
            return entry
        
        if entry.stop + refdist < self.tree_bts[0]:
            return entry

        if entry.start - refdist > self.tree_bts[0]:
            self.tree_bts.pop(0)
            
        if not (size_min <= abs(entry.info["SVLEN"]) <= size_max):
            return entry

        # Don't lookup until we have to
        m_type = None
        candidates = []
        for anno_entry in self.tree.overlap(entry.start - refdist, entry.stop + refdist):
            anno_entry = anno_entry.data
            a_size = truvari.get_vcf_entry_size(anno_entry)
            if not (size_min <= a_size <= size_max):
                continue

            ps, sd = truvari.entry_size_similarity(entry, anno_entry)
            if not ps >= 0.7:
                continue

            mat1 = sv_alt_match.match(anno_entry.alts[0])
            if mat1 is not None:
                a_type = mat1.groupdict()["SVTYPE"]
            else:
                a_type = truvari.get_vcf_variant_type(anno_entry)
            
            # Don't make until we have to, and only do so once
            if m_type is None:
                m_type = truvari.get_vcf_variant_type(entry)

            if not (a_type == m_type
                    or ((a_type == "CN0" or a_type.startswith("DEL")) and m_type == "DEL")
                    or (m_type == "INS" and a_type.startswith("INS"))):
                continue

            # RO doesn't work for INS?
            ro = truvari.entry_reciprocal_overlap(entry, anno_entry)
            if m_type != "INS" and ro < 0.5:
                continue

            candidates.append((ro, ps, anno_entry))

        if candidates:
            truvari.match_sorter(candidates)
            return self.add_info(truvari.copy_entry(entry, self.n_header), candidates[0][-1])
        return entry

    def extract_info(self, annot):
        """MSTART MLEN MEND MEINFO AF EAS_AF EUR_AF AFR_AF AMR_AF SAS_AF ALT"""
        def infoc(key):
            if key in annot.info:
                return key, annot.info[key]
            return None, None

        def altp():
            """reformat the alt seq"""
            ret = []
            for i in annot.alts:
                if i.startswith("<"):
                    ret.append(i[1:-1])
            return "SVTYPE", tuple(ret)

        return [infoc("MSTART"),
                infoc("MLEN"),
                infoc("MEND"),
                infoc("MEINFO"),
                infoc("AF"),
                infoc("EAS_AF"),
                infoc("EUR_AF"),
                infoc("AFR_AF"),
                infoc("AMR_AF"),
                infoc("SAS_AF")]
    def add_info(self, entry, annot):
        """
        Put the relevant info fields into the entry to be annotated
        """
        # Get the annotations out of the annot and add them to the entry
        if not annot:
            return entry
        i = self.extract_info(annot)
        for key, val in self.extract_info(annot):
            if val is not None:
                entry.info["OKG_" + key] = val
        return entry


def parse_args(args):
    setup_logging()
    pass

def main():
    ## Next -- need to figure out how to make an anno per-chromosome
    ## Then figure out how to multi-process it
    ## Also want to accept a --region chrom[:start-end]+ parameter so that I can subset
    ## and, all of that needs to be generalized such that I can reuse it as I expand into other annotations
    #anno = OneKg(pysam.VariantFile("/home/english/science/english/annotation/1kg_phase3/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz"))
    in_vcf = pysam.VariantFile(sys.argv[1])
    logging.info("Loading annotation")
    anno = {}
    n_header = None
    for chrom in in_vcf.header.contigs:
        anno[chrom] = OneKg(pysam.VariantFile(sys.argv[2]), chrom, 0, in_vcf.header.contigs[chrom].length)
        anno[chrom].load_header(in_vcf)
        n_header = anno[chrom].n_header
    out_vcf = pysam.VariantFile("/dev/stdout", 'w', header=n_header)
    logging.info("Annotating Variants")
    for entry in in_vcf:
        out_vcf.write(anno[entry.chrom].annotate(entry))
    out_vcf.close()
    logging.info("Finished")   

if __name__ == '__main__':
    main()
