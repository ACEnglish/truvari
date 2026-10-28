"""
One off unittests
"""
import os
import sys
import unittest

from collections import defaultdict
from intervaltree import IntervalTree
import pysam

# Use the current truvari, not any installed libraries
sys.path.insert(0, os.getcwd())
import truvari
import truvari.collapse
from truvari import phab
from truvari.region_vcf_iter import region_filter_stream, region_filter_fetch

# Assume we're running in truvari root directory


class TestBoundaries(unittest.TestCase):
    def test_boundary_issues(self):
        """
        Test if variants are correctly classified as within or outside a boundary.
        """
        vcf_fn = "repo_utils/test_files/variants/boundary.vcf.gz"
        region_start = 10
        region_end = 20

        vcf = truvari.VariantFile(vcf_fn)
        for entry in vcf:
            # Separate failure reporting per variant
            with self.subTest(entry=str(entry)):
                state = entry.info['include'] == 'in'
                self.assertEqual(
                    state, entry.within(region_start, region_end),
                    f"Bad Boundary {str(entry)}"
                )

    def test_region_filtering(self):
        """
        Test if regions are correctly filtered based on the provided BED file.
        """
        vcf_fn = "repo_utils/test_files/variants/boundary_cpx.vcf.gz"
        bed_fn = "repo_utils/test_files/beds/boundary_cpx.bed"

        tree = defaultdict(IntervalTree)
        with open(bed_fn, 'r') as fh:
            for line in fh:
                data = line.strip().split()
                tree[data[0]].addi(int(data[1]), int(data[2]) + 1)

        vcf = truvari.VariantFile(vcf_fn)

        # Using fetch_regions()
        for entry in vcf.fetch_regions(tree, True, False):
            with self.subTest(method="fetch_regions in", entry=str(entry)):
                self.assertEqual(
                    entry.info['include'], 'in', f"Bad in {str(entry)}")

        vcf.reset()
        for entry in vcf.fetch_regions(tree, False, False):
            with self.subTest(method="fetch_regions out", entry=str(entry)):
                self.assertEqual(
                    entry.info['include'], 'out', f"Bad out {str(entry)}")

        # Using region_filter_stream()
        vcf = truvari.VariantFile(vcf_fn)
        for entry in region_filter_stream(vcf, tree, True, False):
            with self.subTest(method="region_filter_stream in", entry=str(entry)):
                self.assertEqual(
                    entry.info['include'], 'in', f"Bad in {str(entry)}")

        vcf.reset()
        for entry in region_filter_stream(vcf, tree, False, False):
            with self.subTest(method="region_filter_stream out", entry=str(entry)):
                self.assertEqual(
                    entry.info['include'], 'out', f"Bad out {str(entry)}")

        # Using region_filter_fetch()
        vcf = truvari.VariantFile(vcf_fn)
        for entry in region_filter_fetch(vcf, tree, False):
            with self.subTest(method="region_filter_fetch in", entry=str(entry)):
                self.assertEqual(
                    entry.info['include'], 'in', f"Bad in {str(entry)}")


class TestFilteringLogic(unittest.TestCase):
    """
    Filtering logic
    """

    def test_filtering_logic(self):
        p = truvari.VariantParams(sizemin=0, sizefilt=0, passonly=True)
        vcf = truvari.VariantFile(
            "repo_utils/test_files/variants/filter.vcf", params=p)

        for entry in vcf:
            # Ensures all entries are tested separately
            with self.subTest(entry=str(entry)):
                try:
                    self.assertTrue(entry.filter_call(),
                                    f"Didn't filter {str(entry)}")
                except ValueError as e:
                    self.assertTrue(
                        e.args[0].startswith("Cannot compare multi-allelic"),
                        f"Unknown exception {str(entry)}"
                    )


class TestPhab(unittest.TestCase):
    """
    Manually running phab to ensure coverage. This doesn't actually check anything,
    just lets coverage.py see the code being executed. We rely on the phab/refine functional sub_tests 
    for accuracy
    """

    def test_phab(self):
        vcfs = [
            "repo_utils/test_files/variants/phab_base.vcf.gz",
            "repo_utils/test_files/variants/phab_comp2.vcf.gz"
        ]
        vcfhap = phab.VCFtoHaplotypes(
            "repo_utils/test_files/references/phab_ref.fa",
            vcfs, ["HG002", "HG005"]
        )

        reg_fn = truvari.make_temp_filename()
        with open(reg_fn, 'w') as fout:
            fout.write("chr1\t700\t900\n")

        regions = phab.parse_regions(reg_fn)
        vcfhap.set_regions(regions, buff=0)
        haps = vcfhap.build_all()
        key = 'chr1:700-900'
        result = phab.align_wrap(haps[key], phab.run_wfa)
        self.assertFalse(result.startswith("ERROR"), f"WFA failed: {result}")

        result = phab.align_wrap(haps[key], phab.run_poa, dedup=True)
        self.assertFalse(result.startswith("ERROR"), f"POA failed: {result}")

        result = phab.run_mafft(haps[key])
        self.assertTrue(isinstance(result, dict), f"MAFFT failed: {result}")

        with self.assertRaises(TypeError):
            phab.align_wrap("", phab.run_wfa)

        with self.assertRaises(ValueError):
            phab.get_align_method("bob")

        d = {'O':'A', 'ref_':'A'}
        self.assertTrue(phab.deduplicate_haps(d)[0] == {'ref_':'A'}, "Phab dedup failed")

        import pickle
        import multiprocessing.shared_memory as shm
        vcf_info = FakeVCFI()
        data = pickle.dumps(vcf_info)
        shared_info = shm.SharedMemory(create=True, size=len(data))
        shared_info.buf[:len(data)] = data
        mem_vcf_info = (shared_info.name, shared_info.size)
        job = phab.PhabJob("fake", mem_vcf_info)
        self.assertTrue(phab.status_marker(1, {}, "", job) == "ERROR: 'str' object is not callable", "PhabJob")
        phab.cleanup_shared_memory(shared_info)


TEST_VCF_HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr20,length=1000000>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
"""


def write_indexed_vcf(records, samples):
    """
    Write records (lists of the VCF columns, from CHROM to the samples' GT) under TEST_VCF_HEADER
    to a bgzipped, indexed temporary VCF and return its path
    """
    fn = truvari.make_temp_filename(suffix=".vcf")
    with open(fn, 'w') as fout:
        fout.write(TEST_VCF_HEADER)
        fout.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
                             + samples) + "\n")
        for record in records:
            fout.write("\t".join(str(column) for column in record) + "\n")
    return pysam.tabix_index(fn, preset="vcf", force=True)


def run_collapse(input_vcf, *options):
    """
    Run truvari collapse on input_vcf and return (kept record IDs, collapsed record IDs)
    """
    kept_fn = truvari.make_temp_filename(suffix=".vcf")
    collapsed_fn = truvari.make_temp_filename(suffix=".vcf")
    truvari.collapse.collapse_main(["-i", input_vcf, "-o", kept_fn, "-c", collapsed_fn, *options])
    with pysam.VariantFile(kept_fn) as kept, pysam.VariantFile(collapsed_fn) as collapsed:
        return [r.id for r in kept], [r.id for r in collapsed]


class TestCollapseBndSizeSubChunking(unittest.TestCase):
    """
    collapse sub-chunks chunks of 100+ variants by size. BNDs all have size 0, so they must be kept
    in one sub-chunk; otherwise each is alone and none are collapsed.
    """
    N_BNDS = 110

    def make_vcf(self, n_dels=0):
        records = [["chr20", 10000 + i % 50, f"bnd{i}", "N", "N[chr20:50000[", 60, "PASS", "SVTYPE=BND",
                    "GT", "0/1"] for i in range(self.N_BNDS)]
        records += [["chr20", 10010 + i, f"del{i}", "N", "<DEL>", 60, "PASS",
                     f"SVTYPE=DEL;SVLEN=-{100 * (i + 1)};END={10010 + i + 100 * (i + 1)}", "GT", "0/1"]
                    for i in range(n_dels)]
        return write_indexed_vcf(sorted(records, key=lambda r: r[1]), ["sample"])

    def test_bnds_kept_in_one_sub_chunk(self):
        params = truvari.VariantParams()
        entries = list(truvari.VariantFile(self.make_vcf(n_dels=5), params=params))
        sub_chunks = [chunk['base'] for chunk, _ in
                      truvari.collapse.tree_size_chunker(params, [({'base': entries, '__filtered': []}, 1)])]
        bnd_sub_chunks = [chunk for chunk in sub_chunks if any(entry.is_bnd() for entry in chunk)]
        self.assertEqual(len(bnd_sub_chunks), 1, "BNDs were split across sub-chunks")
        self.assertEqual(sum(entry.is_bnd() for entry in bnd_sub_chunks[0]), self.N_BNDS)
        self.assertEqual(sum(len(chunk) for chunk in sub_chunks), self.N_BNDS + 5, "Variants lost")

    def test_bnds_in_large_chunk_collapsed(self):
        kept, collapsed = run_collapse(self.make_vcf())
        self.assertEqual(len(kept), 1, f"Expected 1 kept BND, got {len(kept)}")
        self.assertEqual(len(collapsed), self.N_BNDS - 1)


class FakeVCFI():
    def __init__(self):
        pass
    def get_haplotypes(self, _):
        return "A"

if __name__ == "__main__":
    unittest.main()
