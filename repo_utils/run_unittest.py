"""
One off unittests
"""
import os
import sys
import threading
import unittest

from collections import defaultdict, namedtuple
from unittest.mock import patch

import pandas as pd
from intervaltree import IntervalTree

# Use the current truvari, not any installed libraries
sys.path.insert(0, os.getcwd())
import truvari
from truvari import phab
from truvari import refine as refine_mod
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


FakePoolJob = namedtuple("FakePoolJob", "name")


def _job_that_may_die(job):
    """
    Return the job's name, except "dies", which exits like an OOM-killed worker
    """
    if job.name == "dies":
        os._exit(1)
    return job.name + "\n"


def _job_that_raises(job):
    """
    Return the job's name, except "raises", which raises like a failing aligner
    """
    if job.name == "raises":
        raise ValueError("aligner failed")
    return job.name + "\n"


class TestMonitoredPool(unittest.TestCase):
    def test_survives_dead_worker(self):
        """
        A worker dying mid-task must not stall the other jobs. Without monitoring this hangs,
        so the test fails by timing out
        """
        jobs = [FakePoolJob(n) for n in ("a", "dies", "b", "c")]
        collected = []
        failures = []
        finished = threading.Event()

        def consume():
            try:
                for result in phab.monitored_pool(_job_that_may_die, jobs, 2, failures=failures):
                    collected.append(result)
            finally:
                finished.set()

        # A 1 s window instead of 5 keeps the test quick
        with patch.object(phab, "MONITOR_FAIL_AFTER", 1):
            worker = threading.Thread(target=consume, daemon=True)
            worker.start()
            self.assertTrue(finished.wait(120),
                            "monitored_pool did not return after a worker died")
        self.assertEqual(sorted(collected), ["a\n", "b\n", "c\n"],
                         "every surviving job should still be yielded")
        self.assertEqual(failures, ["dies"],
                         "the lost job should be reported as failed")

    def test_reports_raising_job(self):
        """A job that raises is reported in failures, and the others still yield"""
        jobs = [FakePoolJob(n) for n in ("a", "raises", "b")]
        failures = []
        got = sorted(phab.monitored_pool(_job_that_raises, jobs, 2,
                                         failures=failures))
        self.assertEqual(got, ["a\n", "b\n"])
        self.assertEqual(failures, ["raises"])


class TestRefinedStratify(unittest.TestCase):
    """
    A region phab could not harmonize must keep its unrefined counts
    """

    @staticmethod
    def _regions():
        return pd.DataFrame({
            "chrom": ["chr1", "chr1", "chr2"],
            "start": [100, 500, 100],
            "end": [200, 600, 200],
            "in_tpbase": [1, 3, 5],
            "in_tp": [1, 3, 5],
            "in_fn": [2, 4, 6],
            "in_fp": [2, 4, 6],
            "refined": [True, True, False],
        })

    @staticmethod
    def _zero_counts(*_args, **_kwargs):
        """phab emitted nothing for either refined region"""
        return pd.DataFrame({"tpbase": [0, 0], "tp": [0, 0], "fn": [0, 0], "fp": [0, 0]})

    def test_failed_region_keeps_unrefined_counts(self):
        """The failed region's variants must stay in the comparison."""
        coords = [["chr1", 100, 200], ["chr1", 500, 600]]
        with patch("truvari.benchdir_count_entries", self._zero_counts):
            got = refine_mod.refined_stratify(
                "unused", coords, self._regions(), 1,
                failed_regions=["chr1:100-200"])

        failed = got.iloc[0]
        self.assertFalse(failed["refined"],
                         "a failed region is not a refined region")
        for col in ("tpbase", "tp", "fn", "fp"):
            self.assertEqual(failed[f"out_{col}"], failed[f"in_{col}"],
                             f"failed region lost its {col}")

        # The region phab did handle still takes phab's counts, zeros included.
        self.assertTrue(got.iloc[1]["refined"])
        self.assertEqual(got.iloc[1]["out_fn"], 0)

    def test_no_failures_is_unchanged(self):
        """Without failures the counts come from phab, as before."""
        coords = [["chr1", 100, 200], ["chr1", 500, 600]]
        with patch("truvari.benchdir_count_entries", self._zero_counts):
            got = refine_mod.refined_stratify(
                "unused", coords, self._regions(), 1)
        self.assertEqual(list(got["refined"]), [True, True, False])
        self.assertEqual(list(got["out_fn"]), [0, 0, 6])

    def test_merged_name_demotes_every_row_it_covers(self):
        """phab merges overlapping intervals before naming them, so one name can span rows."""
        regions = self._regions()
        mask = refine_mod.unrefinable_mask(regions, ["chr1:150-550"])
        self.assertEqual(list(mask), [True, True, False])

    def test_mask_is_chromosome_aware(self):
        """A failed name only demotes rows on its own chromosome"""
        regions = self._regions()
        mask = refine_mod.unrefinable_mask(regions, ["chr2:100-200"])
        self.assertEqual(list(mask), [False, False, True])


class FakeVCFI():
    def __init__(self):
        pass
    def get_haplotypes(self, _):
        return "A"

if __name__ == "__main__":
    unittest.main()
