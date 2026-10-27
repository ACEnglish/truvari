"""Check that streaming excluded rows leaves collapse groups and calls unchanged.

Run with the TruVari environment: python repo_utils/test_stream_filtered.py
Only synthetic fixtures in a temporary directory are created.
"""
from collections import Counter
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pysam
import truvari

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##contig=<ID=chr2,length=1000000>
##ALT=<ID=DEL,Description="Deletion">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">
##INFO=<ID=MATEID,Number=.,Type=String,Description="Mate">
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
"""

def row(pos, ident, size=100, chrom="chr1", alt="<DEL>", filt="PASS",
        cells=("0|1:4,7", ".|.:.,.")):
    info = "SVTYPE=BND" if "[" in alt or "]" in alt else f"END={pos+size};SVTYPE=DEL;SVLEN=-{size}"
    return "\t".join(map(str, (chrom, pos, ident, "N", alt, 60, filt,
        info, "GT:AD", *cells)))+"\n"

def write_vcf(directory, rows):
    path=Path(directory)/"input.vcf"
    path.write_text(HEADER+"".join(sorted(rows,key=lambda x:(x.split("\t")[0],int(x.split("\t")[1])))))
    zipped=str(path)+".gz"
    pysam.tabix_compress(str(path),zipped,force=True)
    pysam.tabix_index(zipped,preset="vcf",force=True)
    return zipped

def parameters(**kwargs):
    return truvari.VariantParams(chunksize=500,refdist=500,sizefilt=50,
        sizemin=50,sizemax=50000,short_circuit=True,skip_gt=True,**kwargs)

def grouping(path, params, callback=None):
    with truvari.VariantFile(path,params=params) as source:
        return [({k:[x.id for x in v] for k,v in chunk.items()},number)
                for chunk,number in truvari.chunker(params,("base",iter(source)),
                    filtered_callback=callback)]

def canonical(path):
    with pysam.VariantFile(str(path)) as source:
        records=[]
        for rec in source:
            fields=str(rec).rstrip().split("\t")
            fields[7]=";".join(sorted(fields[7].split(";")))
            tags=fields[8].split(":")
            order=sorted(range(len(tags)),key=lambda i:(tags[i]!="GT",tags[i]))
            fields[8]=":".join(tags[i] for i in order)
            for column in range(9,len(fields)):
                values=fields[column].split(":")
                values.extend(["."]*(len(tags)-len(values)))
                fields[column]=":".join(values[i] for i in order)
            records.append("\t".join(fields))
        return tuple(sorted(records))

class StreamingTests(unittest.TestCase):
    def check_groups(self, rows, **kwargs):
        with tempfile.TemporaryDirectory() as tmp:
            source=write_vcf(tmp,rows)
            baseline=grouping(source,parameters(**kwargs))
            streamed=[]
            actual=grouping(source,parameters(**kwargs),lambda x:streamed.append(x.id))
            expected=[({k:v for k,v in chunk.items() if k!="__filtered"},number)
                      for chunk,number in baseline]
            actual=[({k:v for k,v in chunk.items() if k!="__filtered"},number)
                    for chunk,number in actual]
            self.assertEqual(expected,actual)
            self.assertEqual(Counter(streamed),Counter(x for chunk,_ in baseline
                for x in chunk.get("__filtered",[])))
            return baseline,streamed

    def test_wide_overlap_group_preserved(self):
        baseline,streamed=self.check_groups([row(1000,"left",10000),
            row(2000,"large",100000),row(6000,"right",10000)])
        self.assertEqual(baseline[0][0]["base"],["left","right"])
        self.assertEqual(streamed,["large"])

    def test_small_and_large_bypass(self):
        _,streamed=self.check_groups([row(10,"small",10),row(50,"normal"),
            row(100,"large",100000)])
        self.assertEqual(streamed,["small","large"])

    def test_filter_policy_preserved(self):
        _,streamed=self.check_groups([row(10,"failed",filt="LowQual"),
            row(50,"normal")],passonly=True)
        self.assertEqual(streamed,["failed"])

    def test_bnd_disabled(self):
        _,streamed=self.check_groups([row(10,"bnd",alt="N[chr2:5000["),
            row(50,"normal")],bnddist=-1)
        self.assertEqual(streamed,["bnd"])

    def test_all_excluded(self):
        baseline,streamed=self.check_groups([row(10,"small",10),row(50,"large",100000)])
        self.assertEqual(len(baseline),1)
        self.assertEqual(streamed,["small","large"])

    def test_empty(self):
        self.check_groups([])

    def test_multiple_chromosomes_and_gaps(self):
        self.check_groups([row(10,"small",10),row(1000,"a"),row(5000,"b"),
            row(10,"c",chrom="chr2"),row(100,"large",100000,chrom="chr2")])

    def test_callback_exception_propagates(self):
        with tempfile.TemporaryDirectory() as tmp:
            source=write_vcf(tmp,[row(10,"small",10)])
            def fail(_):
                raise RuntimeError("write failed")
            with self.assertRaisesRegex(RuntimeError,"write failed"):
                grouping(source,parameters(),fail)

    def check_cli(self, rows):
        with tempfile.TemporaryDirectory() as tmp:
            source=write_vcf(tmp,rows)
            environment=os.environ.copy()
            environment.update(PYTHONHASHSEED="0",PYTHONDONTWRITEBYTECODE="1",
                OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1")
            outputs=[]
            for stream in (False,True):
                prefix=Path(tmp)/("stream" if stream else "baseline")
                kept=Path(str(prefix)+"_kept.vcf")
                members=Path(str(prefix)+"_members.vcf")
                command=[sys.executable,"-B","-m","truvari","collapse",
                    "--decompose","-i",source,"-o",str(kept),"-c",str(members)]
                if stream:
                    command.append("--stream-filtered")
                result=subprocess.run(command,cwd=Path(__file__).resolve().parents[1],
                    env=environment,capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stderr)
                outputs.append((canonical(kept),canonical(members)))
            self.assertEqual(outputs[0],outputs[1])
            self.assertEqual(len(outputs[0][0])+len(outputs[0][1]),len(rows))

    def test_cli_records_genotypes_and_clusters(self):
        rows=[row(10,"small",10),row(1000,"a"),row(1001,"b"),
              row(2000,"large",100000,cells=("0/0:20,0","./.:.,.")),
              row(3000,"missing",cells=(".|1:.,4","0:9,0")),
              row(9000,"failed",filt="LowQual")]
        self.check_cli(rows)

    def test_cli_all_excluded(self):
        self.check_cli([row(10,"small",10),row(50,"large",100000)])

    def test_cli_empty(self):
        self.check_cli([])

if __name__=="__main__":
    unittest.main(verbosity=2)
