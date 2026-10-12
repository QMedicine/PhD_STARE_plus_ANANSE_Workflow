"""Duplicate cleanup must preserve distinct loci and sum matching locus annotations."""
import csv,importlib.util,os,tempfile,unittest
from pathlib import Path
ROOT=Path(os.environ.get('PHDSA_TEST_SOURCE',Path(__file__).resolve().parents[1]))
try:
 import pandas
except ImportError:pandas=None

spec=importlib.util.spec_from_file_location('merge',ROOT/'bin/merge_recount.py');merge=importlib.util.module_from_spec(spec)
if pandas is not None:spec.loader.exec_module(merge)
@unittest.skipIf(pandas is None,'Count-cleanup test requires the container pandas runtime')
class Cleanup(unittest.TestCase):
 def test_same_locus_sum_and_distinct_locus_representative(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp)
   (root/'selected.csv').write_text('Geneid,control,treatment\ng2,30,50\ng4,7,11\n')
   annotations=[('g1','1','+','100','Alias',10,20),('g2','1','+','200','Alias',20,30),('g3','2','-','100','Separate',99,88),('g4','3','-','200','Separate',7,11)]
   files=[]
   for sample,column in [('control',5),('treatment',6)]:
    path=root/(sample+'.featureCounts.txt');files.append(path)
    with path.open('w') as f:
     w=csv.writer(f,delimiter='\t');w.writerow(['Geneid','Chr','Start','End','Strand','Length','gene_name','gene_biotype',sample])
     for r in annotations:w.writerow([r[0],r[1],0,int(r[3]),r[2],r[3],r[4],'protein_coding',r[column]])
   previous=Path.cwd();os.chdir(root)
   try:
    merge.merge(root/'selected.csv',files,root/'recount.csv')
    self.assertEqual((root/'recount.csv').read_text(),(root/'selected.csv').read_text())
    self.assertIn('sum_same_chromosome_and_strand',(root/'cleanup_groups.tsv').read_text())
    self.assertIn('retain_checkpoint_representative',(root/'cleanup_groups.tsv').read_text())
   finally:os.chdir(previous)
if __name__=='__main__':unittest.main()
