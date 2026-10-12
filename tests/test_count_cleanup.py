"""Duplicate cleanup must preserve distinct loci and sum matching locus annotations."""
import csv,importlib.util,os,tempfile,unittest
from pathlib import Path
ROOT=Path(os.environ.get('PHDSA_TEST_SOURCE',Path(__file__).resolve().parents[1]))
try:
 import pandas
except ImportError:pandas=None

spec=importlib.util.spec_from_file_location('merge',ROOT/'bin/merge_recount.py');merge=importlib.util.module_from_spec(spec)
if pandas is not None:spec.loader.exec_module(merge)
spec_annotation=importlib.util.spec_from_file_location('annotation',ROOT/'bin/annotation_counts.py');annotation=importlib.util.module_from_spec(spec_annotation)
if pandas is not None:spec_annotation.loader.exec_module(annotation)
@unittest.skipIf(pandas is None,'Count-cleanup test requires the container pandas runtime')
class Cleanup(unittest.TestCase):
 def test_annotation_mode_conserves_counts_and_records_exclusions(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);sheet=root/'samples.tsv';sheet.write_text('sample_name\tassay\ncontrol\tRNA\ntreatment\tRNA\n')
   annotations=[('g1','1','Alpha',100,10,20),('g2','2','ALPHA',200,20,30),('g3','3','Beta',300,7,11),('g4','4','ENSMUSG12345',400,100,200)]
   files=[]
   for sample,column in [('control',4),('treatment',5)]:
    p=root/(sample+'.featureCounts.txt');files.append(p)
    with p.open('w') as f:
     w=csv.writer(f,delimiter='\t');w.writerow(['Geneid','Chr','Start','End','Strand','Length','gene_name','gene_biotype',sample])
     for r in annotations:w.writerow([r[0],r[1],0,r[3],'+',r[3],r[2],'protein_coding',r[column]])
   previous=Path.cwd();os.chdir(root)
   try:
    annotation.build(sheet,files)
    observed=pandas.read_csv('annotation_counts.csv',index_col=0)
    self.assertEqual(observed.to_dict(),{'control':{'g1':30,'g3':7},'treatment':{'g1':50,'g3':11}})
    provenance=__import__('json').loads(Path('annotation_counts_provenance.json').read_text())
    self.assertEqual(provenance['annotated_count_totals'],{'control':37,'treatment':61})
    self.assertEqual(provenance['discarded_unmapped_name_count_totals'],{'control':100,'treatment':200})
    self.assertEqual(provenance['historical_equivalence'],'NOT_CLAIMED')
   finally:os.chdir(previous)
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
