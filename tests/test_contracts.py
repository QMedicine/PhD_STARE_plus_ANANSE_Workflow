import importlib.util,csv,tempfile,unittest,os
from pathlib import Path
ROOT=Path(os.environ.get('PHDSA_TEST_SOURCE',Path(__file__).resolve().parents[1]))
def module(name,path):
 spec=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
validate=module('validate','bin/validate_inputs.py')
assemble=module('assemble','bin/assemble_atac.py')
plot=module('plot','tools/render_historical_plot.py')
prepare=module('prepare','bin/prepare_ananse_de.py')
contacts=module('contacts','bin/split_contacts.py')
report=module('report','bin/build_report.py')
class Contracts(unittest.TestCase):
 def test_both_directions_and_file_routes(self):
  r=validate.validate(ROOT/'templates/samplesheet_Old_and_Young_HSCs_v01.tsv',ROOT/'templates/source2target_Old_and_Young_HSCs_v01.tsv',check_files=False)
  self.assertEqual(r['transitions'],[dict(source='Old_HSC',target='Young_HSC'),dict(source='Young_HSC',target='Old_HSC')]);self.assertEqual(len(r['samples']),14)
 def test_explicit_reverse_does_not_duplicate_fits(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'contrast.tsv';p.write_text('source\ttarget\nOld_HSC\tYoung_HSC\nYoung_HSC\tOld_HSC\n')
   r=validate.validate(ROOT/'templates/samplesheet_Old_and_Young_HSCs_v01.tsv',p,check_files=False)
   self.assertEqual(len(r['transitions']),2);self.assertEqual(len(r['primary_transitions']),1)
   explicit=validate.validate(ROOT/'templates/samplesheet_Old_and_Young_HSCs_v01.tsv',p,reverse=False,check_files=False)
   self.assertEqual(len(explicit['primary_transitions']),2)
 def test_unknown_state_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'contrast.tsv';p.write_text('source\ttarget\nOld_HSC\tYoungHSC\n')
   with self.assertRaises(ValueError):validate.validate(ROOT/'templates/samplesheet_Old_and_Young_HSCs_v01.tsv',p,check_files=False)
 def test_generic_selected_plot_chunk(self):
  code=plot.render(ROOT/'bin/plot_branching_template.R','inf.tsv','diff.tsv','de.tsv','TESTTF','plots','Resting_to_Stimulated')
  self.assertIn('deseq_file <- "de.tsv"',code);self.assertIn('max_steps <- 2',code);self.assertIn('device = cairo_pdf',code)
  self.assertIn('min_edge_weight <- 0.01',code);self.assertIn('source_label <- "Resting"',code)
  for label in ['HSC','Old','Young','JASPAR2024']:self.assertNotIn(label,code)
 def test_generic_cell_states(self):
  with tempfile.TemporaryDirectory() as tmp:
   path=Path(tmp)/'samples.tsv'
   text=(ROOT/'templates/samplesheet_Old_and_Young_HSCs_v01.tsv').read_text().replace('Old_HSC','Resting_Macrophage').replace('Young_HSC','Stimulated_Macrophage').replace('LT_HSC','Macrophage').replace('Old','Resting').replace('Young','Stimulated')
   path.write_text(text)
   contrast=Path(tmp)/'contrast.tsv';contrast.write_text('source\ttarget\nResting_Macrophage\tStimulated_Macrophage\n')
   result=validate.validate(path,contrast,check_files=False)
   self.assertEqual(result['transitions'][1],dict(source='Stimulated_Macrophage',target='Resting_Macrophage'))
 def test_transition_between_cell_types(self):
  with tempfile.TemporaryDirectory() as tmp:
   path=Path(tmp)/'samples.tsv'
   with (ROOT/'templates/samplesheet_Old_and_Young_HSCs_v01.tsv').open() as f:
    reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
   states=list(dict.fromkeys(r['cell_state'] for r in rows))
   for r in rows:r['cell_type']='Macrophage' if r['cell_state']==states[0] else 'Fibroblast'
   with path.open('w') as f:
    writer=csv.DictWriter(f,fields,delimiter='\t');writer.writeheader();writer.writerows(rows)
   result=validate.validate(path,ROOT/'templates/source2target_Old_and_Young_HSCs_v01.tsv',check_files=False)
   self.assertEqual(len(result['transitions']),2)
 def test_unknown_raw_sample_rejected(self):
  with self.assertRaisesRegex(ValueError,'raw_sample_names'):
   validate.validate(ROOT/'templates/samplesheet_RNA_fastq_pilot_v01.tsv',ROOT/'templates/source2target_Old_and_Young_HSCs_v01.tsv',check_files=False,stage='raw_prepare',raw_samples=['unknown'])
 def test_nonmouse_contacts_and_mixed_aliases(self):
  import gzip
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);source=root/'contacts.tsv';source.write_text('2L\t0\t50\t2L\t50\t100\t7\n3R\t0\t50\t3R\t50\t100\t9\n2L\t0\t50\t3R\t50\t100\t11\n')
   contacts.split_contacts(source,root/'split')
   self.assertEqual(sorted(p.name for p in (root/'split').glob('*.gz')),['chr2L.txt.gz','chr3R.txt.gz'])
   with gzip.open(root/'split/chr2L.txt.gz','rt') as f:self.assertEqual(f.read(),'0\t50\t7\n')
   source.write_text(source.read_text()+'chr2L\t0\t50\tchr2L\t50\t100\t1\n')
   with self.assertRaisesRegex(ValueError,'Mixed chromosome aliases'):contacts.split_contacts(source,root/'mixed')
 def test_report_uses_resolved_parameters_and_checkpoint_scope(self):
  context=dict(stage='atac_historical_hic_checkpoint',union_windows='/data/authenticated.bed',chosen_res=50000,summit_flank=100,stare_window=2000000,enhancer_window=3000000,gabc_threshold=.03,min_count=70,min_count_tf=80,min_tpm=.2,min_tpm_tf=.3)
  methods,expression=report.methods_text(context)
  self.assertIn('preserved named summit-window BED',methods);self.assertIn('50000-bp bins',methods);self.assertIn('2000000-bp scoring',methods)
  self.assertIn('≥70',expression);self.assertIn('≥80',expression);self.assertNotIn('100-kb',methods)
 def test_ananse_compatibility_preserves_all_de_values(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);previous=Path.cwd();os.chdir(root)
   try:
    raw=b'gene_id\tlog2FoldChange\tpadj\nGENE_A\t1.123456789012345\t0.002\nGENE_B\t-0.0\tNA\n'
    (root/'de.tsv').write_bytes(raw);prepare.prepare(root/'de.tsv',root/'compat.tsv')
    self.assertEqual((root/'de.tsv').read_bytes(),raw)
    self.assertEqual((root/'compat.tsv').read_bytes().split(b'\n',1)[1],raw.split(b'\n',1)[1])
    self.assertTrue(__import__('json').loads((root/'de_compatibility.json').read_text())['data_bytes_unchanged'])
   finally:os.chdir(previous)
 def test_header_and_data_order_modes(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);samples=root/'samples.tsv'
   names=['O1','O2','Y1','Y2'];values=[27,23,24,26]
   with samples.open('w') as f:
    w=csv.DictWriter(f,['sample_name','assay','condition','replicate'],delimiter='\t');w.writeheader()
    for n in names:w.writerow(dict(sample_name=n,assay='ATAC',condition='Control' if n[0]=='O' else 'Treatment',replicate=n[1]))
   files=[]
   for n,v in zip(names,values):
    p=root/(n+'.window_counts.tsv');p.write_text(f'1\t100\t601\tpeak_1\t{v}\n');files.append(str(p))
   import os
   previous=Path.cwd();os.chdir(root)
   try:
    assemble.assemble(samples,files,'historical_column_order',root/'historic.tsv',['Y1','Y2','O1','O2'],['O1','O2','Y1','Y2'])
    self.assertEqual((root/'historic.tsv').read_text().splitlines(),['peak_id\tY1\tY2\tO1\tO2','peak_1\t27\t23\t24\t26'])
    with self.assertRaises(ValueError):assemble.assemble(samples,files,'historical_column_order',root/'invalid.tsv')
    assemble.assemble(samples,files,'corrected_column_order',root/'corrected.tsv')
    self.assertEqual((root/'corrected.tsv').read_text().splitlines(),['peak_id\tO1\tO2\tY1\tY2','peak_1\t27\t23\t24\t26'])
   finally:os.chdir(previous)
if __name__=='__main__':unittest.main()
