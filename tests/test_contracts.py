import importlib.util,csv,tempfile,unittest,os
from pathlib import Path
ROOT=Path(os.environ.get('PHDSA_TEST_SOURCE',Path(__file__).resolve().parents[1]))
def module(name,path):
 spec=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
validate=module('validate','bin/validate_inputs.py')
assemble=module('assemble','bin/assemble_atac.py')
plot=module('plot','tools/render_historical_plot.py')
class Contracts(unittest.TestCase):
 def test_both_directions_and_file_routes(self):
  r=validate.validate(ROOT/'templates/samplesheet_Old_and_Young_HSCs_v01.tsv',ROOT/'templates/source2target_Old_and_Young_HSCs_v01.tsv',check_files=False)
  self.assertEqual(r['transitions'],[dict(source='Old_HSC',target='Young_HSC'),dict(source='Young_HSC',target='Old_HSC')]);self.assertEqual(len(r['samples']),14)
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
