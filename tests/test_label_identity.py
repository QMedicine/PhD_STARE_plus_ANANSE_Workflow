"""Adversarial identity checks: shuffled producers and ambiguous labels."""
import csv,importlib.util,os,tempfile,unittest
from pathlib import Path
ROOT=Path(os.environ.get('PHDSA_TEST_SOURCE',Path(__file__).resolve().parents[1]))
def load(name):
    spec=importlib.util.spec_from_file_location(name,ROOT/'bin'/f'{name}.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
assemble=load('assemble_atac');validate=load('validate_inputs');de=load('prepare_ananse_de')
try:
    import pandas as pd
except ImportError:pd=None

class LabelIdentity(unittest.TestCase):
    def test_ambiguous_transition_requires_explicit_labels(self):
        spec=importlib.util.spec_from_file_location('render',ROOT/'tools/render_historical_plot.py')
        render=importlib.util.module_from_spec(spec);spec.loader.exec_module(render)
        arguments=[ROOT/'bin/plot_branching_template.R','inf.tsv','diff.tsv','de.tsv','TF','plots','A_to_B_to_C']
        with self.assertRaisesRegex(ValueError,'Ambiguous transition'):render.render(*arguments)
        code=render.render(*arguments,source='A_to_B',target='C')
        self.assertIn('source_label <- "A_to_B"',code);self.assertIn('target_label <- "C"',code)
    def test_atac_shuffled_files_preserve_labels_and_reject_ambiguity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);sheet=root/'samples.tsv';sheet.write_text('sample_name\tassay\nZ\tATAC\nA\tATAC\n')
            a=root/'A.window_counts.tsv';z=root/'Z.window_counts.tsv'
            a.write_text('1\t0\t501\tpeak_1\t17\n1\t502\t1003\tpeak_2\t29\n')
            z.write_text('1\t0\t501\tpeak_1\t113\n1\t502\t1003\tpeak_2\t227\n')
            previous=Path.cwd();os.chdir(root)
            try:
                for files in [[a,z],[z,a]]:
                    assemble.assemble(sheet,files,'corrected_column_order',root/'out.tsv')
                    self.assertEqual((root/'out.tsv').read_text(),'peak_id\tZ\tA\npeak_1\t113\t17\npeak_2\t227\t29\n')
                with self.assertRaisesRegex(ValueError,'Duplicate ATAC'):assemble.assemble(sheet,[a,z,a],'corrected_column_order',root/'out.tsv')
                a.write_text('1\t502\t1003\tpeak_2\t29\n1\t0\t501\tpeak_1\t17\n')
                with self.assertRaisesRegex(ValueError,'identities/order'):assemble.assemble(sheet,[z,a],'corrected_column_order',root/'out.tsv')
                a.write_text('1\t0\t501\tpeak_1\t17\n1\t502\t1003\tpeak_1\t29\n')
                with self.assertRaisesRegex(ValueError,'Duplicate peak'):assemble.assemble(sheet,[a,z],'corrected_column_order',root/'out.tsv')
            finally:os.chdir(previous)

    def test_manifest_reused_bam_and_duplicate_headings_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);original=ROOT/'templates/samplesheet_Old_and_Young_HSCs_v01.tsv'
            with original.open() as f:
                reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
            rna=[r for r in rows if r['assay']=='RNA'];rna[1]['bam']=rna[0]['bam']
            sheet=root/'samples.tsv'
            with sheet.open('w') as f:
                writer=csv.DictWriter(f,fields,delimiter='\t');writer.writeheader();writer.writerows(rows)
            with self.assertRaisesRegex(ValueError,'Input file reused'):validate.validate(sheet,ROOT/'templates/source2target_Old_and_Young_HSCs_v01.tsv',check_files=False)
            sheet.write_text('source\ttarget\tsource\nA\tB\tC\n')
            with self.assertRaisesRegex(ValueError,'duplicate column'):validate.read_tsv(sheet,['source','target'])

    def test_de_rejects_row_number_column_and_duplicate_genes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'in.tsv'
            for contents,pattern in [('row_number\tgene_id\tlog2FoldChange\tpadj\n1\tA\t2\t.01\n','first column'),('gene_id\tlog2FoldChange\tpadj\nA\t2\t.01\nA\t3\t.02\n','Duplicate')]:
                source.write_text(contents)
                with self.assertRaisesRegex(ValueError,pattern):de.prepare(source,root/'out.tsv')

    @unittest.skipIf(pd is None,'Requires the native pandas container')
    def test_rna_count_files_shuffled_and_duplicate_or_extra_columns_rejected(self):
        merge=load('merge_recount')
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);selection=root/'selection.csv';selection.write_text('Geneid,Z,A\ng1,113,17\ng2,227,29\n')
            files=[]
            header='Geneid\tChr\tStart\tEnd\tStrand\tLength\tgene_name\tgene_biotype\t'
            for name,values in [('A',[17,29]),('Z',[113,227])]:
                file=root/f'{name}.featureCounts.txt';files.append(file)
                file.write_text(header+name+'\ng1\t1\t1\t100\t+\t100\tOne\tprotein_coding\t'+str(values[0])+'\ng2\t1\t101\t200\t+\t100\tTwo\tprotein_coding\t'+str(values[1])+'\n')
            previous=Path.cwd();os.chdir(root)
            try:
                for order in [files,list(reversed(files))]:
                    merge.merge(selection,order,root/'out.csv')
                    self.assertEqual((root/'out.csv').read_bytes(),selection.read_bytes())
                with self.assertRaisesRegex(ValueError,'Duplicate RNA'):merge.merge(selection,files+[files[0]],root/'out.csv')
                files[0].write_text(files[0].read_text().replace('\tA\n','\tA\tB\n').replace('\t17\n','\t17\t999\n').replace('\t29\n','\t29\t999\n'))
                with self.assertRaisesRegex(ValueError,'exactly one BAM'):merge.merge(selection,files,root/'out.csv')
            finally:os.chdir(previous)

    @unittest.skipIf(pd is None,'Requires the native pandas container')
    def test_tpm_conflicting_versioned_aliases_rejected(self):
        tpm=load('tpm_filter')
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'map.tsv').write_text('gene_id\tgene_name\ng1.1\tOne\ng1.2\tTwo\n')
            (root/'counts.csv').write_text('Geneid,A\ng1,100\n');(root/'lengths.tsv').write_text('gene\tmerged\ng1\t1000\n')
            with self.assertRaisesRegex(ValueError,'Conflicting gene symbols'):tpm.calculate(root/'counts.csv',root/'lengths.tsv',root/'map.tsv')

if __name__=='__main__':unittest.main()
