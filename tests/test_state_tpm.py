import csv
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

ROOT=Path(os.environ.get('PHDSA_TEST_SOURCE',Path(__file__).resolve().parents[1]))
sys.path.insert(0,str(ROOT/'bin'))
from select_state_tpm import state_vector, select
from prepare_tf_expression import prepare

class StateTPMTests(unittest.TestCase):
    def fixture(self, root):
        samples=root/'samples.tsv'
        # Deliberately interleave states and use names unrelated to their states.
        samples.write_text('sample_name\tassay\tcell_state\nq4\tRNA\tResting\nq2\tRNA\tActivated\nq1\tRNA\tResting\nq3\tRNA\tActivated\na1\tATAC\tResting\n')
        tpm=root/'tpm.tsv'
        tpm.write_text('gene_id\tq3\tq1\tq2\tq4\nTF_A\t9\t1\t11\t3\nTF_B\t2\t7\t4\t9\n')
        return tpm,samples

    def test_state_selection_survives_shuffled_columns_and_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);tpm,samples=self.fixture(root)
            self.assertEqual(state_vector(tpm,samples,'Resting')[3],{'TF_A':2,'TF_B':8})
            self.assertEqual(state_vector(tpm,samples,'Activated')[3],{'TF_A':10,'TF_B':3})
            record=select(tpm,samples,'Resting',root/'selected.tsv',root/'selection.json')
            self.assertEqual(record['selected_rna_samples'],['q4','q1'])
            self.assertEqual((root/'selected.tsv').read_text().splitlines()[0],'gene_id\tq4\tq1')
            self.assertEqual(state_vector(tpm,samples,'Resting','historical_pooled')[3]['TF_A'],6)

    def test_missing_state_sample_and_duplicate_identities_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);tpm,samples=self.fixture(root)
            with self.assertRaisesRegex(ValueError,'No RNA'):state_vector(tpm,samples,'Unknown')
            text=tpm.read_text();tpm.write_text(text.replace('q4','other'))
            with self.assertRaisesRegex(ValueError,'missing from TPM'):state_vector(tpm,samples,'Resting')
            tpm.write_text(text+'TF_A\t1\t1\t1\t1\n')
            with self.assertRaisesRegex(ValueError,'duplicate TPM'):state_vector(tpm,samples,'Resting')
            tpm.write_text(text.replace('TF_A\t9','TF_A\tnan'))
            with self.assertRaisesRegex(ValueError,'Non-finite'):state_vector(tpm,samples,'Activated')

    def test_heatmap_table_is_ranked_and_directional_without_age_labels(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);tpm,samples=self.fixture(root)
            ranks=root/'influence.tsv';de=root/'de.tsv'
            factors=['TF_'+str(i) for i in range(20)]
            ranks.write_text('factor\tinfluence_score\ttarget_score_scaled\tG_score_scaled\n'+''.join(f'{tf}\t1\t.8\t.2\n' for tf in factors))
            tpm.write_text('gene_id\tq3\tq1\tq2\tq4\n'+''.join(f'{tf}\t9\t1\t11\t3\n' for tf in factors))
            de.write_text('gene_id\tlog2FoldChange\tpadj\n'+''.join(f'{tf}\t2\t.001\n' for tf in factors))
            forward=prepare(ranks,de,tpm,samples,'Resting','Activated',root/'forward.tsv',root/'forward.json')
            de.write_text(de.read_text().replace('\t2\t.001','\t-2\t.001'))
            reverse=prepare(ranks,de,tpm,samples,'Activated','Resting',root/'reverse.tsv',root/'reverse.json')
            self.assertEqual([r['factor'] for r in forward],factors)
            self.assertEqual(forward[0]['source_mean_TPM'],reverse[0]['target_mean_TPM'])
            self.assertEqual(forward[0]['target_mean_TPM'],reverse[0]['source_mean_TPM'])
            self.assertEqual(float(forward[0]['log2FoldChange']),-float(reverse[0]['log2FoldChange']))
            self.assertEqual(json.loads((root/'forward.json').read_text())['fc_limits'],json.loads((root/'reverse.json').read_text())['fc_limits'])
            self.assertEqual(forward[0]['source_log2_TPM_plus1'],__import__('math').log2(3))

if __name__=='__main__':unittest.main()
