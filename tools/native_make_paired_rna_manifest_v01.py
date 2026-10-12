#!/usr/bin/env python3
"""Reuse synthetic paired fragments to check RNA trimming/STAR mate handling."""
import argparse,csv,json
from pathlib import Path

def make(source,root,fragment_root):
    source=Path(source);root=Path(root);fragment_root=Path(fragment_root);root.mkdir(parents=True,exist_ok=True)
    params=json.loads((source/'templates/params_RNA_fastq_pilot_v02.json').read_text())
    with (source/'templates/samplesheet_Old_and_Young_HSCs_v01.tsv').open() as f:
        reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
    states=list(dict.fromkeys(r['cell_state'] for r in rows));state_map={state:f'State_{i+1}' for i,state in enumerate(states)}
    selected=next(r for r in rows if r['assay']=='RNA')
    for i,row in enumerate(rows):row.update(sample_name=f'fixture_sample_{i+1:02d}',cell_state=state_map[row['cell_state']],condition=state_map[row['cell_state']],cell_type='Synthetic_geometry_test')
    selected.update(sample_name='paired_rna_fixture',paired='TRUE',fastq_1=str(fragment_root/'paired_atac_fixture_R1.fastq.gz'),fastq_2=str(fragment_root/'paired_atac_fixture_R2.fastq.gz'),bam='',peak='')
    sheet=root/'samplesheet.tsv'
    with sheet.open('w') as f:
        writer=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
    transitions=root/'source2target.tsv';transitions.write_text('source\ttarget\nState_1\tState_2\n')
    params.update(samplesheet=str(sheet),source2target=str(transitions),atac_modes=['corrected_column_order'],historical_header_samples=[],historical_data_samples=[],deseq_state_order=[],plot_settings={},stage='raw_prepare',raw_sample_names=['paired_rna_fixture'],bowtie_index=None,atac_genome_size=None)
    (root/'params.json').write_text(json.dumps(params,indent=2)+'\n')
    identity=json.loads((fragment_root/'fixture_identity.json').read_text());identity.update(purpose='Synthetic mate-handling test for RNA trimming/STAR; no expression/historical-equivalence claim')
    (root/'fixture_identity.json').write_text(json.dumps(identity,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for key in ['source','root','fragment-root']:p.add_argument('--'+key,required=True)
    make(**vars(p.parse_args()))
