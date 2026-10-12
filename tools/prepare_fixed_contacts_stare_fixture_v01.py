#!/usr/bin/env python3
import json
from pathlib import Path
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow');fixture=root/'inputs/fixed_contacts_stare_fixture_v01';fixture.mkdir(exist_ok=True)
params=json.loads((root/'inputs/diffhic_count_fixture_v01/params_fixed_contacts_v01.json').read_text())
gtf=fixture/'synthetic.gtf';lines=[]
for i,start in enumerate([1001,5001,8501],1):
    attributes=f'gene_id "SYNTH_GENE_{i}"; transcript_id "SYNTH_TX_{i}"; gene_name "SYNTH{i}"; gene_biotype "protein_coding";'
    for feature in ['gene','transcript','exon']:lines.append(f'test_contig\tsynthetic\t{feature}\t{start}\t{start+500}\t.\t+\t.\t{attributes}\n')
gtf.write_text(''.join(lines));activity=fixture/'activity.bed'
activity.write_text('#chrom\tstart\tend\tpeak_id\tsynthetic_rep1\tsynthetic_rep2\n'+''.join(f'test_contig\t{start}\t{start+501}\tpeak_{i}\t{20+i}\t{30+i}\n' for i,start in enumerate([750,2750,4250,6250,8750],1)))
base='/vast/projects/Immune_Aging/Itokawa_et_al_2022_nfcore_tfactivity/tfactivity_results/HSC_50kbp_window_H13CORE_motifs/02_motifs/06_psem/motifs.psem'
params.update(fixture_manifest=str(root/'runs/diffhic_fixed_contacts_fixture_v23/results/provenance/validated_inputs.json'),fixture_activity=str(activity),gtf=str(gtf),psem=base,stare_window=20000,enhancer_window=20000,stare_container='/vast/projects/Immune_Aging/apptainer_singularity_cache/stare-abc_1.0.5.sif',outdir=str(root/'runs/fixed_contacts_stare_smoke_v24/results'),trace_suffix='smoke')
(fixture/'params.json').write_text(json.dumps(params,indent=2)+'\n')
