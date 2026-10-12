#!/usr/bin/env python3
"""Reuse the qualified mapper fixture with explicitly relaxed synthetic pruning."""
import csv,json
from pathlib import Path
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow')
fixture=root/'inputs/hic_presplit_fixture_v01'
out=root/'inputs/diffhic_count_fixture_v01';out.mkdir(exist_ok=True)
with (fixture/'samplesheet.tsv').open() as f:
    reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
bam=root/'runs/hic_presplit_fixture_v18/results/raw/HiC/State_1_HiC_1/State_1_HiC_1.name_sorted.bam'
for r in rows:
    if r['assay']=='HiC':r.update(fastq_1='',fastq_2='',bam=str(bam))
with (out/'samplesheet.tsv').open('w') as f:
    writer=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
params=json.loads((fixture/'params.json').read_text())
params.update(stage='hic_diffhic_counts',samplesheet=str(out/'samplesheet.tsv'),raw_sample_names=[],fasta=str(fixture/'reference.fa'),diffhic_container=str(root/'containers/r_diffhic_v01.sif'),resolutions='500,1000',chosen_res=1000,hic_min_inward=0,hic_min_outward=0,hic_max_frag=1000)
(out/'params.json').write_text(json.dumps(params,indent=2)+'\n')
(out/'fixture_identity.json').write_text(json.dumps(dict(status='SYNTHETIC',input='same eight-pair qualified synthetic BAM used twice under distinct test sample labels',pruning='Explicit zero inward/outward distances and 1000-bp max fragment for geometry fixture; historical production defaults remain 1000/1000/800',purpose='Verify restriction assignment, pruning, squareCounts at two resolutions and joint loess offsets; not biological or historical equivalence'),indent=2)+'\n')
