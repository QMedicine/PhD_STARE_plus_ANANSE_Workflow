#!/usr/bin/env python3
"""Two distinguishable library components per synthetic biological pool."""
import csv,json
from pathlib import Path
import pysam
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow');out=root/'inputs/hic_pooling_fixture_v01';out.mkdir(exist_ok=True)
params=json.loads((root/'inputs/diffhic_count_fixture_v01/params_fixed_contacts_v01.json').read_text())
with open(params['samplesheet']) as f:r=csv.DictReader(f,delimiter='\t');fields=r.fieldnames+['pool_group'];old=list(r)
bam=Path(next(r['bam'] for r in old if r['assay']=='HiC'));rows=[]
for row in old:
    if row['assay']!='HiC':rows.append(dict(row,pool_group=''));continue
    for component in ['A','B']:
        name=row['cell_state']+'_HiC_'+component;path=out/(name+'.bam')
        with pysam.AlignmentFile(bam,'rb') as inp,pysam.AlignmentFile(path,'wb',template=inp) as target:
            for read in inp:read.query_name=component+'_'+read.query_name;target.write(read)
        rows.append(dict(row,sample_name=name,bam=str(path),pool_group=row['cell_state']+'_POOL'))
with (out/'samplesheet.tsv').open('w') as f:w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
params.update(samplesheet=str(out/'samplesheet.tsv'),hic_raw_container=str(root/'containers/hic_raw_standard_v01.sif'))
(out/'params.json').write_text(json.dumps(params,indent=2)+'\n')
(out/'fixture_identity.json').write_text(json.dumps(dict(status='SYNTHETIC',components=4,biological_pools=2,qnames='Component-specific prefixes prevent unintended mate merging',expected_retained_pairs_per_pool=14,expected_preserved_library_replicates=2),indent=2)+'\n')
