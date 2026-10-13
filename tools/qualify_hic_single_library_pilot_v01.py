#!/usr/bin/env python3
"""Qualify read conservation and direct routing of the actual six paired prefixes."""
import argparse,csv,json,subprocess
from pathlib import Path
import pysam

def qualify(run,samplesheet,records):
 run=Path(run);result=run/'results'
 with Path(samplesheet).open() as f:reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
 hic=[r for r in rows if r['assay']=='HiC']
 if 'pool_group' in fields:raise ValueError('Expected a samplesheet without the pooling column')
 traces=list((result/'provenance').glob('trace_*.tsv'));trace=[]
 for p in traces:
  with p.open() as f:trace.extend(csv.DictReader(f,delimiter='\t'))
 if any('DIFFHIC_POOL_LIBRARIES' in r['name'] for r in trace):raise ValueError('Unexpected merge task for standalone libraries')
 checks=[]
 for row in hic:
  name=row['sample_name'];bam=result/'raw/HiC'/name/(name+'.name_sorted.bam')
  subprocess.run(['samtools','quickcheck',str(bam)],check=True)
  primary=[0,0];mapped=[0,0]
  with pysam.AlignmentFile(bam,'rb') as f:
   if f.header.to_dict().get('HD',{}).get('SO')!='queryname':raise ValueError('Name-sorted BAM required')
   for read in f:
    if not read.is_paired or read.is_read1==read.is_read2 or read.is_supplementary:raise ValueError('Invalid mate/segment flags')
    if not read.is_secondary:
     mate=0 if read.is_read1 else 1;primary[mate]+=1;mapped[mate]+=not read.is_unmapped
  if primary!=[records,records]:raise ValueError('Primary mate count differs from the input prefix')
  if not (result/'diffhic/pairs'/name/(name+'.trimmed.h5')).is_file():raise ValueError('Per-library retained-pair output missing')
  checks.append(dict(sample_name=name,primary_mates=primary,mapped_primary_mates=mapped))
 contact=json.loads((result/'diffhic/fixed_grid_contacts/fixed_contacts_provenance.json').read_text())
 if set(contact['library_order'])!={r['sample_name'] for r in hic}:raise ValueError('Library identities changed during counting')
 if set(contact['states'])!={r['cell_state'] for r in hic}:raise ValueError('State identities changed')
 if not list((result/'diffhic/counts').glob('interactions_*.rds')):raise ValueError('Rounded count audit missing')
 receipt=dict(status='QUALIFIED_SINGLE_LIBRARY_PAIRED_PREFIX_ROUTE',full_cohort_complete=False,prefix_records_per_library=records,libraries=checks,merge_tasks=0,retained_pair_libraries=len(hic),state_contacts=len(contact['states']),method=contact['method'])
 (result/'qualification.json').write_text(json.dumps(receipt,indent=2)+'\n');return receipt
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--run',required=True);p.add_argument('--samplesheet',required=True);p.add_argument('--records',type=int,required=True)
 print(json.dumps(qualify(**vars(p.parse_args())),indent=2))
