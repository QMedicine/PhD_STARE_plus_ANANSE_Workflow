#!/usr/bin/env python3
"""Generate a bounded, explicitly synthetic paired-end ATAC geometry fixture."""
import argparse,csv,gzip,hashlib,json,subprocess
from pathlib import Path

def make(source,root):
    source=Path(source);root=Path(root)
    if Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow') not in root.resolve().parents:raise ValueError('Fixture must stay on workflow scratch')
    root.mkdir(parents=True,exist_ok=True)
    params=json.loads((source/'templates/params_RNA_fastq_pilot_v02.json').read_text())
    fasta=Path(params['fasta'])
    if not Path(str(fasta)+'.fai').exists():raise ValueError('Read-only reference requires an existing FASTA index')
    # Use the first reference contig, retaining its actual chromosome naming style.
    chrom=Path(str(fasta)+'.fai').read_text().splitlines()[0].split('\t')[0]
    start=30000000;region=f'{chrom}:{start+1}-{start+20000}'
    text=subprocess.check_output(['samtools','faidx',str(fasta),region],text=True)
    sequence=''.join(text.splitlines()[1:]).upper()
    complement=str.maketrans('ACGT','TGCA');read1=[];read2=[];positions=[]
    for i in range(2000):
        local=i%1800+(6000 if i%2 else 0);length=250+i%201
        a=sequence[local:local+75];b=sequence[local+length-75:local+length].translate(complement)[::-1]
        if len(a)!=75 or len(b)!=75 or set(a+b)-set('ACGT'):continue
        name=f'synthetic_fragment_{i:05d}';read1.append(f'@{name}/1\n{a}\n+\n'+75*'I'+'\n');read2.append(f'@{name}/2\n{b}\n+\n'+75*'I'+'\n');positions.append((start+local,start+local+length))
    if len(read1)<1000:raise ValueError('Reference region did not yield enough defined synthetic fragments')
    paths=[]
    for mate,reads in [(1,read1),(2,read2)]:
        path=root/f'paired_atac_fixture_R{mate}.fastq.gz';paths.append(path)
        with path.open('wb') as f:
            with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=1) as compressed:compressed.write(''.join(reads).encode())
    with (source/'templates/samplesheet_Old_and_Young_HSCs_v01.tsv').open() as f:
        reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
    states=list(dict.fromkeys(r['cell_state'] for r in rows));state_map={state:f'State_{i+1}' for i,state in enumerate(states)}
    selected=next(r for r in rows if r['assay']=='ATAC')
    for i,row in enumerate(rows):
        row['sample_name']=f'fixture_sample_{i+1:02d}';row['cell_state']=state_map[row['cell_state']];row['condition']=row['cell_state'];row['cell_type']='Synthetic_geometry_test'
    selected.update(sample_name='paired_atac_fixture',paired='TRUE',fastq_1=str(paths[0]),fastq_2=str(paths[1]),bam='',peak='')
    sheet=root/'samplesheet.tsv'
    with sheet.open('w') as f:
        writer=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
    transitions=root/'source2target.tsv';transitions.write_text('source\ttarget\nState_1\tState_2\n')
    params.update(samplesheet=str(sheet),source2target=str(transitions),atac_modes=['corrected_column_order'],historical_header_samples=[],historical_data_samples=[],deseq_state_order=[],plot_settings={},stage='raw_prepare',atac_samtools='samtools',atac_sort_threads=32,raw_sample_names=['paired_atac_fixture'],rna_adapters=None,star_index=None)
    (root/'params.json').write_text(json.dumps(params,indent=2)+'\n')
    receipt=dict(status='SYNTHETIC_NOT_BIOLOGICAL_DATA',reference=str(fasta),region=region,fragments=len(read1),read_length=75,fragment_length_range=[250,450],paired_fastqs=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths],purpose='Qualify generic labels, paired Bowtie2/Picard/filtering and MACS3 BAMPE; no historical-equivalence claim')
    (root/'fixture_identity.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--root',required=True);a=p.parse_args();make(a.source,a.root)
