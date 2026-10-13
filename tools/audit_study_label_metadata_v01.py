#!/usr/bin/env python3
"""Cross-check study manifests against archived accession and experiment metadata."""
import csv,json,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'evidence/audit/label_identity_audit_v01'
def read(path,sep=','):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter=sep))
def run(value):
    match=re.search(r'SRR[0-9]+',value or '')
    return match.group() if match else None
refs={}
for path in (ROOT/'historical/raw_support').rglob('*sample_information*'):
    if path.suffix not in ['.csv','.tsv']:continue
    assay='ATAC' if 'ATACseq' in str(path) else 'RNA'
    for row in read(path,'\t' if path.suffix=='.tsv' else ','):
        refs[(assay,row['Run'])]=row
cohort=read(ROOT/'templates/samplesheet_Old_and_Young_HSCs_v01.tsv','\t')
checked=[]
for row in cohort:
    if row['assay']=='HiC':continue
    reference=refs[(row['assay'],run(row['sample_name']))]
    expected=(reference['Sample_ID'],reference['Age'],reference['Replicate'].removeprefix('Rep'))
    actual=(row['sample_name'],row['condition'],row['replicate'])
    if actual!=expected:raise ValueError('Accession metadata mismatch '+row['sample_name'])
    for field in ['bam','peak']:
        if row[field] and run(row[field])!=reference['Run']:raise ValueError('Input accession mismatch '+row[field])
    checked.append(dict(sample=row['sample_name'],assay=row['assay'],cell_state=row['cell_state'],archived_age=reference['AGE'],accession_and_replicate_match=True))
hic_path=ROOT/'historical/discovered/20240603_MMP_Hi_Low_HSCs/20240603_Experimental_information_quinnS_Hi-C_YandO_MMP_hi_low_HSC.csv'
hic=read(hic_path);lookup={row['Sample']+'_'+row['Rep']:row for row in hic}
for row in read(ROOT/'templates/samplesheet_HiC_full_raw_cohort_v01.tsv','\t'):
    if row['assay']!='HiC':continue
    original=lookup[row['sample_name']]
    if (Path(row['fastq_1']).name,Path(row['fastq_2']).name,row['condition'],row['replicate'])!=(original['Read1'],original['Read2'],original['Young.or.Old'],original['Rep'].removeprefix('Rep')):raise ValueError('Hi-C experiment metadata mismatch')
    checked.append(dict(sample=row['sample_name'],assay='HiC',pool_group=row['pool_group'],experiment_sheet_mates_condition_replicate_match=True))
# Check the exact historical comma-split parser against a proper CSV parser.
atac_path=ROOT/'historical/raw_support/Itokawa_et_al_2022_ATACseq/sample_information.csv'
lines=atac_path.read_text().splitlines();head=lines[0].replace('"','').split(',');mismatches=[]
for row,line in zip(read(atac_path),lines[1:]):
    for field in ['Read1','Read2','Sample_ID']:
        value=line.split(',')[head.index(field)].strip().strip('"')
        if value!=row[field]:mismatches.append(dict(run=row['Run'],field=field))
if mismatches:raise ValueError('Historical CSV parser assigned different inputs')
csvs=[]
for path in (ROOT/'historical').rglob('*.csv'):
    if not path.name.startswith(('samplesheet','design_matrix','diffbind_samplesheet')):continue
    rows=read(path);sites=0;conflicts=[]
    for row in rows:
        strings=[v for v in row.values() if isinstance(v,str)]
        accessions={run(v) for v in strings if run(v)}
        if len(accessions)>1:conflicts.append(dict(reason='multiple_accessions_in_one_row',accessions=sorted(accessions)))
        if accessions:
            accession=next(iter(accessions));ages={ref['Age'] for (assay,r),ref in refs.items() if r==accession}
            for value in strings:
                ages_in_label=set(re.findall(r'\b(Young|Old)\b',value.replace('_',' ')))
                if ages_in_label and ages_in_label!=ages:conflicts.append(dict(accession=accession,value=value,expected_age=sorted(ages)))
            sites+=1
    if conflicts:raise ValueError('Historical manifest label conflict '+str(path)+' '+str(conflicts))
    csvs.append(dict(path=str(path.relative_to(ROOT)),rows_with_accession=sites,conflicts=conflicts))
notes=[dict(source=str(hic_path.relative_to(ROOT)),sample=r['Sample']+'_'+r['Rep'],field='Mouse.Strain',archived_value=r['Mouse.Strain'],status='UNRESOLVED_SOURCE_METADATA_TYPO_OR_BIOLOGICAL_DIFFERENCE') for r in hic if r['Mouse.Strain']!='C57BL/6']
record=dict(status='PASS_WITH_METADATA_QUALIFICATIONS',sample_checks=checked,historical_csv_checks=csvs,historical_atac_csv_parser=dict(rows=len(read(atac_path)),fields=['Read1','Read2','Sample_ID'],mismatches=mismatches),metadata_qualifications=notes,limits='Archived metadata agrees with declared identities; this does not authenticate physical tube identity or the manual FASTQ merge.')
OUT.mkdir(parents=True,exist_ok=True);(OUT/'study_metadata_checks.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(dict(sample_checks=len(checked),csv_checks=len(csvs),metadata_qualifications=len(notes),status=record['status']),indent=2))
if __name__=='__main__':pass
