"""Verify new identity guards against the completed full historical checkpoint runs."""
import csv,hashlib,importlib.util,json,os,sys
from pathlib import Path
SOURCE,TASK=map(Path,sys.argv[1:]);OUT=Path.cwd()
sys.path.insert(0,str(SOURCE/'bin'))
import assemble_atac,merge_recount,tpm_filter,prepare_ananse_de,validate_inputs

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
checks=[]
def compare(label,new,old):
    record=dict(check=label,observed_sha256=digest(new),reference_sha256=digest(old),byte_identical=digest(new)==digest(old))
    checks.append(record)
    if not record['byte_identical']:raise ValueError('Guarded replay changed '+label)

params=json.loads((SOURCE/'templates/params_rna_recount_v04.json').read_text())
sheet=SOURCE/'templates/samplesheet_Old_and_Young_HSCs_v01.tsv'
atac=TASK/'runs/atac_saved_union_v12/work'
atac_files=list(atac.glob('*/*/*.window_counts.tsv'))
names={r['sample_name'] for r in csv.DictReader(sheet.open(),delimiter='\t') if r['assay']=='ATAC'}
by_name={p.name.removesuffix('.window_counts.tsv'):p for p in atac_files}
if set(by_name)!=names or len(atac_files)!=4:raise ValueError('Actual ATAC producer identity discovery is ambiguous')
for task in atac.glob('*/*'):
    receipt=task/'column_assignment.json'
    if not receipt.is_file():continue
    assignment=json.loads(receipt.read_text());mode=assignment['mode'];folder=OUT/mode;folder.mkdir();os.chdir(folder)
    # Reverse lexical file order intentionally: outputs must depend on explicit identities.
    assemble_atac.assemble(sheet,list(reversed(sorted(atac_files))),mode,'union_counts.tsv',assignment['header_samples'],assignment['data_samples'])
    for name in ['union_counts.tsv','union_base.tsv','column_assignment.json']:compare(mode+'/'+name,folder/name,task/name)
os.chdir(OUT)
rna=TASK/'runs/rna_recount_v08/work/d4/83d493f243e5925e194af9d5840ed2'
rna_files=list((rna/'counts').glob('*.featureCounts.txt'))
folder=OUT/'rna';folder.mkdir();os.chdir(folder)
merge_recount.merge(params['selected_counts'],list(reversed(sorted(rna_files))),'selected_counts.csv')
compare('RNA selected counts',folder/'selected_counts.csv',rna/'selected_counts.csv')
baseline=TASK/'runs/rna_recount_v08/results/rna/tpm'
tpm_filter.calculate('selected_counts.csv',params['gene_lengths'],params['idmap'])
for name in ['counts.clean.tsv','counts.tpm.tsv','counts.counts_filtered.tsv','TFs.counts_filtered.tsv','counts.tpm_filtered.tsv','TFs.tpm_filtered.tsv']:
    old=baseline/name
    if not old.is_file():old=TASK/'runs/rna_recount_v08/work/c3/b18fe2f9d8fbf58872f9091f5d8417'/name
    compare('RNA '+name,folder/name,old)
os.chdir(OUT)
for name in ['samplesheet_Old_and_Young_HSCs_v01.tsv','samplesheet_HiC_full_raw_cohort_v01.tsv','samplesheet_HiC_single_library_practice_v01.tsv']:
    manifest=validate_inputs.validate(SOURCE/'templates'/name,SOURCE/'templates/source2target_Old_and_Young_HSCs_v01.tsv',check_files=True)
    checks.append(dict(check='actual_manifest_'+name,samples=len(manifest['samples']),distinct_input_paths=True,transitions=manifest['transitions']))
(OUT/'guarded_replay.json').write_text(json.dumps(dict(status='PASS',source_snapshot=str(SOURCE),checks=checks),indent=2)+'\n')
print(json.dumps(dict(status='PASS',checks=len(checks)),indent=2))
