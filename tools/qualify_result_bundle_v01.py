#!/usr/bin/env python3
"""Validate producer closure and preserve qualified GRNs/tables with explicit report limits."""
import argparse,csv,hashlib,json,math,shutil
from pathlib import Path
import numpy as np
import pandas as pd

SCRATCH=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow')
PROJECT=Path('/vast/projects/Immune_Aging/PhD_STARE_plus_ANANSE_Workflow')
NORMALIZED=['influence_score','target_score_scaled','G_score','G_score_scaled','direct_targets','factor_fc']

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()

def compare_influence(observed,expected,diff_observed,diff_expected):
    a=pd.read_csv(expected,sep='\t',index_col='factor');b=pd.read_csv(observed,sep='\t',index_col='factor')
    if a.index.has_duplicates or b.index.has_duplicates or set(a.index)!=set(b.index) or list(a.columns)!=list(b.columns):raise ValueError('Ranking identities/components differ')
    aligned=b.loc[a.index]
    exact={c:bool(((a[c]==aligned[c])|(a[c].isna()&aligned[c].isna())).all()) for c in a.columns}
    delta=np.abs(a.to_numpy(float)-aligned.to_numpy(float))
    sequence_equal=np.array_equal(a.influence_score.to_numpy(),b.influence_score.to_numpy())
    x=pd.read_csv(diff_expected,sep='\t').set_index(['source','target']);y=pd.read_csv(diff_observed,sep='\t').set_index(['source','target'])
    edge_ids=not x.index.has_duplicates and not y.index.has_duplicates and set(x.index)==set(y.index)
    edge_weights=edge_ids and bool((x.weight==y.loc[x.index].weight).all())
    record=dict(component_exact=exact,raw_max_absolute_difference=float(np.nanmax(delta)),score_sequence_exact=sequence_equal,top20_set_exact=set(a.head(20).index)==set(b.head(20).index),top20_order_exact=list(a.head(20).index)==list(b.head(20).index),diff_edge_ids_exact=edge_ids,diff_weights_exact=edge_weights,ranking_bytes_exact=sha(expected)==sha(observed),expected_sha256=sha(expected),observed_sha256=sha(observed))
    if not all(exact[c] for c in NORMALIZED) or not sequence_equal or not record['top20_set_exact'] or not edge_weights:raise ValueError('Normalized influence or differential-edge regression failed: '+json.dumps(record))
    return record

def qualify(run,destination,regression,references,software_lock,copy=False):
    run=Path(run).resolve();destination=Path(destination).resolve()
    if SCRATCH not in run.parents or PROJECT not in destination.parents:raise ValueError('Unexpected scratch/project boundary')
    traces=list((run/'results/provenance').glob('trace_*.tsv'))
    if len(traces)!=1:raise ValueError('Choose a run with exactly one producer trace')
    with traces[0].open() as f:tasks=list(csv.DictReader(f,delimiter='\t'))
    if not tasks or any(r['status'] not in ['COMPLETED','CACHED'] or r['exit']!='0' for r in tasks):raise ValueError('Incomplete producer trace')
    params=json.loads((run/'provenance/params.json').read_text())
    if sha(run/'provenance/params.json')!=(run/'provenance/params.sha256').read_text().split()[0]:raise ValueError('Run parameters changed after launch')
    if params['stage']!='atac_historical_hic_checkpoint' or params.get('make_reports',True):raise ValueError('This bundle qualification requires the explicit reports-disabled saved-checkpoint route')
    manifest=json.loads((run/'results/provenance/validated_inputs.json').read_text())
    for process,expected in [('STARE_GABC',len(params['atac_modes'])*len({s for t in manifest['transitions'] for s in [t['source'],t['target']]})),('STARE_TO_NETWORK',len(params['atac_modes'])*len({s for t in manifest['transitions'] for s in [t['source'],t['target']]})),('ANANSE_INFLUENCE',len(params['atac_modes'])*len(manifest['transitions']))]:
        if sum(':'+process+' ' in r['name'] for r in tasks)!=expected:raise ValueError('Producer count does not close for '+process)
    ref_rows=list(csv.DictReader(Path(references).open(),delimiter='\t'))
    affinity=json.loads(Path(regression).read_text())
    historical=[r for r in affinity if r['mode']=='historical_column_order']
    if not historical or any(not r.get('numeric_exact',r.get('network_byte_match',False)) for r in historical):raise ValueError('Historical affinity/GRN regression is not exact')
    files=[];regressions=[];mode_results={};top20_records=[]
    for mode in params['atac_modes']:
        mode_results[mode]={}
        states={s for t in manifest['transitions'] for s in [t['source'],t['target']]}
        for state in sorted(states):
            network=run/'results'/mode/'networks'/state/(state+'_SYMBOLS.tsv')
            n=0
            with network.open() as f:
                reader=csv.DictReader(f,delimiter='\t')
                if reader.fieldnames!=['tf_target','prob']:raise ValueError('Network must use the original two-column hybrid format')
                for row in reader:
                    prob=float(row['prob'])
                    if not row['tf_target'] or not math.isfinite(prob) or not 0<=prob<=1:raise ValueError('Invalid regulatory score')
                    n+=1
            if not n:raise ValueError('Empty GRN')
            recorded=next(r for r in affinity if r['mode']==mode and r['cell_state']==state and 'observed_network_sha256' in r)
            if sha(network)!=recorded['observed_network_sha256']:raise ValueError('Network changed after scientific regression')
            files.append((network,network.relative_to(run/'results')))
            mode_results[mode][state]=dict(edges=n,sha256=sha(network),historical_bytes_exact=recorded['network_byte_match'])
        for t in manifest['transitions']:
            source,target=t['source'],t['target'];folder=run/'results'/mode/(source+'_to_'+target)/'tables'
            table=pd.read_csv(folder/'influence.txt',sep='\t')
            if len(table)<20 or not set(['factor']+NORMALIZED).issubset(table):raise ValueError('Incomplete influence components')
            if table.factor.duplicated().any() or not table.influence_score.is_monotonic_decreasing:raise ValueError('Invalid ranking order')
            if (table.factor_fc<=0).any():raise ValueError('Non-positive directional candidate TF')
            compat=json.loads((folder/'de_compatibility.json').read_text())
            if not compat['data_bytes_unchanged']:raise ValueError('DE body was altered by the compatibility adapter')
            if mode=='historical_column_order':
                matches=[r for r in ref_rows if r['source']==source and r['target']==target]
                if len(matches)!=1:raise ValueError('Missing reference comparison')
                r=matches[0];regressions.append(dict(source=source,target=target,**compare_influence(folder/'influence.txt',r['ranking'],folder/'influence_diffnetwork.tsv',r['diffnetwork'])))
            for path in sorted(folder.iterdir()):
                if path.name not in ['influence.txt','influence_diffnetwork.tsv','directed_de.tsv','de_compatibility.json']:raise ValueError('Unexpected table product '+str(path))
                files.append((path,path.relative_to(run/'results')))
            for position,row in table.head(20).iterrows():top20_records.append(dict(mode=mode,source=source,target=target,position=position+1,**{c:row[c] for c in ['factor']+NORMALIZED}))
        for name in ['column_assignment.json','loessOffsets_qc.tsv','versions.txt']:
            path=run/'results'/mode/'atac'/name;files.append((path,path.relative_to(run/'results')))
    for name in ['params.json','params.sha256','profile_defaults.config','nextflow_version.txt','apptainer_version.txt']:
        path=run/'provenance'/name;files.append((path,Path('provenance')/name))
    files.extend([(traces[0],Path('provenance')/traces[0].name),(run/'results/provenance/validated_inputs.json',Path('provenance/validated_inputs.json'))])
    for name,path in [('samplesheet.tsv',params['samplesheet']),('source2target.tsv',params['source2target']),('software_lock.json',software_lock),('affinity_regression.json',regression),('influence_references.tsv',references)]:files.append((Path(path),Path('provenance')/name))
    plan=[dict(source=str(path),relative_path=str(relative),bytes=path.stat().st_size,sha256=sha(path)) for path,relative in files]
    receipt=dict(status='QUALIFIED_GRNS_AND_TABLES_REPORTS_INCOMPLETE',run=str(run),destination=str(destination),producer_tasks=len(tasks),producer_trace_sha256=sha(traces[0]),modes=mode_results,historical_influence=regressions,report_status='WITHHELD_STRICT_PLOTTING_FAILURE',upstream_scope='Pre-blacklist ATAC BAMs; preserved summit windows, corrected Hi-C H5, RNA TPM and directional H13CORE DE checkpoints',files=plan)
    output=run/'provenance/qualified_bundle_v01.json';output.write_text(json.dumps(receipt,indent=2)+'\n')
    if copy:
        # Idempotent copying is permitted only when all destination bytes agree.
        destination.mkdir(parents=True,exist_ok=True)
        for entry in plan:
            target=destination/entry['relative_path'];target.parent.mkdir(parents=True,exist_ok=True)
            if target.exists():
                if sha(target)!=entry['sha256']:raise ValueError('Refuse to overwrite differing durable product '+str(target))
            else:shutil.copy2(entry['source'],target)
            if target.stat().st_size!=entry['bytes'] or sha(target)!=entry['sha256']:raise ValueError('Durable copy differs')
        top20=pd.DataFrame(top20_records);top20.to_csv(destination/'top20_components_by_mode.tsv',sep='\t',index=False)
        comparisons=[]
        for t in manifest['transitions']:
            subset=top20[(top20.source==t['source'])&(top20.target==t['target'])]
            names={mode:set(subset[subset['mode']==mode].factor) for mode in params['atac_modes']}
            for factor in sorted(set.union(*names.values())):comparisons.append(dict(source=t['source'],target=t['target'],factor=factor,**{mode:factor in values for mode,values in names.items()}))
        pd.DataFrame(comparisons).to_csv(destination/'top20_membership_comparison.tsv',sep='\t',index=False)
        receipt['generated_tables']=[dict(path=name,sha256=sha(destination/name)) for name in ['top20_components_by_mode.tsv','top20_membership_comparison.tsv']]
        receipt['copy_verified']=True;output.write_text(json.dumps(receipt,indent=2)+'\n');shutil.copy2(output,destination/'qualification.json')
    print(json.dumps({k:v for k,v in receipt.items() if k!='files'},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for key in ['run','destination','regression','references','software-lock']:p.add_argument('--'+key,required=True)
    p.add_argument('--copy',action='store_true');qualify(**vars(p.parse_args()))
