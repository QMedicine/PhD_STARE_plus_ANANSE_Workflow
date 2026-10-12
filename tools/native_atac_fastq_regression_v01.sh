#!/usr/bin/env bash
#SBATCH --job-name=phd_atac_raw_regression
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
run="$TASK_ROOT/evidence/atac_fastq_regression_v01";mkdir -p "$run";cd "$run"
image="$TASK_ROOT/containers/raw_tools_v06.sif"
sample=SRR13186609_Young_HSC_B1_Rep1
folder="$TASK_ROOT/runs/raw_ATAC_pilot_v15/results/raw/ATAC/$sample"
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$image" \
 bedtools coverage -a /vast/projects/Immune_Aging/Itokawa_et_al_2022_STARE/union_bed/union_summit500.named.bed \
 -b "$folder/$sample.no_blacklist.bam" -counts > rebuilt.window_counts.tsv
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$image" python - <<'PY'
from pathlib import Path
import hashlib,json,subprocess
root=Path('/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow');name='SRR13186609_Young_HSC_B1_Rep1'
paths=list((root/'runs/atac_saved_union_v12/work').glob('*/*/'+name+'.window_counts.tsv'))
if len(paths)!=1:raise ValueError('Expected one qualified window-count producer')
expected=paths[0];observed=Path('rebuilt.window_counts.tsv');different=0;n=0;max_difference=0
with expected.open() as a,observed.open() as b:
 for x,y in __import__('itertools').zip_longest(a,b):
  if x is None or y is None:raise ValueError('Window counts have different row counts')
  left=x.rstrip('\n').split('\t');right=y.rstrip('\n').split('\t')
  if left[:4]!=right[:4]:raise ValueError('Window identities differ')
  delta=abs(int(left[4])-int(right[4]));different+=delta!=0;max_difference=max(max_difference,delta);n+=1
folder=root/'runs/raw_ATAC_pilot_v15/results/raw/ATAC'/name
expected_bam=Path('/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/filtered')/(name+'_filtered.bam')
total=lambda p:int(subprocess.check_output(['samtools','view','-c',str(p)]))
peaks=folder/'peaks'/(name+'_summits.bed')
if not peaks.is_file() or peaks.stat().st_size==0:raise ValueError('Actual single-end ATAC pilot has no peaks')
record=dict(sample=name,expected_count_producer=str(expected),n_windows=n,different_window_counts=different,max_window_count_difference=max_difference,window_counts_exact=different==0,expected_window_counts_sha256=hashlib.sha256(expected.read_bytes()).hexdigest(),observed_window_counts_sha256=hashlib.sha256(observed.read_bytes()).hexdigest(),expected_pre_blacklist_alignments=total(expected_bam),observed_pre_blacklist_alignments=total(folder/(name+'.filtered.bam')),called_summits=len(peaks.read_text().splitlines()),peak_universe_note='Recalled summit files are distinct from the preserved historical summit-window checkpoint')
Path('atac_fastq_regression.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
if not record['window_counts_exact']:raise SystemExit('ATAC raw count regression differs; preserve and investigate')
PY
