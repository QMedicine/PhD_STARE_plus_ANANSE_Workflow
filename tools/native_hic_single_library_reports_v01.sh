#!/usr/bin/env bash
#SBATCH --job-name=phd_single_hic_reports
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load nextflow/26.04.3
module swap apptainer/1.4.1 apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
PROJECT_RUN=/vast/projects/Immune_Aging/PhD_STARE_plus_ANANSE_Workflow/Old_and_Young_HSCs
SOURCE_DIR="$TASK_ROOT/source_v39"
RUN_DIR="$TASK_ROOT/runs/hic_single_library_reports_v38"
mkdir -p "$RUN_DIR" "$TASK_ROOT/tmp/hic_single_library_reports_v38"
export NXF_HOME="$TASK_ROOT/cache/nextflow" NXF_APPTAINER_CACHEDIR="$TASK_ROOT/containers"
export APPTAINER_CACHEDIR="$TASK_ROOT/cache/apptainer_hic_single_library_reports_v38" APPTAINER_TMPDIR="$TASK_ROOT/tmp/hic_single_library_reports_v38"
export NXF_OPTS='-Xms512m -Xmx6g'
export APPTAINERENV_PREPEND_PATH="$SOURCE_DIR/bin"
cd "$RUN_DIR"
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$TASK_ROOT/containers/ananse_v03.sif" \
 python "$SOURCE_DIR/tools/build_report_replay_manifest_v01.py" \
 --results "$PROJECT_RUN/validated_results/hic_single_library_to_stare_v38" \
 --de-manifest "$PROJECT_RUN/de_Old_and_Young_HSCs_v01.tsv" --output "$RUN_DIR/report_manifest.tsv"
python3 - "$SOURCE_DIR" "$RUN_DIR" <<'PY'
import json,sys
from pathlib import Path
source,run=map(Path,sys.argv[1:]);params=json.loads((source/'templates/params_HiC_single_library_to_STARE_v01.json').read_text())
params['report_manifest']=str(run/'report_manifest.tsv')
(run/'params.json').write_text(json.dumps(params,indent=2)+'\n')
PY
nextflow run "$SOURCE_DIR/report_replay.nf" -c "$SOURCE_DIR/nextflow.config" -profile wehi \
 -params-file "$RUN_DIR/params.json" --trace_suffix "$SLURM_JOB_ID" --outdir "$RUN_DIR/results" -work-dir "$RUN_DIR/work"
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$TASK_ROOT/containers/ananse_v03.sif" \
 python "$SOURCE_DIR/tools/qualify_selected_reports_v01.py" --results "$RUN_DIR/results" --manifest "$RUN_DIR/report_manifest.tsv" \
 --destination "$PROJECT_RUN/validated_results/hic_single_library_reports_v38"
