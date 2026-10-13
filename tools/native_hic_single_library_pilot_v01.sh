#!/usr/bin/env bash
#SBATCH --job-name=phd_hic_single_pilot
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load nextflow/26.04.3
module swap apptainer/1.4.1 apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
SOURCE_DIR="$TASK_ROOT/source_v37"
RUN_DIR="$TASK_ROOT/runs/hic_single_library_pilot_v36"
FIXTURE_DIR="$TASK_ROOT/inputs/hic_single_library_pilot_v36"
mkdir -p "$RUN_DIR" "$TASK_ROOT/tmp/hic_single_library_pilot_v36"
export NXF_HOME="$TASK_ROOT/cache/nextflow" NXF_APPTAINER_CACHEDIR="$TASK_ROOT/containers"
export APPTAINER_CACHEDIR="$TASK_ROOT/cache/apptainer_hic_single_library_pilot_v36" APPTAINER_TMPDIR="$TASK_ROOT/tmp/hic_single_library_pilot_v36"
export NXF_OPTS='-Xms512m -Xmx6g' APPTAINERENV_PREPEND_PATH="$SOURCE_DIR/bin"
cd "$RUN_DIR"
if [[ ${PHDSA_RESUME:-FALSE} != TRUE ]]; then
 apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$TASK_ROOT/containers/hic_raw_standard_v01.sif" \
  python "$SOURCE_DIR/tools/make_hic_single_library_fixture_v01.py" --samplesheet "$SOURCE_DIR/templates/samplesheet_HiC_single_library_practice_v01.tsv" --outdir "$FIXTURE_DIR" --records 5000
 python3 - "$SOURCE_DIR" "$FIXTURE_DIR" "$RUN_DIR" <<'PY'
import json,sys
from pathlib import Path
source,fixture,run=map(Path,sys.argv[1:]);params=json.loads((source/'templates/params_HiC_single_library_full_raw_v01.json').read_text())
params.update(raw_sample_names=[],stage='hic_diffhic_counts',samplesheet=str(fixture/'samplesheet.tsv'),process_fastq=True,hic_contacts_from_diffhic=True,diffhic_container='/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow/containers/r_diffhic_v01.sif',resolutions='100000',chosen_res=100000,hic_count_filter=0)
(run/'params.json').write_text(json.dumps(params,indent=2)+'\n')
(run/'pilot.config').write_text('process { withName: HIC_PRESPLIT_ALIGN { cpus = 4; memory = "16 GB"; time = "1h" } }\n')
PY
fi
cmd=(nextflow run "$SOURCE_DIR/main.nf" -profile wehi -c "$RUN_DIR/pilot.config" -params-file "$RUN_DIR/params.json" --trace_suffix "$SLURM_JOB_ID" --outdir "$RUN_DIR/results" -work-dir "$RUN_DIR/work")
if [[ ${PHDSA_RESUME:-FALSE} == TRUE ]]; then cmd+=(-resume); fi
"${cmd[@]}"
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$TASK_ROOT/containers/hic_raw_standard_v01.sif" \
 python "$SOURCE_DIR/tools/qualify_hic_single_library_pilot_v01.py" --run "$RUN_DIR" --samplesheet "$FIXTURE_DIR/samplesheet.tsv" --records 5000
