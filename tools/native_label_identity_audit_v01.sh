#!/usr/bin/env bash
#SBATCH --job-name=phd_label_identity_audit
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
SOURCE_DIR=${1:?Frozen source directory required}
[[ "$SOURCE_DIR" == "$TASK_ROOT"/source_* ]] || exit 2
run="$TASK_ROOT/runs/label_identity_audit_v01_$SLURM_JOB_ID"
mkdir -p "$run";cd "$run"
export APPTAINER_CACHEDIR="$TASK_ROOT/apptainer_cache"
python_image=/vast/projects/Immune_Aging/apptainer_singularity_cache/community-cr-prod.seqera.io-docker-registry-v2-blobs-sha256-7c-7c256e63e08633ac420692d3ceec1f554fe4fcc794e5bdd331994f743096a46d-data.img
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" --env PHDSA_TEST_SOURCE="$SOURCE_DIR" \
 "$python_image" python -m unittest discover -s "$SOURCE_DIR/tests" -v > contracts.txt 2>&1
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" \
 "$python_image" python "$SOURCE_DIR/tools/native_label_identity_replay_v01.py" "$SOURCE_DIR" "$TASK_ROOT" > replay.log 2>&1
# Reproduce the older archive's incompatible sort/join comparators with GNU coreutils.
printf 'peak_1\t11\npeak_2\t22\npeak_10\t110\n' | LC_ALL=C sort -k1,1V > version_sorted.tsv
set +e
LC_ALL=C join --check-order -t $'\t' version_sorted.tsv version_sorted.tsv > join_output.tsv 2> join_error.txt
join_status=$?
set -e
[[ "$join_status" -ne 0 ]] || { echo 'Expected version/lexical join incompatibility was not detected'; exit 1; }
printf '%s\n' "$join_status" > join_exit_status.txt
# Check the normalization adapter on both actual complete peak matrices.
for task in "$TASK_ROOT"/runs/atac_saved_union_v12/work/*/*; do
 if [[ -f "$task/.command.sh" ]] && grep -q '^normalize_atac.R ' "$task/.command.sh"; then
  mkdir -p "normalization/$(basename "$task")"
  (
   cd "normalization/$(basename "$task")"
   apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$TASK_ROOT/containers/r_plot_v01.sif" \
    Rscript "$SOURCE_DIR/bin/normalize_atac.R" "$task/union_counts.tsv" "$task/union_base.tsv" guarded_activity.tsv > r.log 2>&1
   cmp guarded_activity.tsv "$task/activity.tsv"
   sha256sum guarded_activity.tsv > activity.sha256
  )
 fi
done
printf 'LABEL_IDENTITY_AUDIT_PASS\n'
