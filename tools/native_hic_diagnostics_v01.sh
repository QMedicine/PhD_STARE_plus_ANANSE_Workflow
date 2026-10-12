#!/usr/bin/env bash
#SBATCH --job-name=phd_hic_diagnostics
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
run="$TASK_ROOT/runs/hic_diagnostics_v01_$SLURM_JOB_ID";mkdir -p "$run";cd "$run"
image="$TASK_ROOT/containers/hicexplorer_v01.sif"
exec_hic(){ apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$image" "$@"; }
manifest=${PHDSA_HIC_DIAGNOSTIC_MANIFEST:?Supply the state/cooler TSV}
mapfile -t states < <(awk 'NR>1{print $1}' "$manifest")
mapfile -t matrices < <(awk 'NR>1{print $2}' "$manifest")
[[ ${#states[@]} -eq 2 && ${#matrices[@]} -eq 2 ]] || { echo 'Diagnostic qualification requires exactly two explicit matrices';exit 2; }
for i in 0 1;do
 [[ ${states[$i]} =~ ^[A-Za-z0-9_.-]+$ ]] || exit 2
 exec_hic hicConvertFormat --matrices "${matrices[$i]}" --inputFormat cool --outputFormat h5 --outFileName "${states[$i]}.h5"
done
exec_hic hicNormalize --matrices "${states[0]}.h5" "${states[1]}.h5" --normalize smallest --outFileName "${states[0]}.depth_normalized.h5" "${states[1]}.depth_normalized.h5"
for state in "${states[@]}";do
 exec_hic hicCorrectMatrix diagnostic_plot --matrix "$state.depth_normalized.h5" --plotName "$state.diagnostic.png"
done
cp "$manifest" inputs.tsv
printf '%s\n' 'Exploratory diagnostic before the unresolved historical blacklist step. Do not authenticate estimated thresholds as the lost original thresholds.' > status.txt
sha256sum "$image" *.png inputs.tsv > checksums.sha256
