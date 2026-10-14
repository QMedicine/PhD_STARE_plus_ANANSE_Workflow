#!/usr/bin/env bash
#SBATCH --job-name=phd_state_tpm_smoke
#SBATCH --partition=regular
#SBATCH --account=wehi
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stacpoole.q@wehi.edu.au
set -euo pipefail
module purge
module load apptainer/1.3.5
TASK_ROOT=/vast/scratch/users/stacpoole.q/PhD_STARE_plus_ANANSE_Workflow
SOURCE_DIR=${1:?Frozen source required}
RUN_DIR="$TASK_ROOT/runs/state_tpm_smoke_$SLURM_JOB_ID"
mkdir -p "$RUN_DIR";cd "$RUN_DIR"
export APPTAINERENV_PREPEND_PATH="$SOURCE_DIR/bin"
apptainer exec --cleanenv --no-home --bind "$TASK_ROOT" --env PYTHONPATH="$SOURCE_DIR/bin" \
  "$TASK_ROOT/containers/ananse_v03.sif" python - "$SOURCE_DIR" <<'PY'
import importlib.util,sys
from pathlib import Path
import pandas as pd
from select_state_tpm import state_vector
spec=importlib.util.spec_from_file_location('converter',Path(sys.argv[1])/'historical/scripts/07_network_conversion/01_STARE_TF_Gene_Affinities_to_ananse_network.py')
converter=importlib.util.module_from_spec(spec);spec.loader.exec_module(converter)
Path('samples.tsv').write_text('sample_name\tassay\tcell_state\na\tRNA\tResting\nb\tRNA\tActivated\n')
Path('tpm.tsv').write_text('gene_id\tb\ta\nTF_A\t10\t2\nTF_B\t3\t8\nGENE_X\t5\t4\n')
affinity=pd.DataFrame({'TF_A':[1.,2.,3.],'TF_B':[2.,1.,3.]},index=['TF_A','TF_B','GENE_X'])
values=[]
for state in ['Resting','Activated']:
 tpm=pd.Series(state_vector('tpm.tsv','samples.tsv',state)[3])
 converter.make_network(affinity,tpm,tpm,state)
 network=pd.read_csv(Path(state)/'network.tsv',sep='\t').set_index('tf_target')
 values.append(float(network.loc['TF_A—GENE_X','prob']))
assert abs(values[0]-1/3)<1e-12 and abs(values[1]-9/14)<1e-12,values
print('Verified native multiplicative edge scores after within-state selection:',values)
PY
apptainer exec --cleanenv --no-home --bind /vast/projects/Immune_Aging,"$TASK_ROOT" "$TASK_ROOT/containers/ananse_v03.sif" \
 python - "$SOURCE_DIR" "$TASK_ROOT" <<'PY'
import csv,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(sys.argv[1])/'bin'))
from prepare_tf_expression import prepare
root=Path(sys.argv[2]);params=json.loads((root/'runs/state_specific_tpm_v42/params.json').read_text())
with Path(params['de_manifest']).open() as f:contrast=next(csv.DictReader(f,delimiter='\t'))
ranking=root/'runs/atac_saved_union_v12/results/corrected_column_order'/(contrast['source']+'_to_'+contrast['target'])/'tables/influence.txt'
prepare(ranking,contrast['path'],params['tpm'],params['samplesheet'],contrast['source'],contrast['target'],
 'expression.tsv','settings.json',font_file=params['heatmap_arial_font'])
PY
apptainer exec --cleanenv --no-home --bind "$TASK_ROOT" "$TASK_ROOT/containers/r_plot_v01.sif" \
 Rscript "$SOURCE_DIR/bin/plot_tf_expression.R" expression.tsv settings.json heatmap.pdf
test -s heatmap.pdf
