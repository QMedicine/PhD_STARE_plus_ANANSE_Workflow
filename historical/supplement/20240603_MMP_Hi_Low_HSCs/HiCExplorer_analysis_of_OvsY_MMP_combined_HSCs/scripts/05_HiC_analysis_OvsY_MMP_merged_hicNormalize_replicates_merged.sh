#!/bin/bash
#SBATCH --job-name=normalize_h5
#SBATCH --array=0-8
#SBATCH --output=/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/hicNormalize_logs_merged_reps/hicNormalize_%A_%a.out
#SBATCH --error=/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/hicNormalize_logs_merged_reps/hicNormalize_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G

# Load miniconda and activate the HiCExplorer environment
module load miniconda3
conda activate hicexplorer

# Directories
INPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files_BL_removed/OandY_LTHSC_MMP_combined_Cool_files/merged_replicates_cool_files/hicConvertFormat_h5_files"
OUTPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files_BL_removed/OandY_LTHSC_MMP_combined_Cool_files/merged_replicates_cool_files/hicConvertFormat_h5_files/hicNormalize_merged_replicates_h5_files"

# Ensure output directory exists
mkdir -p $OUTPUT_DIR

# Resolutions (index corresponds to Slurm array task ID)
RESOLUTIONS=(2000000 1000000 500000 200000 100000 50000 30000 20000 10000)
RES=${RESOLUTIONS[$SLURM_ARRAY_TASK_ID]}

# Input files for the current resolution
FILE_1="$INPUT_DIR/O_LTHSC_${RES}_MMP_combined_replicates_merged_BL_removed.h5"
FILE_2="$INPUT_DIR/Y_LTHSC_${RES}_MMP_combined_replicates_merged_BL_removed.h5"

# Output files for the current resolution
OUTPUT_1="$OUTPUT_DIR/O_LTHSC_${RES}_MMP_combined_replicates_merged_BL_removed_hicNormalize.h5"
OUTPUT_2="$OUTPUT_DIR/Y_LTHSC_${RES}_MMP_combined_replicates_merged_BL_removed_hicNormalize.h5"

# Echo input and output details
echo "Processing resolution: ${RES}"
echo "Input files: $FILE_1, $FILE_2"
echo "Output files: $OUTPUT_1, $OUTPUT_2"

# Run hicNormalize
hicNormalize --matrices $FILE_1 $FILE_2 \
             --normalize smallest \
             --outFileName $OUTPUT_1 $OUTPUT_2

# Completion message
echo "Normalization complete for resolution: ${RES}"

# Deactivate conda environment
conda deactivate

# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiC_analysis_OvsY_MMP_merged_hicNormalize_replicates_merged.sh