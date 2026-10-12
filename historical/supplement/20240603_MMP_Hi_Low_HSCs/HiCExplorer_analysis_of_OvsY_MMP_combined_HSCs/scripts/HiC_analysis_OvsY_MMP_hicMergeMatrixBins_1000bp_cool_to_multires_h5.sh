#!/bin/bash
#
#SBATCH --job-name=mergeBins
#SBATCH --output=/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/hicMergeMatrixBins_mergeBins/hicMergeMatrixBins_mergeBins_%A_%a.out
#SBATCH --error=/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/hicMergeMatrixBins_mergeBins/hicMergeMatrixBins_mergeBins_%A_%a.err
#SBATCH --array=0-8
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00

# ------------------------------------------------------------------------
# Load conda and activate HiCExplorer environment
module load miniconda3
conda activate hicexplorer

# ------------------------------------------------------------------------
# Directory containing uncorrected .h5 files at 1 kb resolution
INPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Pairs_files/cool_files_1000bp/OandY_LTHSC_MMP_combined_1000bp_Cool_files/merged_replicates_1000bp_cool_files/merged_replicates_1000bp_cool_files_BL_removed"

# Directory for final .h5 output (merged bins)
OUTPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Pairs_files/cool_files_1000bp/OandY_LTHSC_MMP_combined_1000bp_Cool_files/merged_replicates_1000bp_cool_files/merged_replicates_1000bp_cool_files_BL_removed/merged_replicates_h5_files"

# Input .h5 files (already converted in Script 1)
FILES=("O_LTHSC_1000_MMP_combined_replicates_merged_BL_removed.h5"
       "Y_LTHSC_1000_MMP_combined_replicates_merged_BL_removed.h5")

# Desired resolutions in bp
resolutions=(10000 20000 30000 50000 100000 200000 500000 1000000 2000000)

# Which resolution are we doing in this array job?
res=${resolutions[$SLURM_ARRAY_TASK_ID]}
numBins=$((res / 1000))  # e.g. 10 kb = 10 bins at 1 kb

echo "Job array index: $SLURM_ARRAY_TASK_ID"
echo "Merging bins to resolution = $res bp (numBins = $numBins)"

# Make sure output directory exists
mkdir -p "${OUTPUT_DIR}"

# ------------------------------------------------------------------------
# Merge bins for each base .h5 file
for h5_1kb in "${FILES[@]}"; do
    baseName="${h5_1kb%.h5}"   # e.g. O_LTHSC_1000_MMP_combined_replicates_merged_BL_removed

    # Extract condition & cell type from baseName if needed
    CONDITION="$(echo "${baseName}" | cut -d'_' -f1)"   # e.g. "O" or "Y"
    CELLTYPE="$(echo "${baseName}"   | cut -d'_' -f2)"  # e.g. "LTHSC"

    # Construct output file name
    outputName="${CONDITION}_${CELLTYPE}_${res}_MMP_combined_replicates_merged.h5"
    outPath="${OUTPUT_DIR}/${outputName}"

    echo "---------------------------------------------"
    echo "Input H5:   ${INPUT_DIR}/${h5_1kb}"
    echo "Output H5:  ${outPath}"
    echo "Resolution: ${res} (numBins = ${numBins})"

    hicMergeMatrixBins \
        -m "${INPUT_DIR}/${h5_1kb}" \
        -o "${outPath%.h5}" \
        -nb "${numBins}"

    echo "Done merging for ${h5_1kb} at resolution ${res}."
done

conda deactivate

echo "All done merging!"

# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiC_analysis_OvsY_MMP_hicMergeMatrixBins_1000bp_cool_to_multires_h5.sh