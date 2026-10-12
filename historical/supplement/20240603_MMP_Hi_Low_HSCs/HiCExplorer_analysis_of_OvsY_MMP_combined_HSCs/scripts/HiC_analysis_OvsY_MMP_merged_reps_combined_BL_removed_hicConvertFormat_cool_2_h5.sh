#!/bin/bash
#
#SBATCH --job-name=convertCoolToH5
#SBATCH --output=/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/convertCoolToH5/convertCoolToH5_%j.out
#SBATCH --error=/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/convertCoolToH5/convertCoolToH5_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00

# ------------------------------------------------------------------------
# Load conda and activate HiCExplorer environment
module load miniconda3
conda activate hicexplorer

# ------------------------------------------------------------------------
# Directory containing all .cool files
INPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Pairs_files/cool_files_1000bp/OandY_LTHSC_MMP_combined_1000bp_Cool_files/merged_replicates_1000bp_cool_files/merged_replicates_1000bp_cool_files_BL_removed"

# List all COOL files
FILES=("O_LTHSC_1000_MMP_combined_replicates_merged_BL_removed.cool"
       "Y_LTHSC_1000_MMP_combined_replicates_merged_BL_removed.cool")

echo "Converting the following COOL files to H5:"
printf '%s\n' "${FILES[@]}"
echo "-------------------------------------------"

# Loop through each .cool file
for coolfile in "${FILES[@]}"; do
    baseName="${coolfile%.cool}"
    h5_file="${INPUT_DIR}/${baseName}.h5"

    echo "Processing $coolfile"

    # If the .h5 does NOT exist yet, convert
    if [ ! -f "${h5_file}" ]; then
        echo "Converting ${coolfile} to H5 -> ${h5_file}"
        hicConvertFormat \
            --matrices "${INPUT_DIR}/${coolfile}" \
            --outFileName "${h5_file}" \
            --inputFormat cool \
            --outputFormat h5
        echo "Done converting ${coolfile}"
    else
        echo "H5 file already exists: ${h5_file}"
    fi

    echo "-------------------------------------------"
done

conda deactivate
echo "Conversion complete. Exiting."

# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiC_analysis_OvsY_MMP_merged_reps_combined_BL_removed_hicConvertFormat_cool_2_h5.sh
