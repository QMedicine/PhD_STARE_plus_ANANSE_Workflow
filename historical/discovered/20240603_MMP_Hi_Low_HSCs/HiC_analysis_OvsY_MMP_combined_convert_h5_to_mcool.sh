#!/bin/bash

# Load necessary tools
module load python/3.11.8  # Load Python if using an HPC
module load samtools

# Check if pairtools is installed, if not, install it
if ! command -v pairtools &> /dev/null; then
    echo "pairtools not found. Installing pairtools..."
    pip install pairtools || { echo "Failed to install pairtools. Exiting."; exit 1; }
else
    echo "pairtools is already installed."
fi

# Check if cooler is installed, if not, install it
if ! command -v cooler &> /dev/null; then
    echo "cooler not found. Installing cooler..."
    pip install cooler || { echo "Failed to install cooler. Exiting."; exit 1; }
else
    echo "cooler is already installed."
fi

# Input and output directories
input_dir="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/CO_sorted_bam_files/merged_high_low_bam_files"
temp_dir="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/temp_queryname_sorted_bam_files"
output_dir="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/converted_mcool_files"
chrom_sizes_file="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/mm39.chrom.sizes"

# Create necessary directories
mkdir -p "$temp_dir"
mkdir -p "$output_dir"

# Number of processors
nproc=8

# Loop over all .bam files in the input directory
for bam_file in "$input_dir"/*.bam; do
    # Extract base name
    base_name=$(basename "$bam_file" .bam)
    queryname_sorted_bam="$temp_dir/${base_name}.queryname_sorted.bam"
    pairs_file="$output_dir/${base_name}.pairs"
    sorted_pairs_file="$output_dir/${base_name}.sorted.pairs"
    mcool_file="$output_dir/${base_name}.mcool"

    # Step 0: Sort BAM by query name with parallel processing
    echo "Sorting $bam_file by query name..."
    samtools sort -n -@ "$nproc" -o "$queryname_sorted_bam" "$bam_file" || { echo "Failed to query name sort $bam_file. Exiting."; exit 1; }

    echo "Successfully sorted $bam_file to $queryname_sorted_bam."

    # Step 1: Convert BAM to .pairs using pairtools parse
    echo "Converting $queryname_sorted_bam to .pairs file..."
    pairtools parse \
        --chroms-path "$chrom_sizes_file" \
        --assembly mm39 \
        --min-mapq 30 \
        --output "$pairs_file" \
        "$queryname_sorted_bam" || { echo "Failed at pairtools parse for $queryname_sorted_bam. Exiting."; exit 1; }

    echo "Successfully converted $queryname_sorted_bam to $pairs_file."

    # Step 2: Sort the .pairs file using pairtools sort
    echo "Sorting $pairs_file..."
    pairtools sort \
        --output "$sorted_pairs_file" \
        --nproc "$nproc" \
        "$pairs_file" || { echo "Failed at pairtools sort for $pairs_file. Exiting."; exit 1; }

    echo "Successfully sorted $pairs_file to $sorted_pairs_file."

    # Step 3: Index the sorted .pairs file using pairix
    echo "Indexing $sorted_pairs_file..."
    pairix "$sorted_pairs_file" || { echo "Failed at pairix indexing for $sorted_pairs_file. Exiting."; exit 1; }

    echo "Successfully indexed $sorted_pairs_file."

    # Step 4: Create .cool file using cooler cload
    echo "Creating .cool file from $sorted_pairs_file..."
    cooler cload pairix \
        --assembly mm39 \
        --binsize 10000 \
        "$sorted_pairs_file" \
        "$mcool_file" || { echo "Failed at cooler cload for $sorted_pairs_file. Exiting."; exit 1; }

    echo "Successfully created .cool file: $mcool_file."

    # Step 5: Add multi-resolution levels to .mcool file using cooler zoomify
    echo "Adding multi-resolution levels to $mcool_file..."
    cooler zoomify \
        --nproc "$nproc" \
        --resolutions 10000,25000,50000,100000 \
        "$mcool_file" || { echo "Failed at cooler zoomify for $mcool_file. Exiting."; exit 1; }

    echo "Successfully created multi-resolution .mcool file: $mcool_file."
done

echo "All files processed successfully."


# chmod +x /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiC_analysis_OvsY_MMP_combined_convert_h5_to_mcool.sh
# /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/HiC_analysis_OvsY_MMP_combined_convert_h5_to_mcool.sh
