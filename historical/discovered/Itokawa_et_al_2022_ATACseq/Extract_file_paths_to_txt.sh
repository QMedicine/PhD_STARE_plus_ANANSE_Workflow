#!/bin/bash

# Define directories
fastq_files_dir="/vast/scratch/users/stacpoole.q/Itokawa_et_al_2022_ATACseq/raw_data"
output_file="/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/extracted_file_paths.txt"
csv_file="/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/sample_information.csv"

# Clear the output file if it exists
> "$output_file"

# Get the column indices for Read1 and Read2
header=$(head -n 1 "$csv_file")
read1_col=$(echo "$header" | tr ',' '\n' | grep -n -w "Read1" | cut -d: -f1)
read2_col=$(echo "$header" | tr ',' '\n' | grep -n -w "Read2" | cut -d: -f1)

# Skip the header and process each line
tail -n +2 "$csv_file" | while IFS=, read -r -a columns; do
    # Extract the Read1 and Read2 file names
    Read1="${columns[$((read1_col - 1))]}"
    Read2="${columns[$((read2_col - 1))]}"
    
    # Remove any leading/trailing quotes from Read1 and Read2
    Read1=$(echo "$Read1" | sed 's/^"\|"$//g')
    Read2=$(echo "$Read2" | sed 's/^"\|"$//g')
    
    # Construct full paths
    read1_file="$fastq_files_dir/$Read1"
    
    # If Read2 is not empty or NA, add it to the file paths
    if [[ -n "$Read2" && "$Read2" != "NA" ]]; then
        read2_file="$fastq_files_dir/$Read2"
        echo "$read1_file" >> "$output_file"
        echo "$read2_file" >> "$output_file"
    else
        echo "$read1_file" >> "$output_file"
    fi
done

echo "File paths have been extracted and saved to $output_file"


# save as /vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/Extract_file_paths_to_txt.sh
# use chmod +x /vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/Extract_file_paths_to_txt.sh