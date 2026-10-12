#!/bin/bash

# --- Configuration ---
# Directory where the FastQC .zip reports are located
fastqc_report_dir="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/fastqc_reports/part2_pre_trim"
# Output file to save unique overrepresented sequences
output_adapter_file="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/potential_adapters_part2.txt"
# --- End Configuration ---

echo "Searching for FastQC reports in: $fastqc_report_dir"
echo "Saving potential adapter sequences to: $output_adapter_file"

# Ensure the report directory exists
if [ ! -d "$fastqc_report_dir" ]; then
    echo "Error: FastQC report directory not found: $fastqc_report_dir"
    exit 1
fi

# Temporary file for collecting sequences
tmp_file="${output_adapter_file}.tmp"

# Clear the temporary file if it exists from a previous run
> "$tmp_file"

echo "Extracting overrepresented sequences..."

# Use find with -print0 and while read -d for safe filename handling
find "$fastqc_report_dir" -name '*_fastqc.zip' -print0 | while IFS= read -r -d $'\0' zip_file; do
    # Check if the zip file is readable before processing
    if [ -r "$zip_file" ]; then
        echo "Processing: $zip_file"
        # Use unzip -p to pipe content, extract only fastqc_data.txt
        # Use awk to parse the 'Overrepresented sequences' module
        # Corrected patterns: /^>>Overrepresented sequences/ and /^>>END_MODULE/
        # $1 is the sequence column in fastqc_data.txt for this module
        unzip -p "$zip_file" "*/fastqc_data.txt" 2>/dev/null | \
        awk '
            BEGIN { in_section=0 }
            /^>>Overrepresented sequences/ { in_section=1; next }
            /^>>END_MODULE/ { in_section=0; next } # Stop processing when END_MODULE is hit
            in_section==1 && /^[^#]/ { print $1 }  # Print first field if in section and not header
        ' >> "$tmp_file" || echo "Warning: Failed to process or extract from $zip_file"
    else
      echo "Warning: Cannot read file $zip_file"
    fi
done

echo "Sorting and finding unique sequences..."

# Check if the temporary file has content before sorting
if [ -s "$tmp_file" ]; then # -s checks if file exists and has size > 0
    sort -u "$tmp_file" > "$output_adapter_file"
    rm "$tmp_file"
    found_count=$(wc -l < "$output_adapter_file")
    echo "Found $found_count unique potential adapter sequences."
    if [ "$found_count" -eq 0 ]; then
         echo "Warning: Sorting resulted in 0 unique sequences, check temporary file content if possible before deletion."
    fi
else
    echo "No overrepresented sequences found or temporary file is empty."
    # Optionally remove the empty final file if the temp file was empty/missing
    if [ -f "$output_adapter_file" ]; then
      rm "$output_adapter_file"
    fi
    # Keep the empty temp file for inspection if it exists but has 0 size
    if [ -f "$tmp_file" ] && [ ! -s "$tmp_file" ]; then
        echo "Temporary file $tmp_file exists but is empty."
    elif [ ! -f "$tmp_file" ]; then
         echo "Temporary file $tmp_file was not created."
    fi
fi

echo "Done. Check list at $output_adapter_file (if sequences were found)."

# Example of how to view the file:
# cat /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/potential_adapters_part2.txt