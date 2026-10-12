#!/bin/bash
#SBATCH --job-name=trim_fastqc_array             # Job name
#SBATCH --array=1-4                             # Array range (adjust based on your number of samples)
#SBATCH --cpus-per-task=8                        # Number of CPUs per task (adjust as needed)
#SBATCH --mem=16G                                # Memory per job (adjust as needed)
#SBATCH --time=12:00:00                          # Time limit hrs:min:sec
#SBATCH --output=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/cutadapt_logs/trim_fastqc_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/cutadapt_logs/trim_fastqc_%A_%a.err

# Step 1: Load necessary modules
module purge
module load cutadapt
module load fastqc

# Step 2: Define directories
fastq_files_dir="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/raw_data"
#scratch_output_dir="/vast/scratch/users/stacpoole.q/Itokawa_et_al_2022_RNAseq"
project_output_dir="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq"
tsv_file="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/Itokawa_et_al_2022_RNAseq_sample_information.tsv"

# Create output directories if they don't exist
mkdir -p "$project_output_dir/trimmed"
#mkdir -p "$project_output_dir/fastqc_reports"
#mkdir -p "$project_output_dir/logs"
mkdir -p "$project_output_dir/trimmed/fastqc_reports"

# Step 3: Read the header line from the TSV file
header=$(head -n 1 "$tsv_file")

# Debugging output
echo "Header line: $header"

# Step 4: Use awk to get the column indices for Read1, Read2, and Sample_ID
read1_col=$(echo "$header" | awk -F'\t' '{
    for(i=1;i<=NF;i++) {
        gsub(/^"|"$/, "", $i);  # Remove surrounding quotes
        if($i=="Read1") print i
    }
}')
read2_col=$(echo "$header" | awk -F'\t' '{
    for(i=1;i<=NF;i++) {
        gsub(/^"|"$/, "", $i);
        if($i=="Read2") print i
    }
}')
sample_id_col=$(echo "$header" | awk -F'\t' '{
    for(i=1;i<=NF;i++) {
        gsub(/^"|"$/, "", $i);
        if($i=="Sample_ID") print i
    }
}')

# Debugging outputs for column indices
echo "Read1 column: $read1_col"
echo "Read2 column: $read2_col"
echo "Sample_ID column: $sample_id_col"

# Check if column indices were found
if [[ -z "$read1_col" || -z "$sample_id_col" ]]; then
    echo "Error: Could not find Read1 or Sample_ID columns in header."
    exit 1
fi

# Step 5: Read the specific sample info for this task
sample_line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$tsv_file")

# Debugging output
echo "Sample line: $sample_line"

# Step 6: Extract fields from the sample line
read -r Sample_ID Read1 Read2 <<< $(echo "$sample_line" | awk -F'\t' \
    -v read1_col="$read1_col" -v read2_col="$read2_col" -v sample_id_col="$sample_id_col" '{
    gsub(/^"|"$/, "", $read1_col);      # Remove quotes from Read1
    gsub(/^"|"$/, "", $read2_col);      # Remove quotes from Read2
    gsub(/^"|"$/, "", $sample_id_col);  # Remove quotes from Sample_ID
    printf "%s %s %s", $sample_id_col, $read1_col, $read2_col
}')

# Debugging outputs to check variables
echo "Sample_ID: $Sample_ID"
echo "Read1: $Read1"
echo "Read2: $Read2"

# Step 7: Construct full paths to the reads
read1_file="$fastq_files_dir/$Read1"
read2_file="$fastq_files_dir/$Read2"

# Debugging outputs to check if the paths and variables are set correctly
echo "Read1 file: $read1_file"
echo "Read2 file: $read2_file"

# Check if the Read1 file exists
if [[ ! -f "$read1_file" ]]; then
    echo "Error: Read1 file '$read1_file' does not exist."
    exit 1
fi

# Step 8: Define output file names
trimmed_fastq="$project_output_dir/trimmed/${Sample_ID}_trimmed.fastq.gz"
trimmed_fastq_R1="$project_output_dir/trimmed/${Sample_ID}_R1_trimmed.fastq.gz"
trimmed_fastq_R2="$project_output_dir/trimmed/${Sample_ID}_R2_trimmed.fastq.gz"
fastqc_output_dir="$project_output_dir/trimmed/fastqc_reports"

# Step 9: Define adapter sequences to remove
# These are the overrepresented TruSeq Illumina adapters you identified
adapter1="AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTT"
adapter2="GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACAATCTCGTAT"
adapter3="AGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTTT"
adapter4="GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTAT"
adapter5="GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTATCTCGTAT"
adapter6="GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTATCTCGTAT"

# Create an adapter file
adapter_file="$project_output_dir/adapters.fa"
cat > "$adapter_file" <<EOL
>Adapter1
$adapter1
>Adapter2
$adapter2
>Adapter3
$adapter3
>Adapter4
$adapter4
>Adapter5
$adapter5
>Adapter6
$adapter6
EOL

# Step 10: Run Cutadapt and FastQC
if [[ -z "$Read2" || "$Read2" == "NA" ]]; then
    # Single-end
    echo "Processing single-end reads for sample: $Sample_ID"

    # Step 10a: Run Cutadapt to trim adapters and low-quality bases
    echo "Running Cutadapt for single-end reads..."
    cutadapt -j 8 -q 20 --minimum-length 20 \
        -a file:"$adapter_file" \
        -o "$trimmed_fastq" \
        "$read1_file" > "$project_output_dir/trimmed/${Sample_ID}_cutadapt.log" 2>&1

    # Check if Cutadapt completed successfully
    if [[ $? -ne 0 ]]; then
        echo "Error: Cutadapt failed for sample $Sample_ID"
        exit 1
    fi
    echo "Cutadapt completed for sample $Sample_ID"

    # Step 10b: Run FastQC on trimmed reads
    echo "Running FastQC on trimmed reads..."
    fastqc -t 8 -o "$fastqc_output_dir" "$trimmed_fastq"

    if [[ $? -ne 0 ]]; then
        echo "Error: FastQC failed for sample $Sample_ID"
        exit 1
    fi
    echo "FastQC completed for sample $Sample_ID"

else
    # Paired-end
    echo "Processing paired-end reads for sample: $Sample_ID"

    # Step 10a: Run Cutadapt to trim adapters and low-quality bases
    echo "Running Cutadapt for paired-end reads..."
    cutadapt -j 8 -q 20 --minimum-length 20 \
        -a file:"$adapter_file" \
        -A file:"$adapter_file" \
        -o "$trimmed_fastq_R1" \
        -p "$trimmed_fastq_R2" \
        "$read1_file" "$read2_file" > "$project_output_dir/trimmed/${Sample_ID}_cutadapt.log" 2>&1

    # Check if Cutadapt completed successfully
    if [[ $? -ne 0 ]]; then
        echo "Error: Cutadapt failed for sample $Sample_ID"
        exit 1
    fi
    echo "Cutadapt completed for sample $Sample_ID"

    # Step 10b: Run FastQC on trimmed reads
    echo "Running FastQC on trimmed reads..."
    fastqc -t 8 -o "$fastqc_output_dir" "$trimmed_fastq_R1" "$trimmed_fastq_R2"

    if [[ $? -ne 0 ]]; then
        echo "Error: FastQC failed for sample $Sample_ID"
        exit 1
    fi
    echo "FastQC completed for sample $Sample_ID"

fi

echo "Processing of $Sample_ID completed."

# to run in slurm use this:
# sbatch /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/Itokawa_et_al_2022_RNAseq_job_array_cutadapt.sh