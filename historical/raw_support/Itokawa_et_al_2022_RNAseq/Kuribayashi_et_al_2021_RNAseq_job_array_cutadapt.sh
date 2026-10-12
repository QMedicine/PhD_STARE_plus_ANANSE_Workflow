#!/bin/bash
#SBATCH --job-name=trim_fastqc_array             # Job name
#SBATCH --array=1-9                              # Array range for 9 samples
#SBATCH --cpus-per-task=8                        # Number of CPUs per task (adjust as needed)
#SBATCH --mem=16G                                # Memory per job (adjust as needed)
#SBATCH --time=12:00:00                          # Time limit hrs:min:sec
#SBATCH --output=/vast/projects/Immune_Aging/Kuribayashi_et_al_2021_analysis/cutadapt_logs/trim_fastqc_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/Kuribayashi_et_al_2021_analysis/cutadapt_logs/trim_fastqc_%A_%a.err

# Step 1: Load necessary modules
module purge
module load cutadapt
module load fastqc

# Step 2: Define directories
fastq_files_dir="/vast/projects/Immune_Aging/Kuribayashi_et_al_2021_analysis/RNAseq/RNAseq_raw_files"
scratch_output_dir="/vast/scratch/users/stacpoole.q/Kuribayashi_et_al_2021_analysis"
project_output_dir="/vast/projects/Immune_Aging/Kuribayashi_et_al_2021_analysis"
tsv_file="/vast/projects/Immune_Aging/Kuribayashi_et_al_2021_analysis/Kuribayashi_et_al_2021_experimental_info_RNAseq.tsv"

# Create necessary output directories
# Trimmed FASTQ files go to scratch
mkdir -p "$scratch_output_dir/trimmed"
# FastQC reports and adapter file (and logs) go to the project directory
mkdir -p "$project_output_dir/trimmed/fastqc_reports"
mkdir -p "$project_output_dir/trimmed/logs"

# Create an array of input FASTQ files (with .fastq.bz2 extension)
fastq_files=("$fastq_files_dir"/*.fastq.bz2)
if [ ${#fastq_files[@]} -lt $SLURM_ARRAY_TASK_ID ]; then
    echo "Error: Not enough FASTQ files in $fastq_files_dir"
    exit 1
fi

# Select the input file based on SLURM_ARRAY_TASK_ID (1-indexed)
file="${fastq_files[$SLURM_ARRAY_TASK_ID-1]}"
# Use the basename without the .fastq.bz2 extension as the sample name
Sample=$(basename "$file" .fastq.bz2)

# Define output file names
# Trimmed FASTQ file (large file) is stored on scratch
trimmed_fastq="$scratch_output_dir/trimmed/${Sample}_trimmed.fastq.gz"
# Cutadapt log file (small) goes to the project directory
log_file="$project_output_dir/trimmed/logs/${Sample}_cutadapt.log"
# FastQC reports are stored in the project directory
fastqc_output_dir="$project_output_dir/trimmed/fastqc_reports"
# Adapter file is stored in the project directory
adapter_file="$project_output_dir/adapters.fa"

# Step 3: Define adapter sequence to remove and create the adapter file if it doesn't exist
adapter1="AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTT"
if [ ! -f "$adapter_file" ]; then
    cat > "$adapter_file" <<EOL
>Adapter1
$adapter1
EOL
fi

# Step 4: Decompress the .bz2 file to a temporary file in scratch
temp_file="$scratch_output_dir/${Sample}.fastq"
echo "Decompressing $file to temporary file $temp_file..."
bzip2 -dc "$file" > "$temp_file"
if [ ! -s "$temp_file" ]; then
    echo "Error: Decompressed file is empty."
    exit 1
fi

# Step 5: Run Cutadapt to trim adapters and low-quality bases
echo "Running Cutadapt for sample $Sample..."
cutadapt -j 8 -q 20 --minimum-length 20 \
    -a file:"$adapter_file" \
    -o "$trimmed_fastq" \
    "$temp_file" > "$log_file" 2>&1

if [ $? -ne 0 ]; then
    echo "Error: Cutadapt failed for sample $Sample. Check log: $log_file"
    exit 1
fi
echo "Cutadapt completed for sample $Sample."

# Clean up the temporary decompressed file
rm "$temp_file"

# Step 6: Run FastQC on the trimmed reads
echo "Running FastQC on trimmed reads for sample $Sample..."
fastqc -t 8 -o "$fastqc_output_dir" "$trimmed_fastq"
if [ $? -ne 0 ]; then
    echo "Error: FastQC failed for sample $Sample."
    exit 1
fi
echo "FastQC completed for sample $Sample."

echo "Processing of sample $Sample completed."

# to run in slurm use this:
# sbatch /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/Kuribayashi_et_al_2021_RNAseq_job_array_cutadapt.sh