#!/bin/bash
#SBATCH --job-name=pairtools_process
#SBATCH --output=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/logs/02_pairtools_process_%A_%a.out
#SBATCH --error=/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/logs/02_pairtools_process_%A_%a.err
#SBATCH --array=0-11
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G

# Directories
INPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/filtered_bam_files_BL_removed"
OUTPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/filtered_bam_files_BL_removed/Pairs_files_BL_removed"
TMP_DIR="/vast/scratch/users/stacpoole.q/tmp"
CHROMSIZES="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/mm39.chrom.sizes"

# Ensure output and temp directories exist
mkdir -p $OUTPUT_DIR
mkdir -p $TMP_DIR

# Activate Conda Environment
module load miniconda3
conda activate pairtools_py310

# Get the list of BAM files and select the current one
BAM_FILES=($INPUT_DIR/*_BL_removed.bam)
BAM_FILE=${BAM_FILES[$SLURM_ARRAY_TASK_ID]}
BASENAME=$(basename "$BAM_FILE" _BL_removed.bam)

# Echo input and output details
echo "Processing file: $BAM_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Temporary directory: $TMP_DIR"

# Step 1: Parse BAM to pairs format
PAIRTOOLS_PARSE_OUTPUT="$TMP_DIR/${BASENAME}.pairs.gz"
echo "Step 1: Parsing BAM file to pairs format..."
pairtools parse -o $PAIRTOOLS_PARSE_OUTPUT \
                --min-mapq 30 \
                --walks-policy 5unique \
                --max-inter-align-gap 20 \
                --nproc-in 8 \
                --nproc-out 8 \
                --chroms-path $CHROMSIZES \
                $BAM_FILE
echo "Parsed pairs saved to: $PAIRTOOLS_PARSE_OUTPUT"

# Step 2: Sort pairs by genomic coordinates
PAIRTOOLS_SORT_OUTPUT="$TMP_DIR/${BASENAME}.pairs.sort.gz"
echo "Step 2: Sorting pairs by genomic coordinates..."
pairtools sort -o $PAIRTOOLS_SORT_OUTPUT \
               --tmpdir=$TMP_DIR \
               --nproc 8 \
               $PAIRTOOLS_PARSE_OUTPUT
echo "Sorted pairs saved to: $PAIRTOOLS_SORT_OUTPUT"

# Step 3: Deduplicate pairs
PAIRTOOLS_DEDUP_OUTPUT="$TMP_DIR/${BASENAME}.nodups.pairs.gz"
PAIRTOOLS_DEDUP_STATS="$TMP_DIR/${BASENAME}.dedup.stats"
echo "Step 3: Deduplicating pairs..."
pairtools dedup --max-mismatch 3 \
                --nproc-in 8 \
                --nproc-out 8 \
                --mark-dups \
                --output $PAIRTOOLS_DEDUP_OUTPUT \
                --output-stats $PAIRTOOLS_DEDUP_STATS \
                $PAIRTOOLS_SORT_OUTPUT
echo "Deduplicated pairs saved to: $PAIRTOOLS_DEDUP_OUTPUT"
echo "Deduplication stats saved to: $PAIRTOOLS_DEDUP_STATS"

# Step 4: Filter pairs for meaningful interactions
PAIRTOOLS_SELECT_OUTPUT="$OUTPUT_DIR/${BASENAME}_BL_removed.nodups.pairs.filt.gz"
echo "Step 4: Filtering pairs for meaningful interactions..."
pairtools select "(chrom1 != chrom2) or ((chrom1 == chrom2) and (abs(pos1 - pos2) > 1000))" \
                 $PAIRTOOLS_DEDUP_OUTPUT \
                 -o $PAIRTOOLS_SELECT_OUTPUT
echo "Filtered pairs saved to: $PAIRTOOLS_SELECT_OUTPUT"

# Cleanup intermediate files (optional)
rm -f $PAIRTOOLS_PARSE_OUTPUT $PAIRTOOLS_SORT_OUTPUT $PAIRTOOLS_DEDUP_OUTPUT
echo "Intermediate files removed. Processing complete for: $BAM_FILE"

# sbatch /vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/filtered_bam_files/pairtools_process_bam_BL_removed_to_pairs_BL_removed.sh