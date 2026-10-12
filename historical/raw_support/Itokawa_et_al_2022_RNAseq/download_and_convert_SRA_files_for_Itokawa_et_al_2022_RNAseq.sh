#!/bin/bash

# Directory to save downloaded SRA files
download_dir="/vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/raw_data"
mkdir -p "$download_dir"

# List of SRA accession numbers from the SRR_Acc_List.txt file
sra_ids=$(cat /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/Itokawa_et_al_2022_RNAseq_SRR_Acc_List.txt)

# Loop through each SRA accession ID
for sra_id in $sra_ids
do
    echo "Downloading $sra_id..."
    prefetch --max-size 10G --resume yes --output-directory $download_dir $sra_id

    echo "Converting $sra_id to FASTQ..."
    fasterq-dump --outdir $download_dir --threads 16 $download_dir/$sra_id

    echo "$sra_id processing complete."
done

echo "All files have been downloaded and converted."



#Make the script executable using:
#chmod +x /vast/projects/Immune_Aging/Itokawa_et_al_2022_RNAseq/download_and_convert_SRA_files_for_Itokawa_et_al_2022_RNAseq.sh