---
title: "Untitled"
author: "Quinn Stacpoole"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{bash}
#!/bin/bash

# Directory containing TAD files
INPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files/OandY_LTHSC_MMP_Merged_Cool_files/hicConvertFormat_h5_files/hicNormalize_h5_files/hicCorrectMatrix_h5_files/hicFindTADs_TADs"
OUTPUT_DIR="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files/OandY_LTHSC_MMP_Merged_Cool_files/hicConvertFormat_h5_files/hicNormalize_h5_files/hicCorrectMatrix_h5_files/hicFindTADs_TADs/TAD_analysis_results"
mkdir -p "${OUTPUT_DIR}"

# List all boundary and domain files
BOUNDARY_FILES=($(ls ${INPUT_DIR}/*_boundaries.bed))
DOMAIN_FILES=($(ls ${INPUT_DIR}/*_domains.bed))
SCORE_FILES=($(ls ${INPUT_DIR}/*_score.bedgraph))

# Parse resolutions
RESOLUTIONS=$(ls ${BOUNDARY_FILES[@]} | awk -F'_' '{print $4}' | sort -u)

# Function to extract condition and replicate
extract_condition_replicate() {
    local file=$1
    local condition=$(echo $file | awk -F'_' '{print $1}')
    local replicate=$(echo $file | awk -F'_' '{print $3}')
    echo "${condition}_${replicate}"
}

module load bedtools/2.31.1

# Comparison for each resolution
for RES in ${RESOLUTIONS}; do
    echo "Processing resolution: ${RES}"
    
    # Separate files by condition
    O_BOUNDARIES=($(ls ${INPUT_DIR}/O_*_${RES}_*_boundaries.bed))
    Y_BOUNDARIES=($(ls ${INPUT_DIR}/Y_*_${RES}_*_boundaries.bed))
    O_DOMAINS=($(ls ${INPUT_DIR}/O_*_${RES}_*_domains.bed))
    Y_DOMAINS=($(ls ${INPUT_DIR}/Y_*_${RES}_*_domains.bed))
    O_SCORES=($(ls ${INPUT_DIR}/O_*_${RES}_*_score.bedgraph))
    Y_SCORES=($(ls ${INPUT_DIR}/Y_*_${RES}_*_score.bedgraph))

    # Create combined boundary files for conditions
    O_COMBINED="${OUTPUT_DIR}/O_${RES}_combined_boundaries.bed"
    Y_COMBINED="${OUTPUT_DIR}/Y_${RES}_combined_boundaries.bed"
    cat ${O_BOUNDARIES[@]} > "${O_COMBINED}"
    cat ${Y_BOUNDARIES[@]} > "${Y_COMBINED}"

    # Boundary overlap
    OVERLAP_OUTPUT="${OUTPUT_DIR}/boundary_overlap_${RES}.txt"
    bedtools intersect -a "${O_COMBINED}" -b "${Y_COMBINED}" > "${OUTPUT_DIR}/shared_boundaries_${RES}.bed"
    bedtools intersect -a "${O_COMBINED}" -b "${Y_COMBINED}" -v > "${OUTPUT_DIR}/O_unique_boundaries_${RES}.bed"
    bedtools intersect -a "${Y_COMBINED}" -b "${O_COMBINED}" -v > "${OUTPUT_DIR}/Y_unique_boundaries_${RES}.bed"

    # TAD size comparison
    O_SIZES="${OUTPUT_DIR}/O_${RES}_tad_sizes.txt"
    Y_SIZES="${OUTPUT_DIR}/Y_${RES}_tad_sizes.txt"
    awk '{print $3 - $2}' ${O_DOMAINS[@]} > "${O_SIZES}"
    awk '{print $3 - $2}' ${Y_DOMAINS[@]} > "${Y_SIZES}"
    python3 -c "
import numpy as np
from scipy.stats import mannwhitneyu

O_sizes = np.loadtxt('${O_SIZES}')
Y_sizes = np.loadtxt('${Y_SIZES}')
stat, pval = mannwhitneyu(O_sizes, Y_sizes)
print(f'Resolution: ${RES}, TAD size comparison: Mann-Whitney U={stat}, p-value={pval}')
" > "${OUTPUT_DIR}/tad_size_comparison_${RES}.txt"

    # TAD strength comparison
    O_STRENGTHS="${OUTPUT_DIR}/O_${RES}_tad_scores.txt"
    Y_STRENGTHS="${OUTPUT_DIR}/Y_${RES}_tad_scores.txt"
    awk '{print $4}' ${O_SCORES[@]} > "${O_STRENGTHS}"
    awk '{print $4}' ${Y_SCORES[@]} > "${Y_STRENGTHS}"
    python3 -c "
O_scores = np.loadtxt('${O_STRENGTHS}')
Y_scores = np.loadtxt('${Y_STRENGTHS}')
stat, pval = mannwhitneyu(O_scores, Y_scores)
print(f'Resolution: ${RES}, TAD strength comparison: Mann-Whitney U={stat}, p-value={pval}')
" > "${OUTPUT_DIR}/tad_strength_comparison_${RES}.txt"
done

```

