# List of all chromosomes from the header
chromosomes = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
    "chr18", "chr19", "chrX", "chrY", "chrM", "chr1_GL456210v1_random",
    "chr1_GL456211v1_random", "chr1_GL456212v1_random", "chr1_GL456221v1_random",
    "chr1_MU069434v1_random", "chr1_GL456239v1_random", "chr4_JH584295v1_random",
    "chr5_JH584296v1_random", "chr5_JH584297v1_random", "chr5_JH584298v1_random",
    "chr5_GL456354v1_random", "chr5_JH584299v1_random", "chr7_GL456219v1_random",
    "chrUn_GL456359v1", "chrUn_GL456360v1", "chrUn_GL456366v1", "chrUn_GL456367v1",
    "chrUn_GL456368v1", "chrUn_GL456370v1", "chrUn_GL456372v1", "chrUn_GL456378v1",
    "chrUn_GL456379v1", "chrUn_GL456381v1", "chrUn_GL456382v1", "chrUn_GL456383v1",
    "chrUn_GL456385v1", "chrUn_GL456387v1", "chrUn_GL456389v1", "chrUn_GL456390v1",
    "chrUn_GL456392v1", "chrUn_GL456394v1", "chrUn_GL456396v1", "chrUn_JH584304v1",
    "chrUn_MU069435v1", "chrX_GL456233v2_random", "chrY_JH584300v1_random",
    "chrY_JH584301v1_random", "chrY_JH584302v1_random", "chrY_JH584303v1_random"
]

# Define the region to keep
keep_chrom = "chr5"
keep_start = 10000000
keep_end = 10010000

# Build the condition
conditions = []
for chrom in chromosomes:
    if chrom == keep_chrom:
        conditions.append(f'(chrom1 == "{chrom}" and (pos1 < {keep_start} or pos1 > {keep_end}))')
    else:
        conditions.append(f'(chrom1 == "{chrom}")')

# Combine into a NOT condition
final_condition = "not (\n    " + " and\n    ".join(conditions) + "\n)"

# Save to a file
output_file = "/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Pairs_files/test_condition_not_chr5.txt"
with open(output_file, "w") as f:
    f.write(final_condition)

print(f"Condition saved to {output_file}")
