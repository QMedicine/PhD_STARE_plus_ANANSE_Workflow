OvsY_MMP_combined_DEG_HiC_DI_overlaps <- read.table("/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/OvsY_MMP_combined_DEG_HiC_DI_overlaps.tsv", header=TRUE, sep="\t", stringsAsFactors = FALSE)

# Load necessary libraries
library(dplyr)
library(tidyr)

# Merge duplicated rows
OvsY_MMP_combined_DEG_HiC_DI_overlaps_unique_genes <- OvsY_MMP_combined_DEG_HiC_DI_overlaps %>%
  group_by(gene_symbol, gene_chrom, gene_start, gene_end) %>%
  summarize(
    log2FoldChange = mean(log2FoldChange),
    DI_chrom = paste(unique(DI_chrom), collapse = ","),
    DI_start = paste(unique(DI_start), collapse = ","),
    DI_end = paste(unique(DI_end), collapse = ",")
  ) %>%
  ungroup()

# Display the modified data frame
View (OvsY_MMP_combined_DEG_HiC_DI_overlaps_unique_genes)

GSE162607_A_HSC_vs_Y_HSC_itokawa <- read.table("/vast/projects/Immune_Aging/Itokawa_et_al_2022_analysis/GSE162607_A_HSC-vs-Y-HSC_itokawa.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)

# Ensure the dplyr package is loaded
library(dplyr)

# Extract the list of unique gene symbols from merged_data
gene_list <- OvsY_MMP_combined_DEG_HiC_DI_overlaps_unique_genes$gene_symbol

# Filter the rows in GSE162607_A_HSC_vs_Y_HSC_itokawa based on the gene list
DI_DE_overlap_genes_with_stats <- GSE162607_A_HSC_vs_Y_HSC_itokawa %>%
  filter(ID %in% gene_list)

# Add the prefix "DE_" to all column names except "ID"
DI_DE_overlap_genes_with_stats <- DI_DE_overlap_genes_with_stats %>%
  rename_with(~ paste0("DE_", .), -ID)

# Add a new column DE_gene_symbol
DI_DE_overlap_genes_with_stats <- DI_DE_overlap_genes_with_stats %>%
  mutate(DE_gene_symbol = ID)

# Display the modified data frame to verify the column names
print(head(DI_DE_overlap_genes_with_stats))


OvY_MMP_combined_DI_infor <- readRDS("/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/DI_analysis/DI_results/20240624_OvY_MMP_combined_DI_infor_100kbp.rds")
head

# Split the gene symbols and check for overlaps
extract_rows_with_genes <- function(df, gene_list) {
  df %>%
    filter(
      sapply(strsplit(as.character(anchor1.anno), ","), function(x) any(x %in% gene_list)) |
        sapply(strsplit(as.character(anchor2.anno), ","), function(x) any(x %in% gene_list))
    )
}

# Apply the function to extract relevant rows
DI_overlap_genes <- extract_rows_with_genes(OvY_MMP_combined_DI_infor$OvsY, gene_list)

# Display the filtered data frame
print(head(DI_overlap_genes))

# Perform separate joins for anchor1.anno and anchor2.anno
combined_data_anchor1 <- DI_overlap_genes %>%
  separate_rows(anchor1.anno, sep = ",") %>%
  left_join(DI_DE_overlap_genes_with_stats, by = c("anchor1.anno" = "DE_gene_symbol"))

combined_data_anchor2 <- DI_overlap_genes %>%
  separate_rows(anchor2.anno, sep = ",") %>%
  left_join(DI_DE_overlap_genes_with_stats, by = c("anchor2.anno" = "DE_gene_symbol"))

# Combine the results of the joins
combined_data <- bind_rows(combined_data_anchor1, combined_data_anchor2) %>%
  distinct() %>% # Remove duplicate rows
  filter(!is.na(DE_log2FoldChange)) %>% # Remove rows with NA in DE columns
  arrange(anchor1.anno, anchor2.anno) # Optional: Arrange for readability


# Display the combined data frame
print(head(combined_data))
View((head(combined_data)))
write.table(combined_data, file="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/DI_analysis/DI_results/20240624_OvY_MMP_combined_DI_DE_overlap_with_stats.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# import entire Itokawa gene expression list
GSE162607_RNA_seq_gene_count_normalized_by_DESeq2 <- read.table("/vast/projects/Immune_Aging/Itokawa_et_al_2022_analysis/GSE162607_RNA_seq_gene_count_normalized_by_DESeq2.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)

# Extract the rows that contain the genes in gene_list
filtered_gene_expression <- GSE162607_RNA_seq_gene_count_normalized_by_DESeq2 %>%
  filter(GeneName %in% gene_list)

# Display the filtered data frame
print(head(filtered_gene_expression))

write.table(filtered_gene_expression, file="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/DI_analysis/DI_results/GSE162607_RNA_seq_gene_count_normalized_by_DESeq2_DI_DE_overlap.tsv", sep="\t", quote=FALSE, row.names=FALSE)

library(biomaRt)

# Connect to the Ensembl BioMart database
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# List all available attributes
attributes <- listAttributes(ensembl)

# Extract attributes with the word 'description' in their name
description_attributes <- attributes[grep("description", attributes$description), ]

# Display the attributes with 'description' in their name
print(description_attributes)

View(description_attributes)

# Retrieve functional information for the gene list
biomaRt_gene_annotations <- getBM(
  attributes = c("external_gene_name", "description", "entrezgene_description", "namespace_1003", "name_1006", "go_id", "go_linkage_type", "goslim_goa_description"),
  filters = "external_gene_name",
  values = gene_list,
  mart = ensembl
)

# Display the gene annotations
head(biomaRt_gene_annotations)
View(biomaRt_gene_annotations)

write.table(biomaRt_gene_annotations, file="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/DI_analysis/DI_results/OvY_MMP_combined_DI_DE_overlap_biomaRt_gene_annotations.tsv", sep="\t", quote=FALSE, row.names=FALSE)


library(AnnotationDbi)
library(org.Mm.eg.db)

# Retrieve gene descriptions
AnnotationDbi_gene_descriptions <- select(org.Mm.eg.db, keys = gene_list, columns = "GENENAME", keytype = "SYMBOL")

# Display the gene descriptions
head(AnnotationDbi_gene_descriptions)
write.table(AnnotationDbi_gene_descriptions, file="/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs/DI_analysis/DI_results/OvY_MMP_combined_DI_DE_overlap_AnnotationDbi_gene_descriptions.tsv", sep="\t", quote=FALSE, row.names=FALSE)
