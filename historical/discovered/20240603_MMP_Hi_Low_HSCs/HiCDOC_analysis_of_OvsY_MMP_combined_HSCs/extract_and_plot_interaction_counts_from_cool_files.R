library(HiCExperiment)
library(ggplot2)
library(Cairo)

# Set input directory containing .cool files
input_dir <- "/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files_BL_removed/OandY_LTHSC_MMP_combined_Cool_files"

# Set output PDF file path
output_pdf <- file.path(input_dir, "Interaction_Counts_Histograms.pdf")

# List all .cool files in the input directory
cool_files <- list.files(input_dir, pattern = "\\.cool$", full.names = TRUE)

test_hic_experiment <- HiCExperiment("/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files_BL_removed/OandY_LTHSC_MMP_combined_Cool_files/Y_LTHSC_Rep1_1000000_MMP_combined_BL_removed.cool")
test_interactions_data <- interactions(test_hic_experiment)
test_interactions_data_df <- as.data.frame(test_interactions_data)
head(test_interactions_data_df)
tail(test_interactions_data_df)
View(test_interactions_data_df)
# Filter the data frame to include only rows where seqnames1 is "chrY"
test_interactions_data_df_chrY <- test_interactions_data_df[test_interactions_data_df$seqnames1 == "chrY", ]
View(test_interactions_data_df_chrY)

# Check if there are any .cool files
if (length(cool_files) == 0) {
  stop("No .cool files found in the specified directory.")
}

# Open PDF device using Cairo
CairoPDF(output_pdf, width = 10, height = 8)

# Loop through all .cool files
for (cool_file in cool_files) {
  # Print progress
  cat("Processing file:", cool_file, "\n")
  
  # Import .cool file as HiCExperiment object
  hic_experiment <- import(cool_file, format = "cool")
  
  # Extract interactions and convert to dataframe
  interactions_data <- interactions(hic_experiment)
  interactions_df <- as.data.frame(interactions_data)
  
  # Ensure 'count' column exists and is numeric
  if (!"count" %in% colnames(interactions_df)) {
    warning(paste("No 'count' column found in file:", cool_file, "- skipping."))
    next
  }
  interactions_df$count <- as.numeric(interactions_df$count)
  
  # Plot histogram of interaction counts
  plot_title <- paste("Interaction Counts -", basename(cool_file))
  histogram <- ggplot(interactions_df, aes(x = count)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black") +
    labs(
      title = plot_title,
      x = "Interaction Count",
      y = "Frequency"
    ) +
    theme_minimal()
  
  # Print plot to PDF
  print(histogram)
}

# Close PDF device
dev.off()

cat("Histograms saved to:", output_pdf, "\n")


# cooler dump -t pixels -r chrY /vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files_BL_removed/OandY_LTHSC_MMP_combined_Cool_files/Y_LTHSC_Rep1_100000_MMP_combined_BL_removed.cool > /vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files_BL_removed/OandY_LTHSC_MMP_combined_Cool_files/Y_LTHSC_Rep1_100000_MMP_combined_BL_removed_chrY_interactions.tsv

test_interactions_data_cooler <- read.csv(file="/vast/scratch/users/stacpoole.q/20240603_MMP_Hi_Low_HSCs/Cool_files_BL_removed/OandY_LTHSC_MMP_combined_Cool_files/Y_LTHSC_Rep1_100000_MMP_combined_BL_removed_chrY_interactions.tsv", sep="\t")
