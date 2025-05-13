library(pheatmap)
library(dplyr)
library(viridis)

# Set the directory where the CSV file is located (adjust if needed)
results_dir <- "/home/drc/Downloads/refgenomes_pics_test/12May2025/full_sequence_comparisons/"
csv_file <- file.path(results_dir, "sequence_presence_matrix.csv")

# Read the CSV file
presence_df <- read.csv(csv_file, stringsAsFactors = FALSE, row.names = 1) # IMPORTANT: row.names = 1

# Convert to matrix
presence_matrix <- as.matrix(presence_df)

# Calculate the number of sequences per genome (for sorting rows of the heatmap)
genome_sequence_counts <- colSums(presence_matrix)
sorted_genomes_by_count <- names(sort(genome_sequence_counts, decreasing = TRUE))

# Order the columns of the presence matrix by the genome sequence counts
presence_matrix_ordered <- presence_matrix[, sorted_genomes_by_count]

# --- Annotation for the top N genomes ---
top_n <- 5 # Adjust as needed
top_genomes <- head(sorted_genomes_by_count, top_n)
annotation_col = data.frame(  # Changed to annotation_col
  Group = ifelse(colnames(presence_matrix_ordered) %in% top_genomes, paste("Top", top_n), "Other"),
  row.names = colnames(presence_matrix_ordered)
)

# Define annotation colors
annotation_colors = list(
  Group = c(`Top 5` = "blue", Other = "gray") # Adjust colors as needed
)

# --- Create the heatmap ---
heatmap_plot <- pheatmap(
  mat = presence_matrix_ordered,
  color = colorRampPalette(c("white", "darkgreen"))(2), # Adjust colors for 0 and 1
  border_color = NA,
  cluster_rows = FALSE, # Do not cluster sequences (already sorted by frequency)
  cluster_cols = FALSE, # Do not cluster genomes (already sorted by sequence count)
  show_rownames = TRUE, # Show sequence names
  show_colnames = TRUE,
  annotation_col = annotation_col, # Changed to annotation_col
  annotation_colors = annotation_colors,
  main = "Sequence Presence in Reference Genomes",
  fontsize_row = 6,
  fontsize_col = 8,
  width = 10,
  height = 8,
  silent = TRUE # Add silent = TRUE to prevent double plotting
)

# Save as PNG
png(file.path(results_dir, "sequence_presence_heatmap_r.png"), width = 1200, height = 800, res = 150)
print(heatmap_plot) # Use print() to render the plot within the file device
dev.off()

# Save as SVG (optional)
svg(file.path(results_dir, "sequence_presence_heatmap_r.svg"), width = 10, height = 8)
print(heatmap_plot) # Use print() to render the plot within the file device
dev.off()