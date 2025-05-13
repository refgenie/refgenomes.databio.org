library(pheatmap)
library(dplyr)
library(viridis)
library(stringr) # Load the stringr package

# Set the directory where the CSV file is located (adjust if needed)
results_dir <- "/home/drc/Downloads/refgenomes_pics_test/12May2025/full_sequence_comparisons/"
csv_file <- file.path(results_dir, "sequence_presence_matrix.csv")
annotation_file <- file.path(results_dir, "sample_annotation.csv") # Add path to annotation file

# Read the CSV file
presence_df <- read.csv(csv_file, stringsAsFactors = FALSE, row.names = 1) # IMPORTANT: row.names = 1

# Convert to matrix
presence_matrix <- as.matrix(presence_df)

# Calculate the number of sequences per genome (for sorting rows of the heatmap)
genome_sequence_counts <- colSums(presence_matrix)
sorted_genomes_by_count <- names(sort(genome_sequence_counts, decreasing = TRUE))

# Order the columns of the presence matrix by the genome sequence counts
presence_matrix_ordered <- presence_matrix[, sorted_genomes_by_count]

# --- Annotation for the genomes ---
annotation_df <- read.csv(annotation_file, stringsAsFactors = FALSE) # Read annotation CSV
# Make sure "sample_name" is a column in annotation_df
if (!"sample_name" %in% colnames(annotation_df)) {
  stop("Error: 'sample_name' column is missing in the annotation file.")
}
#check for missing values
if (any(is.na(annotation_df$sample_name))){
  stop("Error: Missing values in sample_name column in annotation file")
}

# Ensure that the annotation file has a 'sample_name' column and a 'group' column
if (!"group" %in% colnames(annotation_df)) {
    stop("Error: 'group' column is missing from the annotation file.")
}

# Create a named vector for the colors
group_colors <- setNames(
  colorRampPalette(c("darkblue", "blue", "lightblue", "pink", "red", "darkred"))(length(unique(annotation_df$group))), # Adjust colors and palette
  unique(annotation_df$group)
)

# Create annotation_col data frame.  We'll use this for column annotation.
# Create a modified sample name in the annotation_df to match the presence_matrix
annotation_df$modified_sample_name <- str_replace_all(annotation_df$sample_name, "-", ".")

annotation_col <- data.frame(
  Group = annotation_df$group[match(colnames(presence_matrix_ordered), annotation_df$modified_sample_name)], # match by modified sample name
  row.names = colnames(presence_matrix_ordered)
)
#if there are NAs after the match
if (any(is.na(annotation_col$Group))){
  warning("Some samples in the presence matrix are missing from the annotation file")
  missing_samples <- colnames(presence_matrix_ordered)[is.na(annotation_col$Group)] #get missing samples
  print("Samples missing from annotation file:")
  print(missing_samples)
  annotation_col <- annotation_col[!is.na(annotation_col$Group), , drop = FALSE]
  presence_matrix_ordered <- presence_matrix_ordered[, rownames(annotation_col)]
}


# --- Create the heatmap ---
heatmap_plot <- pheatmap(
  mat = presence_matrix_ordered,
  color = colorRampPalette(c("white", "darkgreen"))(2), # Adjust colors for 0 and 1
  border_color = NA,
  cluster_rows = FALSE, # Do not cluster sequences (already sorted by frequency)
  cluster_cols = FALSE, # Do not cluster genomes (already sorted by sequence count)
  show_rownames = TRUE, # Show sequence names
  show_colnames = TRUE,
  annotation_col = annotation_col, # Use the annotation data frame
  annotation_colors = list(Group = group_colors), # Pass the named color vector
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