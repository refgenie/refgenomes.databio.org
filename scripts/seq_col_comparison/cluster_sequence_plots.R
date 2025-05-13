library(pheatmap)
library(dplyr)
library(viridis)
library(stringr)
library(RColorBrewer)

# Set the directory where the CSV file is located
results_dir <- "/home/drc/Downloads/refgenomes_pics_test/12May2025/full_sequence_comparisons/"
csv_file <- file.path(results_dir, "sequence_presence_matrix.csv")
annotation_file <- file.path(results_dir, "sample_annotation.csv")

# Read the CSV file
presence_df <- read.csv(csv_file, stringsAsFactors = FALSE, row.names = 1)

# Convert to matrix
presence_matrix <- as.matrix(presence_df)

# Calculate the number of sequences per genome
genome_sequence_counts <- colSums(presence_matrix)
sorted_genomes_by_count <- names(sort(genome_sequence_counts, decreasing = TRUE))

# Order the columns of the presence matrix by the genome sequence counts
presence_matrix_ordered <- presence_matrix[, sorted_genomes_by_count]

# --- Annotation for the genomes ---
annotation_df <- read.csv(annotation_file, stringsAsFactors = FALSE)
if (!"sample_name" %in% colnames(annotation_df)) {
  stop("Error: 'sample_name' column is missing in the annotation file.")
}
if (any(is.na(annotation_df$sample_name))) {
  stop("Error: Missing values in sample_name column in annotation file")
}
if (!"group" %in% colnames(annotation_df)) {
  stop("Error: 'group' column is missing from the annotation file.")
}

# # Create a named vector for the colors
# group_colors <- setNames(
#   c("Decoy" = "#440154", "Primary" = "#fde725", "Other" = "#cccccc"),
#   unique(annotation_df$group)
# )
# Use the Dark2 color palette from RColorBrewer
group_colors <- brewer.pal(n = 8, name = "Dark2")
# Ensure we have enough colors for the unique groups, if not, extend.
unique_groups <- unique(annotation_df$group)
if (length(unique_groups) > length(group_colors)) {
  # You can extend the colors by repeating or interpolating.  A simple repeat:
  group_colors <- rep(group_colors, ceiling(length(unique_groups) / length(group_colors)))[1:length(unique_groups)]
}
# Assign the colors to the groups
group_colors <- setNames(group_colors[1:length(unique_groups)], unique_groups)

# Create annotation_col data frame
annotation_df$modified_sample_name <- str_replace_all(annotation_df$sample_name, "-", ".")
annotation_col <- data.frame(
  Group = annotation_df$group[match(colnames(presence_matrix_ordered), annotation_df$modified_sample_name)],
  row.names = colnames(presence_matrix_ordered)
)

if (any(is.na(annotation_col$Group))) {
  warning("Some samples in the presence matrix are missing from the annotation file")
  missing_samples <- colnames(presence_matrix_ordered)[is.na(annotation_col$Group)]
  print("Samples missing from annotation file:")
  print(missing_samples)
  annotation_col <- annotation_col[!is.na(annotation_col$Group), , drop = FALSE]
  presence_matrix_ordered <- presence_matrix_ordered[, rownames(annotation_col)]
}
title = paste0(" Sequence Presence in Reference Genomes, n=", nrow(presence_matrix_ordered))
# --- Create the heatmap ---
heatmap_plot <- pheatmap(
  mat = t(presence_matrix_ordered), # Transpose the matrix
  color = viridis(100, option = "D"),
  border_color = NA,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE, # Keep row names as they are now columns
  show_colnames = FALSE, #  Do not show column names (originally row names)
  annotation_row = annotation_col, #  Annotation is now on the rows
  annotation_colors = list(Group = group_colors),
  main = title,
  fontsize_row = 8,
  fontsize_col = 2,
  width = 24,
  height = 10,
  silent = TRUE
)

# Save as PNG
png(file.path(results_dir, "sequence_presence_heatmap_r_grouped_authority.png"), width = 2400, height = 1000, res = 150)
print(heatmap_plot)
dev.off()

# Save as SVG (optional)
svg(file.path(results_dir, "sequence_presence_heatmap_r_grouped_authority.svg"), width = 24, height = 10)
print(heatmap_plot)
dev.off()



# BY TYPE

library(pheatmap)
library(dplyr)
library(viridis)
library(stringr)
library(RColorBrewer)

# Set the directory where the CSV file is located
results_dir <- "/home/drc/Downloads/refgenomes_pics_test/12May2025/full_sequence_comparisons/"
csv_file <- file.path(results_dir, "sequence_presence_matrix.csv")
annotation_file <- file.path(results_dir, "sample_annotation_2.csv")

# Read the CSV file
presence_df <- read.csv(csv_file, stringsAsFactors = FALSE, row.names = 1)

# Convert to matrix
presence_matrix <- as.matrix(presence_df)

# Calculate the number of sequences per genome
genome_sequence_counts <- colSums(presence_matrix)
sorted_genomes_by_count <- names(sort(genome_sequence_counts, decreasing = TRUE))

# Order the columns of the presence matrix by the genome sequence counts
presence_matrix_ordered <- presence_matrix[, sorted_genomes_by_count]

# --- Annotation for the genomes ---
annotation_df <- read.csv(annotation_file, stringsAsFactors = FALSE)
if (!"sample_name" %in% colnames(annotation_df)) {
  stop("Error: 'sample_name' column is missing in the annotation file.")
}
if (any(is.na(annotation_df$sample_name))) {
  stop("Error: Missing values in sample_name column in annotation file")
}
if (!"group" %in% colnames(annotation_df)) {
  stop("Error: 'group' column is missing from the annotation file.")
}

# # Create a named vector for the colors
# group_colors <- setNames(
#   c("Decoy" = "#440154", "Primary" = "#fde725", "Other" = "#cccccc"),
#   unique(annotation_df$group)
# )
# Use the Dark2 color palette from RColorBrewer
group_colors <- brewer.pal(n = 8, name = "Dark2")
# Ensure we have enough colors for the unique groups, if not, extend.
unique_groups <- unique(annotation_df$group)
if (length(unique_groups) > length(group_colors)) {
  # You can extend the colors by repeating or interpolating.  A simple repeat:
  group_colors <- rep(group_colors, ceiling(length(unique_groups) / length(group_colors)))[1:length(unique_groups)]
}
# Assign the colors to the groups
group_colors <- setNames(group_colors[1:length(unique_groups)], unique_groups)

# Create annotation_col data frame
annotation_df$modified_sample_name <- str_replace_all(annotation_df$sample_name, "-", ".")
annotation_col <- data.frame(
  Group = annotation_df$group[match(colnames(presence_matrix_ordered), annotation_df$modified_sample_name)],
  row.names = colnames(presence_matrix_ordered)
)

if (any(is.na(annotation_col$Group))) {
  warning("Some samples in the presence matrix are missing from the annotation file")
  missing_samples <- colnames(presence_matrix_ordered)[is.na(annotation_col$Group)]
  print("Samples missing from annotation file:")
  print(missing_samples)
  annotation_col <- annotation_col[!is.na(annotation_col$Group), , drop = FALSE]
  presence_matrix_ordered <- presence_matrix_ordered[, rownames(annotation_col)]
}
title = paste0(" Sequence Presence in Reference Genomes, n=", nrow(presence_matrix_ordered))
# --- Create the heatmap ---
heatmap_plot <- pheatmap(
  mat = t(presence_matrix_ordered), # Transpose the matrix
  color = viridis(100, option = "D"),
  border_color = NA,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE, # Keep row names as they are now columns
  show_colnames = FALSE, #  Do not show column names (originally row names)
  annotation_row = annotation_col, #  Annotation is now on the rows
  annotation_colors = list(Group = group_colors),
  main = title,
  fontsize_row = 8,
  fontsize_col = 2,
  width = 24,
  height = 10,
  silent = TRUE
)

# Save as PNG
png(file.path(results_dir, "sequence_presence_heatmap_r_grouped_type.png"), width = 2400, height = 1000, res = 150)
print(heatmap_plot)
dev.off()

# Save as SVG (optional)
svg(file.path(results_dir, "sequence_presence_heatmap_r_grouped_type.svg"), width = 24, height = 10)
print(heatmap_plot)
dev.off()

