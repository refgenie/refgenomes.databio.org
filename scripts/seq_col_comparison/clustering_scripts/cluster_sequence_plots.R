library(pheatmap)
library(dplyr)
library(viridis)
library(stringr)
library(RColorBrewer)

# Set the directory where the CSV file is located
#results_dir <- "/home/drc/Downloads/refgenomes_pics_test/17june2025/seq_presence/test/"
#results_dir <- "/home/drc/Downloads/refgenomes_pics_test/17june2025/MOUSE/seq_presence/"
#results_dir <- "/home/drc/Downloads/refgenomes_pics_test/09jul2025/seq_presence_testing/"
results_dir <- "/home/drc/Downloads/refgenomes_pics_test/09jul2025/mouse_seq_presence/"

# Define the groups you want to iterate through
#target_groups <- c("hg38", "hg19", "hg18") # Updated target groups
target_groups <- c("m38", "m39", "m37", "mm9", "mm10")
# target_groups <- c("m37")

csv_file <- file.path(results_dir, "sequence_presence_matrix.csv")

# YOU MUST MAKE THIS ANNOTATION FILE MANUALLY, UPDATE IT IF YOU CHANGE THE INPUT SAMPLES.
annotation_file <- file.path(results_dir, "sample_annotation_groups.csv") # Updated annotation file name

# Read the main presence matrix once
presence_df_all <- read.csv(csv_file, stringsAsFactors = FALSE, row.names = 1)
# Explicitly replace hyphens with dots in column names of presence_df_all
# This ensures consistency with how sample_name from annotation_df is handled later
# and resolves the "No matching samples" warning.
colnames(presence_df_all) <- str_replace_all(colnames(presence_df_all), "-", ".")


# Read the annotation file once
annotation_df_all <- read.csv(annotation_file, stringsAsFactors = FALSE)

# Check for essential columns in the annotation file
if (!"sample_name" %in% colnames(annotation_df_all)) {
  stop("Error: 'sample_name' column is missing in the annotation file.")
}
if (any(is.na(annotation_df_all$sample_name))) {
  stop("Error: Missing values in sample_name column in annotation file")
}
if (!"group" %in% colnames(annotation_df_all)) {
  stop("Error: 'group' column is missing from the annotation file.")
}
if (!"authority" %in% colnames(annotation_df_all)) {
  stop("Error: 'authority' column is missing from the annotation file (used for heatmap annotation).")
}

# --- Define the standard color scheme for 0s and 1s using Viridis ---
# Generate a Viridis palette (e.g., 256 colors for fine granularity)
viridis_palette_full <- viridis(256, option = "D") # "D" is the default viridis palette

# Pick a purple/dark color for 0 (from the beginning of the palette)
color_for_zero <- viridis_palette_full[1] # The first color is a dark purple/blue

# Pick a yellow/bright color for 1 (from the end of the palette)
color_for_one <- viridis_palette_full[256] # The last color is bright yellow

binary_colors <- c(color_for_zero, color_for_one)
binary_breaks <- c(-0.5, 0.5, 1.5) # These breaks define the ranges for 0 and 1

# This is height in inches per row for pheatmap
base_row_height_inches <- 0.15 ## Roughly 0.15 to 0.25 inches per row is often good
# This is height in pixels per row for png output
base_row_height_pixels <- 50 

# Loop through each target group to generate separate heatmaps
for (current_group_name in target_groups) {
  message(paste0("Attempting to generate heatmap for group: ", current_group_name))

  # Filter annotation data for the current group
  annotation_df_filtered <- annotation_df_all %>%
    filter(group == current_group_name)

  # Get sample names for the current group
  # Replace hyphens with dots to match column names from presence_matrix (which are now dot-formatted)
  samples_in_current_group <- str_replace_all(annotation_df_filtered$sample_name, "-", ".")

  # Filter presence matrix to include only samples from the current group
  # Ensure only columns that exist in both are kept
  valid_samples_for_matrix <- intersect(samples_in_current_group, colnames(presence_df_all))
  if (length(valid_samples_for_matrix) == 0) {
    warning(paste0("No matching samples found in presence matrix for group: ", current_group_name, ". Skipping heatmap generation for this group."))
    next # Skip to the next group if no samples are found
  }

  presence_df_current_group <- presence_df_all[, valid_samples_for_matrix, drop = FALSE]

  # Convert to matrix
  presence_matrix_current_group <- as.matrix(presence_df_current_group)

  # Check if the matrix is empty after filtering samples
  if (nrow(presence_matrix_current_group) == 0 || ncol(presence_matrix_current_group) == 0) {
    warning(paste0("Filtered matrix for group '", current_group_name, "' is empty after sample filtering. Skipping heatmap generation."))
    next
  }

  # NEW CHECK 1: Skip if all values in the current group's matrix are zero (no sequence presence for any sequence in any sample)
  if (sum(presence_matrix_current_group) == 0) {
    message(paste0("Skipping heatmap for group '", current_group_name, "' as all samples show no sequence presence for any sequence."))
    next
  }

  # NEW CHECK 2: Filter out sequences (rows) that are not present in any sample within the current group
  # Calculate row sums (sum of presence for each sequence across samples in the current group)
  sequence_sums <- rowSums(presence_matrix_current_group)
  # Keep only sequences where the sum is greater than 0
  sequences_to_keep <- names(sequence_sums[sequence_sums > 0])

  if (length(sequences_to_keep) == 0) {
    message(paste0("Skipping heatmap for group '", current_group_name, "' as no sequences are present in any of its samples after initial filtering."))
    next
  }

  presence_matrix_current_group_filtered_seqs <- presence_matrix_current_group[sequences_to_keep, , drop = FALSE]

  # NEW CHECK 3: Filter out samples (columns) that become all zeros after sequence filtering
  # Calculate column sums (sum of presence for each sample across the *filtered* sequences)
  sample_sums <- colSums(presence_matrix_current_group_filtered_seqs)
  # Keep only samples where the sum is greater than 0
  samples_to_keep <- names(sample_sums[sample_sums > 0])

  if (length(samples_to_keep) == 0) {
    message(paste0("Skipping heatmap for group '", current_group_name, "' as all samples are now empty after sequence filtering. No relevant samples left to plot."))
    next
  }

  # This matrix now contains only sequences present and samples that contain at least one of those sequences
  presence_matrix_current_group_filtered_both <- presence_matrix_current_group_filtered_seqs[, samples_to_keep, drop = FALSE]


  # Calculate the number of sequences per genome for the current group (using the doubly filtered matrix)
  genome_sequence_counts_current_group <- colSums(presence_matrix_current_group_filtered_both)
  sorted_genomes_by_count_current_group <- names(sort(genome_sequence_counts_current_group, decreasing = TRUE))

  # Order the columns of the presence matrix by the genome sequence counts
  presence_matrix_ordered_current_group <- presence_matrix_current_group_filtered_both[, sorted_genomes_by_count_current_group, drop = FALSE]

  # --- Annotation for the genomes (based on 'authority' column) ---
  # Create annotation_col data frame for the current group, matching only the *kept* samples
  annotation_col_current_group <- data.frame(
    Group = annotation_df_filtered$authority[match(colnames(presence_matrix_ordered_current_group), str_replace_all(annotation_df_filtered$sample_name, "-", "."))],
    row.names = colnames(presence_matrix_ordered_current_group)
  )

  # Handle missing annotations specifically for the current group's subset (should be less likely now)
  if (any(is.na(annotation_col_current_group$Group))) {
    warning(paste0("Some samples in the presence matrix for group '", current_group_name, "' are missing 'authority' annotation after final filtering."))
    missing_samples_current <- rownames(annotation_col_current_group)[is.na(annotation_col_current_group$Group)]
    message("Samples missing authority annotation:")
    print(missing_samples_current)
    annotation_col_current_group <- annotation_col_current_group[!is.na(annotation_col_current_group$Group), , drop = FALSE]
    presence_matrix_ordered_current_group <- presence_matrix_ordered_current_group[, rownames(annotation_col_current_group), drop = FALSE]
  }

  # If after all cleaning, there are no samples or sequences left to plot, skip
  if (ncol(presence_matrix_ordered_current_group) == 0 || nrow(presence_matrix_ordered_current_group) == 0) {
    warning(paste0("No valid samples or sequences left for group '", current_group_name, "' after all filtering steps. Skipping heatmap generation."))
    next
  }

  # Generate colors for the 'authority' groups based on *all* unique authorities present in the filtered data
  unique_authorities <- unique(annotation_col_current_group$Group)
  group_colors_current <- brewer.pal(n = max(3, length(unique_authorities)), name = "Set3")
  if (length(unique_authorities) > length(group_colors_current)) {
    group_colors_current <- rep(group_colors_current, ceiling(length(unique_authorities) / length(group_colors_current)))[1:length(unique_authorities)]
  }
  group_colors_current <- setNames(group_colors_current, unique_authorities)

  # Ensure annotation_colors list is correctly named for 'Group'
  annotation_colors_list <- list(Group = group_colors_current)


  # --- Calculate dynamic height based on number of samples (rows after transpose) ---
  # After transposition, the number of rows is the number of samples
  num_samples_to_plot <- ncol(presence_matrix_ordered_current_group) # Number of columns in original matrix is rows after transpose
  num_sequences_to_plot <- nrow(presence_matrix_ordered_current_group) # Number of rows in original matrix is cols after transpose

  # Calculate dynamic height for pheatmap (in inches)
  # Add some extra height for title, labels, and color key
  dynamic_height_inches <- num_samples_to_plot * base_row_height_inches + 2 # Adding 2 inches for margins/title/etc.
  dynamic_height_inches <- max(dynamic_height_inches, 2) # Ensure a minimum height if very few samples

  # Calculate dynamic height for png output (in pixels)
  # This should generally align with pheatmap's calculated height * res, plus header/footer space.
  # A simple way to get matching scaling is to determine a reasonable width/height ratio for your cells.
  # If you have fixed `res`, then `width_pixels = width_inches * res` and `height_pixels = height_inches * res`.
  # Let's define a fixed width for `pheatmap` and `png` for consistency, as sequences (columns) can vary widely.
  fixed_width_inches <- 20 # Keep width constant if number of sequences (columns) can vary widely

  # Calculate width for png (pixels)
  # A fixed width ensures that columns (sequences) have consistent width across heatmaps.
  # dynamic_width_pixels <- fixed_width_inches * 150 # (24 inches * 150 dpi)

  # # Calculate height for png (pixels)
  # # This needs to be calculated in the same way as for pheatmap to ensure consistent cell sizes
  # # Add padding for title, axis labels, legend, etc. (roughly 200-300 pixels based on 150 dpi)
  # dynamic_height_pixels <- num_samples_to_plot * base_row_height_pixels + 50 # Adding 300 pixels for top/bottom margins etc.
  # dynamic_height_pixels <- max(dynamic_height_pixels, 300) # Ensure a minimum pixel height

  title <- paste0("Sequence Presence in Reference Genomes (Group: ", current_group_name, "), n=", nrow(presence_matrix_ordered_current_group))

  # --- Create the heatmap ---
  heatmap_plot <- pheatmap(
    mat = t(presence_matrix_ordered_current_group), # Transpose the matrix
    #color = viridis(100, option = "D"),
    color = binary_colors, # Use the Viridis-derived binary colors
    breaks = binary_breaks, # Use the binary breaks
    border_color = NA,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = TRUE, # Keep row names as they are now columns (sample names)
    show_colnames = FALSE, # Do not show column names (originally row names, i.e., sequences)
    annotation_row = annotation_col_current_group, # Annotation is now on the rows (samples)
    annotation_colors = annotation_colors_list,
    main = title,
    fontsize_row = 12,
    fontsize_col = 2,
    width = fixed_width_inches,
    height = dynamic_height_inches,
    silent = TRUE
  )

  # Save as PNG
  # output_png_file <- file.path(results_dir, paste0("sequence_presence_heatmap_r_WITH_GROUPS_", tolower(current_group_name), "_grouped_authority.png"))
  # # png(output_png_file, width = 2400, height = 1000, res = 150)
  # png(output_png_file,
  #     width = dynamic_width_pixels, # Use dynamic width in pixels
  #     height = dynamic_height_pixels, # Use dynamic height in pixels
  #     res = 150
  # )
  # print(heatmap_plot)
  # dev.off()
  # message(paste0("Heatmap for group '", current_group_name, "' exported to: ", output_png_file))

  # # Save as SVG (optional)
  output_svg_file <- file.path(results_dir, paste0("sequence_presence_heatmap_r_", tolower(current_group_name), "_grouped_authority.svg"))
  svg(output_svg_file, width = fixed_width_inches, height = dynamic_height_inches)
  print(heatmap_plot)
  dev.off()
  message(paste0("SVG for group '", current_group_name, "' exported to: ", output_svg_file))
}
