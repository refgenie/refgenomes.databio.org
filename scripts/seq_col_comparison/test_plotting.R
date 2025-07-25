library(pheatmap)
library(dplyr)
library(viridis)
library(stringr)

# Set the directory where the CSV files are located

#results_dir <- "/home/drc/Downloads/refgenomes_pics_test/jun112025/ncbi_patches/"
results_dir <- "/home/drc/Downloads/refgenomes_pics_test/21jul2025/testing_determinism/"

if (!dir.exists(results_dir)) {
  stop(paste("Error: Directory does not exist:", results_dir))
}

# YOUR CSVs MUST HAVE A "sample_name" COLUMN!!!
# Get a list of CSV files in the directory that start with "heatmap_data_"
csv_files <- list.files(path = results_dir, pattern = "heatmap_data_.*\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) {
  stop(paste("Error: No CSV files found in:", results_dir))
}

# Load the sample authority mapping (assuming it's in a CSV file)
# You'll need to create this CSV file from your Python script or have it available
sample_authority_file <- file.path(results_dir, "sample_authority.csv") #  adjust the name/path if needed
if (!file.exists(sample_authority_file)) {
  stop(paste("Error: Sample authority file does not exist:", sample_authority_file))
}
sample_authority <- read.csv(sample_authority_file, stringsAsFactors = FALSE) # Use base R's read.csv, and set stringsAsFactors
# Ensure the file has columns "sample_name" and "authority"

# Loop through each CSV file, read it, and create a heatmap
for (csv_file in csv_files) {
  # Extract the statistic name from the filename
  stat_name <- gsub("heatmap_data_(.*)\\.csv", "\\1", basename(csv_file))

  message(paste("Processing file:", csv_file, "for stat:", stat_name))

  # Read the CSV file into a data frame
  heatmap_data <- tryCatch({
    read.csv(csv_file, stringsAsFactors = FALSE)
  }, error = function(e) {
    stop(paste("Error reading CSV file:", csv_file, "\n", e$message))
  })

  # Convert the data to matrix for pheatmap
  if ("sample_name" %in% colnames(heatmap_data)) {
    if (any(heatmap_data$sample_name == "") || any(is.na(heatmap_data$sample_name))){
      warning(paste("Warning: 'sample_name' column contains missing or empty values in", csv_file))
      heatmap_matrix <- as.matrix(heatmap_data[, -which(colnames(heatmap_data) == "sample_name")])
      rownames(heatmap_matrix) <- paste0("Sample_", 1:nrow(heatmap_matrix))
      colnames(heatmap_matrix) <- paste0("Sample_", 1:ncol(heatmap_matrix)) # Assuming same samples as rows
    }
    else{
      heatmap_matrix <- as.matrix(heatmap_data[, -which(colnames(heatmap_data) == "sample_name")])
      rownames(heatmap_matrix) <- heatmap_data$sample_name
      colnames(heatmap_matrix) <- heatmap_data$sample_name # Assuming same samples as rows
    }
  } else {
    warning(paste("Warning: 'sample_name' column is missing in", csv_file))
    heatmap_matrix <- as.matrix(heatmap_data)
    rownames(heatmap_matrix) <- paste0("Sample_", 1:nrow(heatmap_matrix))
    colnames(heatmap_matrix) <- paste0("Sample_", 1:ncol(heatmap_matrix)) # Assuming same samples as rows
  }

  # Join authority information to heatmap data
  heatmap_data_with_authority <- heatmap_data %>%
    left_join(sample_authority, by = "sample_name")

  # Order the data and matrix by authority
  heatmap_ordered <- heatmap_data_with_authority %>% arrange(authority)
  ordered_samples <- heatmap_ordered$sample_name
  ordered_authorities <- heatmap_ordered$authority

  # Re-index the matrix based on the ordered samples for both rows and columns
  if (all(ordered_samples %in% rownames(heatmap_matrix)) && all(ordered_samples %in% colnames(heatmap_matrix))) {
    heatmap_matrix_ordered <- heatmap_matrix[ordered_samples, ordered_samples]
  } else {
    warning(paste("Warning: Sample names in authority file do not fully match those in", csv_file))
    heatmap_matrix_ordered <- heatmap_matrix # Use original if mismatch
  }

  # Create a named vector for the colors
  authority_colors <- setNames(
    colorRampPalette(c("lightblue", "pink", "darkred"))(length(unique(ordered_authorities))),
    unique(ordered_authorities)
  )

  # Create annotation data frames
  annotation_row = data.frame(Authority = ordered_authorities, row.names = ordered_samples)
  annotation_col = data.frame(Authority = ordered_authorities, row.names = ordered_samples) # Use the same data

  # Define a custom color scale
  my_color <- colorRampPalette(c(viridis(100, option = "D")))
  my_breaks <- seq(0, 1, length.out = 101)

  cleaned_stat_name <- gsub("_", " ", stat_name)
  title_cased_stat_name <- str_to_title(cleaned_stat_name)

  # Save the heatmap (optional)
  # tryCatch({
  # pheatmap(
  #   mat = heatmap_matrix_ordered,
  #   color = my_color(100),
  #   breaks = my_breaks,
  #   border_color = NA,
  #   cluster_rows = FALSE,            # Do not cluster rows
  #   cluster_cols = FALSE,            # Do not cluster columns
  #   order_rows = match(ordered_samples, rownames(heatmap_matrix_ordered)), # Order rows by authority
  #   order_cols = match(ordered_samples, colnames(heatmap_matrix_ordered)), # Order columns by authority
  #   annotation_row = annotation_row,
  #   annotation_col = annotation_col, # Add column annotation
  #   annotation_colors = list(Authority = authority_colors),
  #    #annotation_position = c("row", "bottom"),
  #   show_rownames = TRUE,
  #   show_colnames = TRUE,
  #   main = paste("Heatmap of", title_cased_stat_name),
  #   fontsize = 6,
  #   fontsize_row = 6,
  #   fontsize_col = 6,
  #   cellheight = 20,
  #   cellwidth = 20,
  #   fontfamily = "Arial",
  # )
  #   dev.copy(png, file.path(results_dir, paste0("pheatmap_", stat_name, ".png")), width = 1200, height = 1200, res = 300)
  #   # Create the heatmap with pheatmap, with ordered samples and column annotation
  #   dev.off()
  #   message(paste("Successfully saved plot:", stat_name))
  # }, error = function(e) {
  #   warning(paste("Error saving plot:", stat_name, "\n", e$message))
  # })

  # Save the heatmap as SVG
  tryCatch({
    svg(file.path(results_dir, paste0("pheatmap_", stat_name, ".svg")), width = 9, height = 9) # Adjust width and height as needed
    pheatmap(
      mat = heatmap_matrix_ordered,
      color = my_color(100),
      breaks = my_breaks,
      border_color = NA,
      cluster_rows = FALSE,            # Do not cluster rows
      cluster_cols = FALSE,            # Do not cluster columns
      order_rows = match(ordered_samples, rownames(heatmap_matrix_ordered)), # Order rows by authority
      order_cols = match(ordered_samples, colnames(heatmap_matrix_ordered)), # Order columns by authority
      annotation_row = annotation_row,
      annotation_col = annotation_col, # Add column annotation
      annotation_colors = list(Authority = authority_colors),
       #annotation_position = c("row", "bottom"),
      show_rownames = TRUE,
      show_colnames = TRUE,
      main = paste("Heatmap of", title_cased_stat_name),
      fontsize = 6,
      fontsize_row = 6,
      fontsize_col = 6,
      cellheight = 7,
      cellwidth = 7,
      fontfamily = "Arial",
    )
    dev.off()
    message(paste("Successfully saved SVG plot:", stat_name))
  }, error = function(e) {
    warning(paste("Error saving SVG plot:", stat_name, "\n", e$message))
  })

}
