library(pheatmap)
library(dplyr)

# Set the directory where the CSV files are located
results_dir <- "/home/drc/Downloads/refgenomes_pics_test/08May2025/heatmap_csvs/" # Replace with the actual path
if (!dir.exists(results_dir)) {
  stop(paste("Error: Directory does not exist:", results_dir))
}

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
  
  message(paste("Processing file:", csv_file, "for stat:", stat_name)) # Added message

  # Read the CSV file into a data frame
  heatmap_data <- tryCatch({
    # Use base R read.csv
    read.csv(csv_file, stringsAsFactors = FALSE) # add stringsasfactors
  }, error = function(e) {
    stop(paste("Error reading CSV file:", csv_file, "\n", e$message))
  })

  # Convert the data to matrix for pheatmap
  # Remove the sample_name column before converting to matrix
  if ("sample_name" %in% colnames(heatmap_data)) {
    if (any(heatmap_data$sample_name == "") || any(is.na(heatmap_data$sample_name))){
      warning(paste("Warning: 'sample_name' column contains missing or empty values in", csv_file))
      heatmap_matrix <- as.matrix(heatmap_data[, -which(colnames(heatmap_data) == "sample_name")])
      rownames(heatmap_matrix) <- paste0("Sample_", 1:nrow(heatmap_matrix)) # Create generic row names
    }
    else{
      heatmap_matrix <- as.matrix(heatmap_data[, -which(colnames(heatmap_data) == "sample_name")])
      rownames(heatmap_matrix) <- heatmap_data$sample_name # Set sample names as row names
    }
  } else {
    warning(paste("Warning: 'sample_name' column is missing in", csv_file))
    heatmap_matrix <- as.matrix(heatmap_data)
    rownames(heatmap_matrix) <- paste0("Sample_", 1:nrow(heatmap_matrix)) # Create generic row names
  }

  # Check for non-finite values in the matrix
  if (!all(is.finite(heatmap_matrix))) {
    warning(paste("Warning: heatmap data contains non-finite values (NA, Inf) in", csv_file))
    heatmap_matrix[is.na(heatmap_matrix)] <- 0  # Replace NA with 0, or another appropriate value
    heatmap_matrix[is.infinite(heatmap_matrix)] <- 100 # Replace Inf with a large value
  }
  
  # Join authority information to heatmap_data
  heatmap_data_with_authority <- heatmap_data %>%
    left_join(sample_authority, by = "sample_name")

  # Order the data by authority
  heatmap_ordered <- heatmap_data_with_authority %>% arrange(authority)
  
  # Extract the ordered sample names and authorities
  ordered_samples <- heatmap_ordered$sample_name
  ordered_authorities <- heatmap_ordered$authority
  
  # Create a named vector for the colors
    authority_colors <- setNames(
        rainbow(length(unique(ordered_authorities))),
        unique(ordered_authorities)
    )

  # Create the heatmap with pheatmap, with ordered samples
  pheatmap(
    mat = heatmap_matrix,
    color = colorRampPalette(c("grey", "darkorange"))(100),  # Color scale
    border_color = NA,               # No borders
    cluster_rows = FALSE,            # Do not cluster rows,
    cluster_cols = FALSE,            # Do not cluster columns
    order_rows = match(ordered_samples, rownames(heatmap_matrix)), # Order rows by authority
    order_cols = match(ordered_samples, colnames(heatmap_matrix)), # Order cols by authority
    annotation_row = data.frame(Authority = ordered_authorities, row.names = ordered_samples), # Row annotation
    annotation_colors = list(Authority = authority_colors),
    show_rownames = TRUE,  # Show row names
    show_colnames = TRUE,
    main = paste("Heatmap of", stat_name)
  )

  # Save the heatmap (optional)
  tryCatch({
    dev.copy(png, file.path(results_dir, paste0("pheatmap_", stat_name, ".png")), width = 1800, height = 1800, res = 300)
    dev.off()
    message(paste("Successfully saved plot:", stat_name))
  }, error = function(e) {
    warning(paste("Error saving plot:", stat_name, "\n", e$message))
  })
}
