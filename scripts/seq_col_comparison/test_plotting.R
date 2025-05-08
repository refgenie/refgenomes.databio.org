library(ggplot2)
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

  # Convert the data to long format for ggplot2
  heatmap_long <- heatmap_data %>%
    tidyr::pivot_longer(
      cols = -1, # Keep the first column (sample name) as an identifier
      names_to = "sample2",      # Call the other columns 'sample2'
      values_to = "similarity"     # and the values 'similarity'
    ) %>%
    dplyr::rename(sample1 = 1) # Rename the first column to sample1

  # Print structure of dataframes for debugging
  print("Structure of heatmap_long before join:")
  print(str(heatmap_long))
  print("Structure of sample_authority:")
  print(str(sample_authority))

  # Join authority information
  heatmap_long <- heatmap_long %>%
    left_join(sample_authority, by = c("sample1" = "sample_name")) %>%
    left_join(sample_authority, by = c("sample2" = "sample_name"), suffix = c("_1", "_2"))

  # Get unique authorities and order samples by authority
  all_samples_with_authority <- unique(c(heatmap_long$sample1, heatmap_long$sample2)) %>%
    tibble::enframe(name = NULL, value = "sample_name") %>%
    left_join(sample_authority, by = c("sample_name")) %>%
    arrange(authority)  # Sort by authority

  ordered_samples <- all_samples_with_authority$sample_name
  ordered_authorities <- all_samples_with_authority$authority # Get the ordered authorities.

  # Convert sample names to factors, ensuring the order is preserved.
  heatmap_long$sample1 <- factor(heatmap_long$sample1, levels = ordered_samples)
  heatmap_long$sample2 <- factor(heatmap_long$sample2, levels = ordered_samples)
  
  # Create the heatmap with ggplot2
  heatmap_plot <- ggplot(heatmap_long, aes(x = sample1, y = sample2, fill = similarity)) +
    geom_tile() +
    scale_fill_viridis_c(na.value = "grey", limits = c(0, 1)) + # Ensure 0-1 range
    geom_text(aes(label = sprintf("%.2f", similarity)), size = 2.5, na.rm = TRUE, color="white") + # Add similarity scores as text, increased size.  Reduced size to 3
    labs(
      title = paste("Heatmap of", stat_name),
      x = "Sample 1",
      y = "Sample 2",
      fill = "Similarity"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 10), # Rotate x-axis labels
      axis.text.y = element_text(size = 10),
      legend.title = element_text(size = 10),
      plot.title = element_text(size = 12, hjust = 0.5), # Center the title
      panel.grid.major = element_blank(),  # Remove grid lines
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank() # remove axis ticks
    ) +
    coord_equal() # Ensure the heatmap cells are square
  
  # Add the authority labels as annotations
  # Create a data frame for the annotations
  annotation_data_x <- data.frame(
    sample1 = ordered_samples,
    authority_x = ordered_authorities
  )
  annotation_data_y <- data.frame(
    sample2 = ordered_samples,
    authority_y = ordered_authorities
  )

#   # Add the annotations to the plot - added size parameter
#   heatmap_plot <- heatmap_plot +
#     geom_text(
#       data = annotation_data_x,
#       aes(x = sample1, y = 0, label = authority_x),
#       size = 2,  # Added size parameter here, reduced to 2
#       angle = 0,
#       hjust = 0.5,
#       vjust = 0,
#       position = position_nudge(y = -0.5)
#     ) +
#     geom_text(
#       data = annotation_data_y,
#       aes(x = 0, y = sample2, label = authority_y),
#       size = 2, # Added size parameter here, reduced to 2
#       angle = 0,
#       hjust = 0,
#       vjust = 0.5,
#       position = position_nudge(x = -0.5)
#     )

  # Print the heatmap
  print(heatmap_plot)
  
  # Save the heatmap (optional)
  tryCatch({
    ggsave(
      filename = file.path(results_dir, paste0("heatmap_", stat_name, ".png")),
      plot = heatmap_plot,
      width = 18,  # Adjust as needed
      height = 18, # Adjust as needed
      dpi = 300
    )
    message(paste("Successfully saved plot:", stat_name)) # Added message
  }, error = function(e) {
    warning(paste("Error saving plot:", stat_name, "\n", e$message))
  })
  
}
