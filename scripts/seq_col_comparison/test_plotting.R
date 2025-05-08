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

# Loop through each CSV file, read it, and create a heatmap
for (csv_file in csv_files) {
  # Extract the statistic name from the filename
  stat_name <- gsub("heatmap_data_(.*)\\.csv", "\\1", basename(csv_file))
  
  message(paste("Processing file:", csv_file, "for stat:", stat_name)) # Added message

  # Read the CSV file into a data frame
  heatmap_data <- tryCatch({
      # Use base R read.csv
      read.csv(csv_file)
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

  # Convert sample names to factors, ensuring the order is preserved.
  all_samples <- unique(c(heatmap_long$sample1, heatmap_long$sample2))
  heatmap_long$sample1 <- factor(heatmap_long$sample1, levels = all_samples)
  heatmap_long$sample2 <- factor(heatmap_long$sample2, levels = all_samples)
  
  # Create the heatmap with ggplot2
  heatmap_plot <- ggplot(heatmap_long, aes(x = sample1, y = sample2, fill = similarity)) +
    geom_tile() +
    scale_fill_viridis_c(na.value = "white", limits = c(0, 1)) + # Ensure 0-1 range
    geom_text(aes(label = sprintf("%.2f", similarity)), size = 2, na.rm = TRUE) + # Add similarity scores as text
    labs(
      title = paste("Heatmap of", stat_name),
      x = "Sample 1",
      y = "Sample 2",
      fill = "Similarity"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 8), # Rotate x-axis labels
      axis.text.y = element_text(size = 8),
      legend.title = element_text(size = 10),
      plot.title = element_text(size = 12, hjust = 0.5), # Center the title
      panel.grid.major = element_blank(),  # Remove grid lines
      panel.grid.minor = element_blank()
    ) +
    coord_equal() # Ensure the heatmap cells are square
  
  # Print the heatmap
  print(heatmap_plot)
  
  # Save the heatmap (optional)
  tryCatch({
    ggsave(
      filename = file.path(results_dir, paste0("heatmap_", stat_name, ".png")),
      plot = heatmap_plot,
      width = 8,  # Adjust as needed
      height = 8, # Adjust as needed
      dpi = 300
    )
    message(paste("Successfully saved plot:", stat_name)) # Added message
  }, error = function(e) {
    warning(paste("Error saving plot:", stat_name, "\n", e$message))
  })
  
}


# library(data.table)
# library(ggplot2)

# # Assuming 'dt' data.table from the previous steps exists

# set.seed(123)

# n <- 1000

# dt <- data.table(id = 1:n)
# dt[, x := rnorm(n)]
# dt[, y := 2*x + rnorm(n, sd = 2)]
# dt[, group := sample(c("A", "B", "C"), n, replace = TRUE)]

# # Histogram of the 'x' variable
# histogram_plot <- ggplot(dt, aes(x = x)) +
#   geom_histogram(binwidth = 0.5, fill = "steelblue", color = "black") +
#   labs(title = "Histogram of Random Variable X",
#        x = "Random Variable X",
#        y = "Frequency") +
#   theme_minimal()

# print(histogram_plot)

# # Boxplot of 'y' by 'group'
# boxplot_plot <- ggplot(dt, aes(x = group, y = y, fill = group)) +
#   geom_boxplot() +
#   labs(title = "Boxplot of Random Variable Y by Group",
#        x = "Group",
#        y = "Random Variable Y",
#        fill = "Group") +
#   theme_minimal()

# print(boxplot_plot)