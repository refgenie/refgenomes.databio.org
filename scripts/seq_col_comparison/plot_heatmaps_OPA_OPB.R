library(pheatmap)
library(dplyr)
library(viridis)
library(stringr)

# Set the directory where the CSV files are located
results_dir <- "/home/drc/Downloads/refgenomes_pics_test/28May2025/homo_sapiens/4_comparisons/" 

if (!dir.exists(results_dir)) {
  stop(paste("Error: Directory does not exist:", results_dir))
}

# Load the sample authority mapping (assuming it's in a CSV file)
sample_authority_file <- file.path(results_dir, "sample_authority.csv")
if (!file.exists(sample_authority_file)) {
  stop(paste("Error: Sample authority file does not exist:", sample_authority_file))
}
sample_authority <- read.csv(sample_authority_file, stringsAsFactors = FALSE)

# # THIS DOESNT WORK AS INTENDED
# # --- Special handling for OPA and OPB split triangle heatmap ---
# opa_file_base <- "heatmap_data_opa_name_len.csv"
# opb_file_base <- "heatmap_data_opb_name_len.csv"

# opa_file_path <- file.path(results_dir, opa_file_base)
# opb_file_path <- file.path(results_dir, opb_file_base)

# if (file.exists(opa_file_path) && file.exists(opb_file_path)) {
#   message("Processing combined OPA (lower triangle) and OPB (upper triangle) heatmap...")

#   # Read OPA data
#   opa_data <- tryCatch({
#     read.csv(opa_file_path, stringsAsFactors = FALSE)
#   }, error = function(e) {
#     stop(paste("Error reading OPA CSV file:", opa_file_path, "\n", e$message))
#   })

#   # Read OPB data
#   opb_data <- tryCatch({
#     read.csv(opb_file_path, stringsAsFactors = FALSE)
#   }, error = function(e) {
#     stop(paste("Error reading OPB CSV file:", opb_file_path, "\n", e$message))
#   })

#   # Ensure both dataframes have the 'sample_name' column and are aligned
#   if (!("sample_name" %in% colnames(opa_data)) || !("sample_name" %in% colnames(opb_data))) {
#     stop("Both OPA and OPB CSVs must have a 'sample_name' column.")
#   }

#   # Order both dataframes by sample_name to ensure alignment
#   # This step is crucial for consistent indexing across matrices
#   opa_data <- opa_data %>% arrange(sample_name)
#   opb_data <- opb_data %>% arrange(sample_name)

#   # Convert to matrices, ensuring sample_name is row/col names
#   opa_matrix_raw <- as.matrix(opa_data[, -which(colnames(opa_data) == "sample_name")])
#   rownames(opa_matrix_raw) <- opa_data$sample_name
#   colnames(opa_matrix_raw) <- opa_data$sample_name

#   opb_matrix_raw <- as.matrix(opb_data[, -which(colnames(opb_data) == "sample_name")])
#   rownames(opb_matrix_raw) <- opb_data$sample_name
#   colnames(opb_matrix_raw) <- opb_data$sample_name

#   # Ensure matrices have the same dimensions and sample names
#   if (!all(rownames(opa_matrix_raw) == rownames(opb_matrix_raw)) || !all(colnames(opa_matrix_raw) == colnames(opb_matrix_raw))) {
#     stop("OPA and OPB matrices do not have matching row/column names or dimensions.")
#   }

#   # Join authority information
#   opa_data_with_authority <- opa_data %>%
#     left_join(sample_authority, by = "sample_name")

#   # Order the data and matrix by authority
#   opa_ordered <- opa_data_with_authority %>% arrange(authority)
#   ordered_samples <- opa_ordered$sample_name
#   ordered_authorities <- opa_ordered$authority

#   # Re-index the raw matrices based on the ordered samples for both rows and columns
#   # This is the "canonical" order we will use.
#   opa_matrix_ordered <- opa_matrix_raw[ordered_samples, ordered_samples]
#   opb_matrix_ordered <- opb_matrix_raw[ordered_samples, ordered_samples]

#   # --- Create the combined matrix for coloring and display numbers ---
#   n_samples <- length(ordered_samples)
#   combined_color_matrix <- matrix(NA, nrow = n_samples, ncol = n_samples,
#                                   dimnames = list(ordered_samples, ordered_samples))
#   combined_display_numbers <- matrix("", nrow = n_samples, ncol = n_samples,
#                                      dimnames = list(ordered_samples, ordered_samples))

#   # Populate the combined matrices cell by cell
#   for (i in 1:n_samples) {
#     for (j in 1:n_samples) {
#       if (i >= j) { # Lower triangle and diagonal
#         combined_color_matrix[i, j] <- opa_matrix_ordered[i, j]
#         combined_display_numbers[i, j] <- sprintf("%.2f", opa_matrix_ordered[i, j])
#       } else { # Upper triangle
#         # For the upper triangle cell (i, j), we want the value from OPB's
#         # (j, i) position if OPB's data is only in its lower triangle.
#         # If OPB CSV is also lower triangle, then opb_matrix_ordered[j, i] holds the value.
#         combined_color_matrix[i, j] <- opb_matrix_ordered[j, i]
#         combined_display_numbers[i, j] <- sprintf("%.2f", opb_matrix_ordered[j, i])
#       }
#     }
#   }

#   # Create a named vector for the authority colors
#   authority_colors <- setNames(
#     colorRampPalette(c("lightblue", "pink", "darkred"))(length(unique(ordered_authorities))),
#     unique(ordered_authorities)
#   )

#   # Create annotation data frames
#   annotation_row = data.frame(Authority = ordered_authorities, row.names = ordered_samples)
#   annotation_col = data.frame(Authority = ordered_authorities, row.names = ordered_samples)

#   # Define a custom color scale
#   # Calculate min/max across both matrices to ensure the color scale covers the full range
#   min_val <- min(opa_matrix_ordered, opb_matrix_ordered, na.rm = TRUE)
#   max_val <- max(opa_matrix_ordered, opb_matrix_ordered, na.rm = TRUE)
#   my_color <- colorRampPalette(viridis(100, option = "D"))
#   my_breaks <- seq(min_val, max_val, length.out = 101)

#   # Save the combined heatmap as SVG
#   tryCatch({
#     svg(file.path(results_dir, "pheatmap_OPA_lower_OPB_upper_final_corrected.svg"), width = 10, height = 10)
#     pheatmap(
#       mat = combined_color_matrix,
#       color = my_color(100),
#       breaks = my_breaks,
#       border_color = NA,
#       cluster_rows = FALSE,
#       cluster_cols = FALSE,
#       order_rows = match(ordered_samples, rownames(combined_color_matrix)),
#       order_cols = match(ordered_samples, colnames(combined_color_matrix)),
#       annotation_row = annotation_row,
#       annotation_col = annotation_col,
#       annotation_colors = list(Authority = authority_colors),
#       display_numbers = combined_display_numbers,
#       number_color = "black",
#       fontsize_number = 5,
#       show_rownames = TRUE,
#       show_colnames = TRUE,
#       main = "Heatmap: OPA (Lower Triangle) and OPB (Upper Triangle)",
#       fontsize = 6,
#       fontsize_row = 6,
#       fontsize_col = 6,
#       cellheight = 7,
#       cellwidth = 7,
#       legend_breaks = round(seq(min_val, max_val, length.out = 5), 2),
#       legend_labels = as.character(round(seq(min_val, max_val, length.out = 5), 2))
#     )
#     dev.off()
#     message("Successfully saved final corrected split triangle OPA/OPB SVG plot: pheatmap_OPA_lower_OPB_upper_final_corrected.svg")
#   }, error = function(e) {
#     warning(paste("Error saving final corrected split triangle OPA/OPB SVG plot:\n", e$message))
#   })

# } else {
#   message("OPA and/or OPB files not found. Skipping split triangle heatmap.")
# }


# # -------------------
# # THIS CODE WORKS FOR PLOTTING OPB AS TEXT AND OPA AS COLOR.
# # --- Special handling for OPA and OPB stats ---
# # Define the base names for OPA and OPB files
# opa_file_base <- "heatmap_data_opa_name_len.csv"
# opb_file_base <- "heatmap_data_opb_name_len.csv"

# opa_file_path <- file.path(results_dir, opa_file_base)
# opb_file_path <- file.path(results_dir, opb_file_base)

# if (file.exists(opa_file_path) && file.exists(opb_file_path)) {
#   message("Processing combined OPA and OPB heatmap...")

#   # Read OPA data
#   opa_data <- tryCatch({
#     read.csv(opa_file_path, stringsAsFactors = FALSE)
#   }, error = function(e) {
#     stop(paste("Error reading OPA CSV file:", opa_file_path, "\n", e$message))
#   })

#   # Read OPB data
#   opb_data <- tryCatch({
#     read.csv(opb_file_path, stringsAsFactors = FALSE)
#   }, error = function(e) {
#     stop(paste("Error reading OPB CSV file:", opb_file_path, "\n", e$message))
#   })

#   # Ensure both dataframes have the 'sample_name' column and are aligned
#   if (!("sample_name" %in% colnames(opa_data)) || !("sample_name" %in% colnames(opb_data))) {
#     stop("Both OPA and OPB CSVs must have a 'sample_name' column.")
#   }

#   # Order both dataframes by sample_name to ensure alignment
#   opa_data <- opa_data %>% arrange(sample_name)
#   opb_data <- opb_data %>% arrange(sample_name)

#   # Convert to matrices
#   opa_matrix <- as.matrix(opa_data[, -which(colnames(opa_data) == "sample_name")])
#   rownames(opa_matrix) <- opa_data$sample_name
#   colnames(opa_matrix) <- opa_data$sample_name

#   opb_matrix <- as.matrix(opb_data[, -which(colnames(opb_data) == "sample_name")])
#   rownames(opb_matrix) <- opb_data$sample_name
#   colnames(opb_matrix) <- opb_data$sample_name

#   # Ensure matrices have the same dimensions and sample names
#   if (!all(rownames(opa_matrix) == rownames(opb_matrix)) || !all(colnames(opa_matrix) == colnames(opb_matrix))) {
#     stop("OPA and OPB matrices do not have matching row/column names or dimensions.")
#   }

#   # Join authority information
#   opa_data_with_authority <- opa_data %>%
#     left_join(sample_authority, by = "sample_name")

#   # Order the data and matrix by authority
#   opa_ordered <- opa_data_with_authority %>% arrange(authority)
#   ordered_samples <- opa_ordered$sample_name
#   ordered_authorities <- opa_ordered$authority

#   # Re-index the matrices based on the ordered samples for both rows and columns
#   opa_matrix_ordered <- opa_matrix[ordered_samples, ordered_samples]
#   opb_matrix_ordered <- opb_matrix[ordered_samples, ordered_samples]

#   # Create a named vector for the colors
#   authority_colors <- setNames(
#     colorRampPalette(c("lightblue", "pink", "darkred"))(length(unique(ordered_authorities))),
#     unique(ordered_authorities)
#   )

#   # Create annotation data frames
#   annotation_row = data.frame(Authority = ordered_authorities, row.names = ordered_samples)
#   annotation_col = data.frame(Authority = ordered_authorities, row.names = ordered_samples)

#   # Define a custom color scale for OPA values
#   my_color <- colorRampPalette(c(viridis(100, option = "D")))
#   my_breaks <- seq(0, 1, length.out = 101) # Assuming OPA values are between 0 and 1

#   # Prepare the display_numbers matrix for OPB
#   # Convert OPB values to character for display, you might want to format them
#   opb_display_numbers <- apply(opb_matrix_ordered, 2, function(x) {
#     sprintf("%.2f", x) # Format to 2 decimal places, adjust as needed
#   })

#   # Save the combined heatmap as SVG
#   tryCatch({
#     svg(file.path(results_dir, "pheatmap_OPA_OPB_combined.svg"), width = 10, height = 10)
#     pheatmap(
#       mat = opa_matrix_ordered,
#       color = my_color(100),
#       breaks = my_breaks,
#       border_color = NA,
#       cluster_rows = FALSE,
#       cluster_cols = FALSE,
#       order_rows = match(ordered_samples, rownames(opa_matrix_ordered)),
#       order_cols = match(ordered_samples, colnames(opa_matrix_ordered)),
#       annotation_row = annotation_row,
#       annotation_col = annotation_col,
#       annotation_colors = list(Authority = authority_colors),
#       display_numbers = opb_display_numbers, # Display OPB values as text
#       number_color = "black", # Color of the displayed numbers
#       fontsize_number = 5, # Font size of the displayed numbers
#       show_rownames = TRUE,
#       show_colnames = TRUE,
#       main = "Heatmap of OPA (Color) and OPB (Text) Values",
#       fontsize = 6,
#       fontsize_row = 6,
#       fontsize_col = 6,
#       cellheight = 7,
#       cellwidth = 7
#     )
#     dev.off()
#     message("Successfully saved combined OPA and OPB SVG plot: pheatmap_OPA_OPB_combined.svg")
#   }, error = function(e) {
#     warning(paste("Error saving combined OPA and OPB SVG plot:\n", e$message))
#   })

# } else {
#   message("OPA and/or OPB files not found. Skipping combined heatmap.")
# }
