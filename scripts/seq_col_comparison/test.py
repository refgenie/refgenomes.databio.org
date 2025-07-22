
# FOR DUMPING GPT GEBERATE CODE TO TRY OUT BEFORE INTEGRATION

# import matplotlib.pyplot as plt
# import seaborn as sns
# import numpy as np
# import pandas as pd
# import os

# # --- Dummy Data Generation (replace with your actual pep_df and all_samples) ---
# # This part is just for making the code runnable and demonstrable.
# # In your actual script, pep_df and all_samples would already be defined.
# np.random.seed(42)
# num_samples = 5
# all_samples = [f'sample_{i}' for i in range(num_samples)]

# data_for_df = []
# relevant_stats = ['jaccard_sequences', 'jaccard_name_len', 'jaccard_lengths', 'jaccard_names']

# for i in range(num_samples):
#     for j in range(i + 1, num_samples):
#         sample1 = all_samples[i]
#         sample2 = all_samples[j]
#         row = {'sample_name_1': sample1, 'sample_name_2': sample2}
#         for stat in relevant_stats:
#             row[stat] = np.random.rand() # Random Jaccard scores between 0 and 1
#         data_for_df.append(row)

# pep_df = pd.DataFrame(data_for_df)
# # --- End Dummy Data Generation ---


# # Ensure results_dir exists
# results_dir = 'results' # Or your actual results directory
# os.makedirs(results_dir, exist_ok=True)

# fig_width_mm = 170
# fig_width_inches = fig_width_mm / 25.4
# fig_height_inches = fig_width_inches * (0.4)

# fig, ax = plt.subplots(figsize=(fig_width_inches, fig_height_inches))

# # Convert all relevant stat columns to numeric
# for stat in relevant_stats:
#     pep_df[stat] = pd.to_numeric(pep_df[stat], errors='coerce')

# # Define histogram bins
# # You'll want to choose an appropriate number of bins.
# # For scores from 0 to 1, 10 bins would be 0.1 wide, 20 bins 0.05 wide, etc.
# num_bins = 20 # You can adjust this for more or less detail/granularity
# bins = np.linspace(0, 1, num_bins + 1) # +1 because np.linspace includes start and end

# # Create a list to store histogram data for each statistic
# all_hist_data = {}

# for stat in relevant_stats:
#     current_values = []
#     for i in range(len(all_samples)):
#         for j in range(i + 1, len(all_samples)):
#             sample1 = all_samples[i]
#             sample2 = all_samples[j]

#             comparison = pep_df[
#                 ((pep_df['sample_name_1'] == sample1) & (pep_df['sample_name_2'] == sample2)) |
#                 ((pep_df['sample_name_1'] == sample2) & (pep_df['sample_name_2'] == sample1))
#             ]

#             if not comparison.empty:
#                 current_values.append(comparison[stat].iloc[0])

#     # Filter out NaN values before calculating histogram
#     current_values = np.array([val for val in current_values if not pd.isna(val)])

#     if len(current_values) > 0:
#         # Calculate histogram frequencies
#         # density=False means counts, density=True would normalize to sum to 1
#         hist_counts, _ = np.histogram(current_values, bins=bins, density=False)
#     else:
#         hist_counts = np.zeros(num_bins) # No data, so all counts are zero

#     all_hist_data[stat] = hist_counts

# # Create a DataFrame for the heatmap using the histogram counts
# # The index will represent the center of the bins or the start of the bins.
# # Let's use the start of the bins for clarity.
# bin_labels = [f'{bins[i]:.2f}-{bins[i+1]:.2f}' for i in range(num_bins)]
# heatmap_df = pd.DataFrame(all_hist_data, index=bin_labels).T

# # Plot the heatmap (Jaccard Similarities on Y, Jaccard Score Bins on X)
# # Use a sequential colormap for bright to dark hits (higher counts = darker)
# # For example, 'Blues' will make higher counts darker blue.
# sns.heatmap(heatmap_df, annot=False, cmap="Blues", ax=ax, cbar_kws={'label': 'Frequency (Counts)'})

# # Add labels and title
# ax.set_xlabel('Jaccard Score Bin')
# ax.set_ylabel('Jaccard Similarities')
# ax.set_title('Frequency Distribution of Jaccard Similarities Across Score Bins')

# # Set x-axis ticks and labels to represent the bins clearly
# # Using the bin labels directly makes it clear which range each column represents.
# ax.set_xticks(np.arange(num_bins) + 0.5) # Center ticks in the middle of the bins
# ax.set_xticklabels(bin_labels, rotation=90, ha='center', fontsize=8) # Rotate for readability

# fig.tight_layout()

# # Save the heatmap plot
# output_path = os.path.join(results_dir, 'jaccard_heatmap_histogram.svg') # New filename
# plt.savefig(output_path, dpi=300, bbox_inches='tight')
# plt.close()

# print(f"Heatmap saved to {output_path}")

# # Save the histogram numbers
# hist_file = os.path.join(results_dir, "histogram_counts.csv")
# heatmap_df.to_csv(hist_file, index=True) # index=True to keep bin labels

# print(f"Histogram counts saved to {hist_file}")