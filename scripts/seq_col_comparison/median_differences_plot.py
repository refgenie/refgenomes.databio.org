import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import matplotlib

matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams["text.usetex"] = False

# MUST POINT THIS TO THE SAME DIRECTORY AS THE stats_graphs.py output
# Using the original file path provided by the user
results_dir = "/home/drc/Downloads/refgenomes_pics_test/21jul2025/TESTING_PROVIDER_GROUPING/"
results_dir = "/home/drc/Downloads/refgenomes_pics_test/21jul2025/testing_determinism/"

stats_list = ['jaccard_lengths', 'jaccard_sequences', 'jaccard_names', 'jaccard_name_len']
# List to store DataFrames for each statistic for combined plotting
all_inter_provider_averages = []

for stat in stats_list:
    # Construct the full path for the heatmap data CSV
    temp_string = f'heatmap_data_{stat}.csv'
    heatmap_filepath = os.path.join(results_dir, temp_string)

    # Load the heatmap data
    try:
        jaccard_df = pd.read_csv(heatmap_filepath, index_col='sample_name')
    except FileNotFoundError:
        print(f"Error: Heatmap data file not found at {heatmap_filepath}. Skipping this statistic.")
        continue # Skip to the next statistic if file is not found

    # --- FIX: Make the Jaccard DataFrame symmetrical ---
    # This step is crucial because the upstream code generates lower-triangle only data.
    # Fill NaN values (which are in the upper triangle) with values from the transpose (lower triangle)
    jaccard_df = jaccard_df.fillna(jaccard_df.T)
    # Ensure the diagonal is 1.0 for self-similarity, as Jaccard of an item with itself is 1.
    np.fill_diagonal(jaccard_df.values, 1.0)
    # Ensure all values are numeric after filling, coercing any remaining non-numeric
    jaccard_df = jaccard_df.apply(pd.to_numeric, errors='coerce')
    # --- END FIX ---

    # Load the sample authority mapping
    authority_filepath = os.path.join(results_dir, 'sample_authority.csv')
    try:
        authority_df = pd.read_csv(authority_filepath)
    except FileNotFoundError:
        print(f"Error: Sample authority file not found at {authority_filepath}. Cannot process without it.")
        exit() # Exit if the authority file is missing, as it's critical

    # --- Debugging Sample Name Mismatches (kept for user's diagnostic purposes) ---
    jaccard_samples = set(jaccard_df.index.tolist())
    authority_samples = set(authority_df['sample_name'].tolist())

    common_samples = jaccard_samples.intersection(authority_samples)

    print(f"\n--- Processing Statistic: {stat} ---")
    print(f"Number of samples in Jaccard heatmap: {len(jaccard_samples)}")
    print(f"Number of samples in Authority mapping: {len(authority_samples)}")
    print(f"Number of common samples: {len(common_samples)}")

    # Create a dictionary for sample to authority mapping
    sample_to_authority = dict(zip(authority_df['sample_name'], authority_df['authority']))

    # Initialize dictionaries to store summed jaccard scores and counts for INTER-PROVIDER comparisons
    provider_jaccard_sums_vs_others = {}
    provider_pair_counts_vs_others = {}

    # Get all unique sample names from the heatmap data (index or columns, they are the same)
    # Use only common samples to avoid issues with .get() returning None
    all_samples = list(common_samples)
    print(f"\nProcessing {len(all_samples)} common samples.")

    found_any_inter_matching_pairs = False
    for i in range(len(all_samples)):
        for j in range(i + 1, len(all_samples)):
            sample1 = all_samples[i]
            sample2 = all_samples[j]

            authority1 = sample_to_authority.get(sample1)
            authority2 = sample_to_authority.get(sample2)

            # Check if authorities are DIFFERENT and both authorities are valid (not None)
            if authority1 is not None and authority2 is not None and authority1 != authority2:
                # Retrieve score from the now symmetrical jaccard_df
                score = pd.to_numeric(jaccard_df.loc[sample1, sample2], errors='coerce')

                if not pd.isna(score):
                    found_any_inter_matching_pairs = True

                    # Accumulate for authority1 vs others
                    provider_jaccard_sums_vs_others[authority1] = provider_jaccard_sums_vs_others.get(authority1, 0) + score
                    provider_pair_counts_vs_others[authority1] = provider_pair_counts_vs_others.get(authority1, 0) + 1

                    # Accumulate for authority2 vs others (the comparison is symmetric, so it counts for both)
                    provider_jaccard_sums_vs_others[authority2] = provider_jaccard_sums_vs_others.get(authority2, 0) + score
                    provider_pair_counts_vs_others[authority2] = provider_pair_counts_vs_others.get(authority2, 0) + 1

    # --- Calculate Averages for Inter-Provider Scores ---
    provider_jaccard_averages_vs_others = {}
    for provider in provider_jaccard_sums_vs_others:
        if provider_pair_counts_vs_others.get(provider, 0) > 0: # Avoid division by zero
            provider_jaccard_averages_vs_others[provider] = provider_jaccard_sums_vs_others[provider] / provider_pair_counts_vs_others[provider]
    # --- End Calculate Averages ---

    # Convert the results to a DataFrame for plotting INTER-PROVIDER AVERAGES
    inter_provider_averaged_df = pd.DataFrame(list(provider_jaccard_averages_vs_others.items()), columns=['Provider', 'Average Jaccard Score'])
    inter_provider_averaged_df['Statistic'] = stat # Add a column to identify the statistic

    # Append to the list for combined plotting
    if not inter_provider_averaged_df.empty:
        all_inter_provider_averages.append(inter_provider_averaged_df)
    else:
        print(f"\nNo inter-provider Jaccard averages found for {stat}.")


# --- Combined Plotting of All Inter-Provider Averages as Faceted Median Difference Plots ---
if all_inter_provider_averages:
    combined_df = pd.concat(all_inter_provider_averages)

    # Calculate the median for each statistic separately
    combined_df['Median_Jaccard_Score_Per_Stat'] = combined_df.groupby('Statistic')['Average Jaccard Score'].transform('median')

    # Calculate the difference from the median for each statistic
    combined_df['Median_Difference'] = combined_df['Average Jaccard Score'] - combined_df['Median_Jaccard_Score_Per_Stat']

    # Sort the DataFrame. First by Statistic, then by Median_Difference
    # This ensures that within each facet, the bars will be ordered by Median_Difference
    combined_df_sorted = combined_df.sort_values(by=['Statistic', 'Median_Difference'], ascending=[True, True])

    # Determine the number of unique providers to adjust plot height
    # Use nunique for 'Provider' within the combined_df to get total distinct providers
    num_providers = combined_df['Provider'].nunique()
    # Adjust height based on the number of providers for readability
    plot_height_per_facet = max(2, num_providers * 0.25) # Minimum height of 2, scales with providers

    # Create the faceted plot
    g = sns.catplot(
        x='Median_Difference',
        y='Provider',
        col='Statistic', # Create separate columns for each statistic
        data=combined_df_sorted,
        kind='bar',
        col_wrap=2, # Wrap columns after 2 plots
        height=plot_height_per_facet,
        aspect=1.25, # Adjust aspect ratio for wider bars
        palette='vlag', # Diverging color palette for positive/negative differences
    )

    # Add a vertical line at 0 for each subplot and customize labels
    for ax in g.axes.flat:
        ax.axvline(0, color='grey', linestyle='--', linewidth=0.8)
        ax.set_xlabel(f'Difference from Median Jaccard Score') # Set x-label for each subplot
        ax.set_ylabel('Provider') # Set y-label for each subplot
        
    g.set_titles(col_template="{col_name}") # Set title for each facet to just the statistic name

    plt.suptitle('Median Difference Plots: Each Provider vs. Statistic-Specific Median', y=1.02) # Overall title
    #plt.tight_layout(rect=[0, 0, 1, 0.98]) # Adjust layout to make space for suptitle

    # Save the plot with a distinct name
    save_path_faceted_median_diff_plot = os.path.join(results_dir, 'faceted_median_difference_plots_ALL_STATS.svg')
    plt.savefig(save_path_faceted_median_diff_plot)
    plt.close()
    print(f"\nSaved faceted median difference plot: {save_path_faceted_median_diff_plot}")

else:
    print("\nNo data to generate faceted median difference plots.")

print(f"\nAll processing complete. Results saved to: {results_dir}")