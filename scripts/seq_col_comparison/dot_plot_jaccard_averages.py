import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os


import matplotlib
matplotlib.rcParams["svg.fonttype"] = "none" # do not embed directly, instead the downstream program will view with system fonts
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams["text.usetex"] = False 

# MUST POINT THIS TO THE SAME DIRECTORY AS THEstats_graphs.py output
# Using the original file path provided by the user
results_dir = "/home/drc/Downloads/refgenomes_pics_test/21jul2025/TESTING_PROVIDER_GROUPING/"

#stats_list = ['jaccard_sequences', 'jaccard_names', 'jaccard_lengths', 'jaccard_name_len']
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

    missing_in_authority = jaccard_samples - authority_samples
    missing_in_jaccard = authority_samples - jaccard_samples
    common_samples = jaccard_samples.intersection(authority_samples)

    print(f"\n--- Processing Statistic: {stat} ---")
    print(f"Number of samples in Jaccard heatmap: {len(jaccard_samples)}")
    print(f"Number of samples in Authority mapping: {len(authority_samples)}")
    print(f"Number of common samples: {len(common_samples)}")

    if missing_in_authority:
        print(f"\nSamples in Jaccard heatmap but missing in Authority mapping ({len(missing_in_authority)}):")
        for sample in list(missing_in_authority)[:10]: # Print up to 10 for brevity
            print(f"  - {sample}")

    if missing_in_jaccard:
        print(f"\nSamples in Authority mapping but missing in Jaccard heatmap ({len(missing_in_jaccard)}):")
        for sample in list(missing_in_jaccard)[:10]: # Print up to 10 for brevity
            print(f"  - {sample}")
    # --- End Debugging Sample Name Mismatches ---

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

    print("\nProvider Jaccard Sums (vs. Others):", provider_jaccard_sums_vs_others)
    print("Provider Pair Counts (vs. Others):", provider_pair_counts_vs_others)
    print("Found any inter-group matching pairs:", found_any_inter_matching_pairs)

    # --- Calculate Averages for Inter-Provider Scores ---
    provider_jaccard_averages_vs_others = {}
    for provider in provider_jaccard_sums_vs_others:
        if provider_pair_counts_vs_others.get(provider, 0) > 0: # Avoid division by zero
            provider_jaccard_averages_vs_others[provider] = provider_jaccard_sums_vs_others[provider] / provider_pair_counts_vs_others[provider]
    print("Provider Jaccard Averages (vs. Others):", provider_jaccard_averages_vs_others)
    # --- End Calculate Averages ---

    # Convert the results to a DataFrame for plotting INTER-PROVIDER AVERAGES
    inter_provider_averaged_df = pd.DataFrame(list(provider_jaccard_averages_vs_others.items()), columns=['Provider', 'Average Jaccard Score'])
    inter_provider_averaged_df['Statistic'] = stat # Add a column to identify the statistic

    # Append to the list for combined plotting
    if not inter_provider_averaged_df.empty:
        all_inter_provider_averages.append(inter_provider_averaged_df)
    else:
        print(f"\nNo inter-provider Jaccard averages found for {stat}.")


# --- Combined Plotting of All Inter-Provider Averages as a Dot Plot ---
if all_inter_provider_averages:
    combined_df = pd.concat(all_inter_provider_averages)

    # Ensure the order of statistics on the y-axis is consistent
    combined_df['Statistic'] = pd.Categorical(combined_df['Statistic'], categories=stats_list, ordered=True)
    combined_df = combined_df.sort_values(by=['Statistic', 'Provider'], ascending=[True, True])

    max_jaccard_score = combined_df['Average Jaccard Score'].max()
    x_axis_upper_limit = max_jaccard_score * 1.05 # Add a 5% buffer for aesthetics

    plt.figure(figsize=(8, 4)) # Adjust figure size as needed for a dot plot
    sns.scatterplot(
        x='Average Jaccard Score',
        y='Statistic',
        hue='Provider',
        style='Provider',
        data=combined_df,
        s=125, # Size of the dots
        palette='tab10', # A distinct color palette
        alpha=0.90, # Transparency of the dots
        markers=True

    )
    plt.title('Average Jaccard Scores: Each Provider vs. All Other Providers (Dot Plot)')
    plt.xlabel('Average Jaccard Score')
    plt.ylabel('Jaccard Statistic')

    # # CURRENT DATA IS LESS THAN 0.4!!
    # plt.xlim(0, 1.0) # Set x-axis limit from 0 to 1.0
    plt.xlim(0, x_axis_upper_limit) # Set x-axis limit dynamically
    plt.legend(title='Provider', bbox_to_anchor=(1.05, 1), loc='upper left') # Place legend outside the plot
    plt.grid(True, linestyle='--', alpha=0.6) # Add a grid for readability
    plt.tight_layout()
    # Save the plot with a distinct name
    save_path_combined_avg_dot_plot = os.path.join(results_dir,'dotplot_ALL_STATS_COMBINED_AVERAGES_VS_OTHERS_fixed.svg')
    plt.savefig(save_path_combined_avg_dot_plot)
    plt.close()
    print(f"\nSaved combined dot plot: {save_path_combined_avg_dot_plot}")
else:
    print("\nNo data to generate a combined dot plot.")

print(f"\nAll processing complete. Results saved to: {results_dir}")