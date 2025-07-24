import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

# MUST POINT THIS TO THE SAME DIRECTORY AS THEstats_graphs.py output
results_dir = "/home/drc/Downloads/refgenomes_pics_test/21jul2025/TESTING_PROVIDER_GROUPING/"


stats_list = ['jaccard_sequences', 'jaccard_names', 'jaccard_lengths', 'jaccard_name_len']

#stats_list = ['jaccard_name_len'] # Uncomment this line to run only one statistic

for stat in stats_list:
    # Load the heatmap data
    temp_string = f'heatmap_data_{stat}.csv'
    jaccard_df = pd.read_csv(os.path.join(results_dir,temp_string), index_col='sample_name')

    # Load the sample authority mapping
    authority_df = pd.read_csv(os.path.join(results_dir,'sample_authority.csv'))

    # --- Debugging Sample Name Mismatches ---
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

    # Initialize dictionaries to store summed jaccard scores and counts per provider
    provider_jaccard_sums = {}
    provider_pair_counts = {} # New: To store counts for average calculation

    # Get all unique sample names from the heatmap data (index or columns, they are the same)
    # Use only common samples to avoid issues with .get() returning None
    all_samples = list(common_samples)
    print(f"\nProcessing {len(all_samples)} common samples.")


    # Iterate through unique pairs of samples (i < j to avoid duplicates and self-comparisons)
    found_any_matching_pairs = False
    for i in range(len(all_samples)):
        for j in range(i + 1, len(all_samples)):
            sample1 = all_samples[i]
            sample2 = all_samples[j]

            authority1 = sample_to_authority.get(sample1)
            authority2 = sample_to_authority.get(sample2)

            if authority1 == authority2:
                # Get the Jaccard score, coercing errors to NaN
                score = pd.to_numeric(jaccard_df.loc[sample1, sample2], errors='coerce')

                # Only process if the score is a valid number
                if not pd.isna(score):
                    found_any_matching_pairs = True
                    provider_jaccard_sums[authority1] = provider_jaccard_sums.get(authority1, 0) + score
                    provider_pair_counts[authority1] = provider_pair_counts.get(authority1, 0) + 1 # New: Increment count for average

    print("\nProvider Jaccard Sums:", provider_jaccard_sums)
    print("Provider Pair Counts (for averages):", provider_pair_counts) # New: Print counts
    print("Found any matching pairs:", found_any_matching_pairs)

    # --- New: Calculate Averages ---
    provider_jaccard_averages = {}
    for provider in provider_jaccard_sums:
        if provider_pair_counts.get(provider, 0) > 0: # Avoid division by zero
            provider_jaccard_averages[provider] = provider_jaccard_sums[provider] / provider_pair_counts[provider]
    print("Provider Jaccard Averages:", provider_jaccard_averages)
    # --- End New: Calculate Averages ---


    # Convert the results to a DataFrame for plotting SUMS
    summed_jaccard_df = pd.DataFrame(list(provider_jaccard_sums.items()), columns=['Provider', 'Summed Jaccard Score'])

    # Create the bar graph for sums only if summed_jaccard_df is not empty
    if not summed_jaccard_df.empty:
        # Sort for better visualization
        summed_jaccard_df = summed_jaccard_df.sort_values(by='Summed Jaccard Score', ascending=False)

        plt.figure(figsize=(12, 7))
        sns.barplot(x='Provider', y='Summed Jaccard Score', data=summed_jaccard_df, palette='viridis')
        plt.title(f'Sum of Jaccard Scores Within Provider Groups ({stat})')
        plt.xlabel('Provider')
        plt.ylabel('Summed Jaccard Score')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        save_path_sum = os.path.join(results_dir,f'bargraph_{stat}_SUMS.png' )
        plt.savefig(save_path_sum)
        plt.close()
        print(f"Saved: {save_path_sum}")
    else:
        print(f"\nNo within-provider Jaccard sums found for {stat} to plot a bar graph.")

    # --- New: Plotting Averages ---
    # Convert the results to a DataFrame for plotting AVERAGES
    averaged_jaccard_df = pd.DataFrame(list(provider_jaccard_averages.items()), columns=['Provider', 'Average Jaccard Score'])

    # Create the bar graph for averages only if averaged_jaccard_df is not empty
    if not averaged_jaccard_df.empty:
        # Sort for better visualization
        averaged_jaccard_df = averaged_jaccard_df.sort_values(by='Average Jaccard Score', ascending=False)

        plt.figure(figsize=(12, 7))
        sns.barplot(x='Provider', y='Average Jaccard Score', data=averaged_jaccard_df, palette='viridis')
        plt.title(f'Average of Jaccard Scores Within Provider Groups ({stat})')
        plt.xlabel('Provider')
        plt.ylabel('Average Jaccard Score')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        save_path_avg = os.path.join(results_dir,f'bargraph_{stat}_AVERAGES.png' )
        plt.savefig(save_path_avg)
        plt.close()
        print(f"Saved: {save_path_avg}")
    else:
        print(f"\nNo within-provider Jaccard averages found for {stat} to plot a bar graph.")
    # --- End New: Plotting Averages ---

print(f"\nAll results (sum and average plots) saved to: {results_dir}")