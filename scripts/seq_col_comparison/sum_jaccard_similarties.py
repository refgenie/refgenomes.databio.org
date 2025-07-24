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

    # Initialize dictionaries to store summed jaccard scores and counts for INTER-PROVIDER comparisons
    provider_jaccard_sums_vs_others = {}
    provider_pair_counts_vs_others = {}

    # Get all unique sample names from the heatmap data (index or columns, they are the same)
    # Use only common samples to avoid issues with .get() returning None
    all_samples = list(common_samples)
    print(f"\nProcessing {len(all_samples)} common samples.")


    # Iterate through unique pairs of samples (i < j to avoid duplicates and self-comparisons)
    found_any_inter_matching_pairs = False
    for i in range(len(all_samples)):
        for j in range(i + 1, len(all_samples)):
            sample1 = all_samples[i]
            sample2 = all_samples[j]

            authority1 = sample_to_authority.get(sample1)
            authority2 = sample_to_authority.get(sample2)

            # Check if authorities are DIFFERENT
            if authority1 != authority2:
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


    # --- Plotting Averages for Inter-Provider Scores ---
    # Convert the results to a DataFrame for plotting INTER-PROVIDER AVERAGES
    inter_provider_averaged_df = pd.DataFrame(list(provider_jaccard_averages_vs_others.items()), columns=['Provider', 'Average Jaccard Score vs. Others'])

    # Create the bar graph for inter-provider averages only if dataframe is not empty
    if not inter_provider_averaged_df.empty:
        # Sort for better visualization
        inter_provider_averaged_df = inter_provider_averaged_df.sort_values(by='Average Jaccard Score vs. Others', ascending=False)

        plt.figure(figsize=(12, 7))
        sns.barplot(x='Provider', y='Average Jaccard Score vs. Others', data=inter_provider_averaged_df, palette='viridis')
        plt.title(f'Average Jaccard Scores: Each Provider vs. All Other Providers ({stat})')
        plt.xlabel('Provider')
        plt.ylabel('Average Jaccard Score')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        save_path_avg_vs_others = os.path.join(results_dir,f'bargraph_{stat}_AVERAGES_VS_OTHERS.png' )
        plt.savefig(save_path_avg_vs_others)
        plt.close()
        print(f"Saved: {save_path_avg_vs_others}")
    else:
        print(f"\nNo inter-provider Jaccard averages found for {stat} to plot a bar graph.")
    # --- End Plotting Averages ---

print(f"\nAll inter-provider average plots saved to: {results_dir}")