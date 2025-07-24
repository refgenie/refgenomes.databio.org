import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

# MUST POINT THIS TO THE SAME DIRECTORY AS THEstats_graphs.py output
results_dir = "/home/drc/Downloads/refgenomes_pics_test/21jul2025/TESTING_PROVIDER_GROUPING/"


stats_list = ['jaccard_sequences', 'jaccard_names', 'jaccard_lengths', 'jaccard_name_len']

#stats_list = ['jaccard_name_len']

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

    print(f"\nNumber of samples in Jaccard heatmap: {len(jaccard_samples)}")
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
    #print(sample_to_authority.values())

    # Initialize a dictionary to store summed jaccard scores per provider
    provider_jaccard_sums = {}

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

            authority1 = sample_to_authority.get(sample1) # Should not be None now if using common_samples
            authority2 = sample_to_authority.get(sample2) # Should not be None now

            if authority1 == authority2:
                score = pd.to_numeric(jaccard_df.loc[sample1, sample2], errors='coerce')

                if not pd.isna(score):
                    found_any_matching_pairs = True
                    provider_jaccard_sums[authority1] = provider_jaccard_sums.get(authority1, 0) + score

    print("\nProvider Jaccard Sums:", provider_jaccard_sums)
    print("Found any matching pairs:", found_any_matching_pairs)

    # Convert the results to a DataFrame for plotting
    summed_jaccard_df = pd.DataFrame(list(provider_jaccard_sums.items()), columns=['Provider', 'Summed Jaccard Score'])
    #print(summed_jaccard_df.head(11))

    # Create the bar graph only if summed_jaccard_df is not empty
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
        save_path = os.path.join(results_dir,f'bargraph_{stat}_SUMS.png' )
        plt.savefig(save_path)
        plt.close()
        # plt.show()
        # plt.close()
    else:
        print("\nNo within-provider Jaccard scores found to plot a bar graph.")