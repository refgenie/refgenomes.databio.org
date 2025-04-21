import sys
import os
from pephubclient import PEPHubClient
from refget import fasta_to_digest, fasta_to_seqcol_dict, compare_seqcols, SequenceCollection
from itertools import combinations
from pprint import pprint


#results_pep = "donaldcampbelljr/test_seq_col_results:default"

results_dir = sys.argv[1] # output dir for graphs 
results_pep = sys.argv[2] # input pep for graphing
species_title = sys.argv[3] # additional information for the title


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

phc = PEPHubClient()
pep = phc.load_project(results_pep)
print(pep["_sample_df"])

pep_df = pep["_sample_df"]
print(pep_df)


all_relevant_stats = ['overlap_coefficient_names',	'overlap_coefficient_lengths',	'jaccard_names', 'jaccard_lengths', 'jaccard_similarity_weighted_length', 'jaccard_sequences', 'overlap_coefficient_sequences']

for stat in all_relevant_stats:
    pep_df[stat] = pd.to_numeric(pep_df[stat], errors='coerce')
    # 1. Get all unique sample names
    all_samples = pd.concat([pep_df['sample_name_1'], pep_df['sample_name_2']]).unique()

    # 2. Create an empty DataFrame for the heatmap, initialized with NaN
    heatmap_data = pd.DataFrame(index=all_samples, columns=all_samples)

    # 3. Iterate over all pairs of unique sample names (O(n^2))
    for sample1 in all_samples:
        for sample2 in all_samples:
            # a. If sample1 and sample2 are the same, assign 1 to both (symmetric)
            if sample1 == sample2:
                heatmap_data.loc[sample1, sample2] = 1.0
            else:
                # b. Check if the combination exists in the original dataframe
                comparison = pep_df[
                    ((pep_df['sample_name_1'] == sample1) & (pep_df['sample_name_2'] == sample2)) |
                    ((pep_df['sample_name_1'] == sample2) & (pep_df['sample_name_2'] == sample1))
                ]
                if not comparison.empty:
                    similarity_score = comparison[stat].iloc[0]
                    heatmap_data.loc[sample1, sample2] = similarity_score
                    heatmap_data.loc[sample2, sample1] = similarity_score



    print(heatmap_data)
    heatmap_data = heatmap_data.apply(pd.to_numeric, errors='coerce')

    # 4. Create the heatmap
    plt.figure(figsize=(18, 15))  # Adjust figure size as needed
    sns.heatmap(heatmap_data, annot=True, cmap='viridis', fmt=".2f", linewidths=.5, cbar_kws={'label': stat},annot_kws={"size": 9},vmin=0.0, vmax=1.0)
    plt.title(f'{stat} Heatmap - {species_title}')
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.tight_layout()
    #plt.show()
    output_path = os.path.join(results_dir,stat)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Heatmap saved to: {output_path}")


# Create Histogram
for stat in all_relevant_stats:
    pep_df[stat] = pd.to_numeric(pep_df[stat], errors='coerce')
    similarity_series = pd.Series(pep_df[stat])
    print(similarity_series)
    # 2. Define the bin size
    bin_width = 0.05  # You can adjust this value (e.g., 0.1, 0.02)

    # 3. Create the bins for the histogram
    bins = np.arange(0, 1.0 + bin_width, bin_width)

    # 4. Create the frequency plot (histogram)
    plt.figure(figsize=(10, 6))
    n, bins, patches = plt.hist(similarity_series, bins=bins, edgecolor='black', alpha=0.7)

    for i in range(len(n)):
        count = int(n[i])  # Get the count for the current bin
        if count !=0:
            x_position = (bins[i] + bins[i+1]) / 2  # Center of the bar
            y_position = n[i]  # Top of the bar
            plt.text(x_position, y_position, str(count), ha='center', va='bottom')
    
    # 5. Set the labels and title
    plt.xlabel(f' {stat} (0 to 1.0)')
    plt.ylabel('Count Frequency')
    plt.title(f'Frequency Distribution of {stat} - {species_title}')
    plt.xticks(np.arange(0, 1.0 + bin_width, bin_width))  # Set x-axis ticks
    plt.grid(axis='y', linestyle='--')  # Add a grid for better readability

    # 6. Show the plot
    plt.tight_layout()
    #plt.show()
    output_path = os.path.join(results_dir,stat+'_histogram')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')