from pephubclient import PEPHubClient
from refget import fasta_to_digest, fasta_to_seqcol_dict, compare_seqcols, SequenceCollection
from itertools import combinations
from pprint import pprint

results_pep = "donaldcampbelljr/test_seq_col_results:default"


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

phc = PEPHubClient()
pep = phc.load_project(results_pep)
print(pep["_sample_df"])

pep_df = pep["_sample_df"]
print(pep_df)


all_relevant_stats = ['overlap_coefficient_names',	'overlap_coefficient_lengths',	'jaccard_names', 'jaccard_lengths', 'jaccard_similarity_weighted_length']

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
    plt.figure(figsize=(12, 10))  # Adjust figure size as needed
    sns.heatmap(heatmap_data, annot=True, cmap='viridis', fmt=".2f", linewidths=.5, cbar_kws={'label': stat})
    plt.title(f'{stat} Heatmap')
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.show()