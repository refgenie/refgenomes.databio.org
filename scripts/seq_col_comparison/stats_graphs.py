import sys
import os
import json
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



#---- for running locally only
# json_files = []
# if os.path.isdir("/home/drc/Downloads/jsons_from_rivanna/json/"):
#     for filename in os.listdir("/home/drc/Downloads/jsons_from_rivanna/json/"):
#         if filename.endswith(".json"):
#             full_path = os.path.join("/home/drc/Downloads/jsons_from_rivanna/json/", filename)
#             json_files.append(full_path)

def get_sequence_length(digest):

    json_fp_1 = os.path.join("/home/drc/Downloads/jsons_from_rivanna/json/", digest+".json")
    # json_fp_2 = os.path.join("/home/drc/Downloads/jsons_from_rivanna/json/", digest2+".json")
    with open(json_fp_1, "r") as f:
        reloaded_dict1 = json.load(fp=f)
    
    return len(reloaded_dict1['sorted_sequences'])
    # with open(json_fp_2, "r") as f:
    #     reloaded_dict2 = json.load(fp=f)


# --------

phc = PEPHubClient()
pep = phc.load_project(results_pep)
print(pep["_sample_df"])

pep_df = pep["_sample_df"]
print(pep_df)
#new_df = pep_df.copy()


all_relevant_stats = ['overlap_coefficient_names',	
                      'overlap_coefficient_lengths',
                      'overlap_coefficient_sequences',
                      'overlap_coeff_name_len',	
                      'jaccard_names', 
                      'jaccard_lengths', 
                      'jaccard_similarity_weighted_length', ''
                      'jaccard_sequences',
                      'jaccard_name_len', 
                      'f1_names',
                      'f1_lengths',
                      'f1_sequences',
                      'f1_weighted_lengths',
                      'f1_name_len_pairs'
                      ]
#all_relevant_stats = ['jaccard_names']

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

    #     # 1. Pivot the dataframe to create a matrix
    
    # heatmap_data_asymmetrical = new_df.pivot(index='sample_name_1', columns='sample_name_2', values=stat)

    # # 2. Create the heatmap
    # plt.figure(figsize=(10, 8))  # Adjust figure size as needed
    # sns.heatmap(heatmap_data_asymmetrical, annot=True, cmap='viridis', fmt=".2f", linewidths=.5, cbar_kws={'label': 'Jaccard Similarity'})

    # # Optional: Add a title
    # plt.title('Asymmetrical Jaccard Similarity Heatmap (Direct from Pivot)')

    # # Optional: Rotate x-axis labels for better readability
    # plt.xticks(rotation=90)
    # plt.yticks(rotation=0)

    # # Show the plot
    # plt.tight_layout()  # Adjust layout to prevent labels from being cut off
    # plt.show()

    row_sums = heatmap_data.sum(axis=1)

    # Create row sum graph
    # 2. Display the row sums
    print(f"Sum of {stat} Scores for Each Sample:")
    print(row_sums)
    # 1. Sort the row sums (optional, for better visualization)
    row_sums_sorted = row_sums.sort_values(ascending=False)  # Sort in descending order

    # 2. Create the bar plot
    plt.figure(figsize=(12, 8))
    plt.bar(row_sums_sorted.index, row_sums_sorted.values)

    # 3. Set labels and title
    plt.xlabel("Sample")
    plt.ylabel(f"Sum of {stat}")
    plt.title("Total Similarity of Each Sample")
    plt.xticks(rotation=90, ha="right", fontsize=8)  # Rotate sample names for readability
    plt.tight_layout()

    # 4. Show and save the plot
    #plt.show()
    output_path = os.path.join(results_dir,stat+"_row_sums")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')


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


# Create asymmetrical graph

df = pep["_sample_df"]

stat = 'overlap_coeff_name_len'
similarity_df_sorted = df.sort_values(by=stat, ascending=False)  # Example: Descending order

results_df = pd.DataFrame(columns=['Samples_Superset', 'Samples_Subset', stat])

for index, row in similarity_df_sorted.iterrows():
    sample1 = row['sample_name_1']
    digest1 = row['digest1']
    sample2 = row['sample_name_2']
    digest2 = row['digest2']
    similarity = row[stat]

    if get_sequence_length(digest1) > get_sequence_length(digest2):
        higher_value_sample = sample1
        compared_sample = sample2
    else:
        higher_value_sample = sample2
        compared_sample = sample1

    # Append the result to the new DataFrame
    results_df = pd.concat([results_df, pd.DataFrame({'Samples_Superset': [higher_value_sample],
                                                     'Samples_Subset': [compared_sample],
                                                     stat: [similarity]})], ignore_index=True)

print(results_df)

heatmap_data = results_df.pivot(index='Samples_Superset', columns='Samples_Subset', values=stat)

plt.figure(figsize=(10, 8))  
sns.heatmap(heatmap_data, annot=True, cmap='viridis', fmt=".2f", linewidths=.5, linecolor='black')
plt.title('Heatmap of Result by Samples_Superset and Samples_Subset')
plt.xlabel('Samples_Subset')
plt.ylabel('Samples_Superset')

plt.xticks(rotation=90)
plt.yticks(rotation=0)

plt.tight_layout()  # Adjust layout to prevent labels from being cut off
plt.show()


# Pivot the DataFrame to create the matrix for the heatmap
# heatmap_data = df.pivot(index='sample_name_1', columns='sample_name_2', values='overlap_coefficient_names')
# heatmap_data = df.pivot(index='sample_name_2', columns=['sample_name_1'], values='overlap_coefficient_names')
# # mask = heatmap_data.isnull()


# # Create the heatmap
# plt.figure(figsize=(10, 8))  # Adjust figure size as needed
# sns.heatmap(heatmap_data, annot=True, cmap='viridis', fmt=".2f", linewidths=.5, linecolor='black')
# plt.title('Heatmap of Result by Sample1 and Sample2')
# plt.xlabel('Sample2')
# plt.ylabel('Sample1')
# # Optional: Rotate x-axis labels for better readability
# plt.xticks(rotation=90)
# plt.yticks(rotation=0)

# # Show the plot
# #plt.tight_layout()  # Adjust layout to prevent labels from being cut off
# plt.show()

# 1. Get all unique sample names from both columns
# all_samples = pd.concat([df['sample_name_1'], df['sample_name_2']]).unique()

# # 2. Create a DataFrame with all unique samples as index and columns, initialized with NaN
# heatmap_data_fixed = pd.DataFrame(index=all_samples, columns=all_samples)

# # 3. Populate the heatmap DataFrame with Jaccard similarity values
# for index, row in df.iterrows():
#     sample1 = row['sample_name_1']
#     sample2 = row['sample_name_2']
#     similarity = row['overlap_coefficient_names']
#     if sample1 in heatmap_data_fixed.index and sample2 in heatmap_data_fixed.columns:
#         if similarity:
#             heatmap_data_fixed.loc[sample1, sample2] = similarity

# # 4. Create the heatmap
# plt.figure(figsize=(12, 10))  # Adjust figure size as needed
# sns.heatmap(heatmap_data_fixed, annot=True, cmap='viridis', fmt=".2f", linewidths=.5, cbar_kws={'label': 'Jaccard Similarity'})

# # Customize the plot
# plt.title('Jaccard Similarity Heatmap (All Unique Samples on Both Axes)')
# plt.xlabel('Sample 2')
# plt.ylabel('Sample 1')
# plt.xticks(rotation=90)
# plt.yticks(rotation=0)
# plt.tight_layout()
# plt.show()

# # Save the plot (optional)
# output_path = 'jaccard_similarity_heatmap_fixed_axes.png'
# plt.savefig(output_path, dpi=300, bbox_inches='tight')
# print(f"Heatmap saved to: {output_path}")


