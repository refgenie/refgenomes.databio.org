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

def get_sequence_length_local(digest):

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

desired_order = None

desired_order = [
"GRCh38.p0-fasta-genomic",
"GRCh38.p1-fasta-genomic",
"GRCh38.p2-fasta-genomic",
"GRCh38.p6-fasta-genomic",
"GRCh38.p7-fasta-genomic",
"GRCh38.p8-fasta-genomic",
"GRCh38.p12-fasta-genomic",
"GRCh38.p13-fasta-genomic",
"GRCh38.p14-fasta-genomic",
    
]

#pep_df = pep_df.sort_values(by="ORDER")
#print(pep_df)

# ncbi_target_samples = [  
# "GRCh38.p14-fasta-no-alt-analysis",
# #"GRCh37.p13-fasta-no-alt-analysis",
# #"GRCh37.p13-fasta-genomic",
# "GRCh38.p14-fasta-full-analysis-plus-hs38d1",
# #"GRCh37.p13-fasta-full-analysis",
# "GRCh38.p14-fasta-no-alt-plus-hs38d1",
# "GRCh38.p14-fasta-genomic",
# "GRCh38.p14-fasta-full-analysis",
# ]

# ucsc_target_samples = [
# "hg19-masked-ucsc",
# "hg19-p13-plusMT-masked-ucsc",
# "hg19-p13-no-alt-analysis-ucsc",
# "hg19-p13-full-analysis-ucsc",
# "hg19-initial-ucsc",
# "hg19-p13-plusMT-ucsc",
# ]

# # Pre-filter the DataFrame
# pep_df = pep_df[
#     ((pep_df['sample_name_1'].isin(ucsc_target_samples)) & (pep_df['sample_name_2'].isin(ucsc_target_samples)))
# ]
#new_df = pep_df.copy()


# all_relevant_stats = ['overlap_coefficient_names',	
#                       'overlap_coefficient_lengths',
#                       'overlap_coefficient_sequences',
#                       'overlap_coeff_name_len',	
#                       'jaccard_names', 
#                       'jaccard_lengths', 
#                       'jaccard_similarity_weighted_length', ''
#                       'jaccard_sequences',
#                       'jaccard_name_len', 
#                       'f1_names',
#                       'f1_lengths',
#                       'f1_sequences',
#                       'f1_weighted_lengths',
#                       'f1_name_len_pairs'
#                       ]
#all_relevant_stats = ['jaccard_names']

# for stat in all_relevant_stats:
#     pep_df[stat] = pd.to_numeric(pep_df[stat], errors='coerce')
#     # 1. Get all unique sample names
#     all_samples = pd.concat([pep_df['sample_name_1'], pep_df['sample_name_2']]).unique()

#     # 2. Create an empty DataFrame for the heatmap, initialized with NaN
#     heatmap_data = pd.DataFrame(index=all_samples, columns=all_samples)

#     # 3. Iterate over all pairs of unique sample names (O(n^2))
#     for sample1 in all_samples:
#         for sample2 in all_samples:
#             # a. If sample1 and sample2 are the same, assign 1 to both (symmetric)
#             if sample1 == sample2:
#                 heatmap_data.loc[sample1, sample2] = 1.0
#             else:
#                 # b. Check if the combination exists in the original dataframe
#                 comparison = pep_df[
#                     ((pep_df['sample_name_1'] == sample1) & (pep_df['sample_name_2'] == sample2)) |
#                     ((pep_df['sample_name_1'] == sample2) & (pep_df['sample_name_2'] == sample1))
#                 ]
#                 if not comparison.empty:
#                     similarity_score = comparison[stat].iloc[0]
#                     heatmap_data.loc[sample1, sample2] = similarity_score
#                     heatmap_data.loc[sample2, sample1] = similarity_score



#     print(heatmap_data)
#     heatmap_data = heatmap_data.apply(pd.to_numeric, errors='coerce')

#     # 4. Create the heatmap
#     plt.figure(figsize=(18, 15))  # Adjust figure size as needed
#     sns.heatmap(heatmap_data, annot=True, cmap='viridis', fmt=".2f", linewidths=.5, cbar_kws={'label': stat},annot_kws={"size": 9},vmin=0.0, vmax=1.0)
#     plt.title(f'{stat} Heatmap - {species_title}')
#     plt.xticks(rotation=90)
#     plt.yticks(rotation=0)
#     plt.tight_layout()
#     #plt.show()
#     output_path = os.path.join(results_dir,stat)
#     plt.savefig(output_path, dpi=300, bbox_inches='tight')
#     print(f"Heatmap saved to: {output_path}")

#     #     # 1. Pivot the dataframe to create a matrix
    
#     # heatmap_data_asymmetrical = new_df.pivot(index='sample_name_1', columns='sample_name_2', values=stat)

#     # # 2. Create the heatmap
#     # plt.figure(figsize=(10, 8))  # Adjust figure size as needed
#     # sns.heatmap(heatmap_data_asymmetrical, annot=True, cmap='viridis', fmt=".2f", linewidths=.5, cbar_kws={'label': 'Jaccard Similarity'})

#     # # Optional: Add a title
#     # plt.title('Asymmetrical Jaccard Similarity Heatmap (Direct from Pivot)')

#     # # Optional: Rotate x-axis labels for better readability
#     # plt.xticks(rotation=90)
#     # plt.yticks(rotation=0)

#     # # Show the plot
#     # plt.tight_layout()  # Adjust layout to prevent labels from being cut off
#     # plt.show()

#     row_sums = heatmap_data.sum(axis=1)

#TODO Ensure row sums take direction into consideration
#     # Create row sum graph
#     # 2. Display the row sums
#     print(f"Sum of {stat} Scores for Each Sample:")
#     print(row_sums)
#     # 1. Sort the row sums (optional, for better visualization)
#     row_sums_sorted = row_sums.sort_values(ascending=False)  # Sort in descending order

#     # 2. Create the bar plot
#     plt.figure(figsize=(12, 8))
#     plt.bar(row_sums_sorted.index, row_sums_sorted.values)

#     # 3. Set labels and title
#     plt.xlabel("Sample")
#     plt.ylabel(f"Sum of {stat}")
#     plt.title("Total Similarity of Each Sample")
#     plt.xticks(rotation=90, ha="right", fontsize=8)  # Rotate sample names for readability
#     plt.tight_layout()

#     # 4. Show and save the plot
#     #plt.show()
#     output_path = os.path.join(results_dir,stat+"_row_sums")
#     plt.savefig(output_path, dpi=300, bbox_inches='tight')


# # Create Histogram
# for stat in all_relevant_stats:
#     pep_df[stat] = pd.to_numeric(pep_df[stat], errors='coerce')
#     similarity_series = pd.Series(pep_df[stat])
#     print(similarity_series)
#     # 2. Define the bin size
#     bin_width = 0.05  # You can adjust this value (e.g., 0.1, 0.02)

#     # 3. Create the bins for the histogram
#     bins = np.arange(0, 1.0 + bin_width, bin_width)

#     # 4. Create the frequency plot (histogram)
#     plt.figure(figsize=(10, 6))
#     n, bins, patches = plt.hist(similarity_series, bins=bins, edgecolor='black', alpha=0.7)

#     for i in range(len(n)):
#         count = int(n[i])  # Get the count for the current bin
#         if count !=0:
#             x_position = (bins[i] + bins[i+1]) / 2  # Center of the bar
#             y_position = n[i]  # Top of the bar
#             plt.text(x_position, y_position, str(count), ha='center', va='bottom')
    
#     # 5. Set the labels and title
#     plt.xlabel(f' {stat} (0 to 1.0)')
#     plt.ylabel('Count Frequency')
#     plt.title(f'Frequency Distribution of {stat} - {species_title}')
#     plt.xticks(np.arange(0, 1.0 + bin_width, bin_width))  # Set x-axis ticks
#     plt.grid(axis='y', linestyle='--')  # Add a grid for better readability

#     # 6. Show the plot
#     plt.tight_layout()
#     #plt.show()
#     output_path = os.path.join(results_dir,stat+'_histogram')
#     plt.savefig(output_path, dpi=300, bbox_inches='tight')


# create single column showing subsets based on overlap coefficients
# df = pep["_sample_df"]

# all_relevant_overlap_stats = ['overlap_coefficient_names',	
#                       'overlap_coefficient_lengths',
#                       'overlap_coefficient_sequences',
#                       'overlap_coeff_name_len',]
# #stat = 'overlap_coeff_name_len'

# for stat in all_relevant_overlap_stats[:]:
#     df = pep["_sample_df"].copy()
#     similarity_df_sorted = df.sort_values(by=stat, ascending=False).copy()
#     similarity_df_sorted = similarity_df_sorted[similarity_df_sorted[stat] >= 0.75].copy()
#     similarity_df_sorted = similarity_df_sorted.head(100).copy()
#     comparisons = []
#     highers = []
#     lowers = []
#     for index, row in similarity_df_sorted.iterrows():
#         sample1 = row['sample_name_1']
#         digest1 = row['digest1']
#         sample2 = row['sample_name_2']
#         digest2 = row['digest2']
#         similarity = row[stat]

#         # TODO need to have this not be dependent on local jsons, get from other PEP and pull from HPC storage?
#         if get_sequence_length_local(digest1) > get_sequence_length_local(digest2):
#             higher_value_sample = sample1
#             compared_sample = sample2
#         else:
#             higher_value_sample = sample2
#             compared_sample = sample1

#         comparison_string = compared_sample + ' ⊂ ' + higher_value_sample
#         comparisons.append(comparison_string)
#         highers.append(higher_value_sample)
#         lowers.append(compared_sample)

#     similarity_df_sorted['comparison'] = comparisons
#     similarity_df_sorted['highers'] = highers
#     similarity_df_sorted['lowers'] = lowers

#     heatmap_data_single_col = similarity_df_sorted.set_index('comparison')[[stat]]

#     grid_kws = {"width_ratios": (.8, .05), "wspace": 2.00} # Adjust width ratios and spacing
#     fig, (ax, cbar_ax) = plt.subplots(1, 2, figsize=(12, 15), gridspec_kw=grid_kws)

#     sns.heatmap(heatmap_data_single_col, annot=True, cmap='viridis', fmt=".2f", linewidths=.5, cbar=True, cbar_ax=cbar_ax, ax=ax, cbar_kws={'label': stat})

#     fig.suptitle(f'Overlap Coefficient per Comparison ({stat})', x=0.1, y=0.95, ha='left', va='top',fontweight='bold') # Use fig.suptitle
#     ax.set_ylabel('Subset',rotation=270, labelpad=20)
#     ax.set_xticks(ticks=[])
#     ax.set_yticks(range(len(comparisons)))
#     ax.set_yticklabels(similarity_df_sorted['lowers'].tolist(), rotation=0, ha='left')
#     ax.tick_params(axis='y', which='major', pad=200)
#     ax.text(0.5, 1.05, 'Left Ref ⊂ Right Ref', ha='center', va='top', transform=ax.transAxes,fontweight='bold')

#     ax2 = ax.twinx()  # Create a second y-axis that shares the same x-axis
#     ax2.set_yticks(ax.get_yticks())  # Ensure the ticks are aligned
#     ax2.set_yticklabels(similarity_df_sorted['highers'].tolist(), rotation=0, ha='left')
#     ax2.set_ylabel('Contained by', rotation=270, labelpad=20) # Adjust labelpad
#     #     ax.set_xticks(range(len(all_samples)))
# #     ax.set_yticks(range(len(all_samples)))

#     plt.yticks(rotation=0)
#     output_path = os.path.join(results_dir,stat+'_subset_mapping')
#     plt.savefig(output_path, dpi=300, bbox_inches='tight')


#### PLOT MULTIPLE SUBPLOTS

stats_groups = [['f1_names',
                      'f1_lengths',
                      'f1_sequences',
                      'f1_name_len_pairs',
                      ],['f10_names',
                      'f10_lengths',
                      'f10_sequences',
                      'f10_name_len_pairs',
                      ]]

for all_relevant_stats in stats_groups:
    num_plots = len(all_relevant_stats)
    fig, axes = plt.subplots(1, num_plots, figsize=(6 * num_plots, 6))

    for i, stat in enumerate(all_relevant_stats):
        pep_df[stat] = pd.to_numeric(pep_df[stat], errors='coerce')
        # 1. Get all unique sample names
        all_samples = pd.concat([pep_df['sample_name_1'], pep_df['sample_name_2']]).unique()
        if desired_order:
            all_samples = [sample for sample in desired_order if sample in all_samples]

        # 2. Create an empty DataFrame for the heatmap, initialized with NaN
        heatmap_data = pd.DataFrame(index=all_samples, columns=all_samples)

        for row_idx, sample1 in enumerate(all_samples):
            for col_idx, sample2 in enumerate(all_samples):
                if col_idx >= row_idx:  # Condition to select the upper triangle (including diagonal)
                    if sample1 == sample2:
                        heatmap_data.loc[sample1, sample2] = 1.0
                    else:
                        comparison = pep_df[
                            ((pep_df['sample_name_1'] == sample1) & (pep_df['sample_name_2'] == sample2)) |
                            ((pep_df['sample_name_1'] == sample2) & (pep_df['sample_name_2'] == sample1))
                        ]
                        if not comparison.empty:
                            similarity_score = comparison[stat].iloc[0]
                            heatmap_data.loc[sample1, sample2] = similarity_score
                        else:
                            heatmap_data.loc[sample1, sample2] = np.nan
                else:
                    heatmap_data.loc[sample1, sample2] = np.nan



        #print(heatmap_data)
        heatmap_data = heatmap_data.apply(pd.to_numeric, errors='coerce')
        ax = axes[i]
        # 4. Create the heatmap
        #plt.figure(figsize=(18, 15))  # Adjust figure size as needed
        sns.heatmap(heatmap_data, annot=False, cmap='viridis', fmt=".2f", linewidths=.2, cbar_kws={'label': stat},annot_kws={"size": 3},vmin=0.0, vmax=1.0, ax=ax)
        ax.set_title(f'{stat} Heatmap')
        ax.set_xticks(range(len(all_samples)))
        ax.set_yticks(range(len(all_samples)))
        if i==0:
            ax.set_yticklabels(all_samples, rotation=0, fontsize=2)
        ax.tick_params(axis='both', which='major', labelsize=8)
        ax.set_xticklabels(all_samples, rotation=90, fontsize=2)
        #plt.tight_layout()
        #plt.show()
        # output_path = os.path.join(results_dir,stat)
        # #plt.savefig(output_path, dpi=300, bbox_inches='tight')
        # print(f"Heatmap saved to: {output_path}")
    plt.suptitle(f'Comparison Heatmaps - {species_title}', fontsize=16, y=1.02) # Add a suptitle for the entire figure
    plt.tight_layout(rect=[0, 0, 1, 0.96]) # Adjust layout to make space for suptitle
    #plt.show()
    output_path = os.path.join(results_dir,stat)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')



# # Another way to approach row sums
# for stat in all_relevant_stats:
#     #pep_df[stat] = pd.to_numeric(pep_df[stat], errors='coerce')
#     # 1. Get all unique sample names
#     #all_samples = pd.concat([pep_df['sample_name_1'], pep_df['sample_name_2']]).unique()

#     # 2. Create an empty DataFrame for the heatmap, initialized with NaN
#     #heatmap_data = pd.DataFrame(index=all_samples, columns=all_samples)

#     # 3. Iterate over all pairs of unique sample names (O(n^2))
#     # for sample1 in all_samples:
#     #     for sample2 in all_samples:
#     #         # a. If sample1 and sample2 are the same, assign 1 to both (symmetric)
#     #         if sample1 == sample2:
#     #             heatmap_data.loc[sample1, sample2] = 1.0
#     #         else:
#     #             # b. Check if the combination exists in the original dataframe
#     #             comparison = pep_df[
#     #                 ((pep_df['sample_name_1'] == sample1) & (pep_df['sample_name_2'] == sample2)) |
#     #                 ((pep_df['sample_name_1'] == sample2) & (pep_df['sample_name_2'] == sample1))
#     #             ]
#     #             if not comparison.empty:
#     #                 similarity_score = comparison[stat].iloc[0]
#     #                 heatmap_data.loc[sample1, sample2] = similarity_score
#     #                 heatmap_data.loc[sample2, sample1] = similarity_score
#     # for stat in all_relevant_overlap_stats[:]:
#     df = pep["_sample_df"].copy()
#     similarity_df_sorted = df.sort_values(by=stat, ascending=False).copy()
#     #similarity_df_sorted = similarity_df_sorted[similarity_df_sorted[stat]].copy()
#     comparisons = []
#     highers = []
#     lowers = []
#     for index, row in similarity_df_sorted.iterrows():
#         sample1 = row['sample_name_1']
#         digest1 = row['digest1']
#         sample2 = row['sample_name_2']
#         digest2 = row['digest2']
#         similarity = row[stat]

#         # TODO need to have this not be dependent on local jsons, get from other PEP and pull from HPC storage?
#         if get_sequence_length_local(digest1) > get_sequence_length_local(digest2):
#             higher_value_sample = sample1
#             compared_sample = sample2
#         else:
#             higher_value_sample = sample2
#             compared_sample = sample1

#         comparison_string = compared_sample + ' ⊂ ' + higher_value_sample
#         comparisons.append(comparison_string)
#         highers.append(higher_value_sample)
#         lowers.append(compared_sample)

#     similarity_df_sorted['comparison'] = comparisons
#     similarity_df_sorted['highers'] = highers
#     similarity_df_sorted['lowers'] = lowers

#     heatmap_data = similarity_df_sorted.set_index('highers')[[stat]]




# #     print(heatmap_data)
#     row_sums = heatmap_data.sum(axis=1)

# #TODO Ensure row sums take direction into consideration
#     # Create row sum graph
#     # 2. Display the row sums
#     print(f"Sum of {stat} Scores for Each Sample:")
#     print(row_sums)
#     # 1. Sort the row sums (optional, for better visualization)
#     row_sums_sorted = row_sums.sort_values(ascending=False)  # Sort in descending order

#     # 2. Create the bar plot
#     plt.figure(figsize=(12, 8))
#     plt.bar(row_sums_sorted.index, row_sums_sorted.values)

#     # 3. Set labels and title
#     plt.xlabel("Sample")
#     plt.ylabel(f"Sum of {stat}")
#     plt.title("Total Similarity of Each Sample")
#     plt.xticks(rotation=90, ha="right", fontsize=8)  # Rotate sample names for readability
#     plt.tight_layout()

#     # 4. Show and save the plot
#     #plt.show()
#     output_path = os.path.join(results_dir,stat+"_row_sums_highers")
#     plt.savefig(output_path, dpi=300, bbox_inches='tight')




#CREATE A SPECIAL LOOK AT JACCARD_SEQUENCES when JACCARD_NAME_LEN is perfect

# df = pep["_sample_df"]

# # all_relevant_overlap_stats = ['overlap_coefficient_names',	
# #                       'overlap_coefficient_lengths',
# #                       'overlap_coefficient_sequences',
# #                       'overlap_coeff_name_len',]
# # #stat = 'overlap_coeff_name_len'

# all_relevant_stats = ['jaccard_sequences']

# for stat in all_relevant_stats[:]:
#     df = pep["_sample_df"].copy()
#     similarity_df_sorted = df.sort_values(by=stat, ascending=False).copy()
#     similarity_df_sorted = similarity_df_sorted[similarity_df_sorted['jaccard_name_len'] >= .60].copy()

#     # 1. Get all unique sample names
#     all_samples = pd.concat([pep_df['sample_name_1'], pep_df['sample_name_2']]).unique()

#     # 2. Create an empty DataFrame for the heatmap, initialized with NaN
#     heatmap_data = pd.DataFrame(index=all_samples, columns=all_samples, dtype='float')

#     for row_idx, sample1 in enumerate(all_samples):
#         for col_idx, sample2 in enumerate(all_samples):
#             if col_idx >= row_idx:  # Condition to select the upper triangle (including diagonal)
#                 if sample1 == sample2:
#                     heatmap_data.loc[sample1, sample2] = 1.0
#                 else:
#                     comparison = similarity_df_sorted[
#                         ((similarity_df_sorted['sample_name_1'] == sample1) & (similarity_df_sorted['sample_name_2'] == sample2)) |
#                         ((similarity_df_sorted['sample_name_1'] == sample2) & (similarity_df_sorted['sample_name_2'] == sample1))
#                     ]
#                     if not comparison.empty:
#                         similarity_score = comparison[stat].iloc[0]
#                         heatmap_data.loc[sample1, sample2] = similarity_score
#                     else:
#                         heatmap_data.loc[sample1, sample2] = np.nan
#             else:
#                 heatmap_data.loc[sample1, sample2] = np.nan


#     #heatmap_data_single_col = similarity_df_sorted.set_index('comparison')[[stat]]

#     grid_kws = {"width_ratios": (.8, .05), "wspace": .5} # Adjust width ratios and spacing
#     fig, (ax, cbar_ax) = plt.subplots(1, 2, figsize=(15, 15), gridspec_kw=grid_kws)

#     sns.heatmap(heatmap_data, annot=True, cmap='viridis', fmt=".2f", linewidths=.5, cbar=True, cbar_ax=cbar_ax, ax=ax, cbar_kws={'label': stat},annot_kws={"size": 8},vmin=0.0, vmax=1.0)

#     fig.suptitle(f'({stat}) where jaccard_name_len > 0.60', x=0.5, y=0.95, ha='left', va='top',fontweight='bold') # Use fig.suptitle
#     #ax.set_ylabel('Subset',rotation=270, labelpad=20)
#     # ax.set_xticks(ticks=[])
#     # ax.set_yticks(range(len(similarity_df_sorted)))
#     # ax.set_xticklabels(rotation=0, ha='left')
#     # ax.set_yticklabels(rotation=0, ha='left')
#     # ax.tick_params(axis='y', which='major', pad=200)
#     # ax.text(0.5, 1.05, 'Left Ref ⊂ Right Ref', ha='center', va='top', transform=ax.transAxes,fontweight='bold')
#     # output_path = os.path.join(results_dir,stat+'_vs_name_len1')
#     # plt.savefig(output_path, dpi=300, bbox_inches='tight')
#     plt.show()

# PLOT MOW MEDIAN CHANGES BASED ON JACCARD_NAME_LEN
df = pep_df


stats = ["Mean", "Median"]

for stat in stats:
    # --- Calculate and Plot Median ---
    plt.figure(figsize=(10, 6))  # Adjust figure size as needed

    # Create a list of unique jaccard_name_len values
    #unique_jaccard_lengths = sorted(df['jaccard_name_len'].unique())
    unique_jaccard_lengths = [0.0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0]

    # Calculate and store medians for each jaccard_name_len
    calculated_stats = []
    counts = []  # Store the counts for each bin
    for length in unique_jaccard_lengths:
        bin_data = df[(df['jaccard_name_len'] > (length - 0.05)) & (df['jaccard_name_len'] <= (length + 0.05))]['jaccard_sequences']
        if stat == "Median":
            calculated_stats.append(bin_data.median())
        elif stat == "Mean":
            calculated_stats.append(bin_data.mean())
        counts.append(len(bin_data))  # Count the data points
                           
    # Create the bar plot
    plt.bar(unique_jaccard_lengths, calculated_stats, width=0.08) # Added marker and linestyle
        # Add count labels
    for i, count in enumerate(counts):
        plt.text(unique_jaccard_lengths[i], calculated_stats[i], str(count), ha='center', va='bottom')

    # Add labels and title
    plt.xlabel('Jaccard Similarity of Name Len Pairs')
    plt.ylabel(f'{stat} Jaccard Similarity of Sequences')
    plt.title(f'{stat} Jaccard Sequence Similarity vs. Jaccard Name Length Pair Similarity')
    plt.xticks(unique_jaccard_lengths)  # Ensure all unique lengths are shown on x-axis

    plt.tight_layout()
    output_path = os.path.join(results_dir, 'jaccard_sequence_'+stat)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    #plt.show()

#### PLOT OVERLAP_COEFFICIENT vs OVERLAP_MAX
relevant_stats = [("overlap_coefficient_sequences", "overlap_coefficient_sequences_max")]
#relevant_stats = [("overlap_coefficient_sequences", "overlap_coefficient_sequences_max"),("overlap_coefficient_sequences", "jaccard_sequences")]



for stat in relevant_stats:

    fig, ax = plt.subplots(figsize=(14, 12))  # Create a single subplot

    bottom_stat, top_stat = stat  # Get the stats

    pep_df[bottom_stat] = pd.to_numeric(pep_df[bottom_stat], errors='coerce')
    pep_df[top_stat] = pd.to_numeric(pep_df[top_stat], errors='coerce')

    # 1. Get all unique sample names
    all_samples = pd.concat([pep_df['sample_name_1'], pep_df['sample_name_2']]).unique()

    if desired_order:
        all_samples = [sample for sample in desired_order if sample in all_samples]

    num_samples = len(all_samples)

    # 2. Create an empty DataFrame for the combined heatmap
    combined_heatmap_data = pd.DataFrame(index=all_samples, columns=all_samples)

    for row_idx, sample1 in enumerate(all_samples):
        for col_idx, sample2 in enumerate(all_samples):
            if col_idx >= row_idx:  # Upper triangle (including diagonal)
                if sample1 == sample2:
                    combined_heatmap_data.loc[sample1, sample2] = np.nan  # White diagonal
                else:
                    comparison = pep_df[
                        ((pep_df['sample_name_1'] == sample1) & (pep_df['sample_name_2'] == sample2)) |
                        ((pep_df['sample_name_1'] == sample2) & (pep_df['sample_name_2'] == sample1))
                    ]
                    if not comparison.empty:
                        similarity_score = comparison[top_stat].iloc[0]
                        combined_heatmap_data.loc[sample1, sample2] = similarity_score
                    else:
                        combined_heatmap_data.loc[sample1, sample2] = np.nan
            else:  # Lower triangle
                if sample1 == sample2:
                    combined_heatmap_data.loc[sample1, sample2] = 1.0  # Different diagonal value if needed
                else:
                    comparison = pep_df[
                        ((pep_df['sample_name_1'] == sample1) & (pep_df['sample_name_2'] == sample2)) |
                        ((pep_df['sample_name_1'] == sample2) & (pep_df['sample_name_2'] == sample1))
                    ]
                    if not comparison.empty:
                        similarity_score = comparison[bottom_stat].iloc[0]
                        combined_heatmap_data.loc[sample1, sample2] = similarity_score
                    else:
                        combined_heatmap_data.loc[sample1, sample2] = np.nan

    combined_heatmap_data = combined_heatmap_data.apply(pd.to_numeric, errors='coerce')

    # Create the heatmap
    sns.heatmap(combined_heatmap_data, annot=False, cmap='viridis', fmt=".2f", linewidths=.2,
                cbar_kws={'label': f'{bottom_stat} (Lower), {top_stat} (Upper)'},
                annot_kws={"size": 3}, vmin=0.0, vmax=1.0, ax=ax)
    ax.set_title(f'Combined Comparison Heatmap: {bottom_stat} (Lower), {top_stat} (Upper)')

    # Center the ticks
    ax.set_xticks(np.arange(num_samples) + 0.5)
    ax.set_yticks(np.arange(num_samples) + 0.5)

    # Set the labels to correspond to the new tick positions
    ax.set_xticklabels(all_samples, rotation=90, fontsize=8, ha='center')
    ax.set_yticklabels(all_samples, rotation=0, fontsize=8, va='center')
    ax.tick_params(axis='both', which='major', labelsize=8)

    plt.suptitle(f'Comparison Heatmap - {species_title}', fontsize=16, y=1.02)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    output_path = os.path.join(results_dir, f'combined_overlap_heatmap'+bottom_stat+top_stat)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()

# Plot Hexbin ScatterPlot

# Assuming pep_df, results_dir, and species_title are defined elsewhere

relevant_stats = [("f10_sequences", "f10_name_len_pairs")]
# relevant_stats = [("overlap_coefficient_sequences", "overlap_coefficient_sequences_max"),("overlap_coefficient_sequences", "jaccard_sequences")]

for stat in relevant_stats:
    fig, ax = plt.subplots(figsize=(10, 8))  # Adjust figure size

    bottom_stat, top_stat = stat  # Get the stats

    pep_df[bottom_stat] = pd.to_numeric(pep_df[bottom_stat], errors='coerce')
    pep_df[top_stat] = pd.to_numeric(pep_df[top_stat], errors='coerce')

    # Create lists to store the paired data points
    x_values = []
    y_values = []

    # Iterate through all unique sample pairs (avoiding duplicates and self-comparisons)
    all_samples = pd.concat([pep_df['sample_name_1'], pep_df['sample_name_2']]).unique()

    for i in range(len(all_samples)):
        for j in range(i + 1, len(all_samples)):
            sample1 = all_samples[i]
            sample2 = all_samples[j]

            comparison = pep_df[
                ((pep_df['sample_name_1'] == sample1) & (pep_df['sample_name_2'] == sample2)) |
                ((pep_df['sample_name_1'] == sample2) & (pep_df['sample_name_2'] == sample1))
            ]

            if not comparison.empty:
                x_values.append(comparison[bottom_stat].iloc[0])
                y_values.append(comparison[top_stat].iloc[0])

    # Create the scatter plot
    

    # Create the hexbin plot on top (optional, but shows density)
    hb = ax.hexbin(x_values, y_values, gridsize=10, cmap='plasma', alpha=0.75, label='Density')  # Adjust gridsize and alpha

    ax.scatter(x_values, y_values, alpha=1.0, label='Data Points',color='skyblue')  # Adjust alpha for transparency

    # Add a colorbar for the hexbin plot
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('Count')

    # Set labels and title
    ax.set_xlabel(bottom_stat)
    ax.set_ylabel(top_stat)
    ax.set_title(f'Scatter and Hexbin Plot of {bottom_stat} vs {top_stat}')
    ax.legend()  # Show the legend

    plt.tight_layout()
    output_path = os.path.join(results_dir, f'scatter_hexbin_{bottom_stat}_vs_{top_stat}')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()



# Plot bar graph freq distribution
relevant_stats = [("f10_sequences", "f10_name_len_pairs")]
for stat_pair in relevant_stats:
    fig, ax = plt.subplots(figsize=(10, 6))  # Adjust figure size

    bottom_stat, top_stat = stat_pair

    pep_df[bottom_stat] = pd.to_numeric(pep_df[bottom_stat], errors='coerce')
    pep_df[top_stat] = pd.to_numeric(pep_df[top_stat], errors='coerce')

    # Collect all values for both statistics
    bottom_values = []
    top_values = []
    all_samples = pd.concat([pep_df['sample_name_1'], pep_df['sample_name_2']]).unique()

    for i in range(len(all_samples)):
        for j in range(i + 1, len(all_samples)):
            sample1 = all_samples[i]
            sample2 = all_samples[j]

            comparison = pep_df[
                ((pep_df['sample_name_1'] == sample1) & (pep_df['sample_name_2'] == sample2)) |
                ((pep_df['sample_name_1'] == sample2) & (pep_df['sample_name_2'] == sample1))
            ]

            if not comparison.empty:
                bottom_values.append(comparison[bottom_stat].iloc[0])
                top_values.append(comparison[top_stat].iloc[0])

    # Define bins for the frequency plot
    bins = np.linspace(0, 1, 21)  # Create 20 bins from 0 to 1

    # Calculate frequencies for both statistics
    freq_bottom, _ = np.histogram(bottom_values, bins=bins)
    freq_top, _ = np.histogram(top_values, bins=bins)

    # Set the width of the bars
    width = 0.35

    # Set the positions of the bars on the x-axis
    x = np.arange(len(freq_bottom))

    # Create the bar plot
    rects1 = ax.bar(x - width/2, freq_bottom, width, label=bottom_stat, color='skyblue')
    rects2 = ax.bar(x + width/2, freq_top, width, label=top_stat, color='salmon')

    # Add labels, title, and legend
    ax.set_xlabel('Overlap Coefficient')
    ax.set_ylabel('Frequency')
    ax.set_title(f'Frequency Distribution of {bottom_stat} and {top_stat}')
    ax.set_xticks(x)
    ax.set_xticklabels([f'{b:.2f}-{(b + (bins[1] - bins[0])):.2f}' for b in bins[:-1]], rotation=45, ha='right')
    ax.legend()

    fig.tight_layout()
    output_path = os.path.join(results_dir, f'frequency_plot_{bottom_stat}_vs_{top_stat}')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()