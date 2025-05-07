import sys
import os
import json
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from pephubclient import PEPHubClient
from refget import fasta_to_digest, fasta_to_seqcol_dict, compare_seqcols, SequenceCollection
from itertools import combinations
from pprint import pprint


#results_pep = "donaldcampbelljr/test_seq_col_results:default"

results_dir = sys.argv[1] # output dir for graphs 
results_pep = sys.argv[2] # input pep for graphing
species_title = sys.argv[3] # additional information for the title


LOCAL_JSON_DIRECTORY = "/home/drc/Downloads/jsons_from_rivanna/json/"


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


def get_sequence_length_local(digest):

    # OPENS A LOCAL JSON ONLY
    json_fp_1 = os.path.join(LOCAL_JSON_DIRECTORY, digest+".json")
    with open(json_fp_1, "r") as f:
        reloaded_dict1 = json.load(fp=f)
    
    return len(reloaded_dict1['sorted_sequences'])



# --------

phc = PEPHubClient()
pep = phc.load_project(results_pep)
print(pep["_sample_df"])

pep_df = pep["_sample_df"]


# FOR ORDERING OR SLECTING ONLY SOME SAMPLES
# --------------------------------------------
desired_order = None

# desired_order = [
# "GRCh38.p0-fasta-genomic",
# "GRCh38.p1-fasta-genomic",
# "GRCh38.p2-fasta-genomic",
# "GRCh38.p6-fasta-genomic",
# "GRCh38.p7-fasta-genomic",
# "GRCh38.p8-fasta-genomic",
# "GRCh38.p12-fasta-genomic",
# "GRCh38.p13-fasta-genomic",
# "GRCh38.p14-fasta-genomic",
    
# ]

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

# target_samples = [

# "hg38-primary-113-ensembl",
# "GRCh38-primary-assembly-47-gencode",
# "GRCh38.p14-fasta-genomic",
# "hg38-p14-ucsc",

# ]

# # # Pre-filter the DataFrame
# pep_df = pep_df[
#     ((pep_df['sample_name_1'].isin(target_samples)) & (pep_df['sample_name_2'].isin(target_samples)))
# ]
#new_df = pep_df.copy()
# --------------------------------------------


# BEGIN FIGURES

# --------------------------------------------


all_relevant_stats = [
                    #   'overlap_coefficient_names',	
                    #   'overlap_coefficient_lengths',
                    #   'overlap_coefficient_sequences',
                    #   'overlap_coeff_name_len',	
                      'jaccard_names', 
                      'jaccard_lengths', 
                      'jaccard_sequences',
                      'jaccard_name_len', 
                    #   'f1_names',
                    #   'f1_lengths',
                    #   'f1_sequences',
                    #   'f1_weighted_lengths',
                    #   'f1_name_len_pairs'
                      ]



#### PLOT MULTIPLE SUBPLOTS

# stats_groups = [['f1_names',
#                       'f1_lengths',
#                       'f1_sequences',
#                       'f1_name_len_pairs',
#                       ],['f10_names',
#                       'f10_lengths',
#                       'f10_sequences',
#                       'f10_name_len_pairs',
#                       ]]
stats_groups = [['jaccard_names', 'jaccard_lengths', 'jaccard_sequences', 'jaccard_name_len']]

for all_relevant_stats in stats_groups:
    num_plots = len(all_relevant_stats)
    num_rows = (num_plots + 1) // 2
    fig, axes = plt.subplots(num_rows, 2, figsize=(14, 7 * num_rows), sharex=True, sharey=True)
    axes = np.ravel(axes)

    cbar_ax = fig.add_axes([0.92, 0.15, 0.03, 0.7])  # [left, bottom, width, height] for the colorbar

    for i, stat in enumerate(all_relevant_stats):
        pep_df[stat] = pd.to_numeric(pep_df[stat], errors='coerce')
        all_samples = pd.concat([pep_df['sample_name_1'], pep_df['sample_name_2']]).unique()
        if desired_order:
            all_samples = [sample for sample in desired_order if sample in all_samples]
        heatmap_data = pd.DataFrame(index=all_samples, columns=all_samples)
        for row_idx, sample1 in enumerate(all_samples):
            for col_idx, sample2 in enumerate(all_samples):
                if col_idx >= row_idx:
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
        heatmap_data = heatmap_data.apply(pd.to_numeric, errors='coerce')
        ax = axes[i]
        sns.heatmap(heatmap_data, annot=False, cmap='viridis', fmt=".2f", linewidths=.2, cbar=i == 0, cbar_ax=cbar_ax if i == 0 else None, annot_kws={"size": 3}, vmin=0.0, vmax=1.0, ax=ax)
        ax.set_title(f'{stat} Heatmap')
        ax.set_xticks(np.arange(0.5, len(all_samples), 1))
        ax.set_yticks(np.arange(0.5, len(all_samples), 1))
        ax.set_yticklabels(all_samples, rotation=0, fontsize=8)
        ax.tick_params(axis='both', which='major', labelsize=8)
        ax.set_xticklabels(all_samples, rotation=90, fontsize=8)
        if i == 0:
            cbar_ax.set_ylabel('Similarity Score', fontsize=12) # Set label only once

    if num_plots < num_rows * 2:
        for j in range(num_plots, num_rows * 2):
            fig.delaxes(axes[j])

    plt.suptitle(f'Comparison Heatmaps - {species_title}', fontsize=16, y=1.02)
    #plt.tight_layout(rect=[0, 0, 0.9, 0.96]) # Adjust layout to make space for the colorbar
    output_path = os.path.join(results_dir, f'stacked_heatmap_single_cbar_{"_".join(all_relevant_stats)}.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')


# PLOT OPA AND OPB

stats_groups = [['opa_names', 'opb_names', 'opa_lengths', 'opb_lengths', 'opa_sequences','opb_sequences','opa_name_len', 'opb_name_len']]

for all_relevant_stats in stats_groups:
    num_plots = len(all_relevant_stats)
    num_rows = (num_plots + 1) // 2
    fig, axes = plt.subplots(num_rows, 2, figsize=(14, 7 * num_rows), sharex=True, sharey=True)
    axes = np.ravel(axes)

    cbar_ax = fig.add_axes([0.92, 0.15, 0.03, 0.7])  # [left, bottom, width, height] for the colorbar

    for i, stat in enumerate(all_relevant_stats):
        pep_df[stat] = pd.to_numeric(pep_df[stat], errors='coerce')
        all_samples = pd.concat([pep_df['sample_name_1'], pep_df['sample_name_2']]).unique()
        if desired_order:
            all_samples = [sample for sample in desired_order if sample in all_samples]
        heatmap_data = pd.DataFrame(index=all_samples, columns=all_samples)
        for row_idx, sample1 in enumerate(all_samples):
            for col_idx, sample2 in enumerate(all_samples):
                if col_idx >= row_idx:
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
        heatmap_data = heatmap_data.apply(pd.to_numeric, errors='coerce')
        ax = axes[i]
        sns.heatmap(heatmap_data, annot=False, cmap='viridis', fmt=".2f", linewidths=.2, cbar=i == 0, cbar_ax=cbar_ax if i == 0 else None, annot_kws={"size": 3}, vmin=0.0, vmax=1.0, ax=ax)
        ax.set_title(f'{stat} Heatmap')
        ax.set_xticks(np.arange(0.5, len(all_samples), 1))
        ax.set_yticks(np.arange(0.5, len(all_samples), 1))
        ax.set_yticklabels(all_samples, rotation=0, fontsize=8)
        ax.tick_params(axis='both', which='major', labelsize=8)
        ax.set_xticklabels(all_samples, rotation=90, fontsize=8)
        if i == 0:
            cbar_ax.set_ylabel('Similarity Score', fontsize=12) # Set label only once

    if num_plots < num_rows * 2:
        for j in range(num_plots, num_rows * 2):
            fig.delaxes(axes[j])

    plt.suptitle(f'Comparison Heatmaps - {species_title}', fontsize=16, y=1.02)
    #plt.tight_layout(rect=[0, 0, 0.9, 0.96]) # Adjust layout to make space for the colorbar
    output_path = os.path.join(results_dir, f'stacked_heatmap_single_cbar_{"_".join(all_relevant_stats)}.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')






#CREATE A SPECIAL LOOK AT JACCARD_SEQUENCES when JACCARD_NAME_LEN is perfect

df = pep["_sample_df"].copy()

all_relevant_stats = ['jaccard_sequences']

CUTOFF = .60

for stat in all_relevant_stats[:]:
    df = pep["_sample_df"].copy()
    similarity_df_sorted = df.sort_values(by=stat, ascending=False).copy()
    similarity_df_sorted = similarity_df_sorted[similarity_df_sorted['jaccard_name_len'] >= CUTOFF].copy()

    # 1. Get all unique sample names
    all_samples = pd.concat([pep_df['sample_name_1'], pep_df['sample_name_2']]).unique()

    # 2. Create an empty DataFrame for the heatmap, initialized with NaN
    heatmap_data = pd.DataFrame(index=all_samples, columns=all_samples, dtype='float')

    for row_idx, sample1 in enumerate(all_samples):
        for col_idx, sample2 in enumerate(all_samples):
            if col_idx >= row_idx:  # Condition to select the upper triangle (including diagonal)
                if sample1 == sample2:
                    heatmap_data.loc[sample1, sample2] = 1.0
                else:
                    comparison = similarity_df_sorted[
                        ((similarity_df_sorted['sample_name_1'] == sample1) & (similarity_df_sorted['sample_name_2'] == sample2)) |
                        ((similarity_df_sorted['sample_name_1'] == sample2) & (similarity_df_sorted['sample_name_2'] == sample1))
                    ]
                    if not comparison.empty:
                        similarity_score = comparison[stat].iloc[0]
                        heatmap_data.loc[sample1, sample2] = similarity_score
                    else:
                        heatmap_data.loc[sample1, sample2] = np.nan
            else:
                heatmap_data.loc[sample1, sample2] = np.nan

    grid_kws = {"width_ratios": (.8, .05), "wspace": .5} # Adjust width ratios and spacing
    fig, (ax, cbar_ax) = plt.subplots(1, 2, figsize=(15, 15), gridspec_kw=grid_kws)

    sns.heatmap(heatmap_data, annot=True, cmap='viridis', fmt=".2f", linewidths=.5, cbar=True, cbar_ax=cbar_ax, ax=ax, cbar_kws={'label': stat},annot_kws={"size": 8},vmin=0.0, vmax=1.0)

    fig.suptitle(f'({stat}) where jaccard_name_len > {CUTOFF}', x=0.5, y=0.95, ha='left', va='top',fontweight='bold') # Use fig.suptitle

    output_path = os.path.join(results_dir,stat+'_vs_name_len1')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    #plt.show()

# PLOT MOW MEDIAN CHANGES BASED ON JACCARD_NAME_LEN

# TODO PUT COUNTS ALL AT THE TOP OR THE BOTTOM
df = pep_df.copy()

stats = ["Mean", "Median"]

for stat in stats:
    # --- Calculate and Plot Statistic ---
    plt.figure(figsize=(10, 6))  # Adjust figure size as needed

    unique_jaccard_lengths = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    calculated_stats = []
    counts = []  # Store the counts for each bin

    for length in unique_jaccard_lengths:
        lower_bound = length - 0.05
        upper_bound = length + 0.05
        bin_data = df[(df['jaccard_name_len'] > lower_bound) & (df['jaccard_name_len'] <= upper_bound)]['jaccard_sequences']
        if stat == "Median":
            calculated_stats.append(bin_data.median())
        elif stat == "Mean":
            calculated_stats.append(bin_data.mean())
        counts.append(len(bin_data))  # Count the data points

    # Create the bar plot
    bars = plt.bar(unique_jaccard_lengths, calculated_stats, width=0.08)

    # Add count labels at the bottom of the bars
    for bar, count in zip(bars, counts):
        yval = 0.05  # Position the text at the bottom (y=0)
        plt.text(bar.get_x() + bar.get_width()/2, yval, "n="+str(count), ha='center', va='top',color='white', fontweight='bold')

    # Add labels and title
    plt.xlabel('Jaccard Similarity of Name Len Pairs')
    plt.ylabel(f'{stat} Jaccard Similarity of Sequences')
    plt.title(f'{stat} Jaccard Sequence Similarity vs. Jaccard Name Length Pair Similarity')
    plt.xticks(unique_jaccard_lengths)  # Ensure all unique lengths are shown on x-axis
    plt.ylim(bottom=0) # Ensure the y-axis starts at 0 to accommodate the labels

    plt.tight_layout()
    output_path = os.path.join(results_dir, f'{stat}_counts')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')




# #### PLOT OVERLAP_COEFFICIENT vs OVERLAP_MAX

#TODO CHANGE STATS, PLOT IN SEPARATE FIGURES FOR NOW

# relevant_stats = [("overlap_coefficient_sequences", "overlap_coefficient_sequences_max")]
# #relevant_stats = [("overlap_coefficient_sequences", "overlap_coefficient_sequences_max"),("overlap_coefficient_sequences", "jaccard_sequences")]

# for stat in relevant_stats:

#     fig, ax = plt.subplots(figsize=(14, 12))  # Create a single subplot

#     bottom_stat, top_stat = stat  # Get the stats

#     pep_df[bottom_stat] = pd.to_numeric(pep_df[bottom_stat], errors='coerce')
#     pep_df[top_stat] = pd.to_numeric(pep_df[top_stat], errors='coerce')

#     # 1. Get all unique sample names
#     all_samples = pd.concat([pep_df['sample_name_1'], pep_df['sample_name_2']]).unique()

#     if desired_order:
#         all_samples = [sample for sample in desired_order if sample in all_samples]

#     num_samples = len(all_samples)

#     # 2. Create an empty DataFrame for the combined heatmap
#     combined_heatmap_data = pd.DataFrame(index=all_samples, columns=all_samples)

#     for row_idx, sample1 in enumerate(all_samples):
#         for col_idx, sample2 in enumerate(all_samples):
#             if col_idx >= row_idx:  # Upper triangle (including diagonal)
#                 if sample1 == sample2:
#                     combined_heatmap_data.loc[sample1, sample2] = np.nan  # White diagonal
#                 else:
#                     comparison = pep_df[
#                         ((pep_df['sample_name_1'] == sample1) & (pep_df['sample_name_2'] == sample2)) |
#                         ((pep_df['sample_name_1'] == sample2) & (pep_df['sample_name_2'] == sample1))
#                     ]
#                     if not comparison.empty:
#                         similarity_score = comparison[top_stat].iloc[0]
#                         combined_heatmap_data.loc[sample1, sample2] = similarity_score
#                     else:
#                         combined_heatmap_data.loc[sample1, sample2] = np.nan
#             else:  # Lower triangle
#                 if sample1 == sample2:
#                     combined_heatmap_data.loc[sample1, sample2] = 1.0  # Different diagonal value if needed
#                 else:
#                     comparison = pep_df[
#                         ((pep_df['sample_name_1'] == sample1) & (pep_df['sample_name_2'] == sample2)) |
#                         ((pep_df['sample_name_1'] == sample2) & (pep_df['sample_name_2'] == sample1))
#                     ]
#                     if not comparison.empty:
#                         similarity_score = comparison[bottom_stat].iloc[0]
#                         combined_heatmap_data.loc[sample1, sample2] = similarity_score
#                     else:
#                         combined_heatmap_data.loc[sample1, sample2] = np.nan

#     combined_heatmap_data = combined_heatmap_data.apply(pd.to_numeric, errors='coerce')

#     # Create the heatmap
#     sns.heatmap(combined_heatmap_data, annot=False, cmap='viridis', fmt=".2f", linewidths=.2,
#                 cbar_kws={'label': f'{bottom_stat} (Lower), {top_stat} (Upper)'},
#                 annot_kws={"size": 9}, vmin=0.0, vmax=1.0, ax=ax)
#     ax.set_title(f'Combined Comparison Heatmap: {bottom_stat} (Lower), {top_stat} (Upper)')

#     # Center the ticks
#     ax.set_xticks(np.arange(num_samples) + 0.5)
#     ax.set_yticks(np.arange(num_samples) + 0.5)

#     # Set the labels to correspond to the new tick positions
#     ax.set_xticklabels(all_samples, rotation=90, fontsize=8, ha='center')
#     ax.set_yticklabels(all_samples, rotation=0, fontsize=8, va='center')
#     ax.tick_params(axis='both', which='major', labelsize=8)

#     plt.suptitle(f'Comparison Heatmap - {species_title}', fontsize=16, y=1.02)
#     plt.tight_layout(rect=[0, 0, 1, 0.96])
#     output_path = os.path.join(results_dir, f'combined_overlap_heatmap'+bottom_stat+top_stat)
#     plt.savefig(output_path, dpi=300, bbox_inches='tight')
#     plt.show()

# # Plot Hexbin ScatterPlot

# # Assuming pep_df, results_dir, and species_title are defined elsewhere

#relevant_stats = [('jaccard_sequences', 'jaccard_name_len')]
relevant_stats = ['jaccard_sequences', 'jaccard_name_len', 'jaccard_lengths', 'jaccard_names']
stat_combinations = combinations(iterable=relevant_stats,r=2)

for stat_combination in stat_combinations:
    fig, ax = plt.subplots(figsize=(10, 8))  # Adjust figure size

    bottom_stat, top_stat = stat_combination  # Get the stats

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

    # Define the monochromatic colormap
    colors = ["lightgrey", "darkorange"]
    cmap = LinearSegmentedColormap.from_list("grey_to_orange", colors)

    # Create the hexbin plot with linear scaling first
    hb = ax.hexbin(x_values, y_values, gridsize=25, cmap=cmap, alpha=0.8)
    ax.scatter(x_values, y_values, alpha=.20, label='Data Points',color='darkblue', s=25)
    counts = hb.get_array()

    # Apply logarithmic scaling to the counts, handling zeros
    log_counts = np.log1p(counts)

    # Normalize the logarithmic counts for the colormap
    norm = plt.Normalize(vmin=log_counts.min(), vmax=log_counts.max())
    colored_values = cmap(norm(log_counts))

    # Set the facecolors of the hexbins based on the normalized log counts
    hb.set_array(log_counts)  # Directly set the array to the log-transformed counts
    hb.set_norm(norm)        # Set the normalization for the colormap

    # Add a colorbar for the logarithmic scale
    cb = fig.colorbar(hb, ax=ax, ticks=np.log1p([1, 10, 100, 1000]))  # Example ticks
    cb.set_label('Counts (Log Scale)')
    cb.ax.set_yticklabels([1, 10, 100, 1000])  # Label the ticks

    # Set labels and title
    ax.set_xlabel(bottom_stat)
    ax.set_ylabel(top_stat)
    ax.set_title(f'Hexbin Plot of {bottom_stat} vs {top_stat} (Log Scale)')

    plt.tight_layout()
    output_path = os.path.join(results_dir, f'hexbin_monochromatic_logscale_{bottom_stat}_vs_{top_stat}')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')

# Plot bar graph freq distribution
#TODO change stats
#relevant_stats = [('jaccard_sequences', 'jaccard_name_len')]
relevant_stats = ['jaccard_sequences', 'jaccard_name_len', 'jaccard_lengths', 'jaccard_names']
stat_combinations = combinations(iterable=relevant_stats,r=2)
for stat_combination in stat_combinations:
    fig, ax = plt.subplots(figsize=(10, 6))  # Adjust figure size

    bottom_stat, top_stat = stat_combination
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
    #plt.show()


# PLOT DOT CHART
# TODO also plot bar graphs

def create_comparison_dot_plot(data, labels, title="Comparison Dot Plot", x_limit=10,
                               metric_colors=None, metric_markers=None):
    """
    Creates a comparison dot plot using Matplotlib.

    Args:
        data (dict): A dictionary where keys are sample names (e.g., 'Sample 1')
                      and values are lists of 4 values.
        labels (list): A list of label names (e.g., ['Metric A', 'Metric B', 'Metric C', 'Metric D']).
        title (str, optional): The title of the dot plot. Defaults to "Comparison Dot Plot".
        x_limit (int, optional): The maximum value for the x-axis. Defaults to 10.
        metric_colors (list, optional): A list of colors for each metric. If None, default colors are used.
        metric_markers (list, optional): A list of markers for each metric. If None, default markers are used.

    Returns:
        None: Displays the plot.
    """
    num_vars = len(labels)
    if num_vars != 4:
        raise ValueError("Number of labels must be 4 for this comparison dot plot.")

    fig, ax = plt.subplots(figsize=(8, 10))  # Adjust figure size as needed
    y_positions = np.arange(len(data))  # Create y positions for the samples

    # Default colors and markers
    default_colors = metric_colors = ['#1f77b4', '#ff7f0e', '#9467bd', '#8c564b']
    default_markers = ['o', 's', 'D', '^']
    sample_color = 'k'  # set a default sample color

    # Use provided colors and markers or defaults
    if metric_colors is None:
        metric_colors = default_colors
    if metric_markers is None:
        metric_markers = default_markers

    for i, (sample_name, values) in enumerate(data.items()):
        if len(values) != num_vars:
            raise ValueError(f"Sample '{sample_name}' must have 4 values.")
        ax.plot(values, [i] * num_vars, linestyle='-', color='k', alpha=0.3)  # Connect points with a line
        for j, value in enumerate(values):
            ax.plot(value, i, marker=metric_markers[j],
                    markersize=8, alpha=0.7, label=labels[j], color=metric_colors[j])  # Use color and marker

    ax.set_yticks(y_positions)
    ax.set_yticklabels(list(data.keys()))  # Set sample names as y-axis labels
    ax.set_xlim(0, x_limit)  # Set x-axis limits
    #ax.set_xlabel("Metric Value")
    ax.set_ylabel("Ref Genomes")
    ax.set_title(title, fontsize=14)
    # Create a single legend for the samples
    handles, labels = ax.get_legend_handles_labels()
    unique_labels = list(dict.fromkeys(labels))  # remove duplicate labels
    unique_handles = [handles[labels.index(label)] for label in unique_labels]
    ax.legend(unique_handles, unique_labels, loc='upper right')
    ax.grid(True, axis='x', linestyle='--', alpha=0.6)  # Add grid lines

    plt.tight_layout()
    #plt.show()



primary_samples = ["hg19-initial-ucsc", "GRCh38.p14-fasta-genomic"]

metrics = ['jaccard_names', 'jaccard_lengths', 'jaccard_sequences', 'jaccard_name_len']


for primary_sample in primary_samples:

    new_df = pep_df.copy()

    # Filter the DataFrame
    filtered_df = new_df[(new_df['sample_name_1'] == primary_sample) | (new_df['sample_name_2'] == primary_sample)]

    # Convert the filtered DataFrame to the dictionary format
    data_from_df = {}
    for _, row in filtered_df.iterrows():
        if row['sample_name_1'] == primary_sample:
            sample_name = row['sample_name_2']
        else:
            sample_name = row['sample_name_1']
        if sample_name not in data_from_df:
            data_from_df[sample_name] = []
        data_from_df[sample_name] = [row[metrics[0]], row[metrics[1]], row[metrics[2]], row[metrics[3]]]

    #Use the columns f10_names, f10_lengths, f10_sequences, f10_name_len_pairs
    labels = metrics
    create_comparison_dot_plot(data_from_df, labels, title=f"Comparison of {primary_sample} vs All Other Ref Genomes", x_limit=1.0, metric_markers=['o', 's', 'D', '^'])
    
    
    output_path = os.path.join(results_dir, f'comparison_{primary_sample}_.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    #plt.show()


def create_comparison_bar_graph(data, labels, title="Comparison Bar Graph",
                               metric_colors=None):
    """
    Creates a grouped bar graph using Matplotlib.

    Args:
        data (dict): A dictionary where keys are sample names (e.g., 'Sample 1')
                      and values are lists of 4 values.
        labels (list): A list of label names (e.g., ['Metric A', 'Metric B', 'Metric C', 'Metric D']).
        title (str, optional): The title of the bar graph. Defaults to "Comparison Bar Graph".
        metric_colors (list, optional): A list of colors for each metric. If None, default colors are used.

    Returns:
        None: Displays the plot.
    """
    num_metrics = len(labels)
    if num_metrics != 4:
        raise ValueError("Number of labels must be 4 for this comparison bar graph.")

    if metric_colors is None:
        metric_colors = ['#1f77b4', '#ff7f0e', '#9467bd', '#8c564b']  # Default colors

    fig, ax = plt.subplots(figsize=(10, 8))  # Adjust figure size as needed
    num_samples = len(data)
    bar_width = 0.2  # Width of each individual bar
    group_width = num_metrics * bar_width  # Total width of a group of bars
    index = np.arange(num_samples)  # Positions for the group centers

    for i, (sample_name, values) in enumerate(data.items()):
        if len(values) != num_metrics:
            raise ValueError(f"Sample '{sample_name}' must have {num_metrics} values.")
        for j, value in enumerate(values):
            x_position = index[i] + (j * bar_width) - (group_width / 2) + (bar_width / 2) # Calculate the x position for each bar
            ax.bar(x_position, value, bar_width, color=metric_colors[j], label=labels[j] if i == 0 else "")

    ax.set_xticks(index)
    ax.set_xticklabels(list(data.keys()), rotation=90, ha="right", fontsize=8)  # Sample names as x-axis labels
    ax.set_ylabel("Metric Value")
    ax.set_title(title, fontsize=14)
    ax.legend(loc='upper right')
    ax.grid(True, axis='y', linestyle='--', alpha=0.6)

    plt.tight_layout()
    #plt.show()

# --- Example Usage ---
# Assuming pep_df and results_dir are defined elsewhere

primary_samples = ["hg19-initial-ucsc", "GRCh38.p14-fasta-genomic"]
metrics = ['jaccard_names', 'jaccard_lengths', 'jaccard_sequences', 'jaccard_name_len']

for primary_sample in primary_samples:
    new_df = pep_df.copy()

    # Filter the DataFrame
    filtered_df = new_df[(new_df['sample_name_1'] == primary_sample) | (new_df['sample_name_2'] == primary_sample)]

    # Convert the filtered DataFrame to the dictionary format
    data_from_df = {}
    for _, row in filtered_df.iterrows():
        if row['sample_name_1'] == primary_sample:
            sample_name = row['sample_name_2']
        else:
            sample_name = row['sample_name_1']
        if sample_name not in data_from_df:
            data_from_df[sample_name] = []
        data_from_df[sample_name] = [row[metrics[0]], row[metrics[1]], row[metrics[2]], row[metrics[3]]]

    # Use the columns f10_names, f10_lengths, f10_sequences, f10_name_len_pairs
    labels = metrics
    create_comparison_bar_graph(data_from_df, labels, title=f"Comparison of {primary_sample} vs All Other Ref Genomes",
                                 )

    output_path = os.path.join(results_dir, f'comparison_bar_graph_{primary_sample}.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()