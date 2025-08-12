import sys
import os
import json
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from pephubclient import PEPHubClient
import pipestat
from refget import fasta_to_digest, fasta_to_seqcol_dict, compare_seqcols, SequenceCollection
from itertools import combinations
from pprint import pprint
#import matplotlib.font_manager as fm

# Ensure we default to Arial font, may need to delete font cache json list
# at location given by matplotlib.get_cachedir() so that it rebuilds the cache json
# more general installing arial on Ubuntu -> https://askubuntu.com/questions/1349836/how-to-install-fonts-in-20-04
import matplotlib
matplotlib.rcParams["svg.fonttype"] = "none" # do not embed directly, instead the downstream program will view with system fonts
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams["text.usetex"] = False # this should be default anyway, but ensures we are not ovverriding the above font choices.

#results_pep = "donaldcampbelljr/test_seq_col_results:default"

results_dir = sys.argv[1] # output dir for graphs 
results_pep = sys.argv[2] # input pep for graphing
species_title = sys.argv[3] # additional information for the title
LOCAL_JSON_DIRECTORY  = sys.argv[4] # for lcally stored jsons to pull sequence col data
PEPHUB_PATH_AUTHORITY = sys.argv[5] # a PEP that ties sample name to authority (e.g ncbi, ucsc)

# LOCAL_JSON_DIRECTORY = "/home/drc/Downloads/jsons_from_rivanna/json/"
# PEPHUB_PATH_AUTHORITY = "donaldcampbelljr/human_seqcol_digests:default" # this PEP associates digests/sample_names/authorities together
# PEPHUB_PATH_AUTHORITY = "donaldcampbelljr/ncbi_38_seqcol_digests:default" # this PEP associates digests/sample_names/authorities together

# results_dir = "/home/drc/Downloads/refgenomes_pics_test/jun112025/TEST/"
# #results_pep = "donaldcampbelljr/human_seq_col_results:default"
# results_pep = "donaldcampbelljr/ncbi_38_seqcol_results:default"
# species_title = "human"
# LOCAL_JSON_DIRECTORY = "/home/drc/Downloads/jsons_from_rivanna/json/"
# #PEPHUB_PATH_AUTHORITY = "donaldcampbelljr/human_seqcol_digests:default"
# PEPHUB_PATH_AUTHORITY = "donaldcampbelljr/ncbi_38_seqcol_digests:default"

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
target_samples = None

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

# target_samples = [
# "hg19-p13-plusMT-masked-ucsc",
# "hg19-p13-no-alt-analysis-ucsc",
# "hg19-p13-full-analysis-ucsc",
# "hg19-p13-plusMT-ucsc",
# ]

# target_samples = [

# "hg38-toplevel-113-ensembl",
# "GRCh38-p14-47-gencode",
# "GRCh38.p14-fasta-genomic",
# "hg38-p14-ucsc",

# ]


# target_samples = [

# "hg38-toplevel-113-ensembl",
# "GRCh38-p14-47-gencode",
# "GRCh38.p14-fasta-genomic",
# "hg38-p14-ucsc",
# "hg38-ddbj",
# "GRCh38-ena-29",

# ]

# target_samples =[

# "GRCh38.p14-fasta-no-alt-analysis",
# "GRCh38.p14-fasta-full-analysis-plus-hs38d1",
# "GRCh38.p14-fasta-no-alt-plus-hs38d1",
# "GRCh38.p14-fasta-genomic",
# "GRCh38.p14-fasta-full-analysis",
# "GRCh38.p14-fasta-no-alt-analysis",
# ]

#MOUSE TARGET SAMPLES

# target_samples=[
# "GRCm39-fasta-genomic",
# "GRCm39-toplevel-113-ensembl",
# "mm39-ucsc-initial-soft-masked",
# "GRCm39-all-M36-gencode",
# "GRCm39-ena-09",
# ]     


# # # # Pre-filter the DataFrame
if target_samples:
    pep_df = pep_df[
        ((pep_df['sample_name_1'].isin(target_samples)) & (pep_df['sample_name_2'].isin(target_samples)))
    ]
else:
    pass
#new_df = pep_df.copy()
# --------------------------------------------


# BEGIN FIGURES

# --------------------------------------------


all_relevant_stats = [
                      'jaccard_names', 
                      'jaccard_lengths', 
                      'jaccard_sequences',
                      'jaccard_name_len', 
                      ]



all_samples = pd.concat([pep_df['sample_name_1'], pep_df['sample_name_2']]).unique()
# get authority from another project based on sample_name:
psm = pipestat.PipestatManager(pephub_path=PEPHUB_PATH_AUTHORITY)
results = psm.select_records()

sample_authority = {}
for sample in all_samples:
    for result in results['records']:
        if sample == result['record_identifier']:
            sample_authority[sample] = result['authority']
            break # Assuming one record per sample

# # 1.  Extract sample authorities and create a dictionary
# sample_authority_dict = {}
# for result in results['records']:
#     sample_authority_dict[result['record_identifier']] = result['authority']

# 2. Convert the dictionary to a DataFrame
sample_authority_df = pd.DataFrame(list(sample_authority.items()), columns=['sample_name', 'authority'])

sample_authority_file = os.path.join(results_dir, "sample_authority.csv")

# 4. Save the DataFrame to a CSV file
sample_authority_df.to_csv(sample_authority_file, index=False)  # index=False prevents writing row numbers


#### PLOT MULTIPLE SUBPLOTS

stats_groups = [['jaccard_names', 'jaccard_lengths', 'jaccard_sequences', 'jaccard_name_len']]

for all_relevant_stats in stats_groups:
    num_plots = len(all_relevant_stats)
    num_rows = (num_plots + 1) // 2
    fig, axes = plt.subplots(num_rows, 2, figsize=(14, 7 * num_rows), sharex=True, sharey=True)
    axes = np.ravel(axes)

    cbar_ax = fig.add_axes([0.92, 0.15, 0.03, 0.7])  # [left, bottom, width, height] for the colorbar

    for i, stat in enumerate(all_relevant_stats):
        pep_df[stat] = pd.to_numeric(pep_df[stat], errors='coerce')
        if desired_order:
            all_samples = [sample for sample in desired_order if sample in all_samples]
                # Create a list of (sample, authority) tuples
        sample_authority_list = [(sample, sample_authority.get(sample)) for sample in all_samples]

        # Sort the list based on authority
        sorted_sample_authority = sorted(sample_authority_list, key=lambda item: item[1])

        # Extract the sorted sample names
        all_samples = [item[0] for item in sorted_sample_authority]

        heatmap_data = pd.DataFrame(index=all_samples, columns=all_samples)
        for row_idx, sample1 in enumerate(all_samples):
            for col_idx, sample2 in enumerate(all_samples):
                if col_idx <= row_idx:
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

                # Get the unique authorities in the sorted order of samples
        ordered_authorities = [sample_authority.get(sample) for sample in all_samples]
        ax.set_yticklabels(ordered_authorities, rotation=0, fontsize=8)
        ax.set_xticklabels(ordered_authorities, rotation=90, fontsize=8)

        output_path = os.path.join(results_dir, f'heatmap_data_{stat}.csv')
        heatmap_data.index.name = 'sample_name'  # Set the index name
        heatmap_data.to_csv(output_path,index=True, header=True) # Save the dataframe to a CSV

    if num_plots < num_rows * 2:
        for j in range(num_plots, num_rows * 2):
            fig.delaxes(axes[j])

    plt.suptitle(f'Comparison Heatmaps - {species_title}', fontsize=16, y=1.02)
    #plt.tight_layout(rect=[0, 0, 0.9, 0.96]) # Adjust layout to make space for the colorbar
    #output_path = os.path.join(results_dir, f'stacked_heatmap_single_cbar_{"_".join(all_relevant_stats)}.png')
    output_path = os.path.join(results_dir, f'stacked_heatmap_single_cbar_{"_".join(all_relevant_stats)}.svg')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


# PLOT OPA AND OPB

stats_pairs = [['opa_names', 'opb_names'], ['opa_lengths', 'opb_lengths'], ['opa_sequences', 'opb_sequences'], ['opa_name_len', 'opb_name_len']]
for pair_idx, all_relevant_stats in enumerate(stats_pairs): # Iterate through each pair
        num_plots = len(all_relevant_stats) # Should always be 2 for these pairs
        num_cols_subplot = 2 # Always 2 plots per figure
        num_rows_subplot = 1 # Always 1 row per figure for a pair of plots
        if desired_order:
            all_samples = [sample for sample in desired_order if sample in all_samples]
                # Create a list of (sample, authority) tuples
        sample_authority_list = [(sample, sample_authority.get(sample)) for sample in all_samples]

        # Sort the list based on authority
        sorted_sample_authority = sorted(sample_authority_list, key=lambda item: item[1])

        # Extract the sorted sample names
        all_samples = [item[0] for item in sorted_sample_authority]

        # Create a new figure for each pair of statistics
        fig, axes = plt.subplots(num_rows_subplot, num_cols_subplot, figsize=(14, 7), sharex=True, sharey=True)
        # Flatten axes array for easy iteration
        axes = np.ravel(axes)

        # Define colorbar axis only for the first plot in the pair
        cbar_ax = fig.add_axes([0.92, 0.15, 0.03, 0.7])  # [left, bottom, width, height] for the colorbar

        for i, stat in enumerate(all_relevant_stats):
            # Convert stat column to numeric, coercing errors to NaN
            # This needs to be done for each 'stat' in the current pair
            if pep_df[stat].dtype == 'object': # Only convert if it's not already numeric
                pep_df[stat] = pd.to_numeric(pep_df[stat], errors='coerce')

            # Create an empty DataFrame for the heatmap with sorted samples as index/columns
            heatmap_data = pd.DataFrame(index=all_samples, columns=all_samples, dtype=float)

            # Populate the heatmap data (lower triangle + diagonal)
            for row_idx, sample1 in enumerate(all_samples):
                for col_idx, sample2 in enumerate(all_samples):
                    if col_idx <= row_idx: # Populate only the lower triangle and diagonal
                        if sample1 == sample2:
                            heatmap_data.loc[sample1, sample2] = 1.0 # Diagonal value is 1.0
                        else:
                            # Find the comparison where sample1 and sample2 are involved
                            comparison_forward = pep_df[(pep_df['sample_name_1'] == sample1) & (pep_df['sample_name_2'] == sample2)]
                            comparison_backward = pep_df[(pep_df['sample_name_1'] == sample2) & (pep_df['sample_name_2'] == sample1)]

                            similarity_score = np.nan
                            # Determine the base stat name (e.g., 'names', 'lengths', 'name_len')
                            # Assumes 'opa_' or 'opb_' prefix is always 4 characters
                            base_stat_suffix = stat[4:]

                            if not comparison_forward.empty: # pep_df has (sample1, sample2)
                                if stat.startswith('opa_'):
                                    similarity_score = comparison_forward['opa_' + base_stat_suffix].iloc[0]
                                elif stat.startswith('opb_'):
                                    similarity_score = comparison_forward['opb_' + base_stat_suffix].iloc[0]

                            elif not comparison_backward.empty: # pep_df has (sample2, sample1)
                                if stat.startswith('opa_'):
                                    similarity_score = comparison_backward['opb_' + base_stat_suffix].iloc[0]
                                elif stat.startswith('opb_'):
                                    similarity_score = comparison_backward['opa_' + base_stat_suffix].iloc[0]

                            heatmap_data.loc[sample1, sample2] = similarity_score
                    else:
                        heatmap_data.loc[sample1, sample2] = np.nan # Upper triangle remains NaN

            # Convert to numeric again to ensure NaN handling, especially after assignments
            heatmap_data = heatmap_data.apply(pd.to_numeric, errors='coerce')

            ax = axes[i]
            # Plot the heatmap
            sns.heatmap(heatmap_data, annot=False, cmap='viridis', fmt=".2f", linewidths=.2,
                        cbar=(i == 0), cbar_ax=cbar_ax if i == 0 else None,
                        annot_kws={"size": 3}, vmin=0.0, vmax=1.0, ax=ax)
            ax.set_title(f'{stat} Heatmap')
            ax.set_xticks(np.arange(0.5, len(all_samples), 1)) # Center ticks between cells
            ax.set_yticks(np.arange(0.5, len(all_samples), 1))


            # Get the unique authorities in the sorted order of samples for axis labels
            ordered_authorities = [sample_authority.get(sample, sample) for sample in all_samples]
            ax.set_yticklabels(ordered_authorities, rotation=0, fontsize=8)
            ax.set_xticklabels(ordered_authorities, rotation=90, fontsize=8)
            ax.tick_params(axis='both', which='major', labelsize=8)

            if i == 0: # Only set colorbar label for the first plot in the pair
                cbar_ax.set_ylabel('Similarity Score', fontsize=12)

            # Export the heatmap_data to a CSV file
            output_path_csv = os.path.join(results_dir, f'heatmap_data_{stat}.csv')
            heatmap_data.index.name = 'sample_name'  # Set the index name
            heatmap_data.to_csv(output_path_csv, index=True, header=True)
            print(f"Saved heatmap data for {stat} to {output_path_csv}")

        # --- Final Plot Adjustments and Saving for the current pair ---
        plt.suptitle(f'Comparison Heatmaps ({all_relevant_stats[0].split("_")[0]} vs {all_relevant_stats[1].split("_")[0]}) - {species_title}', fontsize=16, y=1.02)
        plt.subplots_adjust(right=0.9, wspace=0.3, hspace=0.3)

        # Use a more descriptive filename for the combined plot of the pair
        output_filename = f'stacked_heatmap_by_authority_grouped_{all_relevant_stats[0]}_and_{all_relevant_stats[1]}.svg'
        output_path_svg = os.path.join(results_dir, output_filename)
        plt.savefig(output_path_svg, dpi=300, bbox_inches='tight')
        plt.close(fig) # Close the figure to free memory
        print(f"Saved heatmap plot for {all_relevant_stats[0]} and {all_relevant_stats[1]} to {output_path_svg}")


# Heat Map of histogram of counts
from scipy.stats import gaussian_kde
fig_width_mm = 170 # Increased width to accommodate more x-axis points
fig_width_inches = fig_width_mm / 25.4
# Adjust height for heatmap, it can be taller if many bins
fig_height_inches = fig_width_inches * (0.4) # Adjusted for better heatmap display with many points

# Initialize the plot
fig, ax = plt.subplots(figsize=(fig_width_inches, fig_height_inches))

relevant_stats = ['jaccard_sequences', 'jaccard_name_len', 'jaccard_lengths', 'jaccard_names']

# Convert all relevant stat columns to numeric
for stat in relevant_stats:
    pep_df[stat] = pd.to_numeric(pep_df[stat], errors='coerce')

# Define histogram bins
# You'll want to choose an appropriate number of bins.
# For scores from 0 to 1, 10 bins would be 0.1 wide, 20 bins 0.05 wide, etc.
num_bins = 40 # You can adjust this for more or less detail/granularity
bins = np.linspace(0, 1, num_bins + 1) # +1 because np.linspace includes start and end

# Create a list to store histogram data for each statistic
all_hist_data = {}

for stat in relevant_stats:
    current_values = []
    for i in range(len(all_samples)):
        for j in range(i + 1, len(all_samples)):
            sample1 = all_samples[i]
            sample2 = all_samples[j]

            comparison = pep_df[
                ((pep_df['sample_name_1'] == sample1) & (pep_df['sample_name_2'] == sample2)) |
                ((pep_df['sample_name_1'] == sample2) & (pep_df['sample_name_2'] == sample1))
            ]

            if not comparison.empty:
                current_values.append(comparison[stat].iloc[0])

    # Filter out NaN values before calculating histogram
    current_values = np.array([val for val in current_values if not pd.isna(val)])

    if len(current_values) > 0:
        # Calculate histogram frequencies
        # density=False means counts, density=True would normalize to sum to 1
        hist_counts, _ = np.histogram(current_values, bins=bins, density=False)
    else:
        hist_counts = np.zeros(num_bins) # No data, so all counts are zero

    all_hist_data[stat] = hist_counts

# Create a DataFrame for the heatmap using the histogram counts
# The index will represent the center of the bins or the start of the bins.
# Let's use the start of the bins for clarity.
bin_labels = [f'{bins[i]:.2f}-{bins[i+1]:.2f}' for i in range(num_bins)]
heatmap_df = pd.DataFrame(all_hist_data, index=bin_labels).T

global_max_val = heatmap_df.values.max()

# Create a new DataFrame for the normalized data
normalized_heatmap_df = heatmap_df.copy()

if global_max_val > 0:
    # Normalize the entire DataFrame by the global maximum
    normalized_heatmap_df = normalized_heatmap_df / global_max_val
else:
    # Handle the case where all counts are zero
    normalized_heatmap_df = normalized_heatmap_df * 0

# Plot the heatmap using the normalized DataFrame
sns.heatmap(normalized_heatmap_df, annot=False, cmap="Greens", ax=ax, cbar_kws={'label': 'Normalized Frequency'})
#sns.heatmap(heatmap_df, annot=False, cmap="Greens", ax=ax, cbar_kws={'label': 'Frequency (Counts)'})

# Add labels and title
ax.set_xlabel('Jaccard Score Bin')
ax.set_ylabel('Jaccard Similarities')
ax.set_title('Frequency Distribution of Jaccard Similarities Across Score Bins')

# # Set x-axis ticks and labels to represent the bins clearly
# # Using the bin labels directly makes it clear which range each column represents.
# ax.set_xticks(np.arange(num_bins) + 0.5) # Center ticks in the middle of the bins
# ax.set_xticklabels(bin_labels, rotation=90, ha='center', fontsize=8) # Rotate for readability

# Create a list of labels showing only the left edge of each bin
short_bin_labels = [f'{b:.2f}' for b in bins[:-1]]
short_bin_labels[0] = '0.00'  # Explicitly set the first label for clarity

# Define the step for how many ticks to show (e.g., show every 4th tick)
tick_step = 4
num_labels = len(short_bin_labels)

# Select the ticks and labels based on the step, ensuring the first and last are included
# The first tick is always at position 0.5.
# The selected tick positions are at the center of the relevant bins.
tick_positions = np.arange(num_labels)
selected_ticks = tick_positions[::tick_step]
selected_labels = short_bin_labels[::tick_step]

# Set the new ticks and labels
ax.set_xticks(selected_ticks)
ax.set_xticklabels(selected_labels, rotation=45, ha='center', fontsize=8) # Rotate for readability


fig.tight_layout()

# Save the heatmap plot
output_path = os.path.join(results_dir, 'jaccard_heatmap_histogram.svg') # New filename
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"Heatmap saved to {output_path}")

# Save the histogram numbers
hist_file = os.path.join(results_dir, "histogram_counts.csv")
heatmap_df.to_csv(hist_file, index=True) # index=True to keep bin labels

print(f"Histogram counts saved to {hist_file}")

# Heat Map of frequencies using KDE
from scipy.stats import gaussian_kde
fig_width_mm = 170 # Increased width to accommodate more x-axis points
fig_width_inches = fig_width_mm / 25.4
# Adjust height for heatmap, it can be taller if many bins
fig_height_inches = fig_width_inches * (0.4) # Adjusted for better heatmap display with many points

# Initialize the plot
fig, ax = plt.subplots(figsize=(fig_width_inches, fig_height_inches))

relevant_stats = ['jaccard_sequences', 'jaccard_name_len', 'jaccard_lengths', 'jaccard_names']

# Convert all relevant stat columns to numeric
for stat in relevant_stats:
    pep_df[stat] = pd.to_numeric(pep_df[stat], errors='coerce')

# Define a high-resolution spectrum from 0 to 1
x_spectrum = np.linspace(0, 1, 20) 

# Calculate smoothed densities for each statistic
all_densities_data = {}
for stat in relevant_stats:
    current_values = []
    for i in range(len(all_samples)):
        for j in range(i + 1, len(all_samples)):
            sample1 = all_samples[i]
            sample2 = all_samples[j]

            comparison = pep_df[
                ((pep_df['sample_name_1'] == sample1) & (pep_df['sample_name_2'] == sample2)) |
                ((pep_df['sample_name_1'] == sample2) & (pep_df['sample_name_2'] == sample1))
            ]

            if not comparison.empty:
                current_values.append(comparison[stat].iloc[0])

    # Filter out NaN values before calculating KDE
    current_values = np.array([val for val in current_values if not pd.isna(val)])
    actual_count_for_this_stat = len(current_values)
    print(actual_count_for_this_stat)
    #x_spectrum = np.linspace(0, 1, actual_count_for_this_stat) 

    if actual_count_for_this_stat > 1: # KDE requires at least 2 points
        kde = gaussian_kde(current_values, bw_method=0.10)
        densities = kde(x_spectrum)
    else: # Handle cases with insufficient data for KDE
        densities = np.zeros_like(x_spectrum)
    
    all_densities_data[stat] = densities

# Create a DataFrame for the heatmap using the smoothed densities

heatmap_df = pd.DataFrame(all_densities_data, index=x_spectrum).T
output_path_csv = os.path.join(results_dir, f'heatmap_data_kde.csv')
heatmap_df.to_csv(output_path_csv, index=True, header=True)

# Plot the heatmap (Jaccard Similarities on Y, Jaccard Score Spectrum on X)
# annot=False for smoother representation as there are too many points
sns.heatmap(heatmap_df, annot=False, cmap="Greens", ax=ax, cbar_kws={'label': 'Density'})

# Add labels and title
ax.set_xlabel('Jaccard Score')
ax.set_ylabel('Jaccard Similarities')
ax.set_title('Density Distribution of Jaccard Similarities Across Score Spectrum')

# Set x-axis ticks and labels to represent the spectrum clearly
# Show only a few ticks for clarity, e.g., every 0.1 or 0.2
tick_positions = np.arange(0, len(x_spectrum), len(x_spectrum) // 10) # Roughly 10 ticks
tick_labels = [f'{x_spectrum[int(p)]:.1f}' for p in tick_positions]

ax.set_xticks(tick_positions)
ax.set_xticklabels(tick_labels, rotation=45, ha='right')

fig.tight_layout()

# Save the heatmap plot
output_path = os.path.join(results_dir, 'jaccard_heatmap_spectrum.svg') # New filename
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"Heatmap saved to {output_path}")


# # PLOT MOW MEDIAN CHANGES BASED ON JACCARD_NAME_LEN

# # TODO PUT COUNTS ALL AT THE TOP OR THE BOTTOM
# df = pep_df.copy()

# stats = ["Mean", "Median"]

# for stat in stats:
#     # --- Calculate and Plot Statistic ---
#     plt.figure(figsize=(10, 6))  # Adjust figure size as needed

#     unique_jaccard_lengths = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
#     calculated_stats = []
#     counts = []  # Store the counts for each bin

#     for length in unique_jaccard_lengths:
#         lower_bound = length - 0.05
#         upper_bound = length + 0.05
#         bin_data = df[(df['jaccard_name_len'] > lower_bound) & (df['jaccard_name_len'] <= upper_bound)]['jaccard_sequences']
#         if stat == "Median":
#             calculated_stats.append(bin_data.median())
#         elif stat == "Mean":
#             calculated_stats.append(bin_data.mean())
#         counts.append(len(bin_data))  # Count the data points

#     # Create the bar plot
#     bars = plt.bar(unique_jaccard_lengths, calculated_stats, width=0.08)

#     # Add count labels at the bottom of the bars
#     for bar, count in zip(bars, counts):
#         yval = 0.05  # Position the text at the bottom (y=0)
#         plt.text(bar.get_x() + bar.get_width()/2, yval, "n="+str(count), ha='center', va='top',color='white', fontweight='bold')

#     # Add labels and title
#     plt.xlabel('Jaccard Similarity of Name Len Pairs')
#     plt.ylabel(f'{stat} Jaccard Similarity of Sequences')
#     plt.title(f'{stat} Jaccard Sequence Similarity vs. Jaccard Name Length Pair Similarity')
#     plt.xticks(unique_jaccard_lengths)  # Ensure all unique lengths are shown on x-axis
#     plt.ylim(bottom=0) # Ensure the y-axis starts at 0 to accommodate the labels

#     plt.tight_layout()
#     #output_path = os.path.join(results_dir, f'{stat}_counts')
#     output_path = os.path.join(results_dir, f'{stat}_counts.svg')
#     plt.savefig(output_path, dpi=300, bbox_inches='tight')
#     plt.close()



# # # Plot Hexbin ScatterPlot

# # # Assuming pep_df, results_dir, and species_title are defined elsewhere

# #relevant_stats = [('jaccard_sequences', 'jaccard_name_len')]
# relevant_stats = ['jaccard_sequences', 'jaccard_name_len', 'jaccard_lengths', 'jaccard_names']
# stat_combinations = combinations(iterable=relevant_stats,r=2)

# for stat_combination in stat_combinations:
#     fig, ax = plt.subplots(figsize=(10, 8))  # Adjust figure size

#     bottom_stat, top_stat = stat_combination  # Get the stats

#     pep_df[bottom_stat] = pd.to_numeric(pep_df[bottom_stat], errors='coerce')
#     pep_df[top_stat] = pd.to_numeric(pep_df[top_stat], errors='coerce')

#     # Create lists to store the paired data points
#     x_values = []
#     y_values = []


#     for i in range(len(all_samples)):
#         for j in range(i + 1, len(all_samples)):
#             sample1 = all_samples[i]
#             sample2 = all_samples[j]

#             comparison = pep_df[
#                 ((pep_df['sample_name_1'] == sample1) & (pep_df['sample_name_2'] == sample2)) |
#                 ((pep_df['sample_name_1'] == sample2) & (pep_df['sample_name_2'] == sample1))
#             ]

#             if not comparison.empty:
#                 x_values.append(comparison[bottom_stat].iloc[0])
#                 y_values.append(comparison[top_stat].iloc[0])

#     # Define the monochromatic colormap
#     colors = ["lightgrey", "darkorange"]
#     cmap = LinearSegmentedColormap.from_list("grey_to_orange", colors)
#     plot_extent = [0, 1, 0, 1]
#     ax.set_xlim(0, 1)
#     ax.set_ylim(0, 1)

#     # Create the hexbin plot with linear scaling first
#     hb = ax.hexbin(x_values, y_values, gridsize=25, cmap=cmap, alpha=0.8,extent=plot_extent)
#     ax.scatter(x_values, y_values, alpha=.20, label='Data Points',color='darkblue', s=25)
#     counts = hb.get_array()

#     # max_x_value = max(x_values)
#     # max_y_value = max(y_values)
#     # maxval = max(max_x_value,max_y_value)

#     # Apply logarithmic scaling to the counts, handling zeros
#     log_counts = np.log1p(counts)

#     # Normalize the logarithmic counts for the colormap
#     norm = plt.Normalize(vmin=log_counts.min(), vmax=log_counts.max())
#     colored_values = cmap(norm(log_counts))

#     # Set the facecolors of the hexbins based on the normalized log counts
#     hb.set_array(log_counts)  # Directly set the array to the log-transformed counts
#     hb.set_norm(norm)        # Set the normalization for the colormap

#     # Add a colorbar for the logarithmic scale
#     cb = fig.colorbar(hb, ax=ax, ticks=np.log1p([1, 10, 100, 1000]))  # Example ticks
#     cb.set_label('Counts (Log Scale)')
#     cb.ax.set_yticklabels([1, 10, 100, 1000])  # Label the ticks

#     # Set labels and title
#     ax.set_xlabel(bottom_stat)
#     ax.set_ylabel(top_stat)
#     ax.set_title(f'Hexbin Plot of {bottom_stat} vs {top_stat} (Log Scale)')
#     # ax.set_xlim(0, maxval)
#     # ax.set_ylim(0, maxval)

#     plt.tight_layout()
#     #output_path = os.path.join(results_dir, f'hexbin_monochromatic_logscale_{bottom_stat}_vs_{top_stat}')
#     output_path = os.path.join(results_dir, f'hexbin_monochromatic_logscale_{bottom_stat}_vs_{top_stat}.svg')
#     plt.savefig(output_path, dpi=300, bbox_inches='tight')
#     plt.close()

# # Plot bar graph freq distribution
# #TODO change stats
# #relevant_stats = [('jaccard_sequences', 'jaccard_name_len')]
# relevant_stats = ['jaccard_sequences', 'jaccard_name_len', 'jaccard_lengths', 'jaccard_names']
# stat_combinations = combinations(iterable=relevant_stats,r=2)
# for stat_combination in stat_combinations:
#     fig, ax = plt.subplots(figsize=(10, 6))  # Adjust figure size

#     bottom_stat, top_stat = stat_combination
#     pep_df[bottom_stat] = pd.to_numeric(pep_df[bottom_stat], errors='coerce')
#     pep_df[top_stat] = pd.to_numeric(pep_df[top_stat], errors='coerce')

#     # Collect all values for both statistics
#     bottom_values = []
#     top_values = []

#     for i in range(len(all_samples)):
#         for j in range(i + 1, len(all_samples)):
#             sample1 = all_samples[i]
#             sample2 = all_samples[j]

#             comparison = pep_df[
#                 ((pep_df['sample_name_1'] == sample1) & (pep_df['sample_name_2'] == sample2)) |
#                 ((pep_df['sample_name_1'] == sample2) & (pep_df['sample_name_2'] == sample1))
#             ]

#             if not comparison.empty:
#                 bottom_values.append(comparison[bottom_stat].iloc[0])
#                 top_values.append(comparison[top_stat].iloc[0])

#     # Define bins for the frequency plot
#     bins = np.linspace(0, 1, 21)  # Create 20 bins from 0 to 1

#     # Calculate frequencies for both statistics
#     freq_bottom, _ = np.histogram(bottom_values, bins=bins)
#     freq_top, _ = np.histogram(top_values, bins=bins)

#     # Set the width of the bars
#     width = 0.35

#     # Set the positions of the bars on the x-axis
#     x = np.arange(len(freq_bottom))

#     # Create the bar plot
#     rects1 = ax.bar(x - width/2, freq_bottom, width, label=bottom_stat, color='skyblue')
#     rects2 = ax.bar(x + width/2, freq_top, width, label=top_stat, color='salmon')

#     # Add labels, title, and legend
#     ax.set_xlabel('Overlap Coefficient')
#     ax.set_ylabel('Frequency')
#     ax.set_title(f'Frequency Distribution of {bottom_stat} and {top_stat}')
#     ax.set_xticks(x)
#     ax.set_xticklabels([f'{b:.2f}-{(b + (bins[1] - bins[0])):.2f}' for b in bins[:-1]], rotation=45, ha='right')
#     ax.legend()

#     fig.tight_layout()
#     #output_path = os.path.join(results_dir, f'frequency_plot_{bottom_stat}_vs_{top_stat}')
#     output_path = os.path.join(results_dir, f'frequency_plot_{bottom_stat}_vs_{top_stat}.svg')
#     plt.savefig(output_path, dpi=300, bbox_inches='tight')
#     #plt.show()
#     plt.close()


# Attempt frequency plot for all 4 stats in one plot
fig_width_mm = 170
fig_width_inches = fig_width_mm / 25.4
fig_height_inches = fig_width_inches * (0.6) # Adjust this ratio as needed

fig, ax = plt.subplots(figsize=(fig_width_inches, fig_height_inches))
#fig, ax = plt.subplots(figsize=(12, 7)) # Adjust figure size for more bars

relevant_stats = ['jaccard_sequences', 'jaccard_name_len', 'jaccard_lengths', 'jaccard_names']
colors = ['skyblue', 'salmon', 'lightgreen', 'plum'] # Define distinct colors for each stat

# Convert all relevant stat columns to numeric
for stat in relevant_stats:
    pep_df[stat] = pd.to_numeric(pep_df[stat], errors='coerce')

# Define bins for the frequency plot
bins = np.linspace(0, 1, 51)  

# Calculate frequencies for each statistic
all_freqs = {}
for stat in relevant_stats:
    # Collect all values for the current statistic
    current_values = []
    for i in range(len(all_samples)):
        for j in range(i + 1, len(all_samples)):
            sample1 = all_samples[i]
            sample2 = all_samples[j]

            comparison = pep_df[
                ((pep_df['sample_name_1'] == sample1) & (pep_df['sample_name_2'] == sample2)) |
                ((pep_df['sample_name_1'] == sample2) & (pep_df['sample_name_2'] == sample1))
            ]

            if not comparison.empty:
                current_values.append(comparison[stat].iloc[0])
    
    # Filter out NaN values before calculating histogram
    current_values = [val for val in current_values if not pd.isna(val)]
    
    freq, _ = np.histogram(current_values, bins=bins)
    all_freqs[stat] = freq

# Set the width of each individual bar
num_stats = len(relevant_stats)
bar_width = 0.8 / num_stats # Adjust bar width based on number of stats to prevent overlap

# Set the base positions for each group of bars
x_base = np.arange(len(bins) - 1)

# Plot bars for each statistic
rects = []
for i, stat in enumerate(relevant_stats):
    # Calculate the offset for each bar within a group
    offset = (i - (num_stats - 1) / 2) * bar_width
    bars = ax.bar(x_base + offset, all_freqs[stat], bar_width, label=stat, color=colors[i])
    rects.append(bars)

# Add labels, title, and legend
ax.set_xlabel('Overlap Coefficient')
ax.set_ylabel('Frequency')
ax.set_title('Frequency Distribution of All Relevant Statistics')
ax.set_xticks(x_base)
ax.set_xticklabels([f'{b:.2f}-{(b + (bins[1] - bins[0])):.2f}' for b in bins[:-1]], rotation=45, ha='right')
ax.legend()

fig.tight_layout()
output_path = os.path.join(results_dir, 'frequency_plot_all_stats.svg')
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"Frequency plot saved to {output_path}")


# Heat Map of frequencies
# fig_width_mm = 110
# fig_width_inches = fig_width_mm / 25.4
# # Adjust height for heatmap, it can be taller if many bins
# fig_height_inches = fig_width_inches * (0.7) # Adjusted for better heatmap display

# # Initialize the plot
# fig, ax = plt.subplots(figsize=(fig_width_inches, fig_height_inches))

# relevant_stats = ['jaccard_sequences', 'jaccard_name_len', 'jaccard_lengths', 'jaccard_names']

# # Convert all relevant stat columns to numeric
# for stat in relevant_stats:
#     pep_df[stat] = pd.to_numeric(pep_df[stat], errors='coerce')

# # Define bins for the frequency plot (20 bins from 0 to 1)
# bins = np.linspace(0, 1, 21)

# # Calculate frequencies for each statistic and store in a dictionary
# all_freqs_data = {}
# for stat in relevant_stats:
#     current_values = []
#     for i in range(len(all_samples)):
#         for j in range(i + 1, len(all_samples)):
#             sample1 = all_samples[i]
#             sample2 = all_samples[j]

#             comparison = pep_df[
#                 ((pep_df['sample_name_1'] == sample1) & (pep_df['sample_name_2'] == sample2)) |
#                 ((pep_df['sample_name_1'] == sample2) & (pep_df['sample_name_2'] == sample1))
#             ]

#             if not comparison.empty:
#                 current_values.append(comparison[stat].iloc[0])

#     # Filter out NaN values before calculating histogram
#     current_values = [val for val in current_values if not pd.isna(val)]

#     # Calculate frequency for the current stat using the defined bins
#     freq, _ = np.histogram(current_values, bins=bins)
#     all_freqs_data[stat] = freq

# # Create a DataFrame for the heatmap
# # Rows will be the bin labels, columns will be the statistics
# heatmap_df = pd.DataFrame(all_freqs_data)

# # Normalize frequencies to proportions for smoother representation
# # Divide each column by its sum
# heatmap_df_normalized = heatmap_df.apply(lambda x: x / x.sum(), axis=0)


# # Create labels for the x-axis (bins)
# bin_labels = [f'{bins[i]:.2f}-{bins[i+1]:.2f}' for i in range(len(bins) - 1)]
# #bin_labels =["0%","25%", "50%", "75%", "100%"]
# # No need to set index here if we are transposing later, as the original column names become index after transpose

# # Plot the heatmap (swapping x and y axes by transposing the DataFrame)
# # The `T` attribute transposes the DataFrame
# sns.heatmap(heatmap_df_normalized.T, fmt=".2f", cmap="viridis", ax=ax, cbar_kws={'label': 'Proportion'})

# # Add labels and title
# ax.set_xlabel('Jaccard Score')
# ax.set_ylabel('Jaccard Similarities')
# ax.set_title('Proportion Distribution of Jaccard Stats')
# ax.set_xticklabels(bin_labels, rotation=90, ha='right') # Set x-axis labels after transpose

# fig.tight_layout()

# # Save the heatmap plot
# output_path = os.path.join(results_dir, 'jaccard_heatmap_all_stats.svg')
# plt.savefig(output_path, dpi=300, bbox_inches='tight')
# plt.close()

# print(f"Heatmap saved to {output_path}")


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

    #fig_width_mm = 170
    # fig_width_inches = fig_width_mm / 25.4
    # fig_height_inches = fig_width_inches * (0.6) # Adjust this ratio as needed

    # fig, ax = plt.subplots(figsize=(fig_width_inches, fig_height_inches))
    fig, ax = plt.subplots(figsize=(5, 6))  # Adjust figure size as needed
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
    filtered_df = filtered_df.head(10)

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

    labels = metrics
    create_comparison_dot_plot(data_from_df, labels, title=f"Comparison of {primary_sample} vs All Other Ref Genomes", x_limit=1.0, metric_markers=['o', 's', 'D', '^'])
    
    #output_path = os.path.join(results_dir, f'comparison_{primary_sample}_.png')
    output_path = os.path.join(results_dir, f'comparison_{primary_sample}_.svg')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    #plt.show()


# def create_comparison_bar_graph(data, labels, title="Comparison Bar Graph",
#                                metric_colors=None):
#     """
#     Creates a grouped bar graph using Matplotlib.

#     Args:
#         data (dict): A dictionary where keys are sample names (e.g., 'Sample 1')
#                       and values are lists of 4 values.
#         labels (list): A list of label names (e.g., ['Metric A', 'Metric B', 'Metric C', 'Metric D']).
#         title (str, optional): The title of the bar graph. Defaults to "Comparison Bar Graph".
#         metric_colors (list, optional): A list of colors for each metric. If None, default colors are used.

#     Returns:
#         None: Displays the plot.
#     """
#     num_metrics = len(labels)
#     if num_metrics != 4:
#         raise ValueError("Number of labels must be 4 for this comparison bar graph.")

#     if metric_colors is None:
#         metric_colors = ['#1f77b4', '#ff7f0e', '#9467bd', '#8c564b']  # Default colors

#     fig, ax = plt.subplots(figsize=(10, 8))  # Adjust figure size as needed
#     num_samples = len(data)
#     bar_width = 0.2  # Width of each individual bar
#     group_width = num_metrics * bar_width  # Total width of a group of bars
#     index = np.arange(num_samples)  # Positions for the group centers

#     for i, (sample_name, values) in enumerate(data.items()):
#         if len(values) != num_metrics:
#             raise ValueError(f"Sample '{sample_name}' must have {num_metrics} values.")
#         for j, value in enumerate(values):
#             x_position = index[i] + (j * bar_width) - (group_width / 2) + (bar_width / 2) # Calculate the x position for each bar
#             ax.bar(x_position, value, bar_width, color=metric_colors[j], label=labels[j] if i == 0 else "")

#     ax.set_xticks(index)
#     ax.set_xticklabels(list(data.keys()), rotation=90, ha="right", fontsize=8)  # Sample names as x-axis labels
#     ax.set_ylabel("Metric Value")
#     ax.set_title(title, fontsize=14)
#     ax.legend(loc='upper right')
#     ax.grid(True, axis='y', linestyle='--', alpha=0.6)

#     plt.tight_layout()
#     #plt.show()

# # --- Example Usage ---
# # Assuming pep_df and results_dir are defined elsewhere

# primary_samples = ["hg19-initial-ucsc", "GRCh38.p14-fasta-genomic"]
# metrics = ['jaccard_names', 'jaccard_lengths', 'jaccard_sequences', 'jaccard_name_len']

# for primary_sample in primary_samples:
#     new_df = pep_df.copy()

#     # Filter the DataFrame
#     filtered_df = new_df[(new_df['sample_name_1'] == primary_sample) | (new_df['sample_name_2'] == primary_sample)]

#     # Convert the filtered DataFrame to the dictionary format
#     data_from_df = {}
#     for _, row in filtered_df.iterrows():
#         if row['sample_name_1'] == primary_sample:
#             sample_name = row['sample_name_2']
#         else:
#             sample_name = row['sample_name_1']
#         if sample_name not in data_from_df:
#             data_from_df[sample_name] = []
#         data_from_df[sample_name] = [row[metrics[0]], row[metrics[1]], row[metrics[2]], row[metrics[3]]]

#     # Use the columns f10_names, f10_lengths, f10_sequences, f10_name_len_pairs
#     labels = metrics
#     create_comparison_bar_graph(data_from_df, labels, title=f"Comparison of {primary_sample} vs All Other Ref Genomes",
#                                  )

#     #output_path = os.path.join(results_dir, f'comparison_bar_graph_{primary_sample}.png')
#     output_path = os.path.join(results_dir, f'comparison_bar_graph_{primary_sample}.svg')
#     plt.savefig(output_path, dpi=300, bbox_inches='tight')
#     plt.close()