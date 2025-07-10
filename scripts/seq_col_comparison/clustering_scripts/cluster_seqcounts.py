import os
from matplotlib.colors import ListedColormap
import pipestat
import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


# FOR CLUSTERING SEQUENCES
# -------------------------
OUTPUT_PATH = "/home/drc/Downloads/refgenomes_pics_test/09jul2025/seq_presence_testing/"
JSON_DIR = "/home/drc/Downloads/jsons_from_rivanna/json/" # HUMAN JSONS
#JSON_DIR = "/home/drc/Downloads/mouse_jsons_from_rivanna/json/" # MOUSE JSONS

psm = pipestat.PipestatManager(pephub_path="donaldcampbelljr/human_seqcol_digests:default")
#psm = pipestat.PipestatManager(pephub_path="donaldcampbelljr/mouse_seqcol_digests:default")

results = psm.select_records()

digest_samplename = {}
all_digests = []

for result in results['records']:
    digest_samplename.update({result['record_identifier']: result['top_level_digest']})
    all_digests.append(result['top_level_digest'])

json_files = []
if os.path.isdir(JSON_DIR):
    for filename in os.listdir(JSON_DIR):
        if os.path.splitext(filename)[0] in all_digests:
            if filename.endswith(".json"):
                full_path = os.path.join(JSON_DIR, filename)
                json_files.append(full_path)

all_sequences_union = set()
dict_sequences = {}

for fp in json_files:
    with open(fp, "r") as f:
        reloaded_dict1 = json.load(fp=f)
        sequences = reloaded_dict1['sorted_sequences']
        all_sequences_union = all_sequences_union.union(set(sequences))
        dict_sequences.update({fp: sequences})

num_all_seqs = len(all_sequences_union)
sequence_counts = {}
single_sequence_in_fp = {}

for seq in all_sequences_union:
    list_fps = []
    seq_count = 0
    for key, value in dict_sequences.items():
        if seq in value:
            seq_count += 1
            list_fps.append(key)
    sequence_counts.update({seq: seq_count})
    single_sequence_in_fp.update({seq: list_fps})

## Convert the sequence_counts dictionary to a Pandas DataFrame
sequence_counts_df = pd.DataFrame(list(sequence_counts.items()), columns=['Sequence', 'Count'])

print(sequence_counts_df.head())

# Sort the DataFrame by count in descending order
sequence_counts_df_sorted = sequence_counts_df.sort_values(by='Count', ascending=False)
print(sequence_counts_df_sorted.head())

# Get the sorted list of sequences based on their counts
sorted_sequences_by_frequency = sequence_counts_df_sorted['Sequence'].tolist()

print(sorted_sequences_by_frequency[:5])

seq_samples_name = {}
for seq in all_sequences_union:
    sample_names = []
    for key, value in dict_sequences.items():
        if seq in value:
            base_name_with_extension = os.path.basename(key)
            digest = os.path.splitext(base_name_with_extension)[0]
            for k, v in digest_samplename.items():
                if digest in v:
                    sample_names.append(k)
    seq_samples_name.update({seq: sample_names})


# print(seq_samples_name.keys())
sequences = list(seq_samples_name.keys())
#print(sequences[:10])
files = sorted(list(set(file for files in seq_samples_name.values() for file in files)))
print(len(files))
file_hit_counts = {}
for file in files:
    count = 0
    for seq in sequences:
        if file in seq_samples_name[seq]:
            count += 1
    file_hit_counts[file] = count

# Sort files by the number of sequences they contain (row counts)
sorted_files_by_row_count = sorted(file_hit_counts.keys(), key=file_hit_counts.get, reverse=True)

# print(sorted_files_by_row_count.keys()[:10])

# 2. Build the CUSTOM sorted list of files/samples
custom_sorted_files = []
seen_files = set()

# Iterate through sequences from most frequent to least frequent
for seq in sorted_sequences_by_frequency:
    # Get the files associated with the current sequence
    files_for_this_seq = seq_samples_name.get(seq, []) # Use .get() to avoid KeyError if a sequence somehow isn't in seq_samples_name

    # Sort these files by their *overall* hit count (which you calculated earlier)
    # This step ensures that among the files for the current sequence,
    # those that are generally richer in sequences appear earlier.
    # You'll need file_hit_counts from your previous code:
    # file_hit_counts = {}
    # for file in all_sample_names: # You'll need 'all_sample_names' defined from your original code
    #    count = 0
    #    for s in seq_samples_name.keys():
    #        if file in seq_samples_name[s]:
    #            count += 1
    #    file_hit_counts[file] = count

    # Sort files_for_this_seq based on their count (descending)
    # Ensure file_hit_counts is available and correctly populated
    try:
        files_for_this_seq_sorted_by_overall_count = sorted(
            files_for_this_seq,
            key=lambda f: file_hit_counts.get(f, 0), # Use .get(f, 0) to handle potential missing files gracefully
            reverse=True
        )
    except NameError:
        print("Error: 'file_hit_counts' or 'all_sample_names' not defined. Please ensure the preceding code for these variables is run.")
        # If file_hit_counts isn't available, fall back to simple sorting or error out
        files_for_this_seq_sorted_by_overall_count = sorted(files_for_this_seq)


    # Append files to the custom_sorted_files list, only if not already seen
    for file in files_for_this_seq_sorted_by_overall_count:
        if file not in seen_files:
            custom_sorted_files.append(file)
            seen_files.add(file)

print(custom_sorted_files)



data = []
for seq in sorted_sequences_by_frequency:  # Iterate through the sorted sequences
    row = [1 if file in seq_samples_name[seq] else 0 for file in sorted_files_by_row_count]
    data.append(row)

df = pd.DataFrame(data, index=sorted_sequences_by_frequency, columns=sorted_files_by_row_count)
# Export the df DataFrame to a CSV file
# df.to_csv(os.path.join(OUTPUT_PATH, 'sequence_presence_matrix.csv'), index=True)  # IMPORTANT: index=True
# print(f"Sequence presence matrix exported to: {os.path.join(OUTPUT_PATH, 'sequence_presence_matrix.csv')}")
# num_all_seqs = len(all_sequences_union) # Assuming this is defined

# plt.figure(figsize=(24, 20))
# sns.heatmap(df.T, cmap="viridis", cbar=False)  # Transpose the DataFrame
# plt.title("Sequences Present in Reference Genomes (Sorted by Frequency, Rows by Sequence Count)")
# plt.ylabel("Reference Genomes (Sorted by Sequence Count)")  # Updated label
# plt.xlabel(f"Sequences (Sorted by Frequency), n={num_all_seqs}")
# plt.xticks([])
# plt.yticks(rotation=0)
# plt.tight_layout()
# output_path = os.path.join(OUTPUT_PATH, 'sequences_presence_sorted_frequency_rows_by_count')
# plt.savefig(output_path, dpi=300, bbox_inches='tight')

# TRY WITH CUSTOM SORTING
data = []
for seq in sorted_sequences_by_frequency:  # Iterate through the sorted sequences
    row = [1 if file in seq_samples_name[seq] else 0 for file in custom_sorted_files]
    data.append(row)

df = pd.DataFrame(data, index=sorted_sequences_by_frequency, columns=custom_sorted_files)
# Export the df DataFrame to a CSV file
df.to_csv(os.path.join(OUTPUT_PATH, 'sequence_presence_matrix.csv'), index=True)  # IMPORTANT: index=True
print(f"Sequence presence matrix exported to: {os.path.join(OUTPUT_PATH, 'sequence_presence_matrix.csv')}")
# num_all_seqs = len(all_sequences_union) # Assuming this is defined

plt.figure(figsize=(24, 20))
sns.heatmap(df.T, cmap="viridis", cbar=False)  # Transpose the DataFrame
plt.title("Sequences Present in Reference Genomes (Sorted by Frequency, Rows by Sequence Count)")
plt.ylabel("Reference Genomes (Sorted by Sequence Count)")  # Updated label
plt.xlabel(f"Sequences (Sorted by Frequency), n={num_all_seqs}")
plt.xticks([])
plt.yticks(rotation=0)
plt.tight_layout()
output_path = os.path.join(OUTPUT_PATH, 'sequences_presence_custom_sorted_frequency_rows_by_count')
plt.savefig(output_path, dpi=300, bbox_inches='tight')

# df = df.astype(int)
# from scipy.spatial.distance import pdist, squareform
# from scipy.cluster.hierarchy import linkage, leaves_list

# ... (your data loading and initial df construction as before) ...
# Let's assume 'df' is your presence/absence DataFrame with sequences as index and files as columns.
# It's important that df has sequences as index and files as columns if you want to cluster them directly.
# If you used the previous "better way" code, `df` is already in this format.

# # --- Cluster Rows (Sequences) ---
# # Calculate distances between sequences (rows)
# # 'binary' or 'jaccard' are good metrics for binary data
# row_distances = pdist(df.values, metric='jaccard') # Or 'hamming', 'binary'
# row_linkage = linkage(row_distances, method='average') # 'average', 'complete', 'ward'

# # Get the order of leaves from the dendrogram
# row_order = leaves_list(row_linkage)
# sorted_sequences_by_clustering = df.index[row_order].tolist()


# # --- Cluster Columns (Files) ---
# # Calculate distances between files (columns)
# col_distances = pdist(df.T.values, metric='jaccard') # Transpose df for column distances
# col_linkage = linkage(col_distances, method='complete')
# col_order = leaves_list(col_linkage)
# sorted_files_by_clustering = df.columns[col_order].tolist()

# # Reindex the DataFrame according to the clustering orders
# df_clustered = df.reindex(index=sorted_sequences_by_clustering, columns=sorted_files_by_clustering)

# print(df_clustered.head())
# print(df_clustered.tail())
# print(df_clustered.min().min()) # Should be 0
# print(df_clustered.max().max()) # Should be 1
# print(df_clustered.sum().sum()) # Should be > 0 (if there are any 1s)

# # print(f"df_clustered shape: {df_clustered.shape}")
# # calculated_width = max(10, len(df_clustered.columns) * 0.2)
# # calculated_height = max(8, len(df_clustered.index) * 0.2)
# # print(f"Calculated figsize: ({calculated_width}, {calculated_height})")
# plt.figure(figsize=(10, 8))
# # # Generate heatmap
# # plt.figure(figsize=(max(10, len(df_clustered.columns) * 0.2), max(8, len(df_clustered.index) * 0.2)))
# #cmap_binary = ListedColormap(['white', 'darkblue'])
# cmap_binary = ListedColormap(['white', 'black'])

# sns.heatmap(df_clustered.T, cmap="viridis", cbar=False, linecolor='lightgrey', linewidths=0.0)
# plt.title("Sequence Presence Heatmap (Hierarchical Clustering)")
# plt.xlabel("Files (Clustered Order)")
# plt.ylabel("Sequences (Clustered Order)")
# plt.xticks([])
# plt.tight_layout()
# plt.savefig(os.path.join(OUTPUT_PATH, 'sequence_presence_heatmap_clustered.png'))
# plt.show()

# # CLUSTERING BY NAME_LENGTHS
# OUTPUT_PATH = "/home/drc/Downloads/refgenomes_pics_test/17june2025/seq_presence/test/"

# psm = pipestat.PipestatManager(pephub_path="donaldcampbelljr/human_seqcol_digests:default")

# results = psm.select_records()

# digest_samplename = {}
# all_digests = []

# for result in results['records']:
#     digest_samplename.update({result['record_identifier']: result['top_level_digest']})
#     all_digests.append(result['top_level_digest'])

# json_files = []
# if os.path.isdir("/home/drc/Downloads/jsons_from_rivanna/json/"):
#     for filename in os.listdir("/home/drc/Downloads/jsons_from_rivanna/json/"):
#         if os.path.splitext(filename)[0] in all_digests:
#             if filename.endswith(".json"):
#                 full_path = os.path.join("/home/drc/Downloads/jsons_from_rivanna/json/", filename)
#                 json_files.append(full_path)

# all_sequences_union = set()
# dict_sequences = {}

# for fp in json_files:
#     with open(fp, "r") as f:
#         reloaded_dict1 = json.load(fp=f)
#         # CHANGED:  Use name_length_pairs instead of sorted_sequences
#         sequences = reloaded_dict1['name_length_pairs']
#         all_sequences_union = all_sequences_union.union(set(tuple(item.items()) for item in sequences)) #make comparable
#         dict_sequences.update({fp: sequences})
# print("FINISHED CREATING dict_sequuences")
# num_all_seqs = len(all_sequences_union)
# sequence_counts = {}
# single_sequence_in_fp = {}

# for seq in all_sequences_union:
#     list_fps = []
#     seq_count = 0
#     for key, value in dict_sequences.items():
#         # CHANGED: Check for presence of the pair in the list of dictionaries
#         if any(dict(seq).items() <= d.items() for d in value): # changed the logic here
#             seq_count += 1
#             list_fps.append(key)
#     sequence_counts.update({seq: seq_count})
#     single_sequence_in_fp.update({seq: list_fps})

# print("Finsihed sequence counts")

# ## Convert the sequence_counts dictionary to a Pandas DataFrame
# sequence_counts_df = pd.DataFrame(list(sequence_counts.items()), columns=['Sequence', 'Count'])

# # Sort the DataFrame by count in descending order
# sequence_counts_df_sorted = sequence_counts_df.sort_values(by='Count', ascending=False)

# # Get the sorted list of sequences based on their counts
# sorted_sequences_by_frequency = sequence_counts_df_sorted['Sequence'].tolist()

# seq_samples_name = {}
# for seq in all_sequences_union:
#     sample_names = []
#     for key, value in dict_sequences.items():
#         if any(dict(seq).items() <= d.items() for d in value): # changed the logic here
#             base_name_with_extension = os.path.basename(key)
#             digest = os.path.splitext(base_name_with_extension)[0]
#             for k, v in digest_samplename.items():
#                 if digest in v:
#                     sample_names.append(k)
#     seq_samples_name.update({seq: sample_names})

# sequences = sorted(seq_samples_name.keys())
# files = sorted(list(set(file for files in seq_samples_name.values() for file in files)))

# file_hit_counts = {}
# for file in files:
#     count = 0
#     for seq in sequences:
#         if file in seq_samples_name[seq]:
#             count += 1
#     file_hit_counts[file] = count

# # Sort files by the number of sequences they contain (row counts)
# sorted_files_by_row_count = sorted(file_hit_counts.keys(), key=file_hit_counts.get, reverse=True)

# data = []
# for seq in sorted_sequences_by_frequency:  # Iterate through the sorted sequences
#     row = [1 if file in seq_samples_name[seq] else 0 for file in sorted_files_by_row_count]
#     data.append(row)

# df = pd.DataFrame(data, index=sorted_sequences_by_frequency, columns=sorted_files_by_row_count)
# # Export the df DataFrame to a CSV file
# df.to_csv(os.path.join(OUTPUT_PATH, 'sequence_presence_matrix.csv'), index=True)  # IMPORTANT: index=True
# print(f"Sequence presence matrix exported to: {os.path.join(OUTPUT_PATH, 'sequence_presence_matrix.csv')}")
# num_all_seqs = len(all_sequences_union) # Assuming this is defined

# plt.figure(figsize=(24, 10))
# sns.heatmap(df.T, cmap="viridis", cbar=False)  # Transpose the DataFrame
# plt.title("Name-Length Pairs Present in Reference Genomes (Sorted by Frequency, Rows by Sequence Count)") # changed title
# plt.ylabel("Reference Genomes (Sorted by Sequence Count)")
# plt.xlabel(f"Name-Length Pairs (Sorted by Frequency), n={num_all_seqs}") # changed xlabel
# plt.xticks([])
# plt.yticks(rotation=0)
# plt.tight_layout()
# output_path = os.path.join(OUTPUT_PATH, 'sequences_presence_sorted_frequency_rows_by_count')
# plt.savefig(output_path, dpi=300, bbox_inches='tight')