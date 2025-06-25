import os
import pipestat
import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


# FOR CLUSTERING SEQUENCES
# -------------------------
OUTPUT_PATH = "/home/drc/Downloads/refgenomes_pics_test/17june2025/MOUSE/seq_presence/"
#JSON_DIR = "/home/drc/Downloads/jsons_from_rivanna/json/" # HUMAN JSONS
JSON_DIR = "/home/drc/Downloads/mouse_jsons_from_rivanna/json/" # MOUSE JSONS

#psm = pipestat.PipestatManager(pephub_path="donaldcampbelljr/human_seqcol_digests:default")
psm = pipestat.PipestatManager(pephub_path="donaldcampbelljr/mouse_seqcol_digests:default")

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

# Sort the DataFrame by count in descending order
sequence_counts_df_sorted = sequence_counts_df.sort_values(by='Count', ascending=False)

# Get the sorted list of sequences based on their counts
sorted_sequences_by_frequency = sequence_counts_df_sorted['Sequence'].tolist()

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

sequences = sorted(seq_samples_name.keys())
files = sorted(list(set(file for files in seq_samples_name.values() for file in files)))

file_hit_counts = {}
for file in files:
    count = 0
    for seq in sequences:
        if file in seq_samples_name[seq]:
            count += 1
    file_hit_counts[file] = count

# Sort files by the number of sequences they contain (row counts)
sorted_files_by_row_count = sorted(file_hit_counts.keys(), key=file_hit_counts.get, reverse=True)

data = []
for seq in sorted_sequences_by_frequency:  # Iterate through the sorted sequences
    row = [1 if file in seq_samples_name[seq] else 0 for file in sorted_files_by_row_count]
    data.append(row)

df = pd.DataFrame(data, index=sorted_sequences_by_frequency, columns=sorted_files_by_row_count)
# Export the df DataFrame to a CSV file
df.to_csv(os.path.join(OUTPUT_PATH, 'sequence_presence_matrix.csv'), index=True)  # IMPORTANT: index=True
print(f"Sequence presence matrix exported to: {os.path.join(OUTPUT_PATH, 'sequence_presence_matrix.csv')}")
# num_all_seqs = len(all_sequences_union) # Assuming this is defined

# plt.figure(figsize=(24, 10))
# sns.heatmap(df.T, cmap="viridis", cbar=False)  # Transpose the DataFrame
# plt.title("Sequences Present in Reference Genomes (Sorted by Frequency, Rows by Sequence Count)")
# plt.ylabel("Reference Genomes (Sorted by Sequence Count)")  # Updated label
# plt.xlabel(f"Sequences (Sorted by Frequency), n={num_all_seqs}")
# plt.xticks([])
# plt.yticks(rotation=0)
# plt.tight_layout()
# output_path = os.path.join(OUTPUT_PATH, 'sequences_presence_sorted_frequency_rows_by_count')
# plt.savefig(output_path, dpi=300, bbox_inches='tight')


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