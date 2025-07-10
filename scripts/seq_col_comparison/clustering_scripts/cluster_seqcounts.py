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
# OUTPUT_PATH = "/home/drc/Downloads/refgenomes_pics_test/09jul2025/seq_presence_testing/"'
OUTPUT_PATH = "/home/drc/Downloads/refgenomes_pics_test/09jul2025/mouse_seq_presence/"
# JSON_DIR = "/home/drc/Downloads/jsons_from_rivanna/json/" # HUMAN JSONS
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


# TRY WITH CUSTOM SORTING
data = []
for seq in sorted_sequences_by_frequency:  # Iterate through the sorted sequences
    row = [1 if file in seq_samples_name[seq] else 0 for file in custom_sorted_files]
    data.append(row)


df = pd.DataFrame(data, index=sorted_sequences_by_frequency, columns=custom_sorted_files)

df = df.astype(int)

# Get the current row labels (sequence IDs)
current_sequence_labels = df.index.tolist()

sorted_sequence_labels_by_binary_pattern = sorted(
    current_sequence_labels,
    key=lambda seq_label: tuple(df.loc[seq_label, :]),
    reverse=True
)

# Reindex the DataFrame using this new row order
df_final_sorted = df.reindex(index=sorted_sequence_labels_by_binary_pattern)

# Export the df DataFrame to a CSV file
df_final_sorted.to_csv(os.path.join(OUTPUT_PATH, 'sequence_presence_matrix.csv'), index=True)  # IMPORTANT: index=True
print(f"Sequence presence matrix exported to: {os.path.join(OUTPUT_PATH, 'sequence_presence_matrix.csv')}")
# num_all_seqs = len(all_sequences_union) # Assuming this is defined

plt.figure(figsize=(24, 20))
sns.heatmap(df_final_sorted.T, cmap="viridis", cbar=False)  # Transpose the DataFrame
plt.title("Sequences Present in Reference Genomes (Sorted by Frequency, Rows by Sequence Count)")
plt.ylabel("Reference Genomes (Sorted by Sequence Count)")  # Updated label
plt.xlabel(f"Sequences (Sorted by Frequency), n={num_all_seqs}")
plt.xticks([])
plt.yticks(rotation=0)
plt.tight_layout()
output_path = os.path.join(OUTPUT_PATH, 'sequences_presence_custom_sorted_frequency_rows_by_count')
plt.savefig(output_path, dpi=300, bbox_inches='tight')