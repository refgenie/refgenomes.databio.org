# import sys
import os
# import urllib.request
# import urllib.error
import pipestat
from pephubclient import PEPHubClient
from refget import fasta_to_digest, fasta_to_seqcol_dict, compare_seqcols, SequenceCollection, ga4gh_digest
import json

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

OUTPUT_PATH = "/home/drc/Downloads/refgenomes_pics_test/17june2025/seq_presence/"

psm = pipestat.PipestatManager(pephub_path="donaldcampbelljr/human_seqcol_digests_easier_labels:default")

# psm = pipestat.PipestatManager(pephub_path="donaldcampbelljr/ncbi_38_seqcol_digests:default")
# psm = pipestat.PipestatManager(pephub_path="donaldcampbelljr/ncbi_subset_digests:default")
# psm = pipestat.PipestatManager(pephub_path="donaldcampbelljr/ucsc_hg19_subset_digests:default")
# psm = pipestat.PipestatManager(pephub_path="donaldcampbelljr/ncbi_hg38_subset_digests:default")

results = psm.select_records()

digest_samplename = {}
all_digests = []

for result in results['records']:
    digest_samplename.update({result['record_identifier']: result['top_level_digest']})
    all_digests.append(result['top_level_digest'])

json_files = []
if os.path.isdir("/home/drc/Downloads/jsons_from_rivanna/json/"):
    for filename in os.listdir("/home/drc/Downloads/jsons_from_rivanna/json/"):
        if os.path.splitext(filename)[0] in all_digests:
            #print(f"{filename} in all_digest")
            if filename.endswith(".json"):
                full_path = os.path.join("/home/drc/Downloads/jsons_from_rivanna/json/", filename)
                json_files.append(full_path)
        else:
            #print(f"{filename} NOT in all_digest")
            pass

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

sorted_files = sorted(file_hit_counts.keys(), key=file_hit_counts.get, reverse=True)

file_authority={}
for file in sorted_files:
    for result in results['records']:
        if file == result['record_identifier']:
            file_authority[file]=result['authority']

print(file_authority)

# Sort the file_authority dictionary by its values
sorted_file_authority_list = sorted(file_authority.items(), key=lambda item: item[1])
# Create a list of sorted files based on authority
sorted_files_by_authority = [item[0] for item in sorted_file_authority_list]

print(sorted_file_authority_list)
print(f"Sorted files by authority: {sorted_files_by_authority}")

data = []
for seq in sorted_sequences_by_frequency:  # Iterate through the sorted sequences
    row = [1 if file in seq_samples_name[seq] else 0 for file in sorted_files_by_authority]
    data.append(row)

df = pd.DataFrame(data, index=sorted_sequences_by_frequency, columns=sorted_files_by_authority)

num_all_seqs = len(all_sequences_union) # Assuming this is defined

plt.figure(figsize=(24, 10))
sns.heatmap(df.T, cmap="viridis", cbar=False)  # Transpose the DataFrame
plt.title("Sequences Present in Reference Genomes (Sorted by Frequency)")
plt.ylabel("Reference Genomes")  # Swapped labels
plt.xlabel(f"Sequences (Sorted by Frequency), n={num_all_seqs}")  # Swapped labels
plt.xticks([])  # Adjust rotation as needed
plt.yticks(rotation=0)
plt.tight_layout()
output_path = os.path.join(OUTPUT_PATH, 'sequences_presence_sorted_frequency')
#output_path = os.path.join(OUTPUT_PATH, 'sequences_presence_sorted_frequency.svg')
plt.savefig(output_path, dpi=300, bbox_inches='tight')