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
import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from matplotlib.colors import ListedColormap

OUTPUT_PATH = "/home/drc/Downloads/refgenomes_pics_test/06may2025/"

psm = pipestat.PipestatManager(pephub_path="donaldcampbelljr/human_seqcol_digests:default")

# psm = pipestat.PipestatManager(pephub_path="donaldcampbelljr/ncbi_38_seqcol_digests:default")
# psm = pipestat.PipestatManager(pephub_path="donaldcampbelljr/ncbi_subset_digests:default")
# psm = pipestat.PipestatManager(pephub_path="donaldcampbelljr/ucsc_hg19_subset_digests:default")
# psm = pipestat.PipestatManager(pephub_path="donaldcampbelljr/ncbi_hg38_subset_digests:default")

results = psm.select_records()

digest_samplename = {}
samplename_authority = {} # Dictionary to store sample name and its authority
all_digests = []

for result in results['records']:
    record_identifier = result['record_identifier']
    top_level_digest = result['top_level_digest']
    authority = result['authority']
    digest_samplename.update({record_identifier: top_level_digest})
    samplename_authority[record_identifier] = authority # Store sample name and authority
    all_digests.append(top_level_digest)

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
        reloaded_dict1 = json.load(f)
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

# Create a dictionary mapping file (which is the sample name here) to its authority
file_authority = {}
for file in files:
    for result in results['records']:
        if file == result['record_identifier']:
            file_authority[file] = result['authority']
            break

# Create the DataFrame
data = []
for seq in sequences:
    row = [1 if file in seq_samples_name[seq] else 0 for file in files]
    data.append(row)

df = pd.DataFrame(data, index=sequences, columns=files)

# Perform hierarchical clustering on the transposed DataFrame (samples as rows)
linked = linkage(df.T, method='ward')
ordered_sample_indices = leaves_list(linked)
sorted_df = df.iloc[:, ordered_sample_indices]
sorted_files = sorted_df.columns.tolist()

# Create a mapping of the sorted files to their authorities
sorted_file_authorities = {file: file_authority.get(file, 'unknown') for file in sorted_files}
unique_authorities = sorted(list(set(sorted_file_authorities.values())))
authority_colors = plt.cm.viridis(np.linspace(0, 1, len(unique_authorities)))
authority_cmap = ListedColormap(authority_colors)
authority_mapping = {auth: color for auth, color in zip(unique_authorities, authority_colors)}
row_colors = [authority_mapping.get(sorted_file_authorities[file], [0.5, 0.5, 0.5, 1]) for file in sorted_files]

plt.figure(figsize=(12, 10))
clustermap =sns.clustermap(sorted_df.T, cmap="magma", row_cluster=True, col_cluster=False,
               row_linkage=linked, row_colors=row_colors, yticklabels=sorted_files)
plt.title("Sequences Present in Reference Genomes (Clustered by Sample)")
plt.ylabel("Reference Genomes (Clustered)")
plt.xlabel(f"Sequences, n={num_all_seqs}")
plt.yticks(rotation=0)
clustermap.ax_heatmap.set_xticks([]) # Remove x-axis ticks
clustermap.ax_heatmap.set_xticklabels([]) # Remove x-axis tick labels
plt.tight_layout()

# Add a legend for the authority colors
handles = [plt.Rectangle((0, 0), 1, 1, color=authority_mapping[auth]) for auth in unique_authorities]
plt.legend(handles, unique_authorities, title='Authority', bbox_to_anchor=(1.05, 1), loc='upper left')

output_path = os.path.join(OUTPUT_PATH, 'sequences_presence_clustered_authority')
plt.savefig(output_path, dpi=300, bbox_inches='tight')
# plt.show()