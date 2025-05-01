# import sys
import os
# import urllib.request
# import urllib.error
import pipestat
from pephubclient import PEPHubClient
from refget import fasta_to_digest, fasta_to_seqcol_dict, compare_seqcols, SequenceCollection,ga4gh_digest
import json

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np



#psm = pipestat.PipestatManager(pephub_path="donaldcampbelljr/human_seqcol_digests:default")

#psm = pipestat.PipestatManager(pephub_path="donaldcampbelljr/ncbi_38_seqcol_digests:default")

#psm = pipestat.PipestatManager(pephub_path="donaldcampbelljr/ncbi_subset_digests:default")

psm = pipestat.PipestatManager(pephub_path="donaldcampbelljr/ucsc_subset_digests:default")

results = psm.select_records()

#print(results['records'])


digest_samplename = {}

all_digests = []

for result in results['records']:
    #digest_samplename.update({result['top_level_digest']:result['record_identifier']})
    digest_samplename.update({result['record_identifier']:result['top_level_digest']})
    all_digests.append(result['top_level_digest'])


json_files = []
if os.path.isdir("/home/drc/Downloads/jsons_from_rivanna/json/"):
    for filename in os.listdir("/home/drc/Downloads/jsons_from_rivanna/json/"):
        if os.path.splitext(filename)[0] in all_digests:
            print(f"{filename} in all_digest")
            if filename.endswith(".json"):
                full_path = os.path.join("/home/drc/Downloads/jsons_from_rivanna/json/", filename)
                json_files.append(full_path)
        else:
            print(f"{filename} NOT in all_digest")

all_sequences_union = set()

dict_sequences = {}

for fp in json_files:
    with open(fp, "r") as f:
        reloaded_dict1 = json.load(fp=f)
        sequences = reloaded_dict1['sorted_sequences']
        all_sequences_union = all_sequences_union.union(set(sequences))
        #dict_sequences.update({os.path.basename(fp):sequences})
        dict_sequences.update({fp: sequences})

#print(len(all_sequences_union))
num_all_seqs = len(all_sequences_union)

sequence_counts = {}

single_sequence_in_fp={}

for seq in all_sequences_union:
    list_fps = []
    seq_count = 0
    for key, value in dict_sequences.items():
        if seq in value:
            seq_count+=1
            list_fps.append(key)
    sequence_counts.update({seq:seq_count})
    single_sequence_in_fp.update({seq:list_fps})



# Convert the sequence_counts dictionary to a Pandas DataFrame for easier plotting
sequence_counts_df = pd.DataFrame(list(sequence_counts.items()), columns=['Sequence', 'Count'])

# Sort the DataFrame by count in descending order (optional, for better visualization)
sequence_counts_df_sorted = sequence_counts_df.sort_values(by='Count', ascending=False)
#sequence_counts_df_sorted = sequence_counts_df_sorted.head(50).copy()
# --- Plotting Sequence vs Count ---
# plt.figure(figsize=(12, 6))
# plt.bar(sequence_counts_df_sorted['Sequence'], sequence_counts_df_sorted['Count'])
# plt.xlabel("Sequence")
# plt.ylabel("Count (Number of Files Containing Sequence)")
# plt.title("Sequence Frequency Across Ref Genomes")
# plt.xticks(rotation=90, fontsize=8)
# plt.tight_layout()
# plt.show()


# print(single_sequence_in_fp)


# seq_samples_name = {}



# for seq, filelist in single_sequence_in_fp.items():
#     sample_names =[]
#     print(f"SEQUENCE {seq}")
#     for file in filelist:
#         #print(os.path.basename(file))
#         base_name_with_extension = os.path.basename(file)
#         base_name_without_extension = os.path.splitext(base_name_with_extension)[0]

#         for k,v in digest_samplename.items():
#             if base_name_without_extension in v:
#                 print(k)
#                 sample_names.append(k)

#     print("NEXT SEQUENCE")



# sequence_file_map = {
#     "SEQ1": ["file_A.txt", "file_B.txt"],
#     "SEQ2": ["file_B.txt", "file_C.txt", "file_D.txt"],
#     "SEQ3": ["file_A.txt"],
#     "SEQ4": ["file_C.txt", "file_E.txt"],
#     "SEQ5": ["file_B.txt"],
#     "SEQ11": ["file_A.txt", "file_B.txt"],
#     "SEQ12": ["file_B.txt", "file_C.txt", "file_D.txt"],
#     "SEQ13": ["file_A.txt"],
#     "SEQ14": ["file_C.txt", "file_E.txt"],
#     "SEQ15": ["file_B.txt"],
#     "SEQ21": ["file_A.txt", "file_B.txt"],
#     "SEQ22": ["file_B.txt", "file_C.txt", "file_D.txt"],
#     "SEQ23": ["file_A.txt"],
#     "SEQ24": ["file_C.txt", "file_E.txt"],
#     "SEQ25": ["file_B.txt"]
# }

seq_samples_name = {}
for seq in all_sequences_union:
    sample_names = []
    for key, value in dict_sequences.items():
        if seq in value:
            base_name_with_extension = os.path.basename(key)
            digest = os.path.splitext(base_name_with_extension)[0]
            for k,v in digest_samplename.items():
                if digest in v:
                    sample_names.append(k)
    seq_samples_name.update({seq:sample_names})

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


data = []
for seq in sequences:
    row = [1 if file in seq_samples_name[seq] else 0 for file in sorted_files]
    data.append(row)

df = pd.DataFrame(data, index=sequences, columns=sorted_files)

plt.figure(figsize=(10, 8))
sns.heatmap(df, cmap="magma", cbar=False)
plt.title("Sequences Present in Reference Genomes")
plt.xlabel("Reference Genomes")
plt.ylabel(f"Sequences, n={num_all_seqs}")
plt.xticks(rotation=90)
plt.yticks([])
plt.tight_layout()
output_path = "/home/drc/Downloads/refgenomes_pics_test/01may2025/"
output_path = os.path.join(output_path, 'sequences_presence')
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.show()

# from operator import itemgetter
# sorted_dict_ascending_ordered = dict(sorted(sequence_counts.items(), key=itemgetter(1), reverse=True))
# #print("Sorted dictionary (ascending):", sorted_dict_ascending_ordered)

# # Now that the dictionary is sorted

# #print(single_sequence_in_fp['SQ.2LEWMcieZGf9Sx4VpEeWSDcULUVHGm0w'])

# previous_value = None
# dict_with_intersections = {}
# temp_set = None
# orphaned_keys = []
# for key, value in sorted_dict_ascending_ordered.items():
#     if value!=previous_value:
#         if temp_set:
#             dict_with_intersections.update({previous_value:temp_set})
#         previous_value=value
#         temp_set=set(single_sequence_in_fp[key])
#     elif value==previous_value:
#         tmp_temp_set = temp_set.intersection(set(single_sequence_in_fp[key]))
#         #if temp_set:
#         if len(tmp_temp_set)==0:
#             print(f"WARNING, last intersection created 0 length with: {key} {value}")
#             orphaned_keys.append(key)
#             # just skip it
#         else:
#             temp_set = tmp_temp_set

# dict_with_intersections.update({previous_value:temp_set})

# print(dict_with_intersections)

# list_tuples = []

# phc = PEPHubClient()
# pep = phc.load_project("donaldcampbelljr/human_seqcol_digests:default")


# for key, value in dict_with_intersections.items():
#     temp_set = None
#     basename_digests = []
#     for fp in value:
#         with open(fp, "r") as f:
#             reloaded_dict1 = json.load(fp=f)
#             sequences = reloaded_dict1['sorted_sequences']
#             if temp_set is None:
#                 temp_set = set(sequences)
#             else:
#                 temp_set = temp_set.intersection(set(sequences))
#                 if temp_set is None:
#                     print("warning temp_set was set to NONE")
#         filename_with_extension = os.path.basename(fp)
#         basename_without_extension = filename_with_extension.rsplit('.', 1)[0]
#         retrieved_sample_name = None
#         for sample in pep.samples:
#             if sample.top_level_digest == basename_without_extension:
#                 #basename_digests.append(basename_without_extension)
#                 basename_digests.append(sample.sample_name)

#     list_tuples.append((basename_digests,temp_set))
#             # all_sequences_union = all_sequences_union.union(set(sequences))
#             # # dict_sequences.update({os.path.basename(fp):sequences})
#             # dict_sequences.update({fp: sequences})

# #print(list_tuples)


# import networkx as nx
# import matplotlib.pyplot as plt

# G = nx.Graph()

# for strings, sequences in list_tuples:
#     for s in strings:
#         string_node = f"string_{s}"
#         G.add_node(string_node, type="string")  # Label nodes with type

#     if sequences:
#         concat_seq = "".join(sorted(sequences))  # Combine sequences, sorted for consistency
#         joined_seq_digest = ga4gh_digest(concat_seq)  # Calculate digest
#         digest_node_name = f"sequence_{joined_seq_digest}"  # Consistent node naming
#         if digest_node_name not in G.nodes():
#             G.add_node(digest_node_name, type="sequence")
#             print(f"Added sequence node: {digest_node_name}")
#         else:
#             print(f"Sequence node already exists: {digest_node_name}")

#         for s in strings:
#             string_node = f"string_{s}" #redefine string_node
#             G.add_edge(string_node, digest_node_name)  # Connect string to digest node

# # Visualize the graph
# pos = nx.spring_layout(G,k=0.9)  # Define node positions

# # # Color nodes based on their type
# node_color = [
#     "lightblue" if data["type"] == "string" else "lightcoral" for node, data in G.nodes(data=True)
# ]

# # Get node degrees
# node_degrees = dict(G.degree())


# #node_labels = {node: f"{node}\n(Degree: {degree})" for node, degree in node_degrees.items()}
# node_labels = {}
# key = {}
# for node, degree in node_degrees.items():
#     if "string" in node:
#         short_label = node.split("_")[1][:10]  # e.g., "apple"
#         node_labels[node] = f"{short_label}\n(Degree: {degree})"
#         key[node] = f"{node}"
#     elif "sequence" in node:
#         short_label = node.split(":")[1][:10] # First 8 characters of digest.
#         node_labels[node] = f"{short_label}\n(Degree: {degree})"
#         key[node] = f"{node}"

# nx.draw(G, pos,with_labels=True,
#         node_color=node_color, labels=node_labels,bbox=dict(facecolor='gray', alpha=0.1, pad=5, edgecolor='black'))
# plt.title("Network Graph of Strings and Sequences")
# plt.show()


# print(f"HERE ARE ORPHAN COUNTS: {len(orphaned_keys)}")

# # # Plot the adjacency matrix
# nodelist = sorted(G.nodes())
# adj_matrix = nx.to_numpy_array(G, nodelist=nodelist)
# adj_df = pd.DataFrame(adj_matrix, index=nodelist, columns=nodelist)

# plt.figure(figsize=(10, 10))
# sns.heatmap(adj_df, cmap='binary', annot=False, cbar=False, square=True)
# plt.title("Adjacency Matrix of the Network")
# plt.xlabel("Nodes")
# plt.ylabel("Nodes")
# plt.xticks(rotation=90, fontsize=8)
# plt.yticks(rotation=0, fontsize=8)
# plt.tight_layout()
# plt.show()



#Plot adjacency matrix and plot
#-------------------------------------------
# import matplotlib.cm as cm
# import matplotlib.colors as colors
# # Separate lists of string and sequence nodes
# string_nodes = sorted([node for node, data in G.nodes(data=True) if data.get('type') == 'string'])
# sequence_nodes = sorted([node for node, data in G.nodes(data=True) if data.get('type') == 'sequence'])
#
# # Create the nodelist for the adjacency matrix with strings on y and sequences on x
# nodelist_y = string_nodes
# nodelist_x = sequence_nodes
# combined_nodelist = nodelist_y + nodelist_x
#
# # Create the adjacency matrix based on the combined nodelist
# adj_matrix = nx.to_numpy_array(G, nodelist=combined_nodelist)
# adj_df = pd.DataFrame(adj_matrix, index=combined_nodelist, columns=combined_nodelist)
#
# # Select the submatrix corresponding to strings (rows) and sequences (columns)
# adj_df_reordered = adj_df.loc[nodelist_y, nodelist_x]
#
# # Get node degrees for coloring (using the original G and all nodes)
# node_degrees = dict(G.degree())
# max_degree = max(node_degrees.values()) if node_degrees else 1
# min_degree = min(node_degrees.values()) if node_degrees else 0
# norm = colors.Normalize(vmin=min_degree, vmax=max_degree)
# cmap = cm.viridis
#
# plt.figure(figsize=(len(sequence_nodes) * 0.6, len(string_nodes) * 0.6)) # Adjust figure size
#
# sns.heatmap(adj_df_reordered, cmap='binary', annot=False, cbar=False, square=False,
#             linewidths=0.5, linecolor='lightgray', xticklabels=True, yticklabels=True)
#
# # Overlay colored rectangles based on degree
# ax = plt.gca()
# for i, string_node in enumerate(nodelist_y):
#     for j, sequence_node in enumerate(nodelist_x):
#         if adj_df_reordered.loc[string_node, sequence_node] == 1:
#             degree_string = node_degrees.get(string_node, 0)
#             degree_sequence = node_degrees.get(sequence_node, 0)
#             avg_degree = (degree_string + degree_sequence) / 2
#             color = cmap(norm(avg_degree))
#             ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor=color, edgecolor=None))
#
# plt.title("Adjacency Matrix (Strings vs Sequences)")
# plt.xlabel("Sequences")
# plt.ylabel("Strings")
# plt.xticks(rotation=90, fontsize=8)
# plt.yticks(rotation=0, fontsize=8)
# plt.tight_layout()
#
# # Add a colorbar
# sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
# sm.set_array([])
# cbar = plt.colorbar(sm, ax=ax, label='Average Node Degree')
#
# plt.show()

# PLOT BAR CHART
# --------------------

# # Separate nodes into strings and sequences
# string_nodes = [node for node, data in G.nodes(data=True) if data.get('type') == 'string']
# sequence_nodes = [node for node, data in G.nodes(data=True) if data.get('type') == 'sequence']

# # Get degrees for string nodes
# string_degrees = {node: G.degree(node) for node in string_nodes}
# string_df = pd.DataFrame(list(string_degrees.items()), columns=['Node', 'Degree'])
# # Sort string degrees in descending order
# string_df = string_df.sort_values(by='Degree', ascending=False)

# # Get degrees for sequence nodes
# sequence_degrees = {node: G.degree(node) for node in sequence_nodes}
# sequence_df = pd.DataFrame(list(sequence_degrees.items()), columns=['Node', 'Degree'])
# sequence_df = sequence_df.sort_values(by='Degree', ascending=False)

# # --- Plotting String Degrees ---
# plt.figure(figsize=(12, 6))  # Adjust figure size for better readability
# plt.bar(string_df['Node'], string_df['Degree'])
# plt.xlabel("Ref Genome Nodes")
# plt.ylabel("Degree")
# plt.title("Degree Distribution of Ref Genome Nodes")
# plt.xticks(rotation=90, fontsize=8)  # Rotate x-axis labels for readability
# plt.tight_layout()  # Adjust layout to prevent labels from overlapping
# plt.show()

# # --- Plotting Sequence Degrees ---
# plt.figure(figsize=(12, 6))  # Adjust figure size
# plt.bar(sequence_df['Node'], sequence_df['Degree'])
# plt.xlabel("Sequence Nodes")
# plt.ylabel("Degree")
# plt.title("Degree Distribution of Sequence Nodes")
# plt.xticks(rotation=90, fontsize=8)  # Rotate x-axis labels
# plt.tight_layout()
# plt.show()

#
# print(res)

# import networkx as nx
# import matplotlib.pyplot as plt
#
# list_of_tuples = [
#     (["apple", "banana", "cherry"], {"apple_seq", "banana_fruit", "grape_seq"}),
#     (["date", "elderberry"], {"date_fruit", "berry_seq"}),
#     (["fig", "grape"], {"fig_fruit", "grape_vine"}),
# ]
#
# G = nx.Graph()
#
# # Add nodes and edges based on your relationship criteria
# for strings, sequences in list_of_tuples:
#     for s in strings:
#         G.add_node(f"string_{s}", type="string")  # Label nodes with type
#         for seq in sequences:
#             if s in seq:  # Example relationship: string is a substring of sequence
#                 G.add_edge(f"string_{s}", f"sequence_{seq}")
#             # You could add other relationship conditions here
#
#     for seq in sequences:
#         G.add_node(f"sequence_{seq}", type="sequence")
#
# # Visualize the graph
# pos = nx.spring_layout(G)  # Define node positions
#
# # Color nodes based on their type
# node_color = [
#     "lightblue" if data["type"] == "string" else "lightcoral" for node, data in G.nodes(data=True)
# ]
#
# nx.draw(G, pos, with_labels=True, node_color=node_color, node_size=2000, font_size=10, font_weight="bold")
# nx.draw_networkx_edge_labels(G, pos, edge_labels=nx.get_edge_attributes(G, "weight")) # If you have edge weights
#
# plt.title("Network Graph of Strings and Sequences")
# plt.show()