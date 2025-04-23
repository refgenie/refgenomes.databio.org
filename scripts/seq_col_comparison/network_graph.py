import os
from pephubclient import PEPHubClient
from refget import ga4gh_digest
import json
import networkx as nx
import matplotlib.pyplot as plt


json_files = []
if os.path.isdir("/home/drc/Downloads/jsons_from_rivanna/json/"):
    for filename in os.listdir("/home/drc/Downloads/jsons_from_rivanna/json/"):
        if filename.endswith(".json"):
            full_path = os.path.join("/home/drc/Downloads/jsons_from_rivanna/json/", filename)
            json_files.append(full_path)


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

from operator import itemgetter
sorted_dict_ascending_ordered = dict(sorted(sequence_counts.items(), key=itemgetter(1), reverse=True))


previous_value = None
dict_with_intersections = {}
temp_set = None
orphaned_keys = []
for key, value in sorted_dict_ascending_ordered.items():
    if value!=previous_value:
        if temp_set:
            dict_with_intersections.update({previous_value:temp_set})
        previous_value=value
        temp_set=set(single_sequence_in_fp[key])
    elif value==previous_value:
        tmp_temp_set = temp_set.intersection(set(single_sequence_in_fp[key]))
        #if temp_set:
        if len(tmp_temp_set)==0:
            print(f"WARNING, last intersection created 0 length with: {key} {value}")
            orphaned_keys.append(key)
            # just skip it
        else:
            temp_set = tmp_temp_set

dict_with_intersections.update({previous_value:temp_set})

print(dict_with_intersections)

list_tuples = []

phc = PEPHubClient()
pep = phc.load_project("donaldcampbelljr/human_seqcol_digests:default")


for key, value in dict_with_intersections.items():
    temp_set = None
    basename_digests = []
    for fp in value:
        with open(fp, "r") as f:
            reloaded_dict1 = json.load(fp=f)
            sequences = reloaded_dict1['sorted_sequences']
            if temp_set is None:
                temp_set = set(sequences)
            else:
                temp_set = temp_set.intersection(set(sequences))
                if temp_set is None:
                    print("warning temp_set was set to NONE")
        filename_with_extension = os.path.basename(fp)
        basename_without_extension = filename_with_extension.rsplit('.', 1)[0]
        retrieved_sample_name = None
        for sample in pep.samples:
            if sample.top_level_digest == basename_without_extension:
                #basename_digests.append(basename_without_extension)
                basename_digests.append(sample.sample_name)

    list_tuples.append((basename_digests,temp_set))
            # all_sequences_union = all_sequences_union.union(set(sequences))
            # # dict_sequences.update({os.path.basename(fp):sequences})
            # dict_sequences.update({fp: sequences})

#print(list_tuples)

G = nx.Graph()

for strings, sequences in list_tuples:
    for s in strings:
        string_node = f"string_{s}"
        G.add_node(string_node, type="string")  

    if sequences:
        concat_seq = "".join(sorted(sequences)) 
        joined_seq_digest = ga4gh_digest(concat_seq)  
        digest_node_name = f"sequence_{joined_seq_digest}"  
        if digest_node_name not in G.nodes():
            G.add_node(digest_node_name, type="sequence")
            print(f"Added sequence node: {digest_node_name}")
        else:
            print(f"Sequence node already exists: {digest_node_name}")

        for s in strings:
            string_node = f"string_{s}" 
            G.add_edge(string_node, digest_node_name) 


pos = nx.spring_layout(G,k=0.9)  

node_color = [
    "lightblue" if data["type"] == "string" else "lightcoral" for node, data in G.nodes(data=True)
]

node_degrees = dict(G.degree())

#node_labels = {node: f"{node}\n(Degree: {degree})" for node, degree in node_degrees.items()}
node_labels = {}
key = {}
for node, degree in node_degrees.items():
    if "string" in node:
        short_label = node.split("_")[1][:10]  # e.g., "apple"
        node_labels[node] = f"{short_label}\n(Degree: {degree})"
        key[node] = f"{node}"
    elif "sequence" in node:
        short_label = node.split(":")[1][:10] # First 8 characters of digest.
        node_labels[node] = f"{short_label}\n(Degree: {degree})"
        key[node] = f"{node}"

nx.draw(G, pos,with_labels=True,
        node_color=node_color, labels=node_labels,bbox=dict(facecolor='gray', alpha=0.1, pad=5, edgecolor='black'))
plt.title("Network Graph of Strings and Sequences")
plt.show()


print(f"HERE ARE ORPHAN COUNTS: {len(orphaned_keys)}")