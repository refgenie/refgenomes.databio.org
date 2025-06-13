# just print some differences in seqs for illustration purposes
import sys
import os
import urllib.request
import urllib.error
import pipestat
import json
from refget import fasta_to_digest, fasta_to_seqcol_dict,fasta_to_seq_digests,compare_seqcols
import pypiper
import pprint
import peppy
from pephubclient import PEPHubClient
from itertools import combinations
from pprint import pprint
import pandas as pd


#----
def get_dict_seq_col_from_json(digest):

    # OPENS A LOCAL JSON ONLY
    json_fp_1 = os.path.join(LOCAL_JSON_DIRECTORY, digest+".json")
    with open(json_fp_1, "r") as f:
        reloaded_dict1 = json.load(fp=f)
    
    return reloaded_dict1

# Function to create a DataFrame from a list of lists, handling unequal lengths
def create_dataframe(data_list, columns):
    # Pad shorter lists with None to make them equal length for DataFrame creation
    max_len = max(len(sublist) for sublist in data_list)
    padded_data = [sublist + [None] * (max_len - len(sublist)) for sublist in data_list]
    return pd.DataFrame(padded_data).T.set_axis(columns, axis=1)

#----

#looper_config = sys.argv[1]  
#LOCAL_JSON_DIRECTORY  = sys.argv[2] # for lcally stored jsons to pull sequence col data

looper_config = "donaldcampbelljr/human_seqcol_digests:default"
LOCAL_JSON_DIRECTORY = "/home/drc/Downloads/jsons_from_rivanna/json/"
output_folder ="/home/drc/Downloads/refgenomes_pics_test/13june2025/CSV/"

phc = PEPHubClient()
pep = phc.load_project(looper_config)
#print(pep["_sample_df"])

pep_df = pep["_sample_df"]



target_samples = [
"hg38-toplevel-113-ensembl",
"GRCh38-p14-47-gencode",
"GRCh38.p14-fasta-genomic",
"hg38-p14-ucsc",
"hg38-ddbj",
"GRCh38-ena-29",
]

# target_samples = [
# "hg38-toplevel-113-ensembl",
# "GRCh38.p14-fasta-full-analysis",
# "GRCh38.p14-fasta-full-analysis-plus-hs38d1",
# "hg19-p13-full-analysis-ucsc",
# "hg19-p13-plusMT-ucsc",
# "hg38-alt-113-ensembl",
# "GRCh38-full-decoy-hla-ddbj",
# "homo-sapiens-assembly38-ddbj",
# "hg38-ddbj",
# ]

# # # # Pre-filter the DataFrame
pep_df = pep_df[
    ((pep_df['sample_name'].isin(target_samples)))
]

#print(pep_df)

names = []
seqs = []
lengths = []
coords = []
sample_name = []
strings = []

for index, row in pep_df.iterrows():
    reloaded_dict = get_dict_seq_col_from_json(row['top_level_digest'])
    #print(reloaded_dict.keys())
    sample_name.append(row['sample_name'])
    names.append(reloaded_dict['names'])
    seqs.append(reloaded_dict['sequences'])
    lengths.append(reloaded_dict['lengths'])
    coords.append(reloaded_dict['name_length_pairs'])



# for i in range(len(sample_name)):
#     all_string_rows = []
#     for k in range(len(names[i])):
#         new_string = str(names[k]) + str(lengths[k]) + str(seqs[k])+ str(coords[k])
#         all_string_rows.append(new_string)
#     strings.append(all_string_rows)

# print(strings)


names_df = create_dataframe(names, sample_name)
#names_df = names_df.apply(lambda x: x.sort_values().values, axis=0)
seqs_df = create_dataframe(seqs, sample_name)
#seqs_df = seqs_df.apply(lambda x: x.sort_values().values, axis=0)
lengths_df = create_dataframe(lengths, sample_name)
#lengths_df = lengths_df.apply(lambda x: x.sort_values().values, axis=0)
coords_df = create_dataframe(coords, sample_name)

# Save each DataFrame to a CSV file
names_df.to_csv(os.path.join(output_folder, "names.csv"), index=False)
seqs_df.to_csv(os.path.join(output_folder, "sequences.csv"), index=False)
lengths_df.to_csv(os.path.join(output_folder, "lengths.csv"), index=False)
coords_df.to_csv(os.path.join(output_folder, "coords.csv"), index=False)

# print("Names DataFrame:")
# print(names_df)
# print("\nSequences DataFrame:")
# print(seqs_df)
# print("\nLengths DataFrame:")
# print(lengths_df)
# print("\nCoordinates DataFrame:")
# print(coords_df)

# print("FINISHED")


