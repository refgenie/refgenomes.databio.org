import sys
import os
import urllib.request
import urllib.error
import pipestat
from refget import fasta_to_digest, fasta_to_seqcol_dict,fasta_to_seq_digests
import pypiper
import pprint
import peppy
from pephubclient import PEPHubClient
from itertools import combinations

#looper_config = sys.argv[1]  
looper_config = "donaldcampbelljr/human_seqcol_digests_local:default"
# digest = sys.argv[2]
# file_path = sys.argv[3]
# looper_output_dir = sys.argv[4]

# print(f"Here is the file_path: {file_path}")

# exists = os.path.exists(file_path)
# print(f"File exists: {exists}")

# # WE WANT sha512t24u due to robustness wrt collisions

# # for x in fasta_to_seq_digests(file_path):
# #     print(f"{x.id}\t{x.length}\t{x.sha512t24u}\t{x.md5}")

# seq_dict = fasta_to_seqcol_dict(file_path)

# pprint.pprint(seq_dict)

print(f"here is the looper config: {looper_config}")


# initiate pephubclient object
phc = PEPHubClient()
pep = phc.load_project(looper_config)
print(pep)


# GET SEQ_COL_DICTS
# Create a dictionary holding seq collection data for all the samples.
# In the future, should potentially pull this from a database instead of calculating on the fly.
seq_col_all={}
for sample in pep.samples:
    print(f"Here is the sample name: {sample.sample_name} here is top_level_digest: {sample.top_level_digest} Here is brickyard location: {sample.brickyard_location} ")
    #seq_dict = fasta_to_seqcol_dict(sample.brickyard_location)
    seq_dict = {}
    seq_col_all.update({sample.top_level_digest: seq_dict})
    #pprint.pprint(seq_dict)

# Once we have all the seq_col dictionaries computed, we can then do a pair-wise comparison

all_digests_list = list(seq_col_all.keys())

print(all_digests_list)

all_combinations = combinations(iterable=all_digests_list,r=2) # r is for the number of elements in each combination. Just want pairwise for now.

print("combination iterated")
for combination in all_combinations:
    print(combination)