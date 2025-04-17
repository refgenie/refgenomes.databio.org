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

#looper_config = sys.argv[1]  
looper_config = "donaldcampbelljr/human_seqcol_digests_local:default"
# digest = sys.argv[2]
# file_path = sys.argv[3]
# looper_output_dir = sys.argv[4]

print(f"here is the looper config: {looper_config}")


# Calculation Functions
def calc_overlap_coeff(A,B, A_B_intersection):
    minimum_value = min(abs(A),abs(B))
    overlap = abs(A_B_intersection)/minimum_value
    return overlap


def calc_jaccard_similarity():
    pass

def calc_weighted_jaccard():
    pass

# initiate pephubclient object
phc = PEPHubClient()
pep = phc.load_project(looper_config)
print(pep)

# Retrieve all json file paths from the pep
all_jsons = []
key_digest_sample_name = {}

for sample in pep.samples:
    all_jsons.append(sample.brickyard_json_path)
    key_digest_sample_name.update({sample.top_level_digest:sample.sample_name}) # handy later when we need to make this more human readable

print(all_jsons)
print(key_digest_sample_name)

all_combinations = combinations(iterable=all_jsons,r=2)

for combination in all_combinations:
    json_fp_1=combination[0]
    json_fp_2=combination[1]
    with open(json_fp_1, "r") as f:
        reloaded_dict1 = json.load(fp=f)
    with open(json_fp_2, "r") as f:
        reloaded_dict2 = json.load(fp=f)
#
    print(f"COMBINATION: {os.path.basename(json_fp_1)} vs {os.path.basename(json_fp_2)}")
    # print(pprint(compare_seqcols(reloaded_dict1,reloaded_dict2),indent=4))
    comparison = compare_seqcols(reloaded_dict1,reloaded_dict2)

    # overlap for names
    overlap_coefficient_names = calc_overlap_coeff(comparison['array_elements']['a']['names'],comparison['array_elements']['b']['names'],comparison['array_elements']['a_and_b']['names'])
    print(f"Here is the overlap coefficient for names: {overlap_coefficient_names}")




