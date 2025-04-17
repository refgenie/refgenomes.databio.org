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


def calc_jaccard_similarity(A_B_intersection, A_B_union):
    jaccard = abs(A_B_intersection)/abs(A_B_union)
    return jaccard

def calc_weighted_jaccard(list_intersection_lengths, list_union_lengths):
    intersection_total = 0
    union_total = 0

    for length in list_intersection_lengths:
        intersection_total += length
    
    for length in list_union_lengths:
        union_total += length

    return intersection_total/union_total

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

    print(f"Reality Check: {comparison['array_elements']['a_and_b']['names']}  {comparison['array_elements']['a']['names']+ comparison['array_elements']['b']['names'] - comparison['array_elements']['a_and_b']['names']}")
    # overlap for names
    overlap_coefficient_names = calc_overlap_coeff(comparison['array_elements']['a']['names'],comparison['array_elements']['b']['names'],comparison['array_elements']['a_and_b']['names'])
    print(f"Here is the overlap coefficient for names: {overlap_coefficient_names}")

    # overlap for lengths
    overlap_coefficient_lengths = calc_overlap_coeff(comparison['array_elements']['a']['lengths'],comparison['array_elements']['b']['lengths'],comparison['array_elements']['a_and_b']['lengths'])
    print(f"Here is the overlap coefficient for lengths: {overlap_coefficient_names}")

    # jaccard similarity for names, intersection over union where we can find the union -> A+B-ABintersection
    jaccard_names = calc_jaccard_similarity(comparison['array_elements']['a_and_b']['names'],(comparison['array_elements']['a']['names']+comparison['array_elements']['b']['names']-comparison['array_elements']['a_and_b']['names']))
    print(f"Here is the jaccard similarity for names: {jaccard_names}")

    # jaccard similarity for lengths, intersection over union where we can find the union -> A+B-ABintersection
    jaccard_lengths = calc_jaccard_similarity(comparison['array_elements']['a_and_b']['lengths'],(comparison['array_elements']['a']['lengths']+comparison['array_elements']['b']['lengths']-comparison['array_elements']['a_and_b']['lengths']))
    print(f"Here is the jaccard similarity for lengths: {jaccard_lengths}")


    # -------- attempt to get unique lengths vs all shared lengths
    lengths_1 = set(reloaded_dict1['lengths'])
    lengths_2 = set(reloaded_dict2['lengths'])

    lengths_union = lengths_1.union(lengths_2)
    lengths_intersection = lengths_1.intersection(lengths_2)

    jaccard_similarity_weighted_length = calc_weighted_jaccard(lengths_intersection, lengths_union)

    print(f"Here is the weighted jaccard similarity: {jaccard_similarity_weighted_length}")

    # -------- attempt to get unique lengths vs all shared lengths
    # name_length_pairs1 = set(reloaded_dict1['name_length_pairs'])
    # name_length_pairs2 = set(reloaded_dict2['name_length_pairs'])

    # name_length_pairs_intersection_list = list(name_length_pairs1.intersection(name_length_pairs2))
    # name_length_pairs_union_list = list(name_length_pairs1.union(name_length_pairs2))

    # jaccard_similarity_weighted_length = calc_weighted_jaccard(name_length_pairs_intersection_list, name_length_pairs_union_list)

    # print(f"Here is the weighted jaccard similarity: {jaccard_similarity_weighted_length}")

    # build chr1 name : length dicts

    # reloaded_dict1_name_length_dict = {}
    # reloaded_dict2_name_length_dict = {}
    # for i in range(0, len(reloaded_dict1['lengths'])):
    #     reloaded_dict1_name_length_dict.update({reloaded_dict1['names'][i]:reloaded_dict1['lengths'][i]})
    # for i in range(0, len(reloaded_dict2['lengths'])):
    #     reloaded_dict2_name_length_dict.update({reloaded_dict2['names'][i]:reloaded_dict2['lengths'][i]})
    
    # # Build sets for calcs
    # set1 = set(reloaded_dict1['names'])
    # set2 = set(reloaded_dict2['names'])
    # names_intersection = set1.intersection(set2)
    # names_union = set1.union(set2)

    # list_intersection_lengths = []
    # list_union_lengths = []

    # for name in names_intersection:
    #     if reloaded_dict1_name_length_dict.get(name)==reloaded_dict2_name_length_dict.get(name):
    #         list_intersection_lengths.append(reloaded_dict1_name_length_dict.get(name))
    #     else:
    #         pass

    # for name in names_union:
    #     if reloaded_dict1_name_length_dict.get(name)==reloaded_dict2_name_length_dict.get(name):
    #         list_union_lengths.append(reloaded_dict1_name_length_dict.get(name))
    #     else:
    #         list_union_lengths.append(reloaded_dict1_name_length_dict.get(name))
    #         list_union_lengths.append(reloaded_dict2_name_length_dict.get(name))


    # if len(list_union_lengths) != len(names_union):
    #     print(f"The unions were not equal: {len(list_union_lengths)} vs {len(names_union)}") 
    # if names_union:
    #     print(f"names union = {len(names_union)}")
    # if names_intersection:
    #     print(f"names intersection = {len(names_intersection)}")
    

    # pprint(reloaded_dict1_name_length_dict)
    # pprint(reloaded_dict2_name_length_dict)

                                               
