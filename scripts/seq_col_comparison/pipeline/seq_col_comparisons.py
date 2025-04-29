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

# looper_config = sys.argv[1]  
# results_pep = sys.argv[2]
looper_config = "donaldcampbelljr/human_seqcol_digests_local:default"  # input PEP
results_pep = "donaldcampbelljr/test_seq_col_results:default"

print(f"here is the looper config: {looper_config}")


# Calculation Functions
def calc_overlap_coeff(A,B, A_B_intersection):
    minimum_value = min(abs(A),abs(B))
    overlap = abs(A_B_intersection)/minimum_value
    return overlap


def calc_overlap_coeff_MAX(A,B, A_B_intersection):
    minimum_value = max(abs(A),abs(B))
    overlap = abs(A_B_intersection)/minimum_value
    return overlap


def calc_jaccard_similarity(A_B_intersection, A_B_union):
    jaccard = abs(A_B_intersection)/abs(A_B_union)
    return jaccard

# def calc_f1_score(tp, fp, fn):
#     # another way to consider:
#     # tp = a = the number of attrivutes that equal ` for both objects i and j
#     # fp = b = the number of attributes that equal 0 for object i but equal 1 for object j
#     # fn = c = the number of attributes that equal 1 for object i but equal 0 for object j
#     # d = the number of attributes that equal 0 for both objects i and j but we do not need to consider these.
#     return (2*tp)/(2*tp+fp+fn)

def calc_weighted_jaccard(list_intersection_lengths, list_union_lengths):
    intersection_total = 0
    union_total = 0

    for length in list_intersection_lengths:
        if length:
            intersection_total += length
    
    for length in list_union_lengths:
        if length:
            union_total += length

    return intersection_total/union_total

# initiate pephubclient object
phc = PEPHubClient()
pep = phc.load_project(looper_config)
print(pep)

# Retrieve all json file paths from the pep
all_samples = []
key_digest_sample_name = {}
sample_json_path = {}

for sample in pep.samples:
    #all_jsons.append(sample.brickyard_json_path)
    key_digest_sample_name.update({sample.sample_name:sample.top_level_digest}) # handy later when we need to make this more human readable
    sample_json_path.update({sample.sample_name:sample.brickyard_json_path})
    all_samples.append(sample.sample_name)

print(all_samples)
print(key_digest_sample_name)

all_combinations = combinations(iterable=all_samples,r=2)


# psm_input = pipestat.PipestatManager(pephub_path=looper_config)
psm_output = pipestat.PipestatManager(pephub_path=results_pep)

# combos_to_consider = []

combination_count = 0
for combination in all_combinations:
    json_fp_1=sample_json_path[combination[0]]
    json_fp_2=sample_json_path[combination[1]]
    with open(json_fp_1, "r") as f:
        reloaded_dict1 = json.load(fp=f)
    with open(json_fp_2, "r") as f:
        reloaded_dict2 = json.load(fp=f)
 
    #digest1 = os.path.splitext(os.path.basename(json_fp_1))[0]
    digest1 = key_digest_sample_name[combination[0]]
    #digest2 = os.path.splitext(os.path.basename(json_fp_2))[0]
    digest2 = key_digest_sample_name[combination[1]]
    print(f"COMBINATION: samples: {combination[0]} vs {combination[1]} digests: {digest1} vs {digest2}")

    # print(pprint(compare_seqcols(reloaded_dict1,reloaded_dict2),indent=4))
    comparison = compare_seqcols(reloaded_dict1,reloaded_dict2)

    print(f"Reality Check: {comparison['array_elements']['a_and_b']['names']}  {comparison['array_elements']['a']['names']+ comparison['array_elements']['b']['names'] - comparison['array_elements']['a_and_b']['names']}")
    # overlap for names
    overlap_coefficient_names = calc_overlap_coeff(comparison['array_elements']['a']['names'],comparison['array_elements']['b']['names'],comparison['array_elements']['a_and_b']['names'])
    print(f"Here is the overlap coefficient for names: {overlap_coefficient_names}")
    overlap_coefficient_names_max = calc_overlap_coeff_MAX(comparison['array_elements']['a']['names'],comparison['array_elements']['b']['names'],comparison['array_elements']['a_and_b']['names'])

    # overlap for lengths
    overlap_coefficient_lengths = calc_overlap_coeff(comparison['array_elements']['a']['lengths'],comparison['array_elements']['b']['lengths'],comparison['array_elements']['a_and_b']['lengths'])
    print(f"Here is the overlap coefficient for lengths: {overlap_coefficient_names}")
    overlap_coefficient_lengths_max = calc_overlap_coeff_MAX(comparison['array_elements']['a']['lengths'],comparison['array_elements']['b']['lengths'],comparison['array_elements']['a_and_b']['lengths'])

    # jaccard similarity for names, intersection over union where we can find the union -> A+B-ABintersection
    jaccard_names = calc_jaccard_similarity(comparison['array_elements']['a_and_b']['names'],(comparison['array_elements']['a']['names']+comparison['array_elements']['b']['names']-comparison['array_elements']['a_and_b']['names']))
    print(f"Here is the jaccard similarity for names: {jaccard_names}")

    # jaccard similarity for lengths, intersection over union where we can find the union -> A+B-ABintersection
    jaccard_lengths = calc_jaccard_similarity(comparison['array_elements']['a_and_b']['lengths'],(comparison['array_elements']['a']['lengths']+comparison['array_elements']['b']['lengths']-comparison['array_elements']['a_and_b']['lengths']))
    print(f"Here is the jaccard similarity for lengths: {jaccard_lengths}")

    # Get sequences and calculate overlaps for sequences
    set_sequences_1 = set(reloaded_dict1['sorted_sequences'])
    set_sequences_2 = set(reloaded_dict2['sorted_sequences'])

    sequences_intersections = set_sequences_1.intersection(set_sequences_2)
    sequences_intersection_length = len (sequences_intersections)
    sequences_union = set_sequences_1.union(set_sequences_2)
    sequences_union_length = len(sequences_union)


    jaccard_sequences = calc_jaccard_similarity(sequences_intersection_length,sequences_union_length)
    overlap_coeff_sequences = calc_overlap_coeff(len(set_sequences_1),len(set_sequences_2),sequences_intersection_length)
    overlap_coeff_sequences_max =calc_overlap_coeff_MAX(len(set_sequences_1),len(set_sequences_2),sequences_intersection_length)
    print(f"Here is the jaccard similarity for sequences: {jaccard_sequences}")
    print(f"Here is the overlap coeff for sequences: {overlap_coeff_sequences}")

    # Build sets for calculating weighted jaccard score

    reloaded_dict1_name_length_dict = {}
    reloaded_dict2_name_length_dict = {}
    for i in range(0, len(reloaded_dict1['lengths'])):
        reloaded_dict1_name_length_dict.update({reloaded_dict1['names'][i]:reloaded_dict1['lengths'][i]})
    for i in range(0, len(reloaded_dict2['lengths'])):
        reloaded_dict2_name_length_dict.update({reloaded_dict2['names'][i]:reloaded_dict2['lengths'][i]})
    
    set1 = set(reloaded_dict1['names'])
    set2 = set(reloaded_dict2['names'])
    names_intersection = set1.intersection(set2)
    names_union = set1.union(set2)

    print(f"names intersection length: {len(names_intersection)}, names_union_length: {len(names_union)}")

    list_intersection_lengths = []
    list_union_lengths = []

    for name in names_intersection:
        if reloaded_dict1_name_length_dict.get(name)==reloaded_dict2_name_length_dict.get(name):
            list_intersection_lengths.append(reloaded_dict1_name_length_dict.get(name))
        else:
            pass

    for name in names_union:
        list_union_lengths.append(reloaded_dict1_name_length_dict.get(name))
        list_union_lengths.append(reloaded_dict2_name_length_dict.get(name))

    list_union_lengths_new_list = [x for x in list_union_lengths if x is not None]

    print(f"The union lengths vs names union: {len(list_union_lengths_new_list)} vs {len(names_union)}") 
    print(f"The intersection lengths vs names intersection: {len(list_intersection_lengths)} vs {len(names_intersection)}") 
    jaccard_similarity_weighted_length = calc_weighted_jaccard(list_intersection_lengths, list_union_lengths_new_list)
    print(f"Here is the weighted jaccard similarity: {jaccard_similarity_weighted_length}")


    # create set of name_length_pairs
    #print(type(reloaded_dict1['name_length_pairs'][0]))
    set_of_name_len_pairs_1 = {tuple(d.values()) for d in reloaded_dict1['name_length_pairs']}
    set_of_name_len_pairs_2 = {tuple(d.values()) for d in reloaded_dict2['name_length_pairs']}

    name_len_pairs_intersection = set_of_name_len_pairs_1.intersection(set_of_name_len_pairs_2)
    name_len_pairs_union = set_of_name_len_pairs_1.union(set_of_name_len_pairs_2)

    jaccard_name_len = calc_jaccard_similarity(len(name_len_pairs_intersection),len(name_len_pairs_union))
    overlap_coeff_name_len = calc_overlap_coeff(len(set_of_name_len_pairs_1),len(set_of_name_len_pairs_2),len(name_len_pairs_intersection))
    overlap_coeff_name_len_max = calc_overlap_coeff_MAX(len(set_of_name_len_pairs_1),len(set_of_name_len_pairs_2),len(name_len_pairs_intersection))
    print(f"Here is the jaccard similarity for sequences: {jaccard_sequences}")
    print(f"Here is the overlap coeff for sequences: {overlap_coeff_sequences}")

    
    # Calculate F1 scores 

    f1_names = 2*jaccard_names/(jaccard_names+1)
    f1_lengths = 2*jaccard_lengths/(jaccard_lengths+1)
    f1_sequences = 2*jaccard_sequences/(jaccard_sequences+1)
    f1_weighted_lengths = 2*jaccard_similarity_weighted_length/(jaccard_similarity_weighted_length+1)
    f1_name_len_pairs = 2*jaccard_name_len/(jaccard_name_len+1)

    comparison_str = combination[0] +"_vs_" + combination[1]
    psm_output.report(record_identifier=comparison_str, values={"digest1":digest1,"sample_name_1":combination[0],
                                                                "digest2":digest2,"sample_name_2":combination[1], 
                                                                "overlap_coefficient_names":overlap_coefficient_names,
                                                                "overlap_coefficient_lengths":overlap_coefficient_lengths, 
                                                                "jaccard_names":jaccard_names, "jaccard_lengths":jaccard_lengths, 
                                                                "jaccard_similarity_weighted_length":jaccard_similarity_weighted_length, 
                                                                "jaccard_sequences": jaccard_sequences, 
                                                                "overlap_coefficient_sequences": overlap_coeff_sequences,
                                                                "jaccard_name_len":jaccard_name_len,
                                                                "overlap_coeff_name_len":overlap_coeff_name_len,
                                                                "f1_names":f1_names,
                                                                "f1_lengths":f1_lengths,
                                                                "f1_sequences":f1_sequences,
                                                                "f1_weighted_lengths":f1_weighted_lengths,
                                                                "f1_name_len_pairs":f1_name_len_pairs,
                                                                "overlap_coefficient_names_max":overlap_coefficient_names_max,
                                                                "overlap_coefficient_lengths_max":overlap_coefficient_lengths_max,
                                                                "overlap_coefficient_sequences_max": overlap_coeff_sequences_max, 
                                                                "overlap_coeff_name_len_max":overlap_coeff_name_len_max,
                                                                })
    combination_count+=1

    # if len(sequences_intersections) >0:
    #     temp_dict = {"comparison str": comparison_str, "s1":combination[0],"s2":combination[1], "intersection_seq": sequences_intersections, "len_intersection_seq": len(sequences_intersections)}

    # print(f"Here is temp dict: {temp_dict}")

print (f"Finished with {combination_count} combinations processed")
                                               
