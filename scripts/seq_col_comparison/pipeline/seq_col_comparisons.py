import sys
import pipestat
import json
from refget import compare_seqcols, calc_jaccard_similarities
from pephubclient import PEPHubClient
from itertools import combinations

#looper_config = sys.argv[1]  
#results_pep = sys.argv[2]
looper_config = "donaldcampbelljr/human_seqcol_digests_local:default"  # input PEP
results_pep = "donaldcampbelljr/test_seq_col_results:default"

print(f"here is the looper config: {looper_config}")


# Calculation Functions
def calc_jaccard_similarity(A_B_intersection, A_B_union):
    jaccard = abs(A_B_intersection)/abs(A_B_union)
    return jaccard

def overlap_proportion(A_B_intersection, a_or_b):
    # calculate the intersection OVER one of the two sets used to calculate the intersection
    proportion = abs(A_B_intersection)/abs(a_or_b)
    return proportion

# initiate pephubclient object
phc = PEPHubClient()
pep = phc.load_project(looper_config)
print(pep)

# Retrieve all json file paths from the pep
all_samples = []
key_digest_sample_name = {}
sample_json_path = {}

for sample in pep.samples:
    key_digest_sample_name.update({sample.sample_name:sample.top_level_digest}) # handy later when we need to make this more human readable
    sample_json_path.update({sample.sample_name:sample.brickyard_json_path})
    all_samples.append(sample.sample_name)

print(all_samples)
print(key_digest_sample_name)

all_combinations = combinations(iterable=all_samples,r=2)


# psm_input = pipestat.PipestatManager(pephub_path=looper_config)
psm_output = pipestat.PipestatManager(pephub_path=results_pep)

ALL_SEQ_COL_DATA = {}

combination_count = 0
for combination in all_combinations:
    json_fp_1=sample_json_path[combination[0]]
    json_fp_2=sample_json_path[combination[1]]

    if ALL_SEQ_COL_DATA.get(combination[0]):
        reloaded_dict1 = ALL_SEQ_COL_DATA[combination[0]]
    else:
        with open(json_fp_1, "r") as f:
            reloaded_dict1 = json.load(fp=f)
            ALL_SEQ_COL_DATA[combination[0]] = reloaded_dict1
    
    if ALL_SEQ_COL_DATA.get(combination[1]):
        reloaded_dict2 =  ALL_SEQ_COL_DATA[combination[1]]
    else:
        with open(json_fp_2, "r") as f:
            reloaded_dict2 = json.load(fp=f)
            ALL_SEQ_COL_DATA[combination[1]] = reloaded_dict2
 
    digest1 = key_digest_sample_name[combination[0]]
    digest2 = key_digest_sample_name[combination[1]]
    print(f"COMBINATION: samples: {combination[0]} vs {combination[1]} digests: {digest1} vs {digest2}")

    # print(pprint(compare_seqcols(reloaded_dict1,reloaded_dict2),indent=4))
    comparison = compare_seqcols(reloaded_dict1,reloaded_dict2)
    calculated_jaccards = calc_jaccard_similarities(reloaded_dict1,reloaded_dict2)
    #print(f"Reality Check: {comparison['array_elements']['a_and_b']['names']}  {comparison['array_elements']['a']['names']+ comparison['array_elements']['b']['names'] - comparison['array_elements']['a_and_b']['names']}")
    # OPA/OPB Names
    opa_names = overlap_proportion(comparison['array_elements']['a_and_b']['names'], comparison['array_elements']['a']['names'])
    opb_names = overlap_proportion(comparison['array_elements']['a_and_b']['names'], comparison['array_elements']['b']['names'])

    # OPA/OPB Lengths
    opa_lengths = overlap_proportion(comparison['array_elements']['a_and_b']['lengths'], comparison['array_elements']['a']['lengths'])
    opb_lengths = overlap_proportion(comparison['array_elements']['a_and_b']['lengths'], comparison['array_elements']['b']['lengths'])

    # jaccard similarity for names, intersection over union where we can find the union -> A+B-ABintersection
    jaccard_names = calc_jaccard_similarity(comparison['array_elements']['a_and_b']['names'],(comparison['array_elements']['a']['names']+comparison['array_elements']['b']['names']-comparison['array_elements']['a_and_b']['names']))
    #print(f"Here is the jaccard similarity for names: {jaccard_names}")

    # jaccard similarity for lengths, intersection over union where we can find the union -> A+B-ABintersection
    jaccard_lengths = calc_jaccard_similarity(comparison['array_elements']['a_and_b']['lengths'],(comparison['array_elements']['a']['lengths']+comparison['array_elements']['b']['lengths']-comparison['array_elements']['a_and_b']['lengths']))
    #print(f"Here is the jaccard similarity for lengths: {jaccard_lengths}")

    # Get sequences and calculate overlaps for sequences
    set_sequences_1 = set(reloaded_dict1['sorted_sequences'])
    set_sequences_2 = set(reloaded_dict2['sorted_sequences'])

    sequences_intersections = set_sequences_1.intersection(set_sequences_2)
    sequences_intersection_length = len (sequences_intersections)
    sequences_union = set_sequences_1.union(set_sequences_2)
    sequences_union_length = len(sequences_union)


    jaccard_sequences = calc_jaccard_similarity(sequences_intersection_length,sequences_union_length)

    opa_sequences = overlap_proportion(sequences_intersection_length,len(set_sequences_1))
    opb_sequences = overlap_proportion(sequences_intersection_length,len(set_sequences_2))


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

    set1_l = set(reloaded_dict1['lengths'])
    set2_l = set(reloaded_dict2['lengths'])
    lengths_intersection = set1_l.intersection(set2_l)
    lengths_union = set1_l.union(set2_l)


    # create set of name_length_pairs
    set_of_name_len_pairs_1 = {tuple(d.values()) for d in reloaded_dict1['name_length_pairs']}
    set_of_name_len_pairs_2 = {tuple(d.values()) for d in reloaded_dict2['name_length_pairs']}

    name_len_pairs_intersection = set_of_name_len_pairs_1.intersection(set_of_name_len_pairs_2)
    name_len_pairs_union = set_of_name_len_pairs_1.union(set_of_name_len_pairs_2)

    jaccard_name_len = calc_jaccard_similarity(len(name_len_pairs_intersection),len(name_len_pairs_union))

    opa_name_len = overlap_proportion(len(name_len_pairs_intersection),len(set_of_name_len_pairs_1))
    opb_name_len = overlap_proportion(len(name_len_pairs_intersection),len(set_of_name_len_pairs_2))

    print(f"Here is the jaccard similarity for sequences: {jaccard_sequences}")

    comparison_str = combination[0] +"_vs_" + combination[1]
    psm_output.report(record_identifier=comparison_str, values={"digest1":digest1,"sample_name_1":combination[0],
                                                                "digest2":digest2,"sample_name_2":combination[1], 
                                                                "jaccard_names":jaccard_names, 
                                                                "jaccard_lengths":jaccard_lengths, 
                                                                "jaccard_sequences": jaccard_sequences, 
                                                                "jaccard_name_len":jaccard_name_len,
                                                                "opa_names":opa_names,
                                                                "opb_names":opb_names,
                                                                "opa_lengths":opa_lengths,
                                                                "opb_lengths":opb_lengths,
                                                                "opa_sequences":opa_sequences,
                                                                "opb_sequences":opb_sequences,
                                                                "opa_name_len":opa_name_len,
                                                                "opb_name_len":opb_name_len
                                                                })
    combination_count+=1

print (f"Finished with {combination_count} combinations processed")
                                               
