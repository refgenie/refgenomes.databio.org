# A script to troubleshoot direct comparisons

from refget import compare_seqcols
import os
import json
LOCAL_JSON_DIRECTORY = "/home/drc/Downloads/jsons_from_rivanna/json/"
#sample1 = "t-g5olKSC0DsZHICbi3cYR3jF-7i4xaK" # p2
#sample1 = "P33s5fSkktH60MccIpfLNHAA79fJyxlt" # p1
# sample1 = "9VlpW2Uty7rQNSBbYy-nYCNOHqyYoBTt" # p6
# sample1 = "e38absAQBVXM4If9Mp9UhrPrqNGOklZz" # p7
# sample2 = "d2GjQek3_I1_zkq4gLHTUI50BptQsoqZ" # p13

# primary_sample = "d2GjQek3_I1_zkq4gLHTUI50BptQsoqZ" # p13

# samples_to_compare = ["P33s5fSkktH60MccIpfLNHAA79fJyxlt","t-g5olKSC0DsZHICbi3cYR3jF-7i4xaK", "9VlpW2Uty7rQNSBbYy-nYCNOHqyYoBTt", "e38absAQBVXM4If9Mp9UhrPrqNGOklZz", "d2GjQek3_I1_zkq4gLHTUI50BptQsoqZ"]


samples_to_compare = ["XemD97fxYMS4q-FBm_n5CHQgmzh1_67a"] # p0
primary_sample = "u1HyLgIlq8M_XvEwy0oGqAvKGHJMGtxH" # p14

def get_sequence_col_dict(digest):

    # OPENS A LOCAL JSON ONLY
    json_fp_1 = os.path.join(LOCAL_JSON_DIRECTORY, digest+".json")
    with open(json_fp_1, "r") as f:
        reloaded_dict = json.load(fp=f)
    
    return reloaded_dict

def overlap_proportion(A_B_intersection, a_or_b):
    # calculate the intersection OVER one of the two sets used to calculate the intersection
    proportion = abs(A_B_intersection)/abs(a_or_b)
    return proportion


# reloaded_dict1 = get_sequence_col_dict(sample1)
# reloaded_dict2 = get_sequence_col_dict(sample2)
# # print(pprint(compare_seqcols(reloaded_dict1,reloaded_dict2),indent=4))
# comparison = compare_seqcols(reloaded_dict1,reloaded_dict2)

# # create set of name_length_pairs
#     #print(type(reloaded_dict1['name_length_pairs'][0]))
# set_of_name_len_pairs_1 = {tuple(d.values()) for d in reloaded_dict1['name_length_pairs']}
# set_of_name_len_pairs_2 = {tuple(d.values()) for d in reloaded_dict2['name_length_pairs']}

# name_len_pairs_intersection = set_of_name_len_pairs_1.intersection(set_of_name_len_pairs_2)
# name_len_pairs_union = set_of_name_len_pairs_1.union(set_of_name_len_pairs_2)
# opa_name_len = overlap_proportion(len(name_len_pairs_intersection),len(set_of_name_len_pairs_1))
# opb_name_len = overlap_proportion(len(name_len_pairs_intersection),len(set_of_name_len_pairs_2))

# print(f"OPA: {opa_name_len}   OPB: {opb_name_len}")


for sample in samples_to_compare:
    reloaded_dict1 = get_sequence_col_dict(sample)
    reloaded_dict2 = get_sequence_col_dict(primary_sample)
    # print(pprint(compare_seqcols(reloaded_dict1,reloaded_dict2),indent=4))
    comparison = compare_seqcols(reloaded_dict1,reloaded_dict2)

    # create set of name_length_pairs
        #print(type(reloaded_dict1['name_length_pairs'][0]))
    set_of_name_len_pairs_1 = {tuple(d.values()) for d in reloaded_dict1['name_length_pairs']}
    set_of_name_len_pairs_2 = {tuple(d.values()) for d in reloaded_dict2['name_length_pairs']}

    name_len_pairs_intersection = set_of_name_len_pairs_1.intersection(set_of_name_len_pairs_2)
    name_len_pairs_union = set_of_name_len_pairs_1.union(set_of_name_len_pairs_2)
    opa_name_len = overlap_proportion(len(name_len_pairs_intersection),len(set_of_name_len_pairs_1))
    opb_name_len = overlap_proportion(len(name_len_pairs_intersection),len(set_of_name_len_pairs_2))

    print(f"OPA: {opa_name_len}   OPB: {opb_name_len}")


print("Done")