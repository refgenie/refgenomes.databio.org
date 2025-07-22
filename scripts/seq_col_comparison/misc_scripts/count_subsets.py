# COUNT SUBSETS VIA OPA/OPB as well as get percentage of stats that are 1.0 similarity

import sys
from pephubclient import PEPHubClient
import pipestat

RESULTS_PEP= "donaldcampbelljr/human_seq_col_results:default"
#RESULTS_PEP= "donaldcampbelljr/mouse_seq_col_results:default"
COMPARISON = 1.0 


phc = PEPHubClient()
pep = phc.load_project(RESULTS_PEP)

pep_df = pep["_sample_df"]
all_unique_samples = set(pep_df['sample_name_1']).union(set(pep_df['sample_name_2']))

stats = [("opa_name_len","opb_name_len", "name_len"), ("opa_lengths","opb_lengths", "lengths"), ("opa_sequences","opb_sequences", "sequences"), ("opa_names","opb_names", "names")]

# Iterate through the DataFrame rows
for stat1, stat2, stat3 in stats:
    samples_opa_name_len_in_col1 = []
    samples_opb_name_len_in_col2 = []
    for index, row in pep_df.iterrows():
        if row['sample_name_1'] == row['sample_name_2']:
            print("sample names the same, passing")
            pass
        else:
            if row[stat1] == 1.0:
                samples_opa_name_len_in_col1.append(row['sample_name_1'])
            if row[stat2] == 1.0:
                samples_opb_name_len_in_col2.append(row['sample_name_2'])

    # Convert lists to sets to get unique sample names
    unique_samples_col1 = set(samples_opa_name_len_in_col1)
    unique_samples_col2 = set(samples_opb_name_len_in_col2)

    all_col_sample = unique_samples_col1.union(unique_samples_col2)

    print(f"LEN unique_samples_col1: {len(unique_samples_col1)}")
    print(f"LEN unique_samples_col2: {len(unique_samples_col2)}")
    print(f"LEN union: {len(all_col_sample)}")
    print(f"LEN All Unique Samples: {len(all_unique_samples)}")
    print(f"Percentage for OPA/OPB for {stat3} = {(len(all_col_sample)/len(all_unique_samples))*100}")

    # for sample in all_col_sample:
    #     print(sample)


# ------------- Calculating Percentages based on Jaccard SImialrity for Each Stat
all_relevant_stats = [
                      'jaccard_names', 
                      'jaccard_lengths', 
                      'jaccard_sequences',
                      'jaccard_name_len', 
                      ]




for stat in all_relevant_stats:
    print(f"CALCULATING FOR {stat}---------")
    samples_jaccard_name_len = []
    count_all = 0
    count_target = 0
    for index, row in pep_df.iterrows():
        count_all+=1
        if row['sample_name_1'] == row['sample_name_2']:
            print("sample names the same, passing")
            pass
        else:
            if row[stat] >= COMPARISON: # CHECK TO ENSURE THIS IS THE PROPER COMPARISON IF CHANGING NUMBER
                count_target+=1
                samples_jaccard_name_len.append(row['sample_name_1'])
                samples_jaccard_name_len.append(row['sample_name_2'])

    unique_samples_name_len = set(samples_jaccard_name_len)

    # print(f"LEN unique_samples for stat: {STAT}: {len(unique_samples_name_len)}")
    # print(f"LEN All Unique Samples: {len(all_unique_samples)}")
    # print(count_all)
    # print(count_target)

    print(f"here is how many {stat} were equal to {COMPARISON} divided by all comparisons")
    print((count_target/count_all)*100)

    print(f"here is how many samples were equal to {COMPARISON} divided by all total number of samples for stat {stat}")
    print((len(unique_samples_name_len)/len(all_unique_samples))*100)



unique_digest1_values = pep_df['digest1'].unique()
unique_digest2_values = pep_df['digest2'].unique()

all_uniques = set(unique_digest1_values).union(set(unique_digest2_values))
#print(len(all_uniques))
print("Here is unique digests over total samples:")
print((len(all_uniques)/len(all_unique_samples))*100)