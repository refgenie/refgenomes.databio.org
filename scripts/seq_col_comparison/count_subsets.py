import sys
from pephubclient import PEPHubClient
import pipestat



#results_pep = "donaldcampbelljr/test_seq_col_results:default"


#results_pep = sys.argv[1] # input pep for graphing

results_pep = "donaldcampbelljr/human_seq_col_results:default"



phc = PEPHubClient()
pep = phc.load_project(results_pep)
#print(pep["_sample_df"])

pep_df = pep["_sample_df"]

# Initialize lists to store sample names based on 'jaccard_name_len' condition
samples_opa_name_len_in_col1 = []
samples_opb_name_len_in_col2 = []

# opa_name_len
# opb_name_len

# Iterate through the DataFrame rows
for index, row in pep_df.iterrows():
    if row['sample_name_1'] == row['sample_name_2']:
        print("sample names the same, passing")
        pass
    else:
        if row['opa_name_len'] == 1.0:
            samples_opa_name_len_in_col1.append(row['sample_name_1'])
        if row['opb_name_len'] == 1.0:
            samples_opb_name_len_in_col2.append(row['sample_name_2'])

# Convert lists to sets to get unique sample names
unique_samples_col1 = set(samples_opa_name_len_in_col1)
unique_samples_col2 = set(samples_opb_name_len_in_col2)

all_col_sample = unique_samples_col1.union(unique_samples_col2)

# Get all unique sample names from 'sample_name_1' and 'sample_name_2' columns
all_unique_samples = set(pep_df['sample_name_1']).union(set(pep_df['sample_name_2']))

#print(f"Samples where 'opa_name_len' is 1 (from sample_name_1 column): {unique_samples_col1}")
#print(f"Samples where 'opb_name_len' is 1 (from sample_name_2 column): {unique_samples_col2}")
#print(f"All unique samples in the DataFrame: {all_unique_samples}")


print(f"LEN unique_samples_col1: {len(unique_samples_col1)}")
print(f"LEN unique_samples_col2: {len(unique_samples_col2)}")
print(f"LEN union: {len(all_col_sample)}")
print(f"LEN All Unique Samples: {len(all_unique_samples)}")

for sample in all_col_sample:
    print(sample)

# You can then perform comparisons between these sets
# For example, to find samples present in all_unique_samples but not in unique_samples_col1:
# diff_col1 = all_unique_samples - unique_samples_col1
# print(f"Samples in all_unique_samples but not in unique_samples_col1: {diff_col1}")