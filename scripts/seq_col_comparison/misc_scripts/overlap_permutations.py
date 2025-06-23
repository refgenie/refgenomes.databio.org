# import sys
import itertools
import os
import random

# import urllib.request
# import urllib.error
import pipestat
from pephubclient import PEPHubClient
from refget import fasta_to_digest, fasta_to_seqcol_dict, compare_seqcols, SequenceCollection,ga4gh_digest
import json

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


json_files = []
if os.path.isdir("/home/drc/Downloads/jsons_from_rivanna/json/"):
    for filename in os.listdir("/home/drc/Downloads/jsons_from_rivanna/json/"):
        if filename.endswith(".json"):
            full_path = os.path.join("/home/drc/Downloads/jsons_from_rivanna/json/", filename)
            json_files.append(full_path)

all_sequences_union = set()

dict_sequences = {}


bucket_100 = {}
bucket_100_500= {}
bucket_500_1000 = {}
bucket_1000_2000 = {}
bucket_2000_5000 ={}

#all_buckets = [bucket_100,bucket_100_500,bucket_500_1000,bucket_1000_2000,bucket_2000_5000]
all_buckets = [bucket_100]

for fp in json_files:
    with open(fp, "r") as f:
        reloaded_dict1 = json.load(fp=f)
        sequences = reloaded_dict1['sorted_sequences']
        if len(sequences) <= 100:
            bucket_100.update({fp:sequences})
        elif len(sequences) <= 500:
            bucket_100_500.update({fp: sequences})
        elif len(sequences) <= 1000:
            bucket_500_1000.update({fp: sequences})
        elif len(sequences) <= 2000:
            bucket_1000_2000.update({fp: sequences})
        elif len(sequences) <= 5000:
            bucket_2000_5000.update({fp: sequences})

print("done")
#
# for bucket in all_buckets:
#     print("NEW BUCKET")
#     temp_set = None
#     for key, value in bucket.items():
#         print(f"Iterating: {key}")
#         if temp_set is None:
#            temp_set=set(value)
#         else:
#             temp_set = temp_set.intersection(set(value))
#         print(f"new length of set: {len(temp_set)}")


list_keys = list(bucket_100.keys())
num_studies = len(list_keys)
num_permutations_to_try = 1000  # Adjust as needed

sequence_counts = {}

for i in range(num_permutations_to_try):
    random.shuffle(list_keys)  # Create a random permutation of the study keys
    current_intersection = None
    studies_in_permutation = []

    print(f"\n--- Permutation {i + 1} ---")

    for study_key in list_keys:
        if study_key not in studies_in_permutation:
            studies_in_permutation.append(study_key)
            sequences = bucket_100[study_key]
            print(f"Intersecting with study: {study_key}")

            if current_intersection is None:
                current_intersection = set(sequences)
            else:
                previous_intersection = current_intersection
                current_intersection = current_intersection.intersection(set(sequences))
                print(f"New intersection length: {len(current_intersection)}")
                if len(current_intersection) == 0:
                    for seq in previous_intersection:
                        if seq not in sequence_counts:
                            sequence_counts.update({seq:1})
                        else:
                            sequence_counts.update({seq: sequence_counts[seq]+1})

                    print(f"Intersection became empty after including study: {study_key}")

                    break  # Stop if the intersection becomes empty

    if current_intersection is not None:
        print(f"Final intersection for this permutation (studies: {studies_in_permutation}): {current_intersection}")
    else:
        print(f"No intersection calculated for this permutation.")

print("FINAL SEQUENCE COUNTS")

print(sequence_counts)

# --- Plotting ---
# Convert sequence_counts to a Pandas DataFrame for easier handling
sequence_df = pd.DataFrame(list(sequence_counts.items()), columns=['Sequence', 'Count'])

# Sort by count in descending order
sequence_df_sorted = sequence_df.sort_values(by='Count', ascending=False)
sequence_df_sorted = sequence_df_sorted.head(200)
# Plotting the sequence counts
plt.figure(figsize=(12, 6))  # Adjust figure size as needed
plt.bar(sequence_df_sorted['Sequence'], sequence_df_sorted['Count'])
plt.xlabel("Sequences")
plt.ylabel("Number of Times Sequence is Shared")
plt.title(f"Frequency of Sequences Shared across {num_permutations_to_try} Permutations")
plt.xticks(rotation=90, fontsize=8)  # Rotate x-axis labels for readability
plt.tight_layout()
plt.show()
