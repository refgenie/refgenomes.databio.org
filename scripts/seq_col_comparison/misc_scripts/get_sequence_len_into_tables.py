import os
import pipestat
import pephubclient
import pandas as pd


def get_sequence_length_local(digest):

    # OPENS A LOCAL JSON ONLY
    json_fp_1 = os.path.join(LOCAL_JSON_DIRECTORY, digest+".json")
    with open(json_fp_1, "r") as f:
        reloaded_dict1 = json.load(fp=f)
    
    return len(reloaded_dict1['sorted_sequences'])


LOCAL_JSON_DIRECTORY = "/home/drc/Downloads/jsons_from_rivanna/json/"
PEPHUB_PATH_AUTHORITY = "donaldcampbelljr/human_seqcol_digests:default"

results_dir = "/home/drc/Downloads/add_lengths_tables/human/"  # Replace with your directory path

psm = pipestat.PipestatManager(pephub_path=PEPHUB_PATH_AUTHORITY)

# Get all entries (files and directories)
all_entries = os.listdir(results_dir)
print(f"All entries: {all_entries}")

# Filter to get only files
files_only = [entry for entry in all_entries if os.path.isfile(os.path.join(results_dir, entry))]
print(f"Files only: {files_only}")

for file in files_only:
    #results = psm.select_records()
    original_file_name = os.path.basename(file)
    csv_filepath = os.path.join(results_dir, file)
    data_table_df = pd.read_csv(csv_filepath, index_col='sample_name')

    # Iterate over each row, get sample name

    for row in data_table_df.iterrows():
        sample_name =

        # Get digest
        top_level_digest = psm.retrieve_one(record_identifier=sample_name, result_identifier='top_level_digest')

        # Get length by reading the associated json file base don the digest and get length

        seq_len = get_sequence_length_local(digest=top_level_digest)

        # put seq_len and top_level_digest into the datatframe held in memory

        data_table_df['Top Level Digest']=top_level_digest
        data_table_df['Number of Sequences'] = seq_len

    new_file_name = original_file_name+"_modified"+".csv"
    new_file_path = os.path.join(results_dir,new_file_name)
    data_table_df.to_csv(new_file_path, index=True, header=True) 

