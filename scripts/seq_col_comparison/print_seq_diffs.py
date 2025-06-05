# just print some differences in seqs for illustration purposes
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


#----
def get_dict_seq_col_from_json(digest):

    # OPENS A LOCAL JSON ONLY
    json_fp_1 = os.path.join(LOCAL_JSON_DIRECTORY, digest+".json")
    with open(json_fp_1, "r") as f:
        reloaded_dict1 = json.load(fp=f)
    
    return reloaded_dict1



#----

#looper_config = sys.argv[1]  
#LOCAL_JSON_DIRECTORY  = sys.argv[2] # for lcally stored jsons to pull sequence col data

looper_config = "donaldcampbelljr/human_seqcol_digests:default"
LOCAL_JSON_DIRECTORY = "/home/drc/Downloads/jsons_from_rivanna/json/"

phc = PEPHubClient()
pep = phc.load_project(looper_config)
print(pep["_sample_df"])

pep_df = pep["_sample_df"]



target_samples = [
"hg19-p13-plusMT-masked-ucsc",
"hg19-p13-no-alt-analysis-ucsc",
"hg19-p13-full-analysis-ucsc",
"hg19-p13-plusMT-ucsc",
]

# # # # Pre-filter the DataFrame
pep_df = pep_df[
    ((pep_df['sample_name'].isin(target_samples)))
]

#print(pep_df)

for index, row in pep_df.iterrows():
    reloaded_dict = get_dict_seq_col_from_json(row['top_level_digest'])
    print(reloaded_dict.keys())