# takes fasta files, computes json of sequence collection to be used in downstream analysis
# while we are here generate the chrom.sizes files too?

# this will take a pep from pephub and MODIFY it with the new columns.


import pipestat
import pypiper
import os
import sys
import json
from refget import fasta_to_seqcol_dict, SequenceCollection




pep_config= sys.argv[1]   # {looper.pep_config} 
looper_output_dir = sys.argv[2] # {looper.output_dir} 
top_level_digest = sys.argv[3] # {sample.top_level_digest} 
brickyard_location = sys.argv[4] # {sample.brickyard_location} 
json_directory = sys.argv[5] # {sample.json_creation_location} 
chrom_sizes_location = sys.argv[6] # {sample.chrom_sizes_location}
sample_name = sys.argv[7]


pypiper_logs = os.path.join(looper_output_dir, "pipeline_results",top_level_digest)

pm = pypiper.PipelineManager(
    name="SEQ_COL_CREATOR",
    outfolder=pypiper_logs,
    pipestat_record_identifier=top_level_digest,
    recover=True,
)

pm.start_pipeline()

json_creation_path = os.path.join(json_directory,top_level_digest+".json")
chrom_sizes_creation_path = os.path.join(chrom_sizes_location,top_level_digest+".chrom.sizes")

# create seq_col
if not os.path.exists(json_creation_path):
    try:

        res = SequenceCollection.from_dict(
                fasta_to_seqcol_dict(brickyard_location),
                inherent_attrs=["names", "sequences"],
            )

        # write jsons
        level_2 = res.level2()
        with open(json_creation_path, "w") as f:
            f.write(json.dumps(level_2 , indent=2))
    except Exception as e:
        print(f"An unexpected error occurred while creatiing json at path {json_creation_path}: {e}")
else:
    print(f"File exists at: {json_creation_path}")

if not os.path.exists(chrom_sizes_creation_path):
#write chrom.sizes
    try:
        name_length_pairs = level_2['name_length_pairs']
        with open(chrom_sizes_creation_path, 'w') as outfile:
            for item in name_length_pairs:
                line = f"{item['name']}\t{item['length']}\n"
                outfile.write(line)
    except Exception as e:
        print(f"An unexpected error occurred while creatiing chrom.sizes file at path {chrom_sizes_creation_path}: {e}")

else:
    print(f"File exists at: {chrom_sizes_creation_path}")

# report final locations BACK to pephub
psm = pipestat.PipestatManager(pephub_path=pep_config)

psm.report(record_identifier=sample_name, values={"brickyard_json_path": json_creation_path,"brickyard_chrom_sizes_path": chrom_sizes_creation_path})

pm.stop_pipeline()