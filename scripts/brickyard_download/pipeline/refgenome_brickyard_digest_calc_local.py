import sys
import os
import pipestat
from refget import fasta_to_digest
import pypiper


sample_name= sys.argv[1]  
local_fasta_filepath = sys.argv[2]

authority = sys.argv[3]

pephub_path = sys.argv[4]

looper_output_dir = sys.argv[5]


pypiper_logs = os.path.join(looper_output_dir, "pipeline_results",sample_name)

pm = pypiper.PipelineManager(
    name="brickyard_digest_calc",
    outfolder=pypiper_logs,
    pipestat_record_identifier=sample_name,
    recover=True,
)

pm.start_pipeline()


#report final digest and path to a pep on pephub

digest = fasta_to_digest(local_fasta_filepath,inherent_attrs=['names', 'sequences'])

psm = pipestat.PipestatManager(pephub_path=pephub_path)

filename = os.path.basename(local_fasta_filepath)

psm.report(record_identifier=sample_name, values={"top_level_digest":digest, "brickyard_location":local_fasta_filepath, "original_file_name":filename, "authority": authority})
pm.report_result("top_level_digest", digest)
pm.report_result("brickyard_location",local_fasta_filepath)
pm.stop_pipeline()
#psm.set_status(record_identifier=filename, status_identifier='completed')