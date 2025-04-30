import sys
import os
import urllib.request
import urllib.error
import pipestat
from refget import fasta_to_digest, fasta_to_seqcol_dict
import pypiper


sample_name= sys.argv[1]  
ftp_url = sys.argv[2]
species = sys.argv[3]
common_genome_name = sys.argv[4]
authority = sys.argv[5]
file_type = sys.argv[6]
pephub_path = sys.argv[7]
download_location = sys.argv[8]
looper_output_dir = sys.argv[9]


print(f"HERE IS THE FILE PATH:{ftp_url}")

# Construct download location based on pep

# Make digest here and now and download it

pypiper_logs = os.path.join(looper_output_dir, "pipeline_results",sample_name)

pm = pypiper.PipelineManager(
    name="FASTA_DOWLOADER",
    outfolder=pypiper_logs,
    pipestat_record_identifier=sample_name,
    recover=True,
)

pm.start_pipeline()

download_path = os.path.join(download_location, species, authority, common_genome_name, file_type)

try:
    os.makedirs(download_path, exist_ok=True)  
    print(f"Directory created successfully: {download_path}")
except OSError as e:
    print(f"Error creating directory {download_path}: {e}")

#print(f"HERE IS THE DOWNLOAD PATH:{download_path}")

try:
    filename = os.path.basename(ftp_url) 
    filepath = os.path.join(download_path, filename)

    if not os.path.exists(filepath):
        print(f"Downloading {ftp_url} to {filepath}...")
        urllib.request.urlretrieve(ftp_url, filepath)
        print(f"Downloaded {filename} successfully!")
    else:
        print(f"File exists at: {filepath}")

except urllib.error.URLError as e:
    print(f"Error downloading {ftp_url}: {e}")
except Exception as e:
    print(f"An unexpected error occurred while downloading {ftp_url}: {e}")

#report final digest and path to a pep on pephub

digest = fasta_to_digest(filepath,inherent_attrs=['names', 'sequences'])

# seq_col_dict = fasta_to_seqcol_dict(filepath) # 'sorted_name_length_pairs`` is a list of bytes data and is NOT json serializable, so it cannot be uploaded to pephub via pipestat
# print(type(seq_col_dict))

# print(f"Here is the digest: {digest}")
# print(f"Here is the seq_col_dict: {seq_col_dict}")

psm = pipestat.PipestatManager(pephub_path=pephub_path)

psm.report(record_identifier=sample_name, values={"top_level_digest":digest, "brickyard_location":filepath, "original_file_name":filename, "authority": authority})
pm.report_result("top_level_digest", digest)
pm.report_result("brickyard_location",filepath)
pm.stop_pipeline()
#psm.set_status(record_identifier=filename, status_identifier='completed')