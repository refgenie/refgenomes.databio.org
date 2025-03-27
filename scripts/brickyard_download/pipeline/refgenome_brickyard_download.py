import sys
import os
import urllib.request
import urllib.error
import pipestat
from refget import fasta_to_digest


sample_name= sys.argv[1]  
ftp_url = sys.argv[2]
species = sys.argv[3]
common_genome_name = sys.argv[4]
authority = sys.argv[5]
file_type = sys.argv[6]
results_dir = sys.argv[7]
# results_file = sys.argv[3]
# schema_path = sys.argv[4]

print(f"HERE IS THE FILE PATH:{ftp_url}")

# Construct download location based on pep

# Make digest here and now and download it

download_path = os.path.join(results_dir,authority, species, common_genome_name, file_type)

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

print(f"Here is the digest: {digest}")

psm = pipestat.PipestatManager(pephub_path="donaldcampbelljr/human_seqcol_digests:default")

psm.report(record_identifier=filename, values={"top_level_digest":digest, "brickyard_location":filepath})