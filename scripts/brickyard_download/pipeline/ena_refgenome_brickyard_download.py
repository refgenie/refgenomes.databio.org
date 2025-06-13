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


import requests
import gzip
import os
from urllib.parse import urlparse

def download_fasta_stream(url, output_filename, download_dir="downloads"):
    """
    Downloads a FASTA file from a given URL by streaming its content to a local file.
    The file will be saved in the specified download directory with the provided output_filename.
    It checks if the file already exists before downloading.

    Args:
        url (str): The URL of the FASTA file to download.
        output_filename (str): The desired name for the local FASTA file,
                               e.g., "sequence.fa" or "sequence.fa.gz".
        download_dir (str): The directory where the file should be saved.
                            Defaults to "downloads".
    """
    # Create the download directory if it doesn't exist
    os.makedirs(download_dir, exist_ok=True)

    # Construct the full path for the output file
    full_output_path = os.path.join(download_dir, output_filename)

    try:
        # Check if the file already exists
        if os.path.exists(full_output_path):
            print(f"File already exists: {full_output_path}. Skipping download for {url}")
            return full_output_path

        # Initiate the GET request with stream=True to handle large files efficiently
        print(f"Starting download from: {url}")
        with requests.get(url, stream=True) as response:
            response.raise_for_status()  # Raise an HTTPError for bad responses (4xx or 5xx)

            # Informative message if content is gzipped, though no decompression is done here
            if response.headers.get('Content-Encoding') == 'gzip':
                print("Server indicates gzipped content will be downloaded.")
            else:
                print("Server indicates uncompressed content will be downloaded.")

            # Open the local file in binary write mode using the determined full_output_path
            with open(full_output_path, 'wb') as f:
                # Iterate over the response content in chunks and write to file
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:  # Filter out keep-alive new chunks
                        f.write(chunk)
            print(f"Successfully downloaded to: {full_output_path}")
            return full_output_path

    except requests.exceptions.RequestException as e:
        print(f"Network or HTTP error occurred for {url}: {e}")
        return None
    except IOError as e:
        print(f"File system error occurred while writing or reading for {full_output_path}: {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred for {url}: {e}")
        return None

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

print(f"HERE IS THE DOWNLOAD PATH:{download_path}")

parsed_url = urlparse(ftp_url)
base_filename = os.path.basename(parsed_url.path)
base_filename = os.path.splitext(base_filename)[0]

# Determine the desired final filename based on whether 'gzip=true' is in the URL query parameters
if "gzip=true" in parsed_url.query:
    local_fasta_filename = f"{sample_name}_{base_filename}.fa.gz"
else:
    local_fasta_filename = f"{sample_name}_{base_filename}.fa"

# Call the download function with the URL, generated filename, and the specified download directory
filepath = download_fasta_stream(ftp_url, local_fasta_filename, download_dir=download_path)


#report final digest and path to a pep on pephub

if filepath:
    digest = fasta_to_digest(filepath, inherent_attrs=['names', 'sequences'])

    # seq_col_dict = fasta_to_seqcol_dict(filepath) # 'sorted_name_length_pairs`` is a list of bytes data and is NOT json serializable, so it cannot be uploaded to pephub via pipestat
    # print(type(seq_col_dict))

    # print(f"Here is the digest: {digest}")
    # print(f"Here is the seq_col_dict: {seq_col_dict}")

    psm = pipestat.PipestatManager(pephub_path=pephub_path)

    psm.report(record_identifier=sample_name, values={"top_level_digest":digest, "brickyard_location":filepath, "original_file_name":local_fasta_filename, "authority": authority})
    pm.report_result("top_level_digest", digest)
    pm.report_result("brickyard_location",filepath)
else:
    print("No local filepath returned. No digest calcualted nor reported.")
pm.stop_pipeline()
#psm.set_status(record_identifier=filename, status_identifier='completed')

# print(f"HERE IS THE FILE PATH:{ftp_url}")

# # Construct download location based on pep

# # Make digest here and now and download it



# download_path = os.path.join(download_location, species, authority, common_genome_name, file_type)

# try:
#     os.makedirs(download_path, exist_ok=True)  
#     print(f"Directory created successfully: {download_path}")
# except OSError as e:
#     print(f"Error creating directory {download_path}: {e}")

# #print(f"HERE IS THE DOWNLOAD PATH:{download_path}")

# try:
#     filename = os.path.basename(ftp_url) 
#     filepath = os.path.join(download_path, filename)

#     if not os.path.exists(filepath):
#         print(f"Downloading {ftp_url} to {filepath}...")
#         urllib.request.urlretrieve(ftp_url, filepath)
#         print(f"Downloaded {filename} successfully!")
#     else:
#         print(f"File exists at: {filepath}")

# except urllib.error.URLError as e:
#     print(f"Error downloading {ftp_url}: {e}")
# except Exception as e:
#     print(f"An unexpected error occurred while downloading {ftp_url}: {e}")

