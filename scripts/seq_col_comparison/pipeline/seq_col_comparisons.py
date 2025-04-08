import sys
import os
import urllib.request
import urllib.error
import pipestat
from refget import fasta_to_digest, fasta_to_seqcol_dict
import pypiper

sample_name= sys.argv[1]  
digest = sys.argv[2]
file_path = sys.argv[3]
looper_output_dir = sys.argv[4]

print(f"Here is the file_path: {file_path}")

exists = os.path.exists(file_path)
print(f"File exists: {exists}")

