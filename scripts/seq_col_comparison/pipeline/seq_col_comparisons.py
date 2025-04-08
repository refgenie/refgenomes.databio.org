import sys
import os
import urllib.request
import urllib.error
import pipestat
from refget import fasta_to_digest, fasta_to_seqcol_dict,fasta_to_seq_digests
import pypiper
import pprint

sample_name= sys.argv[1]  
digest = sys.argv[2]
file_path = sys.argv[3]
looper_output_dir = sys.argv[4]

print(f"Here is the file_path: {file_path}")

exists = os.path.exists(file_path)
print(f"File exists: {exists}")

# WE WANT sha512t24u due to robustness wrt collisions

# for x in fasta_to_seq_digests(file_path):
#     print(f"{x.id}\t{x.length}\t{x.sha512t24u}\t{x.md5}")

seq_dict = fasta_to_seqcol_dict(file_path)

pprint.pprint(seq_dict)
