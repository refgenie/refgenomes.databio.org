import sys
import os

sample_name= sys.argv[1]  
ftp_file_path = sys.argv[2]
species = sys.argv[3]
common_genome_name = sys.argv[4]
authority = sys.argv[5]
file_type = sys.argv[6]
results_dir = sys.argv[7]
# results_file = sys.argv[3]
# schema_path = sys.argv[4]

print(f"HERE IS THE FILE PATH:{ftp_file_path}")

# Construct download location based on pep

# Make digest here and now and download it

download_path = os.path.join(results_dir, authority, species, common_genome_name, file_type)

print(f"HERE IS THE DOWNLOAD PATH:{download_path}")

#report final digest and path to a pep on pephub