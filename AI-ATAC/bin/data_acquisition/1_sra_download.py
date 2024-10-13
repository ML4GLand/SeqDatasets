import os
import sys
import glob
import subprocess
from tqdm.autonotebook import tqdm
from pysradb.sraweb import SRAweb

# Choose the current dataset we are working with
dataset_name = "AI-ATAC"
srp_id = "SRP110978"

# Set-up directories
base_dir = "/cellar/users/aklie/data/datasets"
cwd = os.path.join(base_dir, dataset_name, "bin", "data_acquisition")
fastq_dir = os.path.join(base_dir,  dataset_name, "fastq", "10Nov23")
metadata_dir = os.path.join(base_dir, dataset_name, "metadata", "10Nov23")

# If fastq_dir does not exist, create it
if not os.path.exists(fastq_dir):
    os.makedirs(fastq_dir)
    
# Connect to SRA
db = SRAweb()

# Grab the metadata for the SRP
metadata = db.sra_metadata(srp_id, detailed=True)

# Save the metadata and the list of srr ids
metadata.to_csv(os.path.join(metadata_dir, f"{srp_id}_metadata.tsv"), index=False, sep="\t")
metadata["run_accession"].to_csv(os.path.join(metadata_dir, f"{srp_id}_srr_ids.txt"), index=False, header=False)

# Download sra files
db.download(df=metadata, out_dir=fastq_dir, skip_confirmation=True)

# Parameters for download
tmp_dir = "/cellar/users/aklie/tmp/fastq-dump"
gzip = True
split_files = True
threads = 64

# Loop through and fastqdump out each SRA download file within the subdirectories of the fastq_dir
for sra_file in glob.glob(os.path.join(fastq_dir, srp_id, "*", "*.sra")):
    sra_dir = os.path.dirname(sra_file)
    if gzip:
        cmd = f"parallel-fastq-dump --threads {threads} --outdir {sra_dir} --split-files --tmpdir {tmp_dir} --gzip -s {sra_file}"
    else:
        cmd = f"parallel-fastq-dump --threads {threads} --outdir {sra_dir} --split-files --tmpdir {tmp_dir} -s {sra_file}"
    print(cmd)
    
    # Check to see if the files have already been downloaded
    if len(glob.glob(os.path.join(sra_dir, "*.fastq*"))) > 0:
        print(f"Files already downloaded for {sra_dir}")
    else:
        subprocess.run(cmd, shell=True)
