#!/bin/bash
#SBATCH --partition carter-compute
#SBATCH --output=/cellar/users/aklie/data/datasets/Horlacher_HepG2_CLIP/bin/data_processing/slurm_logs/%x.%A.%a.out
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=14-00:00:00
#SBATCH --array=1-10%10

#####
# USAGE: 
# sbatch --job-name=preprocess preprocess.sh
#####

# Date and Job ID
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate rbpnet
script_path=/cellar/users/aklie/data/datasets/Horlacher_HepG2_CLIP/bin/data_processing/preprocess_RBP.sh

# Path to given RBP
dir=/cellar/users/aklie/data/datasets/Horlacher_HepG2_CLIP/processed/2023_12_29/encode

# Get all subdirectories in this directory, then choose the one that matches the SLURM array number
sub_dirs=($(ls -d $dir/*/))
sub_dir=${sub_dirs[$SLURM_ARRAY_TASK_ID-1]}  # Assuming SLURM_ARRAY_TASK_ID is set

# Params
rbp=$(basename "$sub_dir")  # Get the basename from this dir
tile_size=100
min_pval=0.01
min_count=15
min_height=5
threads=$SLURM_CPUS_PER_TASK  # Get the number from SLURM

# Files
chrom_sizes=/cellar/users/aklie/data/ref/genomes/hg38/GRCh38_EBV.chrom.sizes
annotation=/cellar/users/aklie/data/ref/genomes/hg38/gencode.v40.basic.annotation.bed

# Check if required variables are set
if [ -z "$rbp" ] || [ -z "$threads" ]; then
  echo "Error: Missing RBP name or SLURM threads configuration."
  exit 1
fi

# Run the command
cmd="bash $script_path $sub_dir $rbp $tile_size $min_pval $min_count $min_height $threads $chrom_sizes $annotation"
echo -e "Running:\n $cmd"
eval $cmd

# End
echo "Done."
date
