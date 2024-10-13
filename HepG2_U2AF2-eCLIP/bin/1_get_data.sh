# Path reference data (TODO: change to your path)
path_ref=/cellar/users/aklie/data/datasets/HepG2_U2AF2-eCLIP/data

# download reference data
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz -O $path_ref/hg38.fa.gz
gunzip $path_ref/hg38.fa.gz

# download reference chromosome sizes 
wget https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv -O $path_ref/hg38.chrom.sizes

# download reference blacklist regions 
wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz -O $path_ref/blacklist.bed.gz
gunzip $path_ref/blacklist.bed.gz

# download HepG2 U2AF2 eCLIP control data for positive and negative strand
wget https://github.com/mhorlacher/rbpnet/raw/main/examples/data/signal/U2AF2_HepG2/control/counts.neg.bw -O $path_ref/U2AF2_HepG2_control_neg.bw
wget https://github.com/mhorlacher/rbpnet/raw/main/examples/data/signal/U2AF2_HepG2/control/counts.pos.bw -O $path_ref/U2AF2_HepG2_control_pos.bw

# download HepG2 U2AF2 eCLIP signal data for positive and negative strand
wget https://github.com/mhorlacher/rbpnet/blob/main/examples/data/signal/U2AF2_HepG2/eCLIP/counts.neg.bw -O $path_ref/U2AF2_HepG2_eCLIP_neg.bw
wget https://github.com/mhorlacher/rbpnet/blob/main/examples/data/signal/U2AF2_HepG2/eCLIP/counts.pos.bw -O $path_ref/U2AF2_HepG2_eCLIP_pos.bw