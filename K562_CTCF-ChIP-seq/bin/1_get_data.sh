# Path reference data (TODO: change to your path)
path_ref=/cellar/users/aklie/data/datasets/K562_CTCF-ChIP-seq/data

# download reference data
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz -O $path_ref/hg38.fa.gz
gunzip $path_ref/hg38.fa.gz

# download reference chromosome sizes 
wget https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv -O $path_ref/hg38.chrom.sizes

# download reference blacklist regions 
wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz -O $path_ref/blacklist.bed.gz
gunzip $path_ref/blacklist.bed.gz

# download K562 ENCODE CTCF ChIP-seq peaks from bpnet-lite examples
wget https://raw.githubusercontent.com/jmschrei/bpnet-lite/master/examples/ENCSR000AKO.bed -O $path_ref/K562_CTCF-ChIP-seq_peaks.bed.gz
gunzip $path_data/K562_CTCF-ChIP-seq_peaks.bed.gz

# download forward and reverse strand counts from bpnet-lite examples
wget https://raw.githubusercontent.com/jmschrei/bpnet-lite/master/examples/ENCSR000AKO_minus.bigWig -O $path_ref/ENCSR000AKO_minus_counts.bw
wget https://raw.githubusercontent.com/jmschrei/bpnet-lite/master/examples/ENCSR000AKO_plus.bigWig -O $path_ref/ENCSR000AKO_plus_counts.bw
