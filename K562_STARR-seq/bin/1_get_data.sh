# Path data (TODO: change to the path you want to store the data in)
path_data=/cellar/users/aklie/data/datasets/K562_STARR-seq/data

# download reference data
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz -O $path_data/hg38.fa.gz
gunzip $path_data/hg38.fa.gz

# download reference chromosome sizes 
wget https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv -O $path_data/hg38.chrom.sizes

# download reference blacklist regions 
wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz -O $path_data/blacklist.bed.gz
gunzip $path_data/blacklist.bed.gz

# download K562 ENCODE STARR-seq bam for reps 1, 2, 3
wget https://www.encodeproject.org/files/ENCFF678CBN/@@download/ENCFF678CBN.bam -O $path_data/ENCSR926NDZ_K562_STARR-seq_rep1.bam
wget https://www.encodeproject.org/files/ENCFF390VZX/@@download/ENCFF390VZX.bam -O $path_data/ENCSR926NDZ_K562_STARR-seq_rep2.bam
wget https://www.encodeproject.org/files/ENCFF405KDE/@@download/ENCFF405KDE.bam -O $path_data/ENCSR926NDZ_K562_STARR-seq_rep3.bam

# download K562 ENCODE STARR-seq peaks
wget https://www.encodeproject.org/files/ENCFF064NPV/@@download/ENCFF064NPV.bed.gz -O $path_data/ENCSR926NDZ_K562_STARR-seq_peaks.bed.gz
gunzip $path_data/ENCSR926NDZ_K562_STARR-seq_peaks.bed.gz
