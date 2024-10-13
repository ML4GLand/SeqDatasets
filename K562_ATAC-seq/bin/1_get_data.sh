# Path data (TODO: change to the path you want to store the data in)
path_data=/cellar/users/aklie/data/datasets/K562_ATAC-seq/data

# download reference data
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz -O $path_data/hg38.fa.gz
gunzip $path_data/hg38.fa.gz

# download reference chromosome sizes 
wget https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv -O $path_data/hg38.chrom.sizes

# download reference blacklist regions 
wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz -O $path_data/blacklist.bed.gz
gunzip $path_data/blacklist.bed.gz

# download K562 ENCODE ATAC-seq peaks
wget https://www.encodeproject.org/files/ENCFF333TAT/@@download/ENCFF333TAT.bed.gz -O $path_data/ENCSR868FGK_K562_ATAC-seq_peaks.bed.gz
gunzip $path_data/ENCSR868FGK_K562_ATAC-seq_peaks.bed.gz

# remove blacklisted regions from peaks
bedtools slop -i blacklist.bed -g hg38.chrom.sizes -b 1057 > temp.bed
bedtools intersect -v -a $path_data/ENCSR868FGK_K562_ATAC-seq_peaks.bed -b $path_data/temp.bed > $path_data/ENCSR868FGK_K562_ATAC-seq_peaks_no_blacklist.bed

# download K562 ENCODE ATAC-seq bam for reps 1, 2, 3
wget https://www.encodeproject.org/files/ENCFF077FBI/@@download/ENCFF077FBI.bam -O $path_data/ENCSR868FGK_K562_ATAC-seq_rep1.bam
wget https://www.encodeproject.org/files/ENCFF128WZG/@@download/ENCFF128WZG.bam -O $path_data/ENCSR868FGK_K562_ATAC-seq_rep2.bam
wget https://www.encodeproject.org/files/ENCFF534DCE/@@download/ENCFF534DCE.bam -O $path_data/ENCSR868FGK_K562_ATAC-seq_rep3.bam
