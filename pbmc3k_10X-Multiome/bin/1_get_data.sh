# Path data (TODO: change to the path you want to store the data in)
path_data=/cellar/users/aklie/data/datasets/pbmc3k_10X-Multiome/data

# download reference data
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz -O $path_data/hg38.fa.gz
gunzip $path_data/hg38.fa.gz

# download reference chromosome sizes 
wget https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv -O $path_data/hg38.chrom.sizes

# download reference blacklist regions 
wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz -O $path_data/blacklist.bed.gz
gunzip $path_data/blacklist.bed.gz

# download count matrices
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_raw_feature_bc_matrix.h5 -P $path_data

# download fragments and peaks from fragments
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_atac_peaks.bed -P $path_data
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz -P $path_data
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz.tbi -P $path_data

# download barcode metadata
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_per_barcode_metrics.csv -P $path_data
