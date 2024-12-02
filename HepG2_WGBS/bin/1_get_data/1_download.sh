# Description: Download the data for the HepG2 WGBS dataset from ENCODE

# Path data (TODO: change to the path you want to store the data in)
path_data=/cellar/users/aklie/data/datasets/SeqDatasets/HepG2_WGBS/processed

# CpG sites coverage as bigwig
wget https://www.encodeproject.org/files/ENCFF782HWK/@@download/ENCFF782HWK.bigWig -O $path_data/ENCFF782HWK_HepG2_WGBS_CpG_sites_coverage_rep1.bigWig
wget https://www.encodeproject.org/files/ENCFF306UZK/@@download/ENCFF306UZK.bigWig -O $path_data/ENCFF306UZK_HepG2_WGBS_CpG_sites_coverage_rep2.bigWig

# CpG methylation state as bedmethyl (see README for details)
wget https://www.encodeproject.org/files/ENCFF820ATI/@@download/ENCFF820ATI.bed.gz -O $path_data/ENCFF820ATI_HepG2_WGBS_CpG_methylation_state_rep1.bed.gz
wget https://www.encodeproject.org/files/ENCFF690FNR/@@download/ENCFF690FNR.bed.gz -O $path_data/ENCFF690FNR_HepG2_WGBS_CpG_methylation_state_rep2.bed.gz

# Total coverage
wget https://www.encodeproject.org/files/ENCFF400QTE/@@download/ENCFF400QTE.bigWig -O $path_data/ENCFF400QTE_HepG2_WGBS_total_coverage_rep1.bigWig
wget https://www.encodeproject.org/files/ENCFF583VWF/@@download/ENCFF583VWF.bigWig -O $path_data/ENCFF583VWF_HepG2_WGBS_total_coverage_rep2.bigWig
