wget -P /cellar/users/aklie/data/datasets/AI-ATAC/analysis/10Nov23/github https://www.dropbox.com/s/r8drj2wxc07bt4j/ImmGenATAC1219.peak_matched.txt
wget -P /cellar/users/aklie/data/datasets/AI-ATAC/analysis/10Nov23/github https://www.dropbox.com/s/7mmd4v760eux755/mouse_peak_heights.csv
cd /cellar/users/aklie/data/datasets/AI-ATAC/analysis/10Nov23/github
cut -f 2,3,4 ImmGenATAC1219.peak_matched.txt > ImmGenATAC1219.peak_matched.bed