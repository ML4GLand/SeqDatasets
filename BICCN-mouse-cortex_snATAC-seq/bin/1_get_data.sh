# URLs
base_url=https://resources.aertslab.org/CREsted/data/mouse_biccn/beds.tar.gz
data_dir=/cellar/users/aklie/projects/ML4GLand/SeqDatasets/BICCN-mouse-cortex_snATAC-seq/data

# Download the data
wget -P $data_dir $base_url
mkdir -p $data_dir/beds
tar -xvzf $data_dir/beds.tar.gz -C $data_dir/beds
