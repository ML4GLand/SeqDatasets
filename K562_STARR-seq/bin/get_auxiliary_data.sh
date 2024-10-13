# Path data (TODO: change to the path you want to store the data in)
path_data=/cellar/users/aklie/data/datasets/K562_STARR-seq/data

# download K562 ENCODE STARR-seq bigWig for reps 1, 2, 3
wget https://www.encodeproject.org/files/ENCFF208RCC/@@download/ENCFF208RCC.bigWig -O $path_data/ENCSR926NDZ_K562_STARR-seq_rep1.bw
wget https://www.encodeproject.org/files/ENCFF381OCU/@@download/ENCFF381OCU.bigWig -O $path_data/ENCSR926NDZ_K562_STARR-seq_rep2.bw
wget https://www.encodeproject.org/files/ENCFF844UAO/@@download/ENCFF844UAO.bigWig -O $path_data/ENCSR926NDZ_K562_STARR-seq_rep3.bw
