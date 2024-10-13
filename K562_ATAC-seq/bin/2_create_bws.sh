# Path data (TODO: change to the path you want to store the data in)
path_data=/cellar/users/aklie/data/datasets/K562_ATAC-seq/data

# Merge and sort bam files
samtools merge -f $path_data/merged_unsorted.bam $path_data/ENCSR868FGK_K562_ATAC-seq_rep1.bam $path_data/ENCSR868FGK_K562_ATAC-seq_rep2.bam $path_data/ENCSR868FGK_K562_ATAC-seq_rep3.bam
samtools sort -@4 $path_data/merged_unsorted.bam -o $path_data/merged.bam
samtools index $path_data/merged.bam

# use chrombpnet to get unstranded counts bigwig with correct shift (+4/-4)
script=/cellar/users/aklie/opt/chrombpnet/chrombpnet/helpers/preprocessing/reads_to_bigwig.py
cmd="python $script \
--genome $path_data/hg38.fa \
--input-bam-file $path_data/merged.bam \
--chrom-sizes $path_data/hg38.chrom.sizes \
--output-prefix $path_data/K562_ATAC-seq_merged \
--data-type ATAC"
echo $cmd
eval $cmd

# clean up intermediate files and rename K562_ATAC-seq_merged_unstranded.bw to K562_ATAC-seq_unstranded_counts.bw
mv $path_data/K562_ATAC-seq_merged_unstranded.bw $path_data/K562_ATAC-seq_unstranded_counts.bw 
rm $path_data/merged_unsorted.bam $path_data/merged.bam
