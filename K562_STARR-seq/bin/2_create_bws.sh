# Path data (TODO: change to the path you specified in bin/1_get_data.sh)
path_data=/cellar/users/aklie/data/datasets/K562_STARR-seq/data

# use samtools to merge the bam files
samtools merge -f $path_data/merged_unsorted.bam $path_data/ENCSR926NDZ_K562_STARR-seq_rep1.bam $path_data/ENCSR926NDZ_K562_STARR-seq_rep2.bam $path_data/ENCSR926NDZ_K562_STARR-seq_rep3.bam
samtools sort -@4 $path_data/merged_unsorted.bam -o $path_data/merged.bam
samtools index $path_data/merged.bam

# get coverage as bedgraph
bedtools genomecov -bg -ibam $path_data/merged.bam -g $path_data/hg38.chrom.sizes | sort -k1,1 -k2,2n > $path_data/merged.bedgraph

# convert bedgraph to bigwig
bedGraphToBigWig $path_data/merged.bedgraph $path_data/hg38.chrom.sizes K562_STARR-seq_unstranded_counts.bw

# remove intermediate files
rm $path_data/merged_unsorted.bam $path_data/merged.bam $path_data/merged.bedgraph
