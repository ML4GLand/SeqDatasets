#!/bin/bash

# USAGE
# preprocess_RBP.sh /cellar/users/aklie/data/datasets/Horlacher_HepG2_CLIP/processed/2023_12_29/encode/U2AF2/chr19 U2AF2 100 0.01 15 5 32

# Directory containing the files
dir="$1"
rbp="$2"  # rbp=$(basename "$(dirname "$dir")")
echo "Processing: $rbp from $dir"

# Tile size, min pval, min count, min height, and threads for enriched windows
tile_size="$3"
min_pval="$4"
min_count="$5"
min_height="$6"
threads="${7}"
chrom_sizes="$8"
annotation="$9"

# Script path
script="/cellar/users/aklie/opt/rbpnet/scripts/enriched-windows.py"


# Change to the directory
cd "$dir" || { echo "Directory not found"; exit 1; }

# Arrays to hold sorted BAM files for later merging
declare -a control_sorted_bams
declare -a eclip_sorted_bams

# Function to merge, sort, and index
merge_sort_index() {
    local merged_bam="$1_merged.bam"
    local sorted_bam="$1_merged.sorted.bam"

    # Merge
    echo "Merging ${1} replicates..."
    samtools merge "$merged_bam" "${@:2}"

    # Sort merged BAM
    echo "Sorting merged ${1} BAM..."
    samtools sort "$merged_bam" -o "$sorted_bam"

    # Index sorted merged BAM
    echo "Indexing sorted merged ${1} BAM..."
    samtools index "$sorted_bam"

    # Clean up unsorted merged BAM
    rm -f "$merged_bam"
}

# Function to extract 5' read-start coverage for plus and minus strands
extract_coverage() {
    local bam_file="$1"
    local prefix="${bam_file%.bam}"
    
    # Extract 5' read-start coverage for plus strand
    bedtools genomecov -ibam "$bam_file" -bg -5 -strand + | sort -k1,1 -k2,2n > "${prefix}.plus.bedGraph"
    
    # Extract 5' read-start coverage for minus strand
    bedtools genomecov -ibam "$bam_file" -bg -5 -strand - | sort -k1,1 -k2,2n > "${prefix}.minus.bedGraph"
    
    # Remove intermediate bam
    rm -f "$bam_file"
    rm -f "$bam_file.bai"
}

# Function to convert bedGraph to bigWig
convert_to_bigwig() {
    local bedgraph_file="$1"
    local bigwig_file="${bedgraph_file%.bedGraph}.bw"
    
    # Rename 'plus' to 'pos' and 'minus' to 'neg' in the output filename
    bigwig_file=$(echo "$bigwig_file" | sed 's/plus/pos/' | sed 's/minus/neg/')
    echo "Converting $bedgraph_file to $bigwig_file"
    bedGraphToBigWig "$bedgraph_file" "$chrom_sizes" "$bigwig_file"
    
    # Clean up bedGraph file
    rm -f "$bedgraph_file"
}

# Function to find enriched windows
find_enriched_windows() {
    local pos_bw="$1"
    local neg_bw="$2"
    local output_file="${pos_bw%.pos.bw}.crosslink.bed"
    
    # Construct and run the command
    cmd="python $script $annotation \
--tile-size $tile_size \
--min-pval $min_pval \
--min-count $min_count \
--min-height $min_height \
--bigWigPlus $pos_bw \
--bigWigMinus $neg_bw \
--output $output_file \
--threads $threads"
    echo "Running: $cmd"
    eval $cmd
}

# Function to rename and select final files
finalize_files() {
    local rbp_prefix="$1"
    local control_prefix="${rbp}_control"
    local eclip_prefix="${rbp}_eCLIP"

    # Check for merged control and eCLIP files and use them if they exist
    if [[ -f "${control_prefix}_merged.sorted.pos.bw" && -f "${control_prefix}_merged.sorted.neg.bw" ]]; then
        cp "${control_prefix}_merged.sorted.pos.bw" control.pos.bw
        cp "${control_prefix}_merged.sorted.neg.bw" control.neg.bw
    elif [[ -f "${control_prefix}_1.R2.sorted.pos.bw" && -f "${control_prefix}_1.R2.sorted.neg.bw" ]]; then
        cp "${control_prefix}_1.R2.sorted.pos.bw" control.pos.bw
        cp "${control_prefix}_1.R2.sorted.neg.bw" control.neg.bw
    fi

    if [[ -f "${eclip_prefix}_merged.sorted.pos.bw" && -f "${eclip_prefix}_merged.sorted.neg.bw" ]]; then
        cp "${eclip_prefix}_merged.sorted.pos.bw" signal.pos.bw
        cp "${eclip_prefix}_merged.sorted.neg.bw" signal.neg.bw
        cp "${eclip_prefix}_merged.sorted.crosslink.bed" peaks.crosslink.bed
    elif [[ -f "${eclip_prefix}_1.R2.sorted.pos.bw" && -f "${eclip_prefix}_1.R2.sorted.neg.bw" ]]; then
        cp "${eclip_prefix}_1.R2.sorted.pos.bw" signal.pos.bw
        cp "${eclip_prefix}_1.R2.sorted.neg.bw" signal.neg.bw
        cp "${eclip_prefix}_1.R2.sorted.crosslink.bed" peaks.crosslink.bed
    fi
}

# Process files for each prefix
for prefix in "${rbp}_control" "${rbp}_eCLIP"; do
    
    # Find files matching the prefix
    files=($(ls ${prefix}*.bam 2> /dev/null))

    # Check if files were found
    if [ ${#files[@]} -eq 0 ]; then
        echo "No files found for ${prefix}"
        continue
    fi

    # Process each BAM file
    for file in "${files[@]}"; do
        
        # Define the output filenames
        r2_bam="${file%.bam}.R2.bam"
        sorted_bam="${file%.bam}.R2.sorted.bam"

        # Extract R2 reads, sort, and index
        echo "Processing $file for R2 reads, sorting, and indexing..."
        samtools view -b -f 0x80 -F 0x4 "$file" > "$r2_bam"
        samtools sort "$r2_bam" -o "$sorted_bam"
        samtools index "$sorted_bam"

        # Clean up intermediate R2 BAM
        rm -f "$r2_bam"

        # Add sorted BAM to the respective array for later merging
        if [[ "$prefix" == "${rbp}_control" ]]; then
            control_sorted_bams+=("$sorted_bam")
        else
            eclip_sorted_bams+=("$sorted_bam")
        fi
    done
done

# Merge, sort, and index control replicates if more than one
if [ ${#control_sorted_bams[@]} -gt 1 ]; then
    merge_sort_index "${rbp}_control" "${control_sorted_bams[@]}"
fi

# Merge, sort, and index eCLIP replicates if more than one
if [ ${#eclip_sorted_bams[@]} -gt 1 ]; then
    merge_sort_index "${rbp}_eCLIP" "${eclip_sorted_bams[@]}"
fi

# Convert bedGraph files to bigWig and extract coverage
for bam_file in "${control_sorted_bams[@]}" "${eclip_sorted_bams[@]}"; do
    if [[ -f "$bam_file" ]]; then
        echo "Extracting 5' read-start coverage for $bam_file"
        extract_coverage "$bam_file"
    fi
done

# Extract coverage for merged control, if it exists
if [[ -f "${rbp}_control_merged.sorted.bam" ]]; then
    echo "Extracting 5' read-start coverage for ${rbp}_control_merged.sorted.bam"
    extract_coverage "${rbp}_control_merged.sorted.bam"
fi

# Extract coverage for merged eCLIP, if it exists
if [[ -f "${rbp}_eCLIP_merged.sorted.bam" ]]; then
    echo "Extracting 5' read-start coverage for ${rbp}_eCLIP_merged.sorted.bam"
    extract_coverage "${rbp}_eCLIP_merged.sorted.bam"
fi

for bedgraph_file in *.bedGraph; do
    if [[ -f "$bedgraph_file" ]]; then
        convert_to_bigwig "$bedgraph_file"
    fi
done

# Run enriched windows finding for eCLIP replicates and merged files
for bw_file in *eCLIP*.pos.bw; do
    neg_bw="${bw_file%.pos.bw}.neg.bw"
    if [[ -f "$neg_bw" ]]; then
        echo "Finding enriched windows for $bw_file and $neg_bw"
        find_enriched_windows "$bw_file" "$neg_bw"
    fi
done

# Rename and select final files
finalize_files "$rbp"
awk '$6 == "+" {print}' peaks.crosslink.bed > peaks.crosslink.pos.bed
awk '$6 == "-" {print}' peaks.crosslink.bed > peaks.crosslink.neg.bed

echo "Pipeline completed. Original BAMs, BigWigs, and peaks.bed files retained."
