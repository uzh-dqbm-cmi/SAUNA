#!/bin/bash

# Check if the input BAM file is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <input.bam>"
    exit 1
fi

# Input BAM file
input_bam="$1"

# Output directory
output_dir="${input_bam%.bam}_by_chromosome_new"
mkdir -p "$output_dir"

# Function to process each chromosome
process_chromosome() {
    chrom=$1

    # Get chromosome length
    chrom_len=$(samtools view -H "$input_bam" | grep "@SQ" | grep "SN:${chrom}" | awk -F 'LN:' '{print $2}')

    # Generate counts for fragment centers, filter by length
    samtools view -F 0x10 "$input_bam" "$chrom" | \
        awk -v chrom="$chrom" '
        BEGIN {OFS="\t"}
        {
            # Calculate fragment length as the absolute value of the template length
            fragment_length = ($9 >= 0) ? $9 : -$9;
            if (fragment_length >= 120 && fragment_length <= 200) {
                # Calculate midpoint
                midpoint = int(($4 + $4 + fragment_length - 1) / 2);
                print midpoint;
            }
        }' | \
        sort -n | \
        uniq -c | \
        awk 'BEGIN {OFS="\t"} {print $2, $1}' > "${output_dir}/${chrom}_counts.tmp"

    # Fill in missing positions with zeros
    awk -v chrom="$chrom" -v chrom_len="$chrom_len" 'BEGIN {OFS="\t"; pos = 1} {
        if (NR == 1) {
            while (pos < $1) {
                print chrom, pos, 0;
                pos++;
            }
        }
        while (pos < $1) {
            print chrom, pos, 0;
            pos++;
        }
        print chrom, $1, $2;
        pos = $1 + 1;
    } END {
        while (pos <= chrom_len) {
            print chrom, pos, 0;
            pos++;
        }
    }' "${output_dir}/${chrom}_counts.tmp" | gzip > "${output_dir}/${chrom}.tsv.gz"

    # Clean up temporary file
    rm "${output_dir}/${chrom}_counts.tmp"

    echo "Processed $chrom"
}

export -f process_chromosome
export input_bam
export output_dir

# Get list of chromosomes from the BAM file
chromosomes=$(samtools idxstats "$input_bam" | cut -f1 | grep -v '*')

# Run the processing in parallel
echo "$chromosomes" | xargs -n 1 -P 4 -I {} bash -c 'process_chromosome "$@"' _ {}

echo "All chromosomes processed. Files are stored in $output_dir"

