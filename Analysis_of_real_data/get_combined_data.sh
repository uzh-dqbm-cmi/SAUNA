#!/bin/bash

# Check if BAM directory and compartment annotation file are provided as arguments
if [ $# -ne 5 ]; then
    echo "Usage: $0 <bam_directory> <AB_compartment_annotation_file> <ChromHMM_compartment_annotation_file> <gene_expression_file> <TSS_annotation_file>"
    exit 1
fi

# Set the input BAM directory and compartment annotation file to the provided arguments
bam_directory="$1"
compartment_annotation_file="$2"
second_compartment_annotation_file="$3"
gene_expression="$4"
TSS_annotation="$5"


# Iterate over each subdirectory in the BAM directory
for subdir in "$bam_directory"/*/; do
    # Define the output directory within the current subdirectory
    output_dir="${subdir}output"
    input_dir="${subdir}input"

    echo "Processing output directory: $output_dir"

    # Iterate over each .bed file in the output directory of the current subdirectory
    for file in "$output_dir"/*.bed; do
        # Check if the file exists
        if [[ -f "$file" ]]; then
            # Extract the filename without extension
            filename=$(basename "$file" .bed)

            # Apply the awk command to the current file and save the output as a .bedgraph file in the output directory
            awk 'BEGIN{OFS="\t"} {midpoint = int(($2 + $3) / 2); print $1, $2,$3, midpoint}' "$file" > "$output_dir/${filename}_processed.bedgraph"

            echo "Processed $file into ${filename}_processed.bedgraph"
        fi
    done

    echo "All .bed files processed in $output_dir."

    # Output file name for concatenated .bedgraph files
    output_file_bedgraph="combined_output.bedgraph"

    # Concatenate all .bedgraph files in the output directory into the output file
    cat "$output_dir"/*.bedgraph >> "$output_dir/$output_file_bedgraph" 2>/dev/null

    echo "All .bedgraph files combined into $output_file_bedgraph in $output_dir."

    # Replace 'chrchr' with 'chr' in the combined .bedgraph file, if it exists
    if [[ -f "$output_dir/$output_file_bedgraph" ]]; then
        sed -i 's/chrchr/chr/' "$output_dir/$output_file_bedgraph"
    fi
    # Run the Python script on the combined_output.bedgraph
    python get_interpeak_distance.py "$output_dir/$output_file_bedgraph"
    echo "Ran Python script on $output_file_bedgraph"
    
    # Add the interpeak_distance.bedgraph file name to the list
    file_name="$output_dir/interpeak_distance.bedgraph"
    file_list+=("$file_name")

    python AB_plots.py "$output_dir/interpeak_distance.bedgraph" "$compartment_annotation_file"
    echo "Ran AB_plots.py on $output_dir/interpeak_distance.bedgraph"
    
    python ChromHMM_plots.py "$output_dir/interpeak_distance.bedgraph" "$second_compartment_annotation_file"
    echo "Ran ChromHMM_plots.py on $output_dir/interpeak_distance.bedgraph"

    python TSS_annotation.py "$input_dir" "$gene_expression" "$TSS_annotation"
    echo "Ran TSS_annotation.py on $input_dir"

    python TSS_annotation_peaks.py "$output_dir" "$gene_expression" "$TSS_annotation"
    echo "Ran TSS_annotation_peaks.py on $output_dir"

done


