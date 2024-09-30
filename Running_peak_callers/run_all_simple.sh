#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <bam_directory>"
    exit 1
fi

# Get the main directory from the first argument
bam_directory="$1"

# Check if the process_files.sh script exists in the current directory
script_dir=$(dirname "$0")
process_script="$script_dir/process_files.sh"

if [ ! -f "$process_script" ]; then
    echo "Error: process_files.sh not found in $script_dir"
    exit 1
fi

# Iterate over each subdirectory in the BAM directory
for subdirectory in "$bam_directory"/*/; do
    # Define the input, output, and archive directories within the subdirectory
    input_directory="${subdirectory}input"
    output_directory="${subdirectory}output"
    archive_directory="${subdirectory}archive"

    # Check if the subdirectory contains the required directories
    if [ -d "$input_directory" ] && [ -d "$output_directory" ] && [ -d "$archive_directory" ]; then
        # Run the process_files.sh script with the appropriate arguments
        echo "Processing subdirectory: $subdirectory"
        "$process_script" "$input_directory" "$output_directory" "$archive_directory"
    else
        echo "Warning: Missing required directories in $subdirectory"
    fi
done

echo "All processing is complete!"
