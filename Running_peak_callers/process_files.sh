#!/bin/bash

# Check if an argument is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <input_folder> <output_folder> <archive_folder>"
    exit 1
fi

# Assign input, output, and archive folder paths provided as arguments
input_folder="$1"
output_folder="$2"
archive_folder="$3"

# Check if the input folder exists
if [ ! -d "$input_folder" ]; then
    echo "Error: Input folder '$input_folder' does not exist."
    exit 1
fi

# Check if the output folder exists, if not create it
if [ ! -d "$output_folder" ]; then
    mkdir -p "$output_folder"
fi

# Check if the archive folder exists, if not create it
if [ ! -d "$archive_folder" ]; then
    mkdir -p "$archive_folder"
fi

# Function to process a single file
process_file() {
    local input_file="$1"  # Local variable declaration
    local output_folder="$2"  # Pass output_folder as an argument
    local archive_folder="$3"
    
    echo "Processing file: $input_file"
    echo "Output folder: $output_folder"  # Ensure output_folder is accessible

    # Run the Python script in the current directory
    python "./MaxFinder.py" "$input_file" "$output_folder"

    # Extract filename without extension
    filename=$(basename "$input_file" .tsv.gz)
    echo "Extracted filename: $filename"

    # Extract the number (or character) from the filename
    extracted=$(echo "$filename" | sed 's/.*_//; s/\.tsv//')
    echo "Extracted number or character: $extracted"
    echo "$output_folder"

    # Check if there are any .tsv files in the output folder with the extracted number (or character) flanked by underscores
    if [ -n "$(find "$output_folder" -maxdepth 1 -type f -name "*_$extracted.*")" ]; then
        echo "Corresponding .tsv file found in the output folder $output_folder, moving $input_file to $archive_folder"
        mv "$input_file" "$archive_folder"
    fi
}

# Export the function to make it available to parallel
export -f process_file

# Process files until there are no .tsv.gz files left in the input folder
while [ "$(find "$input_folder" -maxdepth 1 -type f -name "*.tsv.gz")" ]; do
    # Find all .tsv.gz files in the input folder (excluding subdirectories) and run the process_file function on each file
    find "$input_folder" -maxdepth 1 -type f -name "*.tsv.gz" | parallel -j 24 process_file {} "$output_folder" "$archive_folder"
done
