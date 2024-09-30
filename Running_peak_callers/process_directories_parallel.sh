#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <base_directory>"
    exit 1
fi

# Get the base directory from the command-line argument
base_directory="$1"

# Export the Python script path to be used by parallel
export PYTHON_SCRIPT="smooth_fragments.py"

# Ensure GNU Parallel is installed
if ! command -v parallel &> /dev/null; then
    echo "GNU Parallel is not installed. Please install it to proceed."
    exit 1
fi

# Iterate through all directories in the base directory
find "$base_directory" -mindepth 1 -maxdepth 1 -type d | while read -r dir; do
    # Check if 'input' directory exists in the current directory
    input_dir="$dir/input"
    if [ -d "$input_dir" ]; then
        echo "Processing files in directory: $input_dir"
        
        # Find all .tsv.gz files in the 'input' directory and process them in parallel
        find "$input_dir" -name "*.tsv.gz" | parallel "python $PYTHON_SCRIPT {}"
    else
        echo "No 'input' directory found in $dir"
    fi
done
