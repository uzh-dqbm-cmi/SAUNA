#!/bin/bash

# Check if the input directory and parameter file are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 /path/to/input_directory /path/to/parameter_file"
    exit 1
fi

# Define the python script name
PYTHON_SCRIPT="SAUNA.py"

# Define the input directory containing the tsv.gz files
INPUT_DIR="$1"
PARAM_FILE="$2"


# Check if the directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Directory $INPUT_DIR does not exist."
    exit 1
fi

# Check if the parameter file exists
if [ ! -f "$PARAM_FILE" ]; then
    echo "Error: Parameter file $PARAM_FILE does not exist."
    exit 1
fi

# Loop through each tsv.gz file in the input directory
for GZ_FILE in "$INPUT_DIR"/*.tsv.gz
do
    # Determine the corresponding tsv file name
    BASENAME=$(basename "$GZ_FILE" .tsv.gz)
    TSV_FILE="$INPUT_DIR/MaxFinder_output_${BASENAME}.tsv"

    # Check if the corresponding tsv file exists
    if [ ! -f "$TSV_FILE" ]; then
        echo "Warning: Corresponding TSV file not found for $GZ_FILE. Skipping."
        continue
    fi

    # Print a message to indicate which files are being processed
    echo "Processing $GZ_FILE and $TSV_FILE"

    # Run the python script with both the tsv.gz and tsv files as inputs in the background
    python $PYTHON_SCRIPT "$TSV_FILE" "$GZ_FILE" "$PARAM_FILE" &

    # Wait for 5 minutes before starting the next process (300 seconds)
    sleep 300
done

# Wait for all background processes to finish
wait

# Print a completion message
echo "All files have been processed."
