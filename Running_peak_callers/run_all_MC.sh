#!/bin/bash

# Check if the main directory and parameter file are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 /path/to/main_directory /path/to/parameter_file"
    exit 1
fi

# Define the main directory and parameter file
MAIN_DIR="$1"
PARAM_FILE="$2"


# Check if the main directory exists
if [ ! -d "$MAIN_DIR" ]; then
    echo "Error: Directory $MAIN_DIR does not exist."
    exit 1
fi

# Check if the parameter file exists
if [ ! -f "$PARAM_FILE" ]; then
    echo "Error: Parameter file $PARAM_FILE does not exist."
    exit 1
fi

# Define the path to the run_MC_staggered.sh script within the main directory
script_dir=$(dirname "$0")
SCRIPT_PATH="$script_dir/run_MC_staggered.sh"

# Check if the run_MC_staggered.sh script exists
if [ ! -f "$SCRIPT_PATH" ]; then
    echo "Error: Script $SCRIPT_PATH does not exist."
    exit 1
fi

# Iterate through all subdirectories in the main directory
for SUBDIR in "$MAIN_DIR"/*/
do
    # Define the path to the input directory within the current subdirectory
    INPUT_DIR="${SUBDIR}input"

    # Check if the input directory exists
    if [ -d "$INPUT_DIR" ]; then
        # Print a message indicating the script is being run
        echo "Running run_MC_staggered.sh on $INPUT_DIR with parameters from $PARAM_FILE"

        # Run the run_MC_staggered.sh script on the input directory
        "$SCRIPT_PATH" "$INPUT_DIR" "$PARAM_FILE"

    else
        echo "Warning: Input directory $INPUT_DIR does not exist. Skipping."
    fi
done

# Print a completion message
echo "Script execution complete."
