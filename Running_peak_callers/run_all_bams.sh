#!/bin/bash

# Usage: ./run_analysis.sh /path/to/bam/files
# $1 - Directory containing BAM files

# Check if the correct number of arguments are supplied
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <bam_directory> <params.txt>"
    exit 1
fi

# Assign input argument to a variable
bam_dir="$1"
parameter_file="$2"

# Check if bam_dir exists
if [ ! -d "$bam_dir" ]; then
    echo "Error: BAM directory does not exist!"
    exit 1
fi

# Get the directory where this script is located (and where the analysis script is located)
script_dir=$(dirname "$0")

# Set the path to the analysis script
analysis_script="$script_dir/fragment_center_calculation.sh"

# Check if the analysis script exists and is executable
if [ ! -f "$analysis_script" ] || [ ! -x "$analysis_script" ]; then
    echo "Error: Analysis script does not exist or is not executable!"
    exit 1
fi

# Loop through all BAM files in the directory
for bam_file in "$bam_dir"/*.bam; do
    # Check if bam_file exists to prevent running empty loops
    if [ ! -e "$bam_file" ]; then
        echo "No BAM files found in the directory."
        exit 1
    fi

    # Run the analysis script, passing the BAM file as input
    echo "Running analysis for $bam_file..."
    "$analysis_script" "$bam_file"

    if [ $? -eq 0 ]; then
        echo "Analysis completed for $bam_file"
    else
        echo "Analysis failed for $bam_file"
    fi
done

# Create input, output, archive, and other directories in each BAM file's subdirectory
# Define lists of valid TSV.GZ file names outside the loop
valid_files=("chr"{1..22}".tsv.gz" "chrX.tsv.gz" "chrY.tsv.gz" {1..22}".tsv.gz" "X.tsv.gz" "Y.tsv.gz")

# Create input, output, archive, and other directories in each BAM file's subdirectory
# Define lists of valid TSV.GZ file names outside the loop
valid_files=("chr"{1..22}".tsv.gz" "chrX.tsv.gz" "chrY.tsv.gz" {1..22}".tsv.gz" "X.tsv.gz" "Y.tsv.gz")

# Create input, output, archive, and other directories in each BAM file's subdirectory
for subdir in "$bam_dir"/*; do
    if [ -d "$subdir" ]; then
        echo "Processing directory: $subdir"
        mkdir -p "$subdir/input" "$subdir/output" "$subdir/archive" "$subdir/other"

        # Iterate through each .tsv.gz file in the subdirectory
        for file in "$subdir/"*.tsv.gz; do
            # Check if file exists (avoids issues if no tsv.gz files exist)
            if [ -f "$file" ]; then
                # Get the basename of the file (filename only without the path)
                basename=$(basename "$file")

                # Check if the file matches any in the valid_files list
                if [[ " ${valid_files[*]} " == *" $basename "* ]]; then
                    # Check if the file has at least 10 rows
                    row_count=$(zcat "$file" | wc -l)

                    if [ "$row_count" -ge 10 ]; then
                        # Move file to input directory if it has at least 10 rows
                        mv "$file" "$subdir/input/"
                        echo "Moved $basename to $subdir/input/ (has $row_count rows)"
                    else
                        # Move file to other directory if it has less than 10 rows
                        mv "$file" "$subdir/other/"
                        echo "Moved $basename to $subdir/other/ (only $row_count rows)"
                    fi
                else
                    # Move file to other directory if it doesn't match the valid_files list
                    mv "$file" "$subdir/other/"
                    echo "Moved $basename to $subdir/other/ (invalid name)"
                fi
            fi
        done
    fi
done

# Set the path to the process_directories_parallel.sh script
process_script="$script_dir/process_directories_parallel.sh"

# Check if the process_directories_parallel.sh script exists and is executable
if [ ! -f "$process_script" ] || [ ! -x "$process_script" ]; then
    echo "Error: process_directories_parallel.sh does not exist or is not executable!"
    exit 1
fi

# Run the process_directories_parallel.sh script with the BAM directory as input
echo "Processing all TSV.GZ files in the input directories..."
"$process_script" "$bam_dir"

if [ $? -eq 0 ]; then
    echo "Processing completed successfully!"
else
    echo "Processing failed!"
fi


run_all_script="$script_dir/run_all_simple.sh"
# Finally, run the run_all_simple.sh script with the BAM directory as an argument
echo "Running run_all_simple.sh with the BAM directory..."
"$run_all_script" "$bam_dir"

if [ $? -eq 0 ]; then
    echo "run_all_simple.sh completed successfully!"
else
    echo "run_all_simple.sh failed!"
fi

# Move all files from the archive and output directories to the input directory
echo "Moving files from archive and output directories to input directories..."
for sample_directory in "$bam_dir"/*/; do
    input_directory="$sample_directory/input"
    archive_directory="$sample_directory/archive"
    output_directory="$sample_directory/output"

    # Move files from archive to input directory
    if [ -d "$archive_directory" ]; then
        echo "Moving files from '$archive_directory' to '$input_directory'..."
        mv "$archive_directory/"* "$input_directory/" 2>/dev/null
    else
        echo "Warning: Archive directory '$archive_directory' does not exist."
    fi

    # Move files from output to input directory
    if [ -d "$output_directory" ]; then
        echo "Moving files from '$output_directory' to '$input_directory'..."
        mv "$output_directory/"* "$input_directory/" 2>/dev/null
    else
        echo "Warning: Output directory '$output_directory' does not exist."
    fi
done

echo "All files have been moved to the input directories."


run_all_MC_script="$script_dir/run_all_MC.sh"

# Run the run_all_MC.sh script with the BAM directory as an argument
echo "Running run_all_MC.sh with the BAM directory..."
"$run_all_MC_script" "$bam_dir" "$parameter_file"

if [ $? -eq 0 ]; then
    echo "run_all_MC.sh completed successfully!"
else
    echo "run_all_MC.sh failed!"
fi

# Move all .bed files from the input directory to the output directory
echo "Moving all .bed files from input to output directories..."
for sample_directory in "$bam_dir"/*/; do
    input_directory="$sample_directory/input"
    output_directory="$sample_directory/output"

    # Move .bed files from input to output directory
    if [ -d "$input_directory" ] && [ -d "$output_directory" ]; then
        echo "Moving .bed files from '$input_directory' to '$output_directory'..."
        mv "$input_directory/"*.bed "$output_directory/" 2>/dev/null
    else
        echo "Warning: Input or output directory does not exist in '$sample_directory'."
    fi
done

echo "All .bed files have been moved to the output directories."






