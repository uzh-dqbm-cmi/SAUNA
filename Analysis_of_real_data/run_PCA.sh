
#!/bin/bash

#!/bin/bash

# Check if the directory is passed as an argument
if [ -z "$4" ]; then
    echo "Usage: $0 <directory> <chromosome_lengths> <blacklisted_regions> <genes_file>"
    exit 1
fi

# Set the directory to the first argument passed to the script
dir_to_iterate="$1"
chromosome_lengths="$2"
blacklisted_regions="$3"
genes_file="$4"

# Initialize an empty list to hold the filenames
file_list=""


# Iterate through all files in the specified directory
for file in "$dir_to_iterate"/*; do
    # If it's a file, add it to the list
    if [ -f "$file" ]; then
        # Append the filename to the list, separated by a comma
        file_list+="$file,"
    fi
done

# Remove the trailing comma from the list
file_list=${file_list%,}

# Check if list is not empty
if [ -n "$file_list" ]; then
    # Run the Python script with the list as an argument
    python PCA.py "$chromosome_lengths" "$blacklisted_regions" "$file_list" "$dir_to_iterate" "$genes_file"
else
    echo "No files found in the directory."
fi


