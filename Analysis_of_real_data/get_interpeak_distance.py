import pandas as pd
import os
import sys

def calculate_interpeak_distance(input_file):
    # Read the input bedgraph file
    data = pd.read_csv(input_file, sep='\t', header=None)
    column_names = ["chromosome","start","end", "position"]  # Column names for the data
    data = data.rename(columns=dict(enumerate(column_names)))
    
    # Sort the DataFrame by chromosome and position
    data_sorted = data.sort_values(by=['chromosome', 'position'])
    
    # Calculate interpeak distance for each chromosome separately
    data_sorted['interpeak_distance'] = data_sorted.groupby('chromosome')['position'].diff()
    
    # Prepare output file path
    output_file_path = os.path.join(os.path.dirname(input_file), 'interpeak_distance.bedgraph')
    
    # Save the results to the output file
    data_sorted.to_csv(output_file_path, sep='\t', index=False, header=False)

    print(f"Interpeak distances calculated and saved to: {output_file_path}")

# Example usage
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <path_to_combined_output.bedgraph>")
        sys.exit(1)
    
    input_bedgraph_file = sys.argv[1]  # Get the input file from command line argument
    calculate_interpeak_distance(input_bedgraph_file)
