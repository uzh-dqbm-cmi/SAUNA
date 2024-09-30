#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
import gzip
import os
import sys
from scipy.signal import medfilt

# In[3]:


def running_median(data, window_size):
    return medfilt(data, kernel_size=window_size)

def get_peaks(arr, window, max_gap=35,min_len =30):
    medians = running_median(arr,window)
    above = np.where(arr > medians)[0]
    
    def merge_and_filter_sequences(indices):
        # Convert the list of indices into a NumPy array
        indices_array = np.array(indices)

        # Find the indices where the gaps between elements exceed max_gap
        gap_indices = np.where(np.diff(indices_array) >= max_gap)[0] + 1

        # Split the indices array into subarrays based on the gaps
        sequences = np.split(indices_array, gap_indices)
        return sequences
    
    permissive_windows = merge_and_filter_sequences(above)
    permissive_windows = [sublist for sublist in permissive_windows if len(sublist) > 0]

    if len(permissive_windows) > 0:
    
        # Filter windows based on minimum length (>= 30)
        min_len_windows = [sublist for sublist in permissive_windows if np.max(sublist) - np.min(sublist) + 1 >= min_len]

        # Calculate peaks from the middle of each window        
       
        peaks = [sublist[len(sublist) // 2]  if len(sublist) % 2 != 0 else (sublist[len(sublist) // 2 - 1] + sublist[len(sublist) // 2]) // 2  for sublist in min_len_windows]
        peak_array=np.array(peaks,dtype=int)
        peak_values = arr[peak_array]
        return peak_array, peak_values, medians
    else:
        return [],[],[]


# In[ ]:


def process_file(file_path,output_directory):
    # Extract chromosome information from the file name
    file_name = os.path.basename(file_path)
    chromosome = file_name.split("_")[-1].split(".")[0]
    
    # Read the file using gzip
    with gzip.open(file_path, 'rt') as file:
        df = pd.read_csv(file, sep='\t', header=None)
        
    list2 = list((df[3]).astype(float))
    array_subset = np.array(list2)
    peaks, values, medians = get_peaks(array_subset, 375, 35, 10)  
    pos = np.array(peaks)
    
    # Create DataFrame with chromosome and position information
    df_output = pd.DataFrame({"chromosome": [chromosome] * len(peaks), "position": pos})
    
    # Define the output file path with chromosome information
    output_file_path = os.path.join(output_directory, f"MaxFinder_output_{chromosome}.tsv")
    
    # Write DataFrame to output file
    df_output.to_csv(output_file_path, sep='\t', index=False, header=False)
    print(f"Processed file: {file_path}. Output saved to {output_file_path}")

if __name__ == "__main__":
    # Check if the file path is provided as a command-line argument
    if len(sys.argv) != 3:
        print("Usage: python script.py file_path output_directory")
        sys.exit(1)
    
    # Get the file path from the command-line argument
    file_path = sys.argv[1]
    output_directory = sys.argv[2]
    # Check if the provided path is a file
    if not os.path.isfile(file_path):
        print(f"Error: {file_path} is not a valid file path.")
        sys.exit(1)
    
    # Process the file
    process_file(file_path, output_directory)

