#!/usr/bin/env python
# coding: utf-8

# This code contains code from from the peak-caller developed by Snyder et al. at the Shendure lab, licensed under MIT (see original license in this directory)
# The peak caller's code is Copyright (c) 2016 Jay Shendure's lab at the University of Washington
# For the paper of the project see Snyder, Matthew W. et al. Cell-free DNA Comprises an In Vivo Nucleosome Footprint that Informs Its Tissues-Of-Origin. Cell, Volume 164, Issue 1, 57 - 68, DOI: 10.1016/j.cell.2015.11.050 
#The code was re-written in python3 and refactored 

import pandas as pd
import numpy as np
import gzip
import csv
import statistics
import os
import sys
from scipy.signal import medfilt
from scipy.signal import savgol_filter


# In[120]:


def find_max_sum_contiguous_window_above_median(arr, median):
    above = np.where(arr > median)[0]
    gap_indices = np.where(np.diff(above) > 1)[0] + 1
    sequences = np.split(above, gap_indices)
    max_length_index = max(range(len(sequences)), key=lambda i: len(sequences[i]))
    chosen = sequences[max_length_index]
    start = chosen[0]
    end = chosen[-1]
    center = (start + end)//2
    return start, end, center

# In[89]:


def running_median(data, window_size):
    return medfilt(data, kernel_size=window_size)

def get_peaks_vectorized4(arr, window, max_gap=35,min_len =30, max_len=150):
    medians = running_median(arr,window)
    arr = arr - medians
    arr = savgol_filter(arr, window_length=21, polyorder=2, mode="constant")
    medians=medians-medians
    above = np.where(arr > medians)[0]
    above_values = arr[above]
  
    
    def merge_and_filter_sequences(indices, values):
        # Convert the list of indices into a NumPy array
        indices = np.array(indices)

        # Find the indices where the gaps between elements exceed max_gap
        gap_indices = np.where(np.diff(indices) >= max_gap)[0] + 1

        # Split the indices array into subarrays based on the gaps
        sequences = np.split(indices, gap_indices)
        sequences_values = np.split(values, gap_indices)

        return sequences, sequences_values
    
    permissive_windows, permissive_windows_values = merge_and_filter_sequences(above, above_values)

    if len(permissive_windows) > 0:
    
        # Filter windows based on minimum length (>= 30)
        numbers = [number for number,sublist in enumerate(permissive_windows) if (np.max(sublist) - np.min(sublist) + 1 >= min_len) and (np.max(sublist) - np.min(sublist) + 1 <= 450)]
        min_len_windows = [permissive_windows[number] for number in numbers]
        min_len_windows_values = [permissive_windows_values[number] for number in numbers]
        medians_list = [np.median(arr) for arr in min_len_windows_values]   
        windows_above_medians = [find_max_sum_contiguous_window_above_median(arr, median)[2] for arr, median in zip(min_len_windows_values, medians_list)]
        start = [find_max_sum_contiguous_window_above_median(arr, median)[0] for arr, median in zip(min_len_windows_values, medians_list)]
        start_array=np.array(start,dtype=int)
        
        end = [find_max_sum_contiguous_window_above_median(arr, median)[1] for arr, median in zip(min_len_windows_values, medians_list)]
        end_array=np.array(end,dtype=int)
        distance = end_array - start_array
        selection = distance <= max_len
        
        windows_above_medians_positions = [min_len_windows[number][int(index)] for number, index in enumerate(windows_above_medians)]
        peak_array=np.array(windows_above_medians_positions,dtype=int)
        peak_array = peak_array[selection]
        peak_values = arr[peak_array]
        return peak_array, peak_values, medians,arr
    else:
        return [],[],[],[]

# In[91]:


def process_file(file_path, output_directory):
    # Extract chromosome information from the file name
    file_name = os.path.basename(file_path)
    chromosome = os.path.splitext(file_name)[0]
    chromosome = os.path.splitext(chromosome)[0]    
    # Read the file using gzip
    columns_to_load = [3] 
    with gzip.open(file_path, 'rt') as f:
        dfs = []
        chunk_size = 5000000  # Adjust the chunk size as needed
        for chunk in pd.read_csv(f, sep='\t', chunksize=chunk_size,usecols=columns_to_load):
            dfs.append(chunk)
    df = pd.concat(dfs, ignore_index=True)

    array_subset = np.array((df).astype(float))
    array_subset = array_subset[:,0]
    peaks, values, medians,arr = get_peaks_vectorized4(array_subset, 1001, 5, 50,150)  
    pos = np.array(peaks)
    
    # Create DataFrame with chromosome and position information
    df_output = pd.DataFrame({"chromosome": [chromosome] * len(peaks), "position": pos})
    
    # Define the output file path with chromosome information
    output_file_path = os.path.join(output_directory, f"snyder_{chromosome}.tsv")    
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

