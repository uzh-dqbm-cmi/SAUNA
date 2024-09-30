import pandas as pd
import numpy as np
import scipy.ndimage
from whittaker_eilers import WhittakerSmoother
import gzip
import sys

def smooth_fragment_centers(fragment_center_array):
    whittaker_smoother = WhittakerSmoother(
        lmbda=1000, order=2, data_length=len(fragment_center_array))
    smoothed_fragment_centers = np.array(whittaker_smoother.smooth(fragment_center_array))
    sigma = 30  # Adjust sigma to control the width of the smoothing
    smoothed_fragment_centers = scipy.ndimage.gaussian_filter1d(smoothed_fragment_centers, sigma)
    return smoothed_fragment_centers

def process_file(file_path):
    chunk_size = 5000000  # Adjust chunk size as needed
    # Read the file in chunks without header
    columns_to_load = [0, 1, 2]  # Assume columns are in this order

    # Temporary list to store chunks
    dfs = []

    # Read the file in chunks
    with gzip.open(file_path, 'rt') as f:
        for chunk in pd.read_csv(f, sep='\t', chunksize=chunk_size, header=None, usecols=columns_to_load):
            
            # Extract fragment centers from the third column (index 2)
            fragment_centers = chunk.iloc[:, 2].values  # Third column, which should be fragment_centers
            smoothed_centers = smooth_fragment_centers(fragment_centers)

            # Add the smoothed data as a new column
            chunk['smoothed_fragment_centers'] = smoothed_centers

            # Append the processed chunk to the list
            dfs.append(chunk)

    # Concatenate all chunks into a single DataFrame
    df_output = pd.concat(dfs, ignore_index=True)

    # Write the updated DataFrame back to the same file
    with gzip.open(file_path, 'wt') as f:
        df_output.to_csv(f, sep='\t', index=False, header=False)

def main():
    if len(sys.argv) != 2:
        print("Usage: python smooth_fragments.py <file_path>")
        sys.exit(1)
    
    file_path = sys.argv[1]
    
    process_file(file_path)

if __name__ == '__main__':
    main()
