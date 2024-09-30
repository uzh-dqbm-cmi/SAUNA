#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import random
from whittaker_eilers import WhittakerSmoother
import pandas as pd
import gzip
import scipy.ndimage

# In[ ]:


def save_data_as_tsv_gz(array, filename):
    # Create DataFrame
    df = pd.DataFrame({
        'Column1': [1] * len(array),
        'Column2': range(1, len(array) + 1),
        'Column3': [0] * len(array),
        'Column4': array
    })

    # Save DataFrame as gzipped TSV
    with gzip.open(filename, 'wt', encoding='utf-8') as f:
        df.to_csv(f, sep='\t', index=False, header=False)


# In[ ]:


def get_nucleosome_positions_phase_shift(length, num_regular_peaks, num_phase_shifted_peaks):    
    regular_peak_positions = np.linspace(0, length-200, num_regular_peaks, dtype=int)
    regular_peak_positions +=100
    phase_shift = 100
    shifted_peak_positions = regular_peak_positions[0:num_phase_shifted_peaks]+phase_shift

    nucleosome_positions = np.concatenate((regular_peak_positions, shifted_peak_positions))
    nucleosome_positions = np.sort(nucleosome_positions)

    return nucleosome_positions


# In[ ]:


def generate_random_nucleosome_positions(length,num_positions):
    nucleosome_positions = []
    for _ in range(num_positions):
        position = random.randint(0, length-200)
        position +=100
        nucleosome_positions.append(position)
    return np.array(nucleosome_positions)


# In[ ]:


def get_nucleosome_positions_with_random(length, num_regular_peaks, num_phase_shifted_peaks):
    regular_peak_positions = np.linspace(0, length-200, num_regular_peaks, dtype=int)
    regular_peak_positions +=100
    random_indices = np.random.choice(num_regular_peaks, num_phase_shifted_peaks, replace=False)
    chosen_regular_peak_positions = regular_peak_positions[random_indices]
    num_extra_shifts = num_phase_shifted_peaks * 2  # Generate twice as many to ensure enough non-zero shifts
    shifts = np.random.randint(-150, 150, num_extra_shifts)
    # Filter out 0s from the shifts
    shifts = shifts[shifts != 0]
    # Take only the required number of non-zero shifts
    shifts = shifts[:num_phase_shifted_peaks]    
    phase_shifted_peaks = chosen_regular_peak_positions + shifts
    # Ensure all phase-shifted peaks are within the range of 0 and 30000
    phase_shifted_peak_positions = phase_shifted_peaks[(phase_shifted_peaks >= 0) & (phase_shifted_peaks <= 20000)]
    # Ensure we have exactly num_phase_shifted_peaks phase-shifted peaks
    phase_shifted_peak_positions = phase_shifted_peak_positions[:num_phase_shifted_peaks]
    # Concatenate regular and phase-shifted peaks
    nucleosome_positions = np.concatenate((regular_peak_positions, phase_shifted_peak_positions))
    # If the total number of nucleosomes is less than 150, randomly add more
    num_additional_nucleosomes = 150 - len(nucleosome_positions)
    if num_additional_nucleosomes > 0:
        additional_positions = np.random.randint(0, length, num_additional_nucleosomes)
        nucleosome_positions = np.concatenate((nucleosome_positions, additional_positions))
    nucleosome_positions = np.sort(nucleosome_positions)
    return nucleosome_positions


# In[ ]:


def get_nucleosome_positions_with_remove(length, num_regular_peaks, num_phase_shifted_peaks,num_remove):
    regular_peak_positions = np.linspace(0, length-200, num_regular_peaks, dtype=int)
    regular_peak_positions += 100
    random_indices = np.random.choice(num_regular_peaks, num_phase_shifted_peaks, replace=False)
    chosen_regular_peak_positions = regular_peak_positions[random_indices]
    #randomly remove number og elements
    indices_to_remove = random.sample(range(len(regular_peak_positions)), num_remove)
    regular_peak_positions = np.delete(regular_peak_positions, indices_to_remove)
   
    num_extra_shifts = num_phase_shifted_peaks * 2  # Generate twice as many to ensure enough non-zero shifts
    shifts = np.random.randint(-150, 150, num_extra_shifts)
    # Filter out 0s from the shifts
    shifts = shifts[shifts != 0]
    # Take only the required number of non-zero shifts
    shifts = shifts[:num_phase_shifted_peaks]    
    phase_shifted_peaks = chosen_regular_peak_positions + shifts
    # Ensure all phase-shifted peaks are within the range of 0 and 30000
    phase_shifted_peak_positions = phase_shifted_peaks[(phase_shifted_peaks >= 0) & (phase_shifted_peaks <= 20000)]
    # Ensure we have exactly num_phase_shifted_peaks phase-shifted peaks
    phase_shifted_peak_positions = phase_shifted_peak_positions[:num_phase_shifted_peaks]
    # Concatenate regular and phase-shifted peaks
    nucleosome_positions = np.concatenate((regular_peak_positions, phase_shifted_peak_positions))
    # If the total number of nucleosomes is less than 150, randomly add more
    num_additional_nucleosomes = num_regular_peaks+num_phase_shifted_peaks-num_remove - len(nucleosome_positions)
    if num_additional_nucleosomes > 0:
        additional_positions = np.random.randint(0, length, num_additional_nucleosomes)
        nucleosome_positions = np.concatenate((nucleosome_positions, additional_positions))
    nucleosome_positions = np.sort(nucleosome_positions)

    return nucleosome_positions


# In[ ]:


def get_probabilities(nucleosome_positions, length):
    probabilities = np.zeros(length+200)
    for i in nucleosome_positions:

        probabilities[(i-20):(i+21)] += 30
    
    
    return probabilities


# In[ ]:


def smooth_probabilities(probabilities):
    whittaker_smoother = WhittakerSmoother(
        lmbda=1000, order=1, data_length=len(probabilities)
    )

    smoothed_probabilities = np.array(whittaker_smoother.smooth(probabilities))
    return smoothed_probabilities


# In[ ]:


def smooth_fragment_centers(fragment_center_array):
    whittaker_smoother = WhittakerSmoother(
    lmbda=1000, order=2, data_length=len(fragment_center_array))
    smoothed_fragment_centers = np.array(whittaker_smoother.smooth(fragment_center_array))
    sigma = 30  # Adjust sigma to control the width of the smoothing
    smoothed_fragment_centers = scipy.ndimage.gaussian_filter1d(smoothed_fragment_centers, sigma)
    return smoothed_fragment_centers


# In[ ]:


def normalize_probabilities(probs):
    min_val = np.min(probs)
    max_val = np.max(probs)
    normalized_probs = (probs - min_val) / (max_val - min_val)
    return normalized_probs


# In[ ]:


def get_positions(normalized_probabilities, number_of_fragments):
    new_normalized_probalities = normalized_probabilities / np.sum(normalized_probabilities)  # Normalize

    # Create the cumulative distribution function (CDF)
    cdf = np.cumsum(new_normalized_probalities)

    # Generate uniform random numbers
    random_numbers = np.random.rand(number_of_fragments)

    # Use the CDF to get positions
    positions = np.searchsorted(cdf, random_numbers)
    return positions


# In[ ]:


def create_fragment_center_array(centers, length_of_sequence):
    counts = np.bincount(centers, minlength=length_of_sequence)
    return counts


# In[ ]:


def create_fragment_array(centers, length_of_sequence, min_fragment_length, max_fragment_length):
    num_fragments = len(centers)
    
    # Create a zero array with dimensions (num_fragments, length_of_sequence)
    fragment_array = np.zeros((num_fragments, length_of_sequence+200), dtype=int)
    
    # Iterate over each fragment center position
    for i, center in enumerate(centers):
        # Generate a random fragment length within the specified range
        fragment_length = np.random.randint(min_fragment_length, max_fragment_length + 1)
        
        # Calculate the start and end positions of the fragment
        start_pos = max(center - fragment_length // 2, 0)
        end_pos = min(center + fragment_length // 2 + 1, length_of_sequence+200)
        
        # Fill in the positions of the fragment in the 2D array
        fragment_array[i, start_pos:end_pos] = 1
    
    return fragment_array


# In[ ]:


def calculate_wps_optimized(fragment_array, length_of_sequence, window_size=60):
    num_fragments = fragment_array.shape[0]
    window_span = 2 * window_size +1
    
    # Calculate cumulative sums along the sequence dimension
    cumulative_sum = np.cumsum(fragment_array, axis=1)
    
    # Get the sum of windows using cumulative sums
    padded_cumsum = np.pad(cumulative_sum, ((0, 0), (window_size, window_size)), mode='constant')
    window_sums = padded_cumsum[:, window_span:] - padded_cumsum[:, :-window_span]
    
    # Calculate spanning counts
    spanning_counts = np.sum(window_sums == window_span, axis=0)
    
    # Calculate ending counts
    ending_counts = np.sum(window_sums > 0, axis=0) - spanning_counts
    
    # Calculate the WPS for each position
    wps_array = spanning_counts - ending_counts
    wps_array = wps_array[98:length_of_sequence+98]
    return wps_array



# In[ ]:


def get_data_phase_shift(length, num_regular_peaks, num_phase_shifted_peaks, min_fragment_length, max_fragment_length, number_of_fragments):
    nucleosome_positions = get_nucleosome_positions_phase_shift(length, num_regular_peaks, num_phase_shifted_peaks)    
    probabilities = get_probabilities(nucleosome_positions,length)
    smoothed_probabilities = smooth_probabilities(probabilities)
    normalized_probabilities = normalize_probabilities(smoothed_probabilities)
    positions = get_positions(normalized_probabilities, number_of_fragments)
    wps_array = create_fragment_center_array(positions, length)
    wps_array = smooth_fragment_centers(wps_array)
    wps_array = wps_array[:length]

    link = "smoothed_array_phase_shift.tsv.gz"
    link2 = "nucleosome_positions_phase_shift.tsv"
    np.savetxt(link2, nucleosome_positions, delimiter='\t', fmt='%d')
    save_data_as_tsv_gz(wps_array,link )    

    


# In[ ]:


def get_data_only_regular(length, num_regular_peaks, min_fragment_length, max_fragment_length,number_of_fragments):
    nucleosome_positions = get_nucleosome_positions_phase_shift(length, num_regular_peaks, 0)
    probabilities = get_probabilities(nucleosome_positions,length)
    smoothed_probabilities = smooth_probabilities(probabilities)
    normalized_probabilities = normalize_probabilities(smoothed_probabilities)
    positions = get_positions(normalized_probabilities, number_of_fragments)

    wps_array = create_fragment_center_array(positions, length)
    wps_array = smooth_fragment_centers(wps_array)
    wps_array = wps_array[:length]

    
    link = "smoothed_array_only_regular.tsv.gz"
    link2 = "nucleosome_positions_only_regular.tsv"
    np.savetxt(link2, nucleosome_positions, delimiter='\t', fmt='%d')
    save_data_as_tsv_gz(wps_array,link )


# In[ ]:


def get_data_phase_only_random(length, num_peaks,min_fragment_length,max_fragment_length,number_of_fragments):
    nucleosome_positions = generate_random_nucleosome_positions(length,num_peaks)
    probabilities = get_probabilities(nucleosome_positions,length)
    smoothed_probabilities = smooth_probabilities(probabilities)
    normalized_probabilities = normalize_probabilities(smoothed_probabilities)
    positions = get_positions(normalized_probabilities, number_of_fragments)

    wps_array = create_fragment_center_array(positions, length)
    wps_array = smooth_fragment_centers(wps_array)
    wps_array = wps_array[:length]

    link = "smoothed_array_only_random.tsv.gz"
    link2 = "nucleosome_positions_only_random.tsv"
    np.savetxt(link2, nucleosome_positions, delimiter='\t', fmt='%d')
    save_data_as_tsv_gz(wps_array,link )


# In[ ]:


def get_data_with_random(length, num_regular_peaks, num_phase_shifted_peaks, min_fragment_length, max_fragment_length,number_of_fragments):
    nucleosome_positions = get_nucleosome_positions_with_random(length, num_regular_peaks, num_phase_shifted_peaks)
    probabilities = get_probabilities(nucleosome_positions,length)
    smoothed_probabilities = smooth_probabilities(probabilities)

    normalized_probabilities = normalize_probabilities(smoothed_probabilities)
    positions = get_positions(normalized_probabilities, number_of_fragments)

    wps_array = create_fragment_center_array(positions, length)
    wps_array = smooth_fragment_centers(wps_array)
    wps_array = wps_array[:length]

    
    link = "smoothed_array_with_random.tsv.gz"
    link2 = "nucleosome_positions_with_random.tsv"
    np.savetxt(link2, nucleosome_positions, delimiter='\t', fmt='%d')
    save_data_as_tsv_gz(wps_array,link )
    


# In[ ]:


def get_data_phase_with_remove(length, num_regular_peaks, num_phase_shifted_peaks, min_fragment_length, max_fragment_length,num_remove,number_of_fragments):
    nucleosome_positions = get_nucleosome_positions_with_remove(length, num_regular_peaks, num_phase_shifted_peaks,num_remove)
    probabilities = get_probabilities(nucleosome_positions,length)
    smoothed_probabilities = smooth_probabilities(probabilities)
    normalized_probabilities = normalize_probabilities(smoothed_probabilities)
    positions = get_positions(normalized_probabilities, number_of_fragments)

    wps_array = create_fragment_center_array(positions, length)
    wps_array = smooth_fragment_centers(wps_array)
    wps_array = wps_array[:length]


    link = "smoothed_array_with_remove.tsv.gz"
    link2 = "nucleosome_positions_with_remove.tsv"
    np.savetxt(link2, nucleosome_positions, delimiter='\t', fmt='%d')
    save_data_as_tsv_gz(wps_array,link )


# In[ ]:


min_fragment_length = 120
max_fragment_length = 180
length = 20_000
coverage = 100
number_of_fragments = int(coverage * length / 167)



get_data_phase_with_remove(length, 60, 50, min_fragment_length, max_fragment_length,10,number_of_fragments)
get_data_phase_only_random(length, 100,min_fragment_length,max_fragment_length,number_of_fragments)
get_data_with_random(length, 66, 34, min_fragment_length, max_fragment_length,number_of_fragments)
get_data_only_regular(length, 100, min_fragment_length, max_fragment_length,number_of_fragments)
get_data_phase_shift(length, 66, 34, min_fragment_length, max_fragment_length,number_of_fragments)
        
        

