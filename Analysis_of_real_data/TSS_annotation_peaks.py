#!/usr/bin/env python
# coding: utf-8

# In[125]:


import pandas as pd
import numpy as np
import gzip
import argparse
import os
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt

# In[3]:


def main(output_dir, gene_expression_file, tss_annotation_file):
    # Check if the output directory exists
    if not os.path.isdir(output_dir):
        raise FileNotFoundError(f"The output directory {output_dir} does not exist.")

    # Check if the gene expression file exists
    if not os.path.isfile(gene_expression_file):
        raise FileNotFoundError(f"The gene expression file {gene_expression_file} does not exist.")
    
    # Check if the TSS annotation file exists
    if not os.path.isfile(tss_annotation_file):
        raise FileNotFoundError(f"The TSS annotation file {tss_annotation_file} does not exist.")

    # Placeholder for further processing logic
    print(f"Output Directory: {output_dir}")
    print(f"Gene Expression File: {gene_expression_file}")
    print(f"TSS Annotation File: {tss_annotation_file}")

if __name__ == "__main__":
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Process gene expression and TSS annotation data.")
    parser.add_argument("output_dir", type=str, help="Path to the output directory.")
    parser.add_argument("gene_expression_file", type=str, help="Path to the gene expression file.")
    parser.add_argument("tss_annotation_file", type=str, help="Path to the TSS annotation file.")

    # Parse the arguments
    args = parser.parse_args()
    output_dir = args.output_dir
    gene_expression_file = args.gene_expression_file
    tss_annotation_file = args.tss_annotation_file
    # Call the main function with the provided arguments
    main(output_dir, gene_expression_file, tss_annotation_file)


expression = pd.read_csv(gene_expression_file,delimiter='\t')


# In[9]:


expression_subset = expression[["bone_marrow","GeneID"]]

# In[13]:


TSS = pd.read_csv(tss_annotation_file,delimiter='\t')


# In[65]:


chromosomes_to_keep = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
       'chr9', 'chrX', 'chrY', 'chr10', 'chr11', 'chr12', 'chr13',
       'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
       'chr21', 'chr22']


# In[66]:


TSS = TSS[TSS["chrom"].isin(chromosomes_to_keep)]


# In[70]:


TSS = TSS[TSS["strand"] == '+']


# In[73]:


TSS = TSS.drop_duplicates(subset='name2', keep='first').reset_index(drop=True)


# In[81]:


TSS.columns.values[0] = 'name'



# In[48]:

file_name = 'interpeak_distance.bedgraph'
full_path = os.path.join(output_dir, file_name)

BH01_simulation = pd.read_csv(full_path,delimiter='\t', header=None)
new_column_names = ['chromosome', 'start',"end", 'position',"interpeak_distance"]

BH01_simulation.columns = new_column_names
BH01_simulation = BH01_simulation.iloc[:, [0, 3, 4]]  # Assuming indices start from 0


# In[91]:


def save_compressed_tsv(data, path):
    with gzip.open(path, 'wt') as f:
        np.savetxt(f, data, fmt='%s', delimiter='\t')


# In[133]:


def get_heatmap_matrix(data_df, tss_df, gene_expression_df,path ):
    # Convert pandas DataFrames to NumPy arrays
    data_array = data_df.to_numpy()
    tss_array = tss_df.to_numpy()
    gene_expression_array = gene_expression_df.to_numpy()

    # Extract chromosome and position columns from data_array
    chromosomes = data_array[:, 0]
    positions = data_array[:, 1]

    # Create an array to store the heatmap
    heatmap = np.zeros((len(tss_array), 2001), dtype=int)  # 2001 positions (-1000 to 1000)

    # Iterate over TSS positions
    for i, (_, chrom, _, tss_position, _) in enumerate(tss_array):
        # Find positions within +/- 1000 bp of the TSS
        within_range = np.logical_and(chromosomes == chrom,
                                      np.logical_and(positions >= tss_position - 1000,
                                                     positions <= tss_position + 1000))
        # Get positions relative to the TSS
        relative_positions = positions[within_range] - (tss_position - 1000)

        # Convert relative_positions to integer type
        relative_positions = relative_positions.astype(int)

        # Increment the corresponding positions in the heatmap
        heatmap[i, relative_positions] = 1

    # Calculate gene expression levels and sort heatmap by expression levels
    gene_expression_levels = gene_expression_array[:, 0]  
    gene_names = gene_expression_array[:, 1]  
    sorted_indices = np.argsort(gene_expression_levels)[::-1]  # Sort indices in descending order

    # Sort heatmap and gene_names by sorted_indices
    sorted_heatmap = heatmap[sorted_indices]
    sorted_gene_names = gene_names[sorted_indices]
    combined_data = np.column_stack((sorted_gene_names, sorted_heatmap))
    
    # Save combined_data as TSV
    save_compressed_tsv(combined_data, path)

# In[90]:

full_path = os.path.join(output_dir, "heatmap_with_gene_names_nucleosome_centers.tsv.gz")
get_heatmap_matrix(BH01_simulation, TSS, expression_subset,path = full_path)


# %%
def read_matrix(path):
    BH01_simple = pd.read_csv(path,delimiter='\t',compression = 'gzip', header=None)
    BH01_simple = BH01_simple.to_numpy()
    BH01_simple_values = BH01_simple[:,1:].astype(int)

    BH01_simple_gene_names = BH01_simple[:,0]
    BH01_simple_values_added = np.sum(BH01_simple_values[:2000, :], axis=0)
    BH01_simple_values_added_2D = BH01_simple_values_added.reshape(1, -1)

    return BH01_simple_values, BH01_simple_gene_names, BH01_simple_values_added_2D, BH01_simple_values_added


# %%
BH01_centers, BH01_centers_gene_names, BH01_centers_values_added2D,BH01_centers_values_added = read_matrix(full_path)

# Create the X-axis values representing positions for BH01
positions = np.arange(len(BH01_centers_values_added))

# Smooth the data using a Savitzky-Golay filter
smoothed_data = savgol_filter(BH01_centers_values_added, window_length=100, polyorder=2)

# Create the plot
fig, ax = plt.subplots(figsize=(10, 3), dpi=250)  # Adjust the figure size as needed
plt.plot(positions, smoothed_data, color='b', label='BH01')

# Set the x-axis ticks and labels
ax.set_xticks(np.arange(0, 2001, 500))  # Adjust the tick positions
ax.set_xticklabels(np.arange(-1000, 1001, 500))  # Label the ticks

# Set labels and title
title = 'Peak counts for 2000 most expressed genes in BH01'
plt.xlabel('Position relative to TSS', fontsize=13, labelpad=10)
plt.ylabel('Peak Counts', fontsize=13, labelpad=12)
plt.tick_params(axis='both', which='major', labelsize=12)

# Customize the plot appearance
plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.2)
ax.spines['top'].set_visible(False)  # Hide the top spine
ax.spines['right'].set_visible(False)  # Hide the right spine
ax.spines['left'].set_linewidth(0.5)  # Customize the left spine
ax.spines['bottom'].set_linewidth(0.5)  # Customize the bottom spine

# Add legend and save the plot
plt.legend(fontsize=13)

file_name = 'TSS_annotation_peaks_plot.svg'
full_path = os.path.join(output_dir, file_name)
plt.savefig(full_path, format='svg')
# Show the plot
plt.show()
