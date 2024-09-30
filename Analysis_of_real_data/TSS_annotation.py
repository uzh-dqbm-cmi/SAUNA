#!/usr/bin/env python
# coding: utf-8

# In[15]:


import argparse
import os
import pandas as pd
import gzip
import numpy as np
import matplotlib.pyplot as plt

# In[ ]:


def main(input_dir, gene_expression_file, tss_annotation_file):
    # Check if the input directory exists
    if not os.path.isdir(input_dir):
        raise FileNotFoundError(f"The input directory {input_dir} does not exist.")

    # Check if the gene expression file exists
    if not os.path.isfile(gene_expression_file):
        raise FileNotFoundError(f"The gene expression file {gene_expression_file} does not exist.")
    
    # Check if the TSS annotation file exists
    if not os.path.isfile(tss_annotation_file):
        raise FileNotFoundError(f"The TSS annotation file {tss_annotation_file} does not exist.")

    # Placeholder for further processing logic
    print(f"Input Directory: {input_dir}")
    print(f"Gene Expression File: {gene_expression_file}")
    print(f"TSS Annotation File: {tss_annotation_file}")

if __name__ == "__main__":
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Process gene expression and TSS annotation data.")
    parser.add_argument("input_dir", type=str, help="Path to the input directory.")
    parser.add_argument("gene_expression_file", type=str, help="Path to the gene expression file.")
    parser.add_argument("tss_annotation_file", type=str, help="Path to the TSS annotation file.")

    # Parse the arguments
    args = parser.parse_args()
    input_dir = args.input_dir
    gene_expression_file = args.gene_expression_file
    tss_annotation_file = args.tss_annotation_file
    # Call the main function with the provided arguments
    main(input_dir, gene_expression_file, tss_annotation_file)
 


# In[4]:


expression = pd.read_csv(gene_expression_file,delimiter='\t')
expression_subset = expression[["bone_marrow","GeneID"]]
del expression


# In[7]:


TSS = pd.read_csv(tss_annotation_file,delimiter='\t')
TSS = TSS[TSS['strand'] == "+"]
chromosomes_to_keep = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
       'chr9', 'chr10','chr11', 'chr12', 'chr13',
       'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
       'chr21', 'chr22']
TSS = TSS[TSS["chrom"].isin(chromosomes_to_keep)]
TSS.columns.values[0] = 'name'
TSS = TSS.drop_duplicates(subset='name2', keep='first').reset_index(drop=True)


# In[9]:


import glob

def center_scores_matrix(df_genes, df_expression, part_of_path):
    # Initialize variables
    window_size = 1000
    result = []

    # Group genes by chromosome for faster lookup
    genes_by_chrom = df_genes.groupby('chrom')

    # List all chromosome files
    chrom_files = sorted(glob.glob(part_of_path+'*.tsv.gz'))

    # Define a function to process each chromosome file
    def process_chrom_file(chrom_file, genes_by_chrom, window_size):
        chrom_result = []
        print("Processing: ",chrom_file)
        # Read relevant columns from the chromosome file

        columns_to_load = [0,1,3] 
        with gzip.open(chrom_file, 'rt') as f:
            dfs = []
            chunk_size = 5000000  # Adjust the chunk size as needed
            for chunk in pd.read_csv(f, sep='\t', chunksize=chunk_size,usecols=columns_to_load):
                dfs.append(chunk)
        chunk = pd.concat(dfs, ignore_index=True)
        chunk.columns = ['chromosome', 'position', 'wps score']
        print("File has been loaded.")

        # Extract chromosome name and add "chr" prefix
        chrom_name = 'chr' + str(chunk['chromosome'].iloc[0])
        print("Chromosome name extracted: ",chrom_name)

        if chrom_name not in genes_by_chrom.groups:
            return []

        genes_on_chrom = genes_by_chrom.get_group(chrom_name)

        for _, gene in genes_on_chrom.iterrows():
            tx_start = gene['txStart']
            name2 = gene['name2']

            # Calculate the start and end indices for the window around txStart
            start_idx = max(tx_start - window_size - 1, 0)  # -1 because row numbers start from 0
            end_idx = min(tx_start + window_size - 1, len(chunk) - 1)  # -1 because row numbers start from 0

            # Ensure indices are within the bounds of the data
            if start_idx < 0 or end_idx >= len(chunk):
                print("Gene indices not within bound of data: ", start_idx, end_idx)
                continue

            # Extract the WPS scores directly using row indices
            wps_scores = chunk.iloc[start_idx:end_idx + 1]['wps score'].values

            # Add the gene name and WPS scores to the result list
            chrom_result.append([name2] + list(wps_scores))

        print("Chromosome ", chrom_name, "done.")
        return chrom_result

    # Process each chromosome file
    for chrom_file in chrom_files:
        chrom_result = process_chrom_file(chrom_file, genes_by_chrom, window_size)
        result.extend(chrom_result)

    # Convert the result to a NumPy array
    columns = ['Gene'] + [f'pos_{i}' for i in range(-window_size, window_size + 1)]
    df_wps_matrix = pd.DataFrame(result, columns=columns)

    # Merge with expression data
    df_wps_matrix = df_wps_matrix.merge(df_expression, left_on='Gene', right_on='GeneID')

    # Sort by expression values
    df_wps_matrix = df_wps_matrix.sort_values(by='bone_marrow', ascending=False)

    # Drop unnecessary columns and convert to NumPy array
    wps_matrix_np = df_wps_matrix.drop(columns=['GeneID', 'bone_marrow']).to_numpy()

    return wps_matrix_np


# In[12]:


center_matrix_np_BH01 = center_scores_matrix(TSS, expression_subset,input_dir )


# In[ ]:


file_name = 'TSS_annotation_data.tsv.gz'
full_path = os.path.join(input_dir, file_name)

# Chunk size (number of rows per chunk)
chunk_size = 10000

# Open the gzip file for writing
with gzip.open(full_path, 'wt', encoding='utf-8') as f:
    for i in range(0, len(center_matrix_np_BH01), chunk_size):
        chunk = center_matrix_np_BH01[i:i + chunk_size]
        # Convert the chunk to tab-separated string
        chunk_str = '\n'.join(['\t'.join(map(str, row)) for row in chunk])
        f.write(chunk_str + '\n')


# In[16]:


BH01_wps = center_matrix_np_BH01[:,1:]
BH01_gene_names = center_matrix_np_BH01[:,0]
BH01_wps_added = np.sum(BH01_wps[:2000, :], axis=0)
BH01_wps_added_2D = BH01_wps_added.reshape(1, -1)
BH01_wps_added_last = np.sum(BH01_wps[-2000:, :], axis=0)
BH01_wps_added_last_2D = BH01_wps_added_last.reshape(1, -1)


# In[19]:


# Ensure BH01_wps_added_2D is numeric
BH01_wps_added_2D_numeric = pd.DataFrame(BH01_wps_added_2D).apply(pd.to_numeric, errors='coerce')

BH01_wps_added_2D_numeric = BH01_wps_added_2D_numeric.fillna(0)

# Now plot the heatmap
fig, ax = plt.subplots(figsize=(15, 3), dpi=250)

# Plot the heatmap for BH01
im = ax.imshow(BH01_wps_added_2D_numeric, aspect='auto', cmap='viridis')

# Set x-axis and y-axis labels and ticks
ax.set_xticks([])  # No x-axis ticks
ax.set_yticks([0])
ax.set_yticklabels(['BH01'], fontsize=25)  # BH01 label
ax.tick_params(axis='y', which='both', length=0)  # Remove tick lines but keep labels

# Add the color bar
cax = fig.add_axes([0.97, 0.14, 0.02, 0.7])
cbar = fig.colorbar(im, cax=cax)
cbar.ax.set_title('Fragment Centers', fontsize=18, pad=12)
cbar.ax.yaxis.set_ticks_position('left')
cbar.ax.tick_params(labelsize=12)

# Save the figure
file_name = 'TSS_annotation_plot.svg'
full_path = os.path.join(input_dir, file_name)
plt.savefig(full_path, format='svg')

plt.show()

