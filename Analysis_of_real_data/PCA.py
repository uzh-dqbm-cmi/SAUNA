#!/usr/bin/env python
# coding: utf-8

# In[6]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import argparse
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import os

# In[ ]:


# Main function
def main():
    global chromosome_lengths, blacklisted_regions, sample_files, output_dir, genes_file
    # Set up argument parsing

    parser = argparse.ArgumentParser(description="Process chromosome lengths, blacklisted regions, and sample file names.")
    
    # Add arguments for chromosome lengths, blacklisted regions, and sample files
    parser.add_argument('chromosome_lengths', type=str, help="Path to the chromosome lengths file.")
    parser.add_argument('blacklisted_regions', type=str, help="Path to the blacklisted regions file.")
    parser.add_argument('sample_files', type=str, help="Path to the file containing list of sample file names.")
    parser.add_argument('output_dir', type=str, help="Path to the output directory.")
    parser.add_argument('genes_file', type=str, help="Tsv.gz file with gene names and positions.")

    # Parse the arguments
    args = parser.parse_args()
    
    # Load the input files
    chromosome_lengths = args.chromosome_lengths
    blacklisted_regions = args.blacklisted_regions
    sample_files = args.sample_files
    output_dir = args.output_dir
    genes_file = args.genes_file


if __name__ == "__main__":
    main()


# In[ ]:


chromosome_information = pd.read_csv(chromosome_lengths, sep='\t', header = None)
chromosome_information = chromosome_information[chromosome_information[0]!='chrY']
chromosome_information = chromosome_information[chromosome_information[0]!='chrX']
chromosomes = np.array(chromosome_information[0])
chromosome_lengths = np.array(chromosome_information[1])


# In[ ]:


def read_bed_file(bed_file):
    """
    Reads a BED file containing blacklisted regions.

    Args:
    bed_file (str): Path to the BED file.

    Returns:
    pd.DataFrame: DataFrame with columns ['chromosome', 'start', 'end'].
    """
    # Read the BED file with no headers and provide column names
    blacklisted_regions = pd.read_csv(bed_file, sep='\t', header=None, names=['chromosome', 'start', 'end'])
    return blacklisted_regions

blacklisted_regions = read_bed_file(blacklisted_regions)


# In[ ]:


# Step 1: Concatenate all DataFrames
def z_scale(input_dfs):

    normalized_dfs = []
    for df in input_dfs:
        mean_interpeak = df['interpeak_distance'].mean()
        std_interpeak = df['interpeak_distance'].std()
        df['z_scaled_interpeak_distance'] = (df['interpeak_distance'] - mean_interpeak) / std_interpeak
        normalized_dfs.append(df)

    return normalized_dfs


# In[ ]:


def get_mean_interpeak_distances(sample_dfs):
    mean_interpeak_distances = []
    for df in sample_dfs:
        mean_distances = df.groupby('bin')['z_scaled_interpeak_distance'].mean()
        mean_interpeak_distances.append(mean_distances)
    return mean_interpeak_distances


# In[ ]:


def get_combined_df(unique_sample_names, mean_interpeak_distances):
    combined_data = []  # Create an empty list to collect DataFrames
    
    # Loop through each Series (df) in mean_interpeak_distances
    for i, df in enumerate(mean_interpeak_distances):
        sample_name = unique_sample_names[i]  # Get the sample name
        
        # Convert the Series to a DataFrame using its index as the 'bin' column
        mean_interpeak_distances[i] = pd.DataFrame({
            'sample_name': [sample_name] * len(df),  # Repeat the sample name for all rows
            'mean_interpeak_distance': df.values,  # Extract the values from the Series
            'bin': df.index  # Use the index as the 'bin' numbers
        })
        
        # Append the resulting DataFrame to the list
        combined_data.append(mean_interpeak_distances[i])
    
    # Concatenate all the DataFrames into one combined DataFrame
    combined_df = pd.concat(combined_data, ignore_index=True)
    
    return combined_df


# In[ ]:

def filter_blacklisted_regions(dfs, blacklisted_regions):
    """
    Removes rows from each DataFrame in dfs where the position falls within blacklisted regions.
    
    Args:
    dfs (list of pd.DataFrame): List of pandas DataFrames with columns ['bin', 'chromosome', 'start', 'end', 'position', 'interpeak_distance'].
    blacklisted_regions (pd.DataFrame): DataFrame of blacklisted regions with columns ['chromosome', 'start', 'end'].
    
    Returns:
    list of pd.DataFrame: List of DataFrames with rows filtered based on blacklisted regions.
    """
    # Create an IntervalIndex for blacklisted regions per chromosome
    blacklisted_intervals = {}
    for chrom in blacklisted_regions['chromosome'].unique():
        chrom_regions = blacklisted_regions[blacklisted_regions['chromosome'] == chrom]
        intervals = pd.IntervalIndex.from_tuples(list(zip(chrom_regions['start'], chrom_regions['end'])))
        blacklisted_intervals[chrom] = intervals
    
    filtered_dfs = []
    
    # Iterate over each DataFrame
    count = 0
    for df in dfs:
        count+=1
        print(count)
        # Filter the rows based on blacklist intervals
        mask = pd.Series([True] * len(df))
        
        for chrom in df['chromosome'].unique():
            if chrom in blacklisted_intervals:
                intervals = blacklisted_intervals[chrom]
                chrom_mask = (df['chromosome'] == chrom)
                positions = df.loc[chrom_mask, 'position']
                
                # Check if positions fall in any blacklist interval
                mask.loc[chrom_mask] = ~positions.apply(lambda pos: intervals.contains(pos).any())
        
        filtered_df = df[mask].copy()
        filtered_dfs.append(filtered_df)
    
    return filtered_dfs


# In[ ]:


def filter_interpeak_distance(dfs, threshold):
    """
    Filters rows in each DataFrame where interpeak_distance is above a certain threshold.

    Args:
    dfs (list of pd.DataFrame): List of pandas DataFrames.
    threshold (float): Threshold value for interpeak_distance.

    Returns:
    list of pd.DataFrame: List of DataFrames with rows filtered based on interpeak_distance.
    """
    filtered_dfs = []
    for df in dfs:
        # Remove rows where interpeak_distance is above the threshold
        filtered_df = df[df['interpeak_distance'] <= threshold].copy()
        filtered_dfs.append(filtered_df)
    
    return filtered_dfs


# In[ ]:


def bin_genome(chromosome_names, chromosome_lengths, num_bins=3000):
    bin_edges = []
    bin_size = sum(chromosome_lengths) / num_bins
    
    current_edge = 0
    
    # Iterate through chromosomes and calculate bin edges
    for chrom_name, chrom_length in zip(chromosome_names, chromosome_lengths):
        if chrom_name not in ['chrX', 'chrY']:
            num_bins_chrom = int(np.ceil(chrom_length / bin_size))
            for _ in range(num_bins_chrom):
                bin_edges.append(current_edge)
                current_edge += bin_size
    
    # Ensure the last bin edge covers the full range
    if len(bin_edges) < num_bins:
        bin_edges.append(current_edge)
    
    # Ensure bin edges are unique and sorted
    bin_edges = sorted(set(bin_edges))
    
    return bin_edges


bin_edges = bin_genome(chromosomes, chromosome_lengths, num_bins=3000)


# In[ ]:

sample_files = sample_files.split(',')
sample_dfs_SeCT = []

for file in sample_files:
    df = pd.read_csv(file, sep='\t', header=None, names=['chromosome','start','end', 'position', 'interpeak_distance'])
    sample_dfs_SeCT.append(df)


# In[ ]:


sample_dfs_SeCT = [df[~df['chromosome'].isin(['chrX', 'chrY'])] for df in sample_dfs_SeCT]


# In[ ]:


# Create a cumulative length mapping
cumulative_lengths = np.cumsum(chromosome_lengths)
chromosome_to_cumulative = dict(zip(chromosomes, np.concatenate(([0], cumulative_lengths[:-1]))))


# In[ ]:


def convert_to_global_position(df, cumulative_lengths_mapping):
    # Add a new column with the global position
    df['global_position'] = df.apply(lambda row: row['position'] + cumulative_lengths_mapping[row['chromosome']], axis=1)
    return df


# In[ ]:


def assign_bins(df, bin_edges):
    # Assign bins based on global coordinates
    df['bin'] = pd.cut(df['global_position'], bins=bin_edges, labels=False, include_lowest=True)
    return df



# In[ ]:


for df in sample_dfs_SeCT:
    df = convert_to_global_position(df, chromosome_to_cumulative)
    df = assign_bins(df, bin_edges)


# In[ ]:


without_blacklisted_sample_dfs_SeCT = filter_blacklisted_regions(sample_dfs_SeCT, blacklisted_regions)
filtered_sample_dfs_SeCT = filter_interpeak_distance(without_blacklisted_sample_dfs_SeCT, 350)
normalized_sample_dfs_SeCT = z_scale(filtered_sample_dfs_SeCT)
mean_interpeak_distances_SeCT = get_mean_interpeak_distances(normalized_sample_dfs_SeCT)
unique_sample_names_SeCT = ['SeCT-20','SeCT-22','SeCT-26_t1','SeCT-26_t2', 'SeCT-26_t3','SeCT-47','SeCT-58','SeCT-61','SeCT-82']
combined_df_SeCT = get_combined_df(unique_sample_names_SeCT, mean_interpeak_distances_SeCT)


# In[ ]:


name = ""

# Define the sample categories
high_tumor_fraction = ['SeCT-26_t1', 'SeCT-26_t2', 'SeCT-26_t3']
benign = ['SeCT-58', 'SeCT-82']
low_tumor_fraction = ['SeCT-20', 'SeCT-22', 'SeCT-47', 'SeCT-61']

# Create a mapping of sample names to categories
sample_category = {}
for sample in high_tumor_fraction:
    sample_category[sample] = 'High Tumor Fraction'
for sample in benign:
    sample_category[sample] = 'Benign'
for sample in low_tumor_fraction:
    sample_category[sample] = 'Low Tumor Fraction'

# Pivot the DataFrame to have a row for each sample and a column for each bin
pivoted_df = combined_df_SeCT.pivot(index='sample_name', columns='bin', values='mean_interpeak_distance').fillna(0)

# Standardize the pivoted DataFrame
scaler = StandardScaler()
scaled_df = scaler.fit_transform(pivoted_df)

# Perform PCA
pca = PCA()
pca_result = pca.fit_transform(scaled_df)

# Create a DataFrame containing the PCA results
pca_df = pd.DataFrame(data=pca_result, columns=[f'PC{i+1}' for i in range(pca_result.shape[1])], index=pivoted_df.index)

# Add category information to pca_df
pca_df['category'] = pca_df.index.map(sample_category)

# Define colors for each category
colors = {
    'High Tumor Fraction': 'red',
    'Low Tumor Fraction': 'orange',
    'Benign': 'blue'
}

# Plot the PCA results
plt.figure(figsize=(4, 4), dpi=200)

# Plot each category with a different color
for category, color in colors.items():
    category_df = pca_df[pca_df['category'] == category]
    plt.scatter(category_df['PC1'], category_df['PC2'], c=color, label=category, alpha=0.5)

# title = "PCA of Mean Interpeak Distances " + name
plt.xlabel('Principal Component 1', labelpad=15, fontsize=12)
plt.ylabel('Principal Component 2', fontsize=12)

# Add sample names
for i, sample in enumerate(pca_df.index):
    plt.text(pca_df.loc[sample, 'PC1'], pca_df.loc[sample, 'PC2'], sample, fontsize=10)

plt.legend(title='', loc='upper left', bbox_to_anchor=(1, 1))  # Position legend outside the plot
plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.2)
ax = plt.gca()  # Get current axis

ax.spines['top'].set_visible(False)  # Hide the top spine
ax.spines['right'].set_visible(False)  # Hide the right spine
ax.spines['left'].set_linewidth(0.5)  # Customize the left spine
ax.spines['bottom'].set_linewidth(0.5)  # Customize the bottom spine

file_name = 'PCA_plot.svg'
full_path = os.path.join(output_dir, file_name)
plt.savefig(full_path, format='svg')
plt.show()


# In[3]:


def get_genes_df(path,chromosomes,chromosome_lengths):
    # Read the compressed TSV file into a DataFrame
    df = pd.read_csv(path, sep='\t', compression='gzip')

    # Extract relevant columns
    genes_df = df[['#name', 'chrom', 'txStart', 'txEnd']]
    genes_df.columns = ['gene_name', 'chromosome', 'start', 'end']

    # Filter the DataFrame to include only the chromosomes in the list
    filtered_genes_df = genes_df[genes_df['chromosome'].isin(chromosomes)]

    # Create a cumulative length mapping (from earlier)
    cumulative_lengths = np.cumsum(chromosome_lengths)
    chromosome_to_cumulative = dict(zip(chromosomes, np.concatenate(([0], cumulative_lengths[:-1]))))

    # Function to convert start and end to global positions
    def convert_to_global_positions(df, cumulative_lengths_mapping):
        df['global_start'] = df.apply(lambda row: row['start'] + cumulative_lengths_mapping[row['chromosome']], axis=1)
        df['global_end'] = df.apply(lambda row: row['end'] + cumulative_lengths_mapping[row['chromosome']], axis=1)
        return df

    # Convert start and end to global positions
    filtered_genes_df = convert_to_global_positions(filtered_genes_df, chromosome_to_cumulative)

    return filtered_genes_df


# In[ ]:


genes_df = get_genes_df(genes_file, chromosomes, chromosome_lengths)


# In[7]:


def get_genes(top_positive, top_negative, bin_edges_SeCT, genes_df):
    def assign_bins(positions, bin_edges):
        # Find the bin index for each position
        return np.digitize(positions, bin_edges) - 1  # Subtract 1 to get zero-based index

    # Assign bin indices to start and end positions
    genes_df['start_bin'] = assign_bins(genes_df['global_start'].values, bin_edges_SeCT)
    genes_df['end_bin'] = assign_bins(genes_df['global_end'].values, bin_edges_SeCT)

    genes_fully_in_bins = genes_df[genes_df['start_bin'] == genes_df['end_bin']]

    # Find genes that are fully contained in the top 100 positive bins
    positive_genes_full = genes_fully_in_bins[genes_fully_in_bins['start_bin'].isin(top_positive)]['gene_name'].tolist()

    # Find genes that are fully contained in the top 100 negative bins
    negative_genes_full = genes_fully_in_bins[genes_fully_in_bins['start_bin'].isin(top_negative)]['gene_name'].tolist()
    return positive_genes_full, negative_genes_full


# In[5]:


import mygene

def get_gene_symbols(positive_genes, negative_genes):
    # Initialize MyGene.info query service
    mg = mygene.MyGeneInfo()

    # Step 1: Strip the version number from ENST IDs (e.g., removing '.4' at the end)
    positive_genes_clean = [gene.split('.')[0] for gene in positive_genes]
    negative_genes_clean = [gene.split('.')[0] for gene in negative_genes]

    # Step 2: Convert ENST IDs to standard gene symbols using MyGene
    positive_gene_info = mg.querymany(positive_genes_clean, scopes='ensembl.transcript', fields='symbol', species='human')
    negative_gene_info = mg.querymany(negative_genes_clean, scopes='ensembl.transcript', fields='symbol', species='human')

    # Step 3: Extract the gene symbols from the MyGene results
    positive_gene_symbols = [entry['symbol'] for entry in positive_gene_info if 'symbol' in entry]
    negative_gene_symbols = [entry['symbol'] for entry in negative_gene_info if 'symbol' in entry]
    
    return positive_gene_symbols, negative_gene_symbols


# In[ ]:


cancer_samples = ['SeCT-26_t1', 'SeCT-26_t2', 'SeCT-26_t3']
benign_samples = ['SeCT-58', 'SeCT-82','SeCT-20', 'SeCT-22', 'SeCT-47', 'SeCT-61']

# Separate cancer and benign samples
cancer_df = combined_df_SeCT[combined_df_SeCT['sample_name'].isin(cancer_samples)]
benign_df = combined_df_SeCT[combined_df_SeCT['sample_name'].isin(benign_samples)]

# Group by 'bin' and aggregate mean_interpeak_distance for each group
grouped_cancer = cancer_df.groupby('bin')['mean_interpeak_distance'].apply(list)
grouped_benign = benign_df.groupby('bin')['mean_interpeak_distance'].apply(list)

# Initialize lists to store results
cancer_lower_bins = []
benign_lower_bins = []

# Perform t-tests for each bin
for bin_id in grouped_cancer.index:
    if bin_id in grouped_benign.index:
        cancer_values = grouped_cancer[bin_id]
        benign_values = grouped_benign[bin_id]
        
        # Perform t-test
        t_stat, p_value = ttest_ind(cancer_values, benign_values, equal_var=False)
        
        # Check if significant (you can adjust the p-value threshold)
        if p_value <= 0.01:
            # Compare means to determine direction of difference
            if pd.Series(cancer_values).mean() < pd.Series(benign_values).mean():
                cancer_lower_bins.append(bin_id)
            else:
                benign_lower_bins.append(bin_id)

# cancer_lower_bins contains bins where cancer mean_interpeak_distance is lower
# benign_lower_bins contains bins where benign mean_interpeak_distance is lower


# In[ ]:


positive_transcripts, negative_transcripts = get_genes(cancer_lower_bins,benign_lower_bins, bin_edges, genes_df)    


# In[ ]:


positive_genes, negative_genes = get_gene_symbols(positive_transcripts, negative_transcripts)


# In[ ]:


from gprofiler import GProfiler

gp = GProfiler(return_dataframe=True)

positive_enrichment_fdr = gp.profile(organism='hsapiens',significance_threshold_method='fdr', query=positive_genes, sources=[ "REAC"])
negative_enrichment_fdr = gp.profile(organism='hsapiens', significance_threshold_method='fdr',query=negative_genes, sources=["REAC"])


# In[ ]:


file_name = 'lower_in_cancer_enrichment_terms.tsv'
full_path = os.path.join(output_dir, file_name)

positive_enrichment_fdr.to_csv(full_path, sep='\t', index=False)

file_name = 'lower_in_benign_enrichment_terms.tsv'
full_path = os.path.join(output_dir, file_name)

negative_enrichment_fdr.to_csv(full_path, sep='\t', index=False)

