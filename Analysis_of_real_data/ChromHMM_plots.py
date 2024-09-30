#!/usr/bin/env python
# coding: utf-8

# In[16]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os
import numpy as np


# In[18]:


def normalize_AUC_new(df):
    # Step 1: Calculate the total number of counts for each compartment
    total_counts = {}
    for compartment in ['enhancer', 'promoter', 'repressed', 'transcription']:
        total_counts[compartment] = df[compartment].sum()
    # Step 2: Normalize counts to frequencies (i.e., make each compartment's counts sum to 1)
    for compartment in ['enhancer', 'promoter', 'repressed', 'transcription']:
        df[compartment] = df[compartment] / total_counts[compartment]

    # Step 3: Calculate the AUC for each compartment (should be equal after normalization)
    auc_values = {}
    for compartment in ['enhancer', 'promoter', 'repressed', 'transcription']:
        interpeak_distance = df['interpeak_distance'].values
        counts = df[compartment].values
        auc_values[compartment] = np.trapz(counts, interpeak_distance)

    # Step 4: Determine the target AUC (average AUC)
    target_auc = sum(auc_values.values()) / len(auc_values)

    # Step 5: Scale frequencies to equalize AUCs
    for compartment, auc in auc_values.items():
        norm_factor = target_auc / auc  # Scale to make AUC equal to target AUC
        df[compartment] *= norm_factor

    return df


# In[2]:


# Check if the correct number of arguments is provided
if len(sys.argv) != 3:
    print("Usage: python your_script.py <path_to_bedgraph> <path_to_compartment_annotation>")
    sys.exit(1)

# Get the file paths from the command-line arguments
file_path = sys.argv[1]
compartment_annotation_path = sys.argv[2]

# Read the data from the specified file path
BH01_centers = pd.read_csv(file_path, delimiter='\t', header=None)

# Define new column names
new_column_names = ['chromosome', 'start', 'end', 'position', 'interpeak_distance']

# Assign new column names
BH01_centers.columns = new_column_names

# Select the relevant columns
BH01_centers = BH01_centers.iloc[:, [0, 3, 4]]  


# In[3]:


compartments = pd.read_csv(compartment_annotation_path,delimiter='\t', header=None,usecols=[0, 1, 2, 3])
new_column_names2 = ['chromosome', 'start',"end", 'compartment']

compartments.columns = new_column_names2
compartments = compartments.sort_values(by=['chromosome', 'start'])


# In[8]:


# List of specific values to filter by
specific_values = ['12_Repressed','4_Strong_Enhancer', '5_Strong_Enhancer','1_Active_Promoter','11_Weak_Txn']

# Get the subset of the DataFrame
subset_df = compartments[compartments['compartment'].isin(specific_values)]


# In[9]:


# Define the substring to search for and the new value
def change_names(subset_df, substring, new_value):
    subset_df.loc[subset_df['compartment'].str.contains(substring), 'compartment'] = new_value
    return subset_df


# In[10]:


subset_df = change_names(subset_df, 'Enhancer','enhancer')
subset_df = change_names(subset_df, 'Promoter','promoter')
subset_df = change_names(subset_df, 'Repressed','repressed')
subset_df = change_names(subset_df, 'Txn','transcription')


# In[12]:


def compartment_annotation(data, compartment):
    data = data.dropna() 
    compartment = compartment.dropna()

    # Perform a merge_asof operation
    merged_df = pd.merge_asof(data.sort_values(["position","chromosome"]), compartment.sort_values(["start","chromosome"]), by='chromosome', left_on='position', right_on='start', direction='backward')

    # Filter rows where position is within the start-end range
    merged_df = merged_df[(merged_df['position'] >= merged_df['start']) & (merged_df['position'] <= merged_df['end'])]

    return merged_df


# In[19]:


df = BH01_centers
compartment = subset_df
start = 120
end = 260
# title = 'BH01'
df.sort_values(['chromosome', 'position'], inplace=True)
excluded_chromosomes = ['chrX', 'chrY']
df_filtered = df[~df['chromosome'].isin(excluded_chromosomes)]
df_filtered = df_filtered[(df_filtered['interpeak_distance']>=start) & (df_filtered['interpeak_distance']<=end) ]
df_comp = compartment_annotation(df_filtered, compartment)
   
bin_size = 100_000
df_comp['bin'] = df_comp.apply(lambda row: row['position'] // bin_size, axis=1)

median_interpeak_distance = df_comp.groupby(['chromosome', 'bin', 'compartment'])['interpeak_distance'].median().reset_index()
median_interpeak_distance["interpeak_distance"]=median_interpeak_distance["interpeak_distance"].astype(int)

counts = median_interpeak_distance.groupby(['interpeak_distance', 'compartment']).size().reset_index(name='count')
    
pivot_counts = counts.pivot(index='interpeak_distance', columns='compartment', values='count')

    # Fill missing values with 0
pivot_counts = pivot_counts.fillna(0)
pivot_counts = pivot_counts.reset_index()
pivot_counts = normalize_AUC_new(pivot_counts)
#     Plot the data using matplotlib
plt.figure(figsize=(7, 6), dpi=250)  # Adjust the figure size if needed

    # Plot each compartment separately
for compartment in pivot_counts.columns[1:]:
    plt.plot(pivot_counts["interpeak_distance"].tolist(), pivot_counts[compartment].tolist(), label=compartment, linewidth = 4)

plt.title('BH01', fontsize =25,pad = 20)
plt.xlabel('Interpeak Distance', fontsize=18, labelpad= 15)
plt.ylabel('Frequency', fontsize=18, labelpad = 15)
# plt.legend(title='Compartment', fontsize=16, title_fontsize=18,loc='upper center', bbox_to_anchor=(0.45, -0.25), ncol=2, borderaxespad=0.)
plt.tick_params(axis='both', which='major', labelsize=16)

plt.xlim(150,225)
plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.2)
ax = plt.gca()  # Get current axis

ax.spines['top'].set_visible(False)  # Hide the top spine
ax.spines['right'].set_visible(False)  # Hide the right spine
ax.spines['left'].set_linewidth(0.5)  # Customize the left spine
ax.spines['bottom'].set_linewidth(0.5)  # Customize the bottom spine

output_directory = os.path.dirname(file_path)

output_file_path = os.path.join(output_directory, 'ChromHMM_compartment_plot.svg')
plt.savefig(output_file_path, format='svg', bbox_inches='tight')
plt.show()

