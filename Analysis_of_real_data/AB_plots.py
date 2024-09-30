#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os


# In[9]:


def compartment_annotation(data, compartment):
    data = data.dropna() 
    compartment = compartment.dropna()
    # Perform a merge_asof operation
    merged_df = pd.merge_asof(data.sort_values(["position","chromosome"]), compartment.sort_values(["start","chromosome"]), by='chromosome', left_on='position', right_on='start', direction='backward')

    # Filter rows where position is within the start-end range
    merged_df = merged_df[(merged_df['position'] >= merged_df['start']) & (merged_df['position'] <= merged_df['end'])]

    return merged_df


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


# In[ ]:


compartments = pd.read_csv(compartment_annotation_path,delimiter='\t', header=None)
new_column_names2 = ['chromosome', 'start',"end", 'compartment']

compartments.columns = new_column_names2


# In[10]:


df = BH01_centers
compartment = compartments
df.sort_values(['chromosome', 'position'], inplace=True)
excluded_chromosomes = ['chrX', 'chrY']
df_filtered = df[~df['chromosome'].isin(excluded_chromosomes)]
df_comp = compartment_annotation(df_filtered, compartment)
df_comp['compartment_type'] = df_comp['compartment'].str[0]
binned_df = df_comp.groupby(['start','compartment',"chromosome"])["interpeak_distance"].median().reset_index(name='interpeak_distance')
binned_df["interpeak_distance"]=binned_df["interpeak_distance"].astype(int)
median_distances_df = binned_df.groupby(['interpeak_distance',"compartment"]).size().reset_index(name="counts")


# In[ ]:


# Get the directory of the bedgraph file
output_directory = os.path.dirname(file_path)

# Define your color mapping
color_mapping = {
    'A1': '#FF6347',    # Tomato
    'A2': '#F08080',    # LightCoral
    'B3': '#4682B4',    # Steel Blue
    'B4': '#87CEEB',    # Sky Blue
    'B1': '#1E90FF',    # Dodger Blue
    'B2': '#6495ED'     # Cornflower Blue
}

# Create the figure
plt.figure(figsize=(7, 6), dpi=250)  # Adjust the figure size if needed
sns.lineplot(x='interpeak_distance', y='counts', hue="compartment", data=median_distances_df, palette=color_mapping, linewidth=4)

# Set title and labels
plt.title('BH01', fontsize=25, pad=20)
plt.xlabel('Interpeak Distance', fontsize=18, labelpad=12)
plt.ylabel('Count', fontsize=18)
plt.xlim(185, 215)

# Customize legend
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [1, 2, 4, 0, 3, 5]  # New order of legend items
#plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], title='Compartment', fontsize=16, title_fontsize=18, loc='upper right', bbox_to_anchor=(1.4, 1), borderaxespad=0.)
plt.legend()
plt.tick_params(axis='both', which='major', labelsize=16)

# Save the figure in the same directory as the bedgraph file
output_file_path = os.path.join(output_directory, 'AB_compartment_plot.svg')
plt.savefig(output_file_path, format='svg', bbox_inches='tight')

# Show the figure
plt.show()

