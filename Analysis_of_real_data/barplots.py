#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import textwrap
import sys
import os

# In[ ]:


# Check if the correct number of arguments is provided
if len(sys.argv) != 3:
    print("Usage: python your_script.py <path_to_enrichment_results> <output_dir>")
    sys.exit(1)

# Get the file paths from the command-line arguments
file_path = sys.argv[1]
output_dir = sys.argv[2]

enrichment_results = pd.read_csv(file_path,usecols=[ 1, 5], sep=',')



# In[ ]:


# Function to wrap long text
def wrap_labels(label, width=30):
    return '\n'.join(textwrap.wrap(label, width))

# Apply text wrapping to the term names
enrichment_results['wrapped_term_name'] = enrichment_results['term_name'].apply(lambda x: wrap_labels(x, width=70))

# Sort the DataFrame by p-value for better visual representation
enrichment_results = enrichment_results.sort_values('negative_log10_of_adjusted_p_value', ascending=False)

# Set the figure size
enrichment_results_first_20 = enrichment_results.head(5)
plt.figure(figsize=(16, len(enrichment_results_first_20) *0.75), dpi=200)

enrichment_results_first_20 = enrichment_results_first_20.sort_values('negative_log10_of_adjusted_p_value', ascending=True)

# Create the horizontal bar plot
bar_width = 0.7  # Adjust this value for thinner bars
plt.barh(enrichment_results_first_20['wrapped_term_name'], 
         enrichment_results_first_20['negative_log10_of_adjusted_p_value'], 
         alpha=0.7, 
         color='midnightblue', 
         height=bar_width)  # Set the height for thinner bars

# Adjust y-limits to decrease space between bars
plt.ylim(-0.5, len(enrichment_results_first_20) - 0.5)  # Adjust as needed to control spacing

# Add labels and title
plt.xlabel('-log10(p-value)', fontsize=21,labelpad= 15)
ax = plt.gca()  # Get current axis

ax.spines['top'].set_visible(False)  # Hide the top spine
ax.spines['right'].set_visible(False)  # Hide the right spine
ax.spines['left'].set_linewidth(0.5)  # Customize the left spine
ax.spines['bottom'].set_linewidth(0.5)  # Customize the bottom spine
plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.2)
ax.tick_params(axis='both', which='major', labelsize=19)  # Adjust label size here
output_file_path = os.path.join(output_dir, 'enrichment_barplot.svg')
plt.savefig(output_file_path, format='svg', bbox_inches='tight')

# Display the plot
plt.tight_layout()
plt.show()


# In[3]:


data = {
    'sample_name': ['SeCT-20', 'SeCT-22', 'SeCT-26_t1','SeCT-26_t2','SeCT-26_t3','SeCT-47', 'SeCT-58', 'SeCT-61', 'SeCT-82' ],
    'tumor_fraction': [1.7, 0.7, 29, 32.7, 40.3, 0.5, 0.6, 2.5, 0.7],
    'sample_type': ['Low Tumor Fraction', 'Low Tumor Fraction', 'High Tumor Fraction', 'High Tumor Fraction', 'High Tumor Fraction', 'Low Tumor Fraction', 'Benign', 'Low Tumor Fraction','Benign']
}

df = pd.DataFrame(data)
# Set the aesthetic style of the plots

# Define the color palette
palette = {'High Tumor Fraction': 'red', 'Low Tumor Fraction': 'orange', 'Benign': 'blue'}
legend_order = ['High Tumor Fraction', 'Low Tumor Fraction', 'Benign']

# Create the bar plot with a narrower figure
plt.figure(figsize=(10, 7), dpi=250)
bar_plot = sns.barplot(
    x='sample_name', 
    y='tumor_fraction', 
    hue='sample_type', 
    data=df, 
    palette=palette,
    alpha=0.85
)

# Customize the title and axis labels with larger font sizes
bar_plot.set_xlabel('Sample Name', fontsize=22, labelpad=20)
bar_plot.set_ylabel('Tumor Fraction [%]', fontsize=22, labelpad=15)
bar_plot.tick_params(axis='x', labelsize=18)
bar_plot.tick_params(axis='y', labelsize=18)

# Customize the legend
handles, labels = plt.gca().get_legend_handles_labels()

# Create a mapping of labels to handles
label_to_handle = dict(zip(labels, handles))

# Sort handles and labels based on the defined order
sorted_handles = [label_to_handle[label] for label in legend_order]
sorted_labels = legend_order

# Add the legend with sorted handles and labels
plt.legend(
    sorted_handles, sorted_labels,
    fontsize='18',
    title='', title_fontsize='18',
    loc='upper center', bbox_to_anchor=(0.5, -0.45),
    ncol=3
)

# Rotate x-tick labels
for label in bar_plot.get_xticklabels():
    label.set_rotation(45)
    label.set_ha('right')  # Align labels to the right
bar_plot.grid(False)
ax = plt.gca()  # Get current axis

ax.spines['top'].set_visible(False)  # Hide the top spine
ax.spines['right'].set_visible(False)  # Hide the right spine
ax.spines['left'].set_linewidth(0.5)  # Customize the left spine
ax.spines['bottom'].set_linewidth(0.5)  # Customize the bottom spine

output_file_path = os.path.join(output_dir, 'tumor_fractions_SeCT.svg')
plt.savefig(output_file_path, format='svg', bbox_inches='tight')# Display the plot
plt.show()

