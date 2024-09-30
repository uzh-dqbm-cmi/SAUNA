#!/usr/bin/env python
# coding: utf-8

# In[3]:


import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import os
import shutil
import sys
import pandas as pd


# In[ ]:


def process_file(input_file, output_directory):
    # Check if the input file exists
    if not os.path.isfile(input_file):
        print(f"Error: The file '{input_file}' does not exist.")
        return

    # Check if the output directory exists, create it if it doesn't
    if not os.path.exists(output_directory):
        print(f"The output directory '{output_directory}' does not exist. Creating it.")
        os.makedirs(output_directory)

if __name__ == "__main__":
    # Check if the correct number of arguments are provided
    if len(sys.argv) != 3:
        print("Usage: python process_file.py <input_file> <output_directory>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_directory = sys.argv[2]

    # Call the process_file function
    process_file(input_file, output_directory)


# In[ ]:


mean_scores = pd.read_csv(input_file, sep='\t', header=None)
column_names = ["type", "coverage", "MaxFinder", "Snyder","SAUNA"]
mean_scores.columns = column_names


# In[ ]:


# Define the linestyles for peak callers
line_styles = {
    "SAUNA": '-',
    "MaxFinder": ':',
    "Snyder": '--',
}

# Nucleosome configuration: Regular
nuc_config = "Regular"

# Create a plot for "Regular" nucleosome configuration
fig, (ax, ax_legend) = plt.subplots(2, 1, figsize=(6,4), gridspec_kw={'height_ratios': [4, 1]})

# Plot each peak caller with a unique linestyle
for peak_caller, linestyle in line_styles.items():
    subset = mean_scores[(mean_scores['type'] == nuc_config)]
    if peak_caller in subset.columns:
        ax.plot(subset['coverage'], subset[peak_caller], linestyle=linestyle, color='black', lw=3, 
                label=f'{peak_caller}')



# Create custom legend entries for scenario types (linestyles)
linestyle_handles = [mlines.Line2D([], [], color='black', linestyle=linestyle, label=linestyle_name) for linestyle_name, linestyle in line_styles.items()]


# Add the legends to the legend subplot

ax_legend.axis('off')  # Hide the axis for the legend subplot

ax.set_xlabel('Coverage', fontsize=17, labelpad = 10)
ax.set_ylabel('F1 Score', fontsize=17, labelpad = 10)
ax.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.2)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(0.5)
ax.spines['bottom'].set_linewidth(0.5)
ax.tick_params(axis='both', which='major', labelsize=14)  # Adjust label size here


# Place the legend outside of the plot
ax.legend(fontsize=11, title="Peak Caller", title_fontsize='12', frameon=False, loc='center left', bbox_to_anchor=(1, 0.5))

# Adjust layout and save the figure
plt.tight_layout(rect=[0, 0, 0.85, 1])
file_name = "F1_scores_regular.svg"
output_file_path = os.path.join(output_directory, file_name)

plt.savefig(output_file_path, format='svg', bbox_inches='tight')

plt.show()


# In[ ]:


# Define the linestyles for peak callers
line_styles = {
    "SAUNA": '-',
    "MaxFinder": ':',
    "Snyder": '--',
}


# Nucleosome configuration: Regular
nuc_config = "Regular + Phase-Shifted"

# Create a plot for "Regular" nucleosome configuration
fig, (ax, ax_legend) = plt.subplots(2, 1, figsize=(6,4), gridspec_kw={'height_ratios': [4, 1]})

# Plot each peak caller with a unique linestyle
for peak_caller, linestyle in line_styles.items():
    subset = mean_scores[(mean_scores['type'] == nuc_config)]
    if peak_caller in subset.columns:
        ax.plot(subset['coverage'], subset[peak_caller], linestyle=linestyle, color='black', lw=3, 
                label=f'{peak_caller}')


ax_legend.axis('off')  # Hide the axis for the legend subplot

ax.set_xlabel('Coverage', fontsize=17, labelpad = 10)
ax.set_ylabel('F1 Score', fontsize=17, labelpad = 10)
ax.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.2)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(0.5)
ax.spines['bottom'].set_linewidth(0.5)
ax.tick_params(axis='both', which='major', labelsize=14)  # Adjust label size here

# Place the legend outside of the plot
ax.legend(fontsize=11, title="Peak Caller", title_fontsize='12', frameon=False, loc='center left', bbox_to_anchor=(1, 0.5))

# Adjust layout and save the figure
plt.tight_layout(rect=[0, 0, 0.85, 1])
file_name = "F1_scores_phase_shifted.svg"
output_file_path = os.path.join(output_directory, file_name)

plt.savefig(output_file_path, format='svg', bbox_inches='tight')
plt.show()


# In[ ]:


# Define the linestyles for peak callers
line_styles = {
    "SAUNA": '-',
    "MaxFinder": ':',
    "Snyder": '--',
}


# Nucleosome configuration: Regular
nuc_config = "Random"

# Create a plot for "Regular" nucleosome configuration
fig, (ax, ax_legend) = plt.subplots(2, 1, figsize=(6,4), gridspec_kw={'height_ratios': [4, 1]})

# Plot each peak caller with a unique linestyle
for peak_caller, linestyle in line_styles.items():
    subset = mean_scores[(mean_scores['type'] == nuc_config)]
    if peak_caller in subset.columns:
        ax.plot(subset['coverage'], subset[peak_caller], linestyle=linestyle, color='black', lw=3, 
                label=f'{peak_caller}')



ax_legend.axis('off')  # Hide the axis for the legend subplot

ax.set_xlabel('Coverage', fontsize=17, labelpad = 10)
ax.set_ylabel('F1 Score', fontsize=17, labelpad = 10)
ax.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.2)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(0.5)
ax.spines['bottom'].set_linewidth(0.5)
ax.tick_params(axis='both', which='major', labelsize=14)  # Adjust label size here


# Place the legend outside of the plot
ax.legend(fontsize=11, title="Peak Caller", title_fontsize='12', frameon=False, loc='center left', bbox_to_anchor=(1, 0.5))

# Adjust layout and save the figure
plt.tight_layout(rect=[0, 0, 0.85, 1])
file_name = "F1_scores_random.svg"
output_file_path = os.path.join(output_directory, file_name)

plt.savefig(output_file_path, format='svg', bbox_inches='tight')
plt.show()

