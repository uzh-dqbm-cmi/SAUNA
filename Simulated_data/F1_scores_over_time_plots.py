#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

# In[ ]:


def process_file(temperature_file,F1_score_file, output_directory):
    # Check if the input file exists
    if not os.path.isfile(temperature_file):
        print(f"Error: The file '{temperature_file}' does not exist.")
        return
    if not os.path.isfile(F1_score_file):
        print(f"Error: The file '{F1_score_file}' does not exist.")
        return

    # Check if the output directory exists, create it if it doesn't
    if not os.path.exists(output_directory):
        print(f"The output directory '{output_directory}' does not exist. Creating it.")
        os.makedirs(output_directory)

if __name__ == "__main__":
    # Check if the correct number of arguments are provided
    if len(sys.argv) != 4:
        print("Usage: python script.py <temperature_file> <F1_score_file> <output_directory>")
        sys.exit(1)

    temperature_file = sys.argv[1]
    F1_score_file = sys.argv[2]
    output_directory = sys.argv[3]

    # Call the process_file function
    process_file(temperature_file,F1_score_file, output_directory)


# In[ ]:


temperature_df = pd.read_csv(temperature_file, sep='\t', header=None)
temperature = temperature_df[0].to_numpy()  
annealing_steps = temperature_df[1].to_numpy()  


# In[ ]:


plt.figure(figsize=(7, 5), dpi = 250)
x = annealing_steps
y = temperature
plt.plot(x,y,c="blue", zorder=1)

plt.xlabel("Annealing Steps", fontsize = 18, labelpad = 10)
plt.ylabel("Temperature [K]", fontsize =18, labelpad =10)
plt.tick_params(axis='both', which='major', labelsize=13)
plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.2)
ax = plt.gca()  # Get current axis

ax.spines['top'].set_visible(False)  # Hide the top spine
ax.spines['right'].set_visible(False)  # Hide the right spine
ax.spines['left'].set_linewidth(0.5)  # Customize the left spine
ax.spines['bottom'].set_linewidth(0.5)  # Customize the bottom spine
file_name = "temperature_over_time.svg"
output_file_path = os.path.join(output_directory, file_name)

plt.savefig(output_file_path, format='svg', bbox_inches='tight')
plt.show()


# In[ ]:


F1_df = pd.read_csv(F1_score_file, sep='\t', header=None)
regular = F1_df[0].to_numpy()  
phase_shifted = F1_df[1].to_numpy()  
random = F1_df[2].to_numpy() 


# In[ ]:


plt.figure(figsize=(6, 4), dpi = 250)
x = annealing_steps
plt.plot(x,regular,c="indigo", zorder=1, label='Regular')
plt.plot(x,phase_shifted,c="darkgreen", zorder=1, label='Regular + Phase-Shifted')
plt.plot(x,random,c="orange", zorder=1, label='Random')

plt.xlabel("Annealing Steps", fontsize = 14, labelpad = 10)
plt.ylabel("F1 Score", fontsize =14, labelpad =10)
plt.legend()
plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.2)
ax = plt.gca()  # Get current axis

ax.spines['top'].set_visible(False)  # Hide the top spine
ax.spines['right'].set_visible(False)  # Hide the right spine
ax.spines['left'].set_linewidth(0.5)  # Customize the left spine
ax.spines['bottom'].set_linewidth(0.5)  # Customize the bottom spine

file_name = "F1_scores_over_time.svg"
output_file_path = os.path.join(output_directory, file_name)

plt.savefig(output_file_path, format='svg', bbox_inches='tight')
plt.show()

