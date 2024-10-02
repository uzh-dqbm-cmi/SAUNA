#!/usr/bin/env python
# coding: utf-8

# In[8]:


import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os


# In[ ]:


# Main function
def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="")
    
    # Add arguments for chromosome lengths, blacklisted regions, and sample files
    parser.add_argument('different_iterations_file', type=str)
    parser.add_argument('different_peak_callers_file', type=str)
    parser.add_argument('F1_score_file', type=str)
    parser.add_argument('output_dir', type=str, help="Path to the output directory.")

    # Parse the arguments
    args = parser.parse_args()
    global different_iterations_file, different_peak_callers_file, F1_score_file,output_dir  # Declare all variables as global    
    # Load the input files
    different_iterations_file = args.different_iterations_file
    different_peak_callers_file = args.different_peak_callers_file
    F1_score_file = args.F1_score_file
    output_dir = args.output_dir


if __name__ == "__main__":
    main()


# In[ ]:


different_iterations = pd.read_csv(different_iterations_file, sep='\t', header=None)
column_names = ["type", "1x (With Input)", "1x", "10x","100x","1000x"]
different_iterations.columns = column_names


# In[ ]:


new_column_names = ['type','MaxFinder', 'Snyder', 'SAUNA1', 'SAUNA2','SAUNA3']
different_peak_callers = pd.read_csv(different_peak_callers_file, sep='\t', header=None)
different_peak_callers.columns = new_column_names


# In[14]:


F1_df = pd.read_csv('/mnt/DATA3/nucleosome/testing_scripts/simulated_data/F1_score_file_supplementary.tsv', sep='\t', header=None)
annealing_steps = F1_df[0].to_numpy()  

regular = F1_df[1].to_numpy()  
phase_shifted = F1_df[2].to_numpy()  
random = F1_df[3].to_numpy() 
with_random = F1_df[4].to_numpy() 
with_remove = F1_df[5].to_numpy() 


# In[ ]:


F1_df = pd.read_csv(F1_score_file, sep='\t', header=None)
annealing_steps = F1_df[0].to_numpy()  

regular = F1_df[1].to_numpy()  
phase_shifted = F1_df[2].to_numpy()  
random = F1_df[3].to_numpy() 
with_random = F1_df[4].to_numpy() 
with_remove = F1_df[5].to_numpy() 


# In[10]:


filtered_df = different_iterations
# Assign new column names
    
mean_scores = filtered_df.groupby('type').mean()
desired_row_order = [0, 1, 2,3, 4]  # Replace ... with remaining indices

# Reorder rows of the DataFrame
mean_scores = mean_scores.iloc[desired_row_order]
# Fill NaN values with a placeholder, e.g., 0
mean_scores.fillna(0, inplace=True)

# Create the heatmap using matplotlib
plt.figure(figsize=(12, 9), dpi =250)

plt.imshow(mean_scores, aspect='auto', cmap='coolwarm')

# Add color bar
cbar = plt.colorbar()
cbar.set_label('Mean F1 Score', rotation=270, labelpad=30,fontsize=23)
cbar.ax.yaxis.set_tick_params(labelsize=20) 

# Set axis labels
plt.xlabel('Number of Iterations',fontsize=30,labelpad =30)
plt.ylabel('Nucleosome Configuration', fontsize =30,labelpad=22)
# plt.title('Mean Scores Heatmap',fontsize=20, pad=40)

plt.xticks(np.arange(mean_scores.shape[1]), mean_scores.columns, rotation=45, ha='left', fontsize=23)
plt.yticks(np.arange(mean_scores.shape[0]), mean_scores.index, fontsize=23)

# Annotate each cell with the numeric value
for i in range(mean_scores.shape[0]):
    for j in range(mean_scores.shape[1]):
        plt.text(j, i, f'{mean_scores.iloc[i, j]:.2f}', ha='center', va='center', color='black', fontsize =20)

# Move x-axis labels to the top
plt.gca().xaxis.tick_top()
plt.gca().xaxis.set_label_position('top')
file_name = 'number_of_iteratios_heatmap.svg'
full_path = os.path.join(output_dir, file_name)
plt.savefig(full_path, format='svg')

plt.show()


# In[12]:


filtered_df = different_peak_callers
  
mean_scores = filtered_df.groupby('type').mean()
mean_scores.fillna(0, inplace=True)

plt.figure(figsize=(10, 9), dpi =250)

plt.imshow(mean_scores, aspect='auto', cmap='coolwarm')

cbar = plt.colorbar()
cbar.set_label('Mean F1 Score', rotation=270, labelpad=30,fontsize=23)
cbar.ax.yaxis.set_tick_params(labelsize=20) 

plt.xlabel('Peak Caller',fontsize=30,labelpad =22)
plt.ylabel('Nucleosome Configuration', fontsize =30,labelpad=22)

plt.xticks(np.arange(mean_scores.shape[1]), mean_scores.columns, rotation=45, ha='left', fontsize=23)
plt.yticks(np.arange(mean_scores.shape[0]), mean_scores.index, rotation=45, fontsize=23)  # Rotate y-axis labels

for i in range(mean_scores.shape[0]):
    for j in range(mean_scores.shape[1]):
        plt.text(j, i, f'{mean_scores.iloc[i, j]:.2f}', ha='center', va='center', color='black', fontsize =23)

plt.gca().xaxis.tick_top()
plt.gca().xaxis.set_label_position('top')
file_name = 'different_peak_callers_heatmap.svg'
full_path = os.path.join(output_dir, file_name)
plt.savefig(full_path, format='svg')

plt.show()


# In[19]:


plt.figure(figsize=(6, 4), dpi = 250)
x = annealing_steps
plt.plot(x,regular,c="indigo", zorder=1, label='Regular')
plt.plot(x,phase_shifted,c="darkgreen", zorder=1, label='Regular + Phase-Shifted')
plt.plot(x,random,c="orange", zorder=1, label='Random')
plt.plot(x,with_random,c="red", zorder=1, label='Regular + Random')
plt.plot(x,with_remove,c="turquoise", zorder=1, label='Regular + Random - Regular')

plt.xlabel("Annealing Steps", fontsize = 14, labelpad = 10)
plt.ylabel("F1 Score", fontsize =14, labelpad =10)
plt.legend()
plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.2)
ax = plt.gca()  # Get current axis

ax.spines['top'].set_visible(False)  # Hide the top spine
ax.spines['right'].set_visible(False)  # Hide the right spine
ax.spines['left'].set_linewidth(0.5)  # Customize the left spine
ax.spines['bottom'].set_linewidth(0.5)  # Customize the bottom spine

file_name = "F1_scores_over_time_supplementary.svg"
output_file_path = os.path.join(output_dir, file_name)

plt.savefig(output_file_path, format='svg', bbox_inches='tight')
plt.show()

