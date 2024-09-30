# Analysis Scripts for Different Nucleosome Configurations on Simulated Data

This folder contains scripts for simulating nucleosome configurations and analyzing their corresponding F1 scores under various conditions. The scripts included are designed to facilitate the generation of simulated data and visualization of F1 scores across different scenarios and coverages.

## Scripts Overview

### 1. `data_simulation.py`
This script creates 5 `tsv.gz` files with simulated data for 5 different nucleosome configurations. It does not take any input parameters. However, it contains a parameter called `coverage`, which can be adjusted to simulate fragment center data for different coverages.

### 2. `F1_scores_over_time.py`
- `<temperature_file>`: A TSV file with two columns. The first column contains the temperatures, and the second column contains the annealing steps.
- `<F1_score_file>`: A TSV file with three columns. The first column contains the mean F1 scores for each annealing step for the scenario with only regular nucleosomes. The second column contains scores for the scenario with phase-shifted nucleosomes, and the third column contains scores for random nucleosomes.
- `<output_directory>`: The directory where the output plot will be saved.

**Output**: The script generates a plot showing mean F1 scores for each nucleosome configuration scenario over time (for each annealing step) and plots the temperature over time (for different annealing steps). The number of annealing steps (rows) in both the temperature file and the F1 score file must be the same.

### 3. `Mean_F1_different_coverages_plots.py`
- `<input_file>`: A TSV file containing the mean F1 scores for different scenarios, coverages, and peak callers. It should have 5 columns:
  - Column 1: Scenario (values should be 'Random', 'Regular', or 'Regular + Phase-Shifted')
  - Column 2: Coverage
  - Column 3: Mean F1 scores for the MaxFinder peak caller
  - Column 4: Mean F1 scores for the Snyder peak caller
  - Column 5: Mean F1 scores for the SAUNA peak caller
- `<output_directory>`: The directory where the output plot will be saved.

**Output**: The script creates a plot showing the mean F1 scores of different peak callers over different coverages, which is saved in the specified output directory.
