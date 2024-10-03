# Analysis Scripts for Different Nucleosome Configurations on Simulated Data

This folder contains scripts for simulating nucleosome configurations and analyzing their corresponding F1 scores under various conditions. The scripts included are designed to facilitate the generation of simulated data and visualization of F1 scores across different scenarios and coverages.

## Input File Creation

### Simulated Data Generation
- The `data_simulation.py` script generates 5 `tsv.gz` files containing simulated nucleosome configuration data. The `coverage` parameter in the script can be adjusted to simulate fragment center data for various coverages.

### F1 Score and Temperature Over Time
  - **Temperature File**: The temperature for each annealing step is saved in a list during the simulation, and at the end of the simulation, both the temperature and corresponding annealing steps are saved as a TSV file.
  - **F1 Score File**: The creation of the F1 score file involves multiple steps:
    1. The `data_simulation.py` script is run 5 times, each time generating 5 datasets (one for each nucleosome configuration). 
    2. The SAUNA peak caller is run on each of the 5 generated datasets. For each simulation run, the F1 score is calculated during each annealing step and is saved in a list during the simulation. At the end of the simulation the F1 scores over time are saved as a TSV file. We end up with one TSV file for each nucleosome configuration.
    3. After the 5 runs, the mean F1 scores over time for each scenario are calculated, and saved in a final TSV file with three columns: mean F1 scores for regular nucleosomes, phase-shifted nucleosomes, and random nucleosomes (The other scenarios are not included here).

### Mean F1 Score Input
- The input file for `Mean_F1_different_coverages_plots.py` is created through multiple steps:
  1. The `data_simulation.py` script is run 1000 times, generating datasets for different scenarios each time.
  2. For each iteration the MaxFinder, Snyder, and SAUNA peak callers are run on these datasets, and the F1 scores for each peak caller and scenario are added to a TSV file.
  3. After 1000 iterations, the mean F1 scores for each peak caller and scenario are calculated.
  4. The resulting TSV file is filtered to retain only information for three specific scenarios: Regular, Phase-Shifted, and Random.
  5. This process is repeated for different coverages, and the TSV files are combined into one file, resulting in a final TSV file with 5 columns:
     - Column 1: Scenario (values should be 'Random', 'Regular', or 'Regular + Phase-Shifted')
     - Column 2: Coverage
     - Column 3: Mean F1 scores for the MaxFinder peak caller
     - Column 4: Mean F1 scores for the Snyder peak caller
     - Column 5: Mean F1 scores for the SAUNA peak caller

## Scripts Overview

### 1. `data_simulation.py`
This script creates 5 `tsv.gz` files with simulated data for 5 different nucleosome configurations. It does not take any input parameters. However, it contains a parameter called `coverage`, which can be adjusted to simulate fragment center data for different coverages.

### 2. `F1_scores_over_time_plots.py`
**Inputs**:
- `<temperature_file>`: A TSV file with two columns. The first column contains the temperatures, and the second column contains the annealing steps.
- `<F1_score_file>`: A TSV file with three columns. The first column contains the mean F1 scores for each annealing step for the scenario with only regular nucleosomes. The second column contains scores for the scenario with phase-shifted nucleosomes, and the third column contains scores for random nucleosomes.
- `<output_directory>`: The directory where the output plot will be saved.

**Usage**:
```shell
python F1_scores_over_time_plots.py <temperature_file> <F1_score_file> <output_directory>
```

**Output**: The script generates a plot showing mean F1 scores for each nucleosome configuration scenario over time and a plot showing the temperature over time. The number of annealing steps (rows) in both the temperature file and the F1 score file must be the same.

### 3. `Mean_F1_different_coverages_plots.py`
- `<input_file>`: A TSV file containing the mean F1 scores for different scenarios, coverages, and peak callers. It should have 5 columns:
  - Column 1: Scenario (values should be 'Random', 'Regular', or 'Regular + Phase-Shifted')
  - Column 2: Coverage
  - Column 3: Mean F1 scores for the MaxFinder peak caller
  - Column 4: Mean F1 scores for the Snyder peak caller
  - Column 5: Mean F1 scores for the SAUNA peak caller
- `<output_directory>`: The directory where the output plot will be saved.

**Usage**:
```shell
python Mean_F1_different_coverages_plots.py <input_file> <output_directory>
```

**Output**: The script creates a plot showing the mean F1 scores of different peak callers over different coverages, which is saved in the specified output directory.
