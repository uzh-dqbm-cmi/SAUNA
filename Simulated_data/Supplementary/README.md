# Supplementary Figures

This folder contains the `Supplementary_figures.py` script, which generates supplementary figures.


## Inputs

The script accepts the following inputs:

1. **`different_iterations_file`**:  
   A file with 6 columns:
   - **Column 1**: Nucleosome configuration information.
   - **Column 2**: F1 scores of the SAUNA peak caller.
   - **Column 3**: F1 scores of the SAUNA3 (SAUNA with no input) peak caller with the regular number of iterations.
   - **Column 4**: F1 scores of SAUNA3 with 10x the regular number of iterations.
   - **Column 5**: F1 scores of SAUNA3 with 100x the regular number of iterations.
   - **Column 6**: F1 scores of SAUNA3 with 1000x the regular number of iterations.

2. **`different_peak_callers_file`**:  
   A file with 6 columns:
   - **Column 1**: Nucleosome configuration information.
   - **Column 2**: F1 scores of MaxFinder.
   - **Column 3**: F1 scores of the Snyder peak caller.
   - **Column 4**: F1 scores of SAUNA.
   - **Column 5**: F1 scores of SAUNA2 (SAUNA with Snyder input).
   - **Column 6**: F1 scores of SAUNA3.

3. **`F1_score_file`**:  
   A file with 6 columns:
   - **Column 1**: Annealing steps.
   - **Column 2**: Mean F1 score for the annealing step for regularly spaced nucleosomes.
   - **Column 3**: Mean F1 score for regularly spaced nucleosomes with phase-shifted nucleosomes.
   - **Column 4**: Mean F1 score for random nucleosomes.
   - **Column 5**: Mean F1 score for regularly spaced nucleosomes with random ones.
   - **Column 6**: Mean F1 score for regularly spaced nucleosomes with random where some of the regular nucleosomes were removed.

4. **`output_dir`**:  
   The path to the output directory where the generated figures will be saved.

## Usage

To run the script, use the following command:

```shell
python Supplementary_figures.py <different_iterations_file> <different_peak_callers_file> <F1_score_file> <output_dir>
```
