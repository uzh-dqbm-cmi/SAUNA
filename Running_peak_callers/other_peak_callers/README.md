# Using Other Peak Callers: SAUNA_no_input and Snyder

This folder contains two peak callers designed for nucleosome position detection in cfDNA sequencing data.

## Table of Contents

- [Overview](#overview)
- [Usage](#usage)
  - [SAUNA_no_input](#sauna_no_input)
  - [Snyder](#snyder)
- [File Requirements](#file-requirements)

## Overview

This folder includes two distinct peak callers:

1. **SAUNA_no_input**: This peak caller executes nucleosome position detection using the SAUNA algorithm without requiring input from another peak caller.
2. **Snyder**: This peak caller, designed by Snyder et al., is used to analyze nucleosome positions based on specified input files.

## Usage

### SAUNA_no_input

To use the SAUNA peak caller, you can run the following command:
```shell
python SAUNA_no_input.py <nucleosome_center_data.tsv.gz> <params.txt> [output-path]
```

- `<nucleosome_center_data.tsv.gz>`: This is a tsv.gz input file that contains nucleosome center data. It should have the smoothed fragment center data in the 4th column, and all data should be from the same chromosome.
- `<params.txt>`: A parameter file specifying the necessary settings for the peak caller.
- `[output-path]`: This is an optional argument where you can specify an alternative output directory.

### Snyder

To use the Snyder peak caller, execute the script with the following command:
```shell
python Snyder.py <nucleosome_center_data.tsv.gz> <output_directory>
```

- `<nucleosome_center_data.tsv.gz>`: This should be the path to the tsv.gz file, which should have the same format as described above for the SAUNA peak caller.
- `<output_directory>`: This is the directory where the output will be stored.

## File Requirements

Both peak callers require a tsv.gz file where the 4th column contains the smoothed fragment center data. The data should be consistent, coming from the same chromosome for proper analysis.

