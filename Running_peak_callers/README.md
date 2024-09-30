# Running SAUNA and MaxFinder

This project provides scripts for running SAUNA and MaxFinder, tools designed to find nucleosome positions in cfDNA sequencing data.

## Table of Contents

- [Overview](#overview)
- [Usage](#usage)
- [Features](#features)

## Overview

To get started, download all files in the `Running_peak_callers` directory within this repository and place them in the same directory on your local machine. This will ensure that all necessary scripts and resources are available for execution.

The `run_everything` script is designed to analyze BAM files containing cfDNA sequencing data in the specified directory. These BAM files must be sorted by coordinates and should all be located in one folder, which is the `bam_directory`. The script creates a folder for each BAM file with the following subdirectories:
- **input**: Contains outputs from the MaxFinder peak caller (TSV files) and TSV.GZ files with fragment center data.
- **output**: Contains outputs from the SAUNA peak caller (BED files).
- **archive**: Used temporarily during script execution but should be empty upon completion.
- **other**: Contains unused intermediate files that are generated during execution.

## Usage

The main script for running the analysis is `run_everything`. This script orchestrates the execution of all other scripts in the 'Running_peak_callers' folder. 

To use the script, run the following command in your terminal:
./run_everything <bam_directory> <params.txt>


- `<bam_directory>`: The path to the directory containing all BAM files to be analyzed.
- `<params.txt>`: A text file containing parameters that can be adjusted for the analysis.

## Features

- **Nucleosome Position Detection**: Efficiently identifies nucleosome positions in cfDNA sequencing data using the SAUNA and MaxFinder algorithms.
- **Automated Workflow**: The `run_everything` script streamlines the process by automatically organizing and managing input and output files for each BAM file, reducing the need for manual file handling.
- **Customizable Parameters**: Users can easily adjust analysis parameters in the `params.txt` file, allowing for flexibility in the analysis based on specific experimental conditions.


