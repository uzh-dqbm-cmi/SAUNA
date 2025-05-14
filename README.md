# Compliance update

This repository contains code (in the SAUNA folder) adapted from https://doi.org/10.1093/bioinformatics/btt404 by Schöpflin et al. 2013, and also code (in the Running_peak_callers/other_peak_callers) folder adapted from https://doi.org/10.1016/j.cell.2015.11.050 by Snyder et al. 2016. The documentation and licensing of the repository has been updated between 8th and 14th May 2025 to comply with the GPL-3.0 license of the original code and to attribute code to the original authors and describe changes. The original licenses can be found in the repository's main folder. We apologize for the oversight.
Do not use information you obtained from this repository prior to 14th May 2025!

# Repository Overview

This repository contains scripts and tools for analyzing nucleosome positions in cfDNA sequencing data, simulating nucleosome configurations, and performing downstream analyses.

## System Dependencies

In addition to the Python packages listed in `requirements.txt`, the following system dependencies are required to run the scripts:

- **GNU Parallel**: Install using your package manager. For example:
  - On Ubuntu:
    ```shell
    sudo apt-get install parallel
    ```
  - On macOS:
    ```shell
    brew install parallel
    ```
  
- **samtools**: A suite of programs for interacting with high-throughput sequencing data. Install using your package manager:
  - On Ubuntu:
    ```shell
    sudo apt-get install samtools
    ```
  - On macOS:
    ```shell
    brew install samtools
    ```
    
Make sure to install these dependencies before running the scripts.


## Additional Information

For more detailed information about specific scripts and analyses, please refer to the README files located in the different folders of the repository.
