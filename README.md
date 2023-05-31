# GHOST-MAGNET
## Excuse the mess
Untill the completion of the project the repository will be messy with undocumented files and folders. A thorough clean up of files and structures will be commenced towards the end of the project.



![workflow_final2](https://github.com/ndreey/ghost-magnet/assets/72344166/764b5951-d6f6-4fd3-b9e5-f56642bf4cf0)


# Evaluation of Host Contamination Removal Post Sequencing

This repository contains the code and data for my BSc thesis project in Molecular Bioinformatics at the University of Skövde. The project evaluated the removal of host contamination post sequencing of non-model organisms.

## Project Aim

The aim of this project was two-fold:

1. To determine how different levels of host contamination affects assembly and genome binning, and whether a significant correlation exists.
2. To evaluate if removal of reads mapping to host bins is a valid method for host contamination depletion post sequencing.


![image](https://github.com/ndreey/ghost-magnet/assets/72344166/234b0e6f-b582-441b-99f0-809a2b30d2e7)


## Study Design

This project is an exploratory study that aims to determine if the method of host contamination removal has validity and should be explored further to generate a package/program that classifies and removes host bins, and re-runs metagenomic analysis.

## Repository Contents

This repository contains the following files:

- `README.md`: This file, which provides an overview of the project.
- `code/`: This directory contains the code used in the project.
- `data/`: This directory contains the data used in the project.

### Subdirectories

Within the `code/` directory, the following subdirectories serve specific functions:

- `scripts/`: This directory contains all the necessary scripts and programs for processing and cleaning the data.

Within the `data/` directory, the following subdirectories are further categorized:

- `raw/`: This directory contains the raw data that is utilized throughout the project.
- `processed/`: This directory contains the processed data generated during the project.
- `references/`: This directory contains reference data used for analysis.
- `training/`: This directory contains training data used for analysis.

Additionally, the `exploratory/` subdirectory is designated for storing files, notes, papers, and other various files that can be of interest to the project. In contrast, files determined relevant to the final report are placed in the `submission/` subdirectory.
 
## Using GitHub for Scientific Data Analysis

The development of software is a complex process that requires managing various versions of code, tracking progress, and ensuring that the final product works as intended. These criteria are similarly applicable and helpful for scientific data analysis, where reproducibility is of great importance. GitHub is an online platform that allows developers to create repositories for storing and managing code changes. 

The platform enables changes made to the code to be tracked, with the possibility to revert to previous versions if necessary. It also provides the potential for managing different versions of the code for different environments. With the intention of ensuring reproducibility and documentation, the main idea was to create a GitHub repository for each step of the data analysis. For the mock data simulation, “01_clean_mock_data” was created. 

With the GitHub repository, the project will also align with the FAIR principles described by (Wilkinson et al., 2016), which aim to make data and code findable, accessible, interoperable, and reusable. Further, when working with high-performance computing (HPC), it is essential to carefully optimize and test the programs to ensure a clean run on the server. By utilizing the local Git repository, thorough tests of scripts, programs, and program configurations can be determined and adjusted for using a subset of the data before the highly resourceful run on the server. The use of a local Git repository also facilitates the uploading of the correct file structure for the proper execution of programs without significant modifications to paths. 

Finally, with the use of the issue feature on GitHub, one can track, report and resolve problems or tasks related to the project. The feature is helpful for collaborative development, where multiple contributors are involved in different aspects. Nonetheless



