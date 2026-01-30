# Behavioral Analysis – BRAIN Finger Tapping Task

This folder contains **MATLAB code for analysing behavioural metrics**
derived from the BRAIN finger-tapping task.

The scripts implement the full behavioural analysis pipeline used in the
associated study, including metric computation, data restructuring, and
pre/post comparisons.

All analyses can be reproduced by running a single MATLAB script.
---

## Quick start (recommended)

1. Open MATLAB
2. Set the current folder to this directory
3. Run:```matlab
run_analysis

The script will execute the full analysis pipeline automatically.

---

##  What this script does

The run_analysis.m script performs the following steps:

This script contains all functions and analysis steps internally and:

1. Loads the required data files

2. Computes behavioral metrics (e.g. velocity, dwell time, travel time)

3. Merges behavioral data with demographic information

4. Reconstruct data structure 

5. Performs pre/post comparisons

6. Saves all summary tables and results in .csv 

Key analysis parameters can be modified within the script if required.


---


## Main script

run_analysis.m
Entry point for the entire behavioural analysis pipeline.
All functions and analysis steps are called internally from this script.


## Input data

The analysis scripts require locally provided behavioural and demographic
input data generated during earlier stages of the analysis pipeline.
These input files are not included in this repository.

Intermediate results are generated and stored internally as MATLAB `.mat`
files during the analysis but are not included in this repository.



> Note:  Do not rename or move files, as the scripts rely on the current file
names and directory structure.

Note: Raw participant-level keystroke data generated during task execution
are not included in this repository.
---


## Output

Running run_analysis.m will generate:

- Updated behavioural summary tables

- Combined behavioural–demographic datasets

- MATLAB .mat files containing intermediate and final results (generated locally)

All outputs are saved directly into this folder.

---

## Requirements

- MATLAB (R2021a or newer recommended)

- No additional toolboxes required

---

## Notes

- The analysis assumes that the folder structure and filenames remain unchanged.

- If errors occur, ensure MATLAB’s current working directory is set to this folder.
 