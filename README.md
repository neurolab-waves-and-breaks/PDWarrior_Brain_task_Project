# BRAIN Finger-Tapping Task  
## Data Collection and Analysis Pipeline

This repository contains the complete data collection and analysis pipeline
for the Bradykinesia Akinesia Incoordination (BRAIN) finger-tapping task.
The code was developed to support behavioural analyses reported in the associated study.

The workflow consists of three main stages:
1. Task execution and data collection
2. Data preprocessing and metric computation
3. Statistical analysis and visualisation

---

## Repository structure

brain-task/
├── brain-task/ # Pygame task for data collection
├── preprocessing_code/ # MATLAB preprocessing and metric computation
├── stats_analysis_code/ # Python statistical analysis and plotting
└── README.md


---

## 1. Task execution (`brain-task/`)

This folder contains the **Pygame implementation of the BRAIN finger-tapping task**.

### Purpose
- Run the finger-tapping task with auditory start (`Go.wav`) and end (`Trumpet.wav`) cues
- Record key press events over a predefined task duration
- Export raw keystroke timing data for each recording

### Output
Each run produces a `.csv` file containing:
- `strokeOnset` – timestamp of key press onset  
- `dwellTime` – key press duration  
- `letters` – key identity  

Files are named using the following convention:

<date><participantID><pre/post><hand>.csv


Example:
14_01_2026_P5_preExerc_right.csv


### Requirements
- Python (≥ 3.8 recommended)
- pygame
- pandas
- numpy

---

## 2. Preprocessing and metric computation (`preprocessing_code/`)

This folder contains **MATLAB scripts for preprocessing the raw task output
and computing behavioural metrics**.

### Purpose
- Load raw keystroke `.csv` files
- Reconstruct data at the participant level
- Compute behavioural metrics (e.g. velocity, dwell time, travel time)
- Generate pre/post summary tables
- Merge behavioural metrics with demographic information

### Output
- Processed `.mat` files
- Summary `.csv` tables
- Combined behavioural–demographic datasets

These outputs are used as inputs for the statistical analyses.

### Requirements
- MATLAB (R2021a or newer recommended)

---

## 3. Statistical analysis and visualisation (`stats_analysis_code/`)

This folder contains **Python scripts for statistical testing and figure generation**.

### Purpose
- Load preprocessed datasets generated in MATLAB
- Perform statistical analyses (e.g. pre- vs post-condition comparisons)
- Generate publication-ready figures and summary statistics

### Output
- Statistical result tables
- Figures for manuscripts or reports

### Requirements
- Python (≥ 3.9 recommended)
- numpy
- pandas
- scipy
- matplotlib
- seaborn

---

## Typical workflow

1. Run the finger-tapping task (`brain-task/`) to collect raw data
2. Preprocess data and compute behavioural metrics using MATLAB (`preprocessing_code/`)
3. Perform statistical analyses and generate figures in Python (`stats_analysis_code/`)

---

## Notes

- Folder names and file naming conventions should not be changed, as later
  stages rely on a consistent directory structure.
- Each stage of the pipeline can be run independently once the required inputs exist.
- Raw participant data are not included in this repository.
- MATLAB scripts require a licensed version of MATLAB.

---

## Contact

Rui  
University of Bristol

---

## License

This project is licensed under the MIT License.