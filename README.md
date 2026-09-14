# gc_peak_identification_tool
This tool automatically identifies peaks in chromatograms using:
1. Retention-time similarity against reference_chromatogram.txt.
2. An optional compound library in CSV format for peaks not identified in the reference chromatogram.

## Requirements
It is advised to create a virtual environment using Anaconda Prompt: conda create --name NAME python=3.13\
Once created, install the libraries in the requirements.txt file.

## Required input files
All files must be in the same folder as peak_identification_base.py.
### Reference chromatogram
The default name is reference_chromatogram.txt.\
It must be tab separated and contain: Peak #, Sample, Area S1, Name, 1st Dimension Time (min).
### Input files
Input files must be tab-separated .txt files with the same columns.\
The reference file and files ending in _identified.txt are excluded automatically.
### Optional CSV library
A CSV file is optional. It must contain: Name,Min_RT_min,Max_RT_min.\
Example:\
Name,Min_RT_min,Max_RT_min\
Compound A,10.00,10.40\
Compound B,12.10,12.70\
The script calculates:\
Mean_RT_min: midpoint between minimum and maximum RT.\
Std_RT_min: half the interval width.\
The CSV library is applied only to peaks still unidentified after reference chromatogram matching. If no valid CSV file is found, this step is skipped automatically.

## Usage
Place the files in the script folder and run: peak_identification_base.py (_e.g._, using VS code).\
The script creates files with the suffix: _identified.txt.\
For example:\
A1.txt\
A1_identified.txt

## Processing steps
The script:
- Loads and cleans the input files.
- Removes duplicate peak numbers.
- Matches unknown peaks against the reference chromatogram.
- Matches remaining unknown peaks against the optional CSV library.
- Removes low-area peaks.
- Removes early unidentified peaks.
- Removes configured unwanted compounds.
- Reports the percentage of remaining unidentified peak area.
- Prints warnings for peak area masks.

## Configuration
The main settings can be changed at the top of the script, including:\
REFERENCE_FILE\
OUTPUT_FOLDER\
FILE_PATTERN\
RT_TOLERANCE\
LOW_AREA_THRESHOLD\
EARLY_RT_THRESHOLD\
COMPOUNDS_TO_REMOVE
