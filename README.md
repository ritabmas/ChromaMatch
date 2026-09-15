# ChromaMatch
`ChromaMatch` is a Python tool for automatically identifying peaks in chromatograms using based on Retention-Time (RT) matching.\
The toll identified compounds using:
1. RT similarity against a reference chromatogram.
2. An optional compound library in CSV format for peaks that remain unindetified after chromatogram matching.\
The resulting identified peak tables are saved as new `.txt`files with the suffix `_identified`.

## Requirements
- Python 3.13
- pandas 3.0.5\
A dedicated Python environment is recommended.\
For example, uding _Anaconda Prompt_):
```
conda create --name NAME python=3.13
conda activate NAME
```
Then install the required Python packages:
```
pip install -r requirements.txt
```
## Repository structure
ChromaMatch/ 
│ 
├── peak_identification_base.py 
├── requirements.txt 
├── README.md 
├── LICENSE 
│ 
└── examples/ 
    ├── reference_chromatogram.txt 
    ├── example_input.txt 
    └── example_library.csv
The `examples/` folder contains input files that can be used to understand the required file formats and test the script.

## Required input files
All files must be in the same folder as `peak_identification_base.py`.
### Reference chromatogram
The reference chromatogram contains peaks that have already been assigned compound names and is used as the primary source for peak identification.\
The default name is `reference_chromatogram.txt`.\
The file must be tab-separated and contain at least the following columns: 
```
Peak #
Sample
Area S1
Name
1st Dimension Time (min)
```
An example reference chromatogram is provided in `examples/reference_chromatogram.txt`.
### Input files
Input chromatogram files must be tab-separated `.txt` files containing at least the same columns as the reference chromatogram:
```
Peak #
Sample
Area S1
Name
1st Dimension Time (min)
```
The reference chromatogram and files ending in `_identified.txt` are automatically excluded from processing.\
For example:
```
A1.txt
A2.txt
reference_chromatogram.txt
```
Only the appropriate inpt files will be processed.\
An example input file is provided in `examples/example_input.txt`.
### Optional compound library
An optional compound library can be provided as a `.csv` file.
The library must contain the following columns:
```
Name
Min_RT_min
Max_RT_min
```
For example:
```
Name,Min_RT_min,Max_RT_min
Compound A,10.00,10.40
Compound B,12.10,12.70
```
For each compound, the script calculates:
- Mean_RT_min: midpoint between minimum and maximum RT.
- Std_RT_min: half the interval width.
The compound library is applied only to peaks that remain unidentified after matching against the reference chromatogram.\
If no valid CSV library is found, this step is skipped automatically.\
An example compound library is provided in `examples/example_library.csv`

## Usage
Place the required input files in the same folder as `peak_identification_base.py`.\
Then run:
```
python peak_identification_base.py
```
The script can also be run directly from an IDE such as _Visual Studio Code_ (recommended).\
For example, if the input file is:
```
A1.txt
```
The script creates:
```
A1_identified.txt
```
The output file contains the processed peak information and the compound identifications obtained from the reference chromatogram and, when available, the optional CSV compound library.


## Processing steps
The script performs the following main steps:
- Loads and cleans the input files.
- Removes duplicate peak numbers.
- Matches unknown peaks against the reference chromatogram based on RT.
- Matches remaining unknown peaks against the optional CSV compound library.
- Removes peak areas below the configured area threshold.
- Removes early unidentified peaks according to the configured RT threshold.
- Removes configured unwanted compounds.
- Calculates the percentage of remaining unidentified peak area.
- Prints warnings for configured peak-area masks.
- Saves the processed data as `_identified.txt` files.

## Configuration
The main processing parameters can be modified at the top of `peak_identification_base.py`.
- REFERENCE_FILE: name of the reference chromatogram.
- OUTPUT_SUFFIX: suffix added to the processed input file names.
- RT_TOLERANCE: maximum RT difference allowed when matching peaks.
- LOW_AREA_THRESHOLD: area theshold used to remove low-area peaks.
- EARLY_RT_THRESHOLD: RT threshold used to remove early unidentifiable peaks.
- WARNING_MEDIUM_MIN_AREA, WARNING_MEDIUM_MAX_AREA, WARNING_HIGH_AREA: area threshold for warnings.
- COMPOUNDS_TO_REMOVE: list of compounds that should be excluded from the final results.
Changing these parameters allows the processing criteria to be adapted to different chromatographic datasets.

## Examples
The `examples/` folder provides files illustrating the expected input formats:
- `reference_chromatogram.txt`: example reference chromatogram.
- `example_input.txt`: example chromatogram to be identified.
- `example_library.csv`: optional compound library.

## Output
For each processed input file, the script creates a new file with the suffix:
```
_identified.txt
```
For example:
```
A1.txt
A1_identified.txt
```
The original input files are not overwritten.

## License
This project is licensed under the MIT License. See the `LICENSE` file for the full license text.
