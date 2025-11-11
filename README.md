# Brain Connectivity Analysis in MATLAB
Code to compute brain connectivity features from EEG data.

## Overview

### `conn.m`
Main script that computes connectivity features from patients EEG following the classical steps: pre-processing, nodes and edges computations, and features extraction.

This is where you need to specify (1) all parameters concerning your connectivity analyses:
- *f_bands*: list of frequency bands of interest
- *channels*: list of EEG channels of interest
- *ref*: original reference
- *fs*: sampling frequency
- *reref*: list of re-referencing schemes of interest
- *select_meas*: list of association measures of interest
- *window*: window  size (in seconds)
- *nb_features* = number of connectivity features to extract

(2) the specifics regarding your database:
- *files*: path to your files
- *save*: path to the folder in which you want to save your analyses
- *start*: analysis start time (in seconds)
- *duration*: length of analysis (in seconds)
- *patients*: number of patients you want to process
- *diff_pat*: list of the patients you wish to re-reference to the *ref* specified above and/or to resample to the *fs* specified above

The data of each patient should be stored in a MAT structure with 2 fields:
- channels: structure of size 1 x *number of channels* with 2 fields: Name and Fs
- data: double of size *length of recordings* [samples] x *number of channels* containing the EEG voltages (in microVolts)

### `data.m`
Class handling all processing on the raw data, including data preparation, pre-processing, and window splitting.

### `measures.m`
Class handling all association measures implementations (definition and parametrization).

### `associations.m`
Class handling the computations of the associations between nodes to generate the adjacency matrix.

### `network.m`
Class handling the extraction of network features from the adjacency matrix.

## Installation and dependencies
1. Ensure you have MATLAB installed (recommended version: R2024b or later).
2. Clone or download this repository.
3. Add the project folder to your MATLAB path.

**Dependencies**

MATLAB TOOLBOXES:
- Signal Processing Toolbox
- Parallel Computing Toolbox

OTHER TOOLBOXES:
- Brain Connectivity Toolbox. This is a freely available and open source Matlab toolbox for complex brain-network analysis. You can download it [here](brain-connectivity-toolbox.net). For more details, please refer to: *Complex network measures of brain connectivity: Uses and interpretations. Rubinov M, Sporns O (2010) NeuroImage 52:1059-69.*


## Usage
1. Open MATLAB and navigate to the project directory.
2. Update parameters and database information in the `conn.m` file.
3. Run `conn.m` to compute brain connectivity features for every file you specified.

## Contributors
- Lise Cottin

## License
This project is licensed under MIT License.
