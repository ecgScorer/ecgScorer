# ecgScorer: An open source MATLAB Toolbox for ECG Signal Quality Assessment

ecgScorer is a MATLAB toolbox that provides tools for assessing the quality of electrocardiogram (ECG) signals and extracting signal quality indices (SQIs). It is designed to be modular, flexible, and adaptable to various research workflows. It can operate through a standalone graphical interface or using provided functions.

## Features

- Load ECG signal databases from various file formats
- Compute quality scores based on a novel robust SQA method [1]
- Export analysis results with generated reports
- Extract up to 37 SQIs 

## Installation
1. Prerequisites:

MATLAB (tested with R2022b) with signal processing toolbox,  Statistics and Machine Learning Toolbox, and Parallel Computing Toolbox.

2. install

To install ecgScorer, download or clone this repository and add it to your MATLAB path. You can do this by running the following commands in MATLAB:

```matlab
% Download the repository
websave('ecgScorer.zip', 'https://github.com/ecgScorer/ecgScorer/archive/refs/heads/main.zip');
% Unzip the repository
unzip('ecgScorer.zip');
% Add the repository to the MATLAB path
addpath(genpath('ecgScorer-main'));
```

Alternatively, you can manually download the repository from [here](https://github.com/ecgScorer/ecgScorer/archive/refs/heads/main.zip) and unzip it to a folder of your choice. Then, add the folder and its subfolders to your MATLAB path (or simply use `instal_ecgScorer.m` function)

## Usage

To use ecgScorer, you can either launch the graphical interface or use the functions.

### Graphical Interface

To launch the graphical interface, run the following command in MATLAB:

```matlab
ecgScorer
```

This will open a window with two tabs: ECG Quality Assessment and Signal Quality Indices Extraction.

#### ECG Quality Assessment

In this tab, you can load an ECG signal database, input the sampling frequency, and compute the quality scores. the result is automatically saved to the `base` workspace, but you can choose to export them to a file.

#### SQIs Extraction
In this tab, you can load an ECG signal database, input the sampling frequency, and compute the up to 37 SQIs. the result is also saved to the `base` workspace or export to a file.

## Functions
#### ECG Quality Assessment

use the function `scorer12.m`

#### SQIs Extraction
The following functions are used :
- `statisticSQI.m`   for statistics-based SQIs extraction.
- `frequencySQI.m`   for frequency domain based SQIs extraction.
- `nonLinearSQI.m`   for non-linear tools based SQIs extraction.
- `qrsDetectorSQI.m` for QRS detectors based SQIs extraction.

These functions are localized in the **compute_SQIs** folder.
The **aScorer** folder contains all files for automatically perform the ECG Quality assessment : 
- scorer12.m
- simSQI.m 
- ecgScorer.mlapp

The **SQIs Functions** folder hold :
- avecorr.m
- encodingLZC.m
- FuzzyEn2.m
- fuzzyMEn.m
- hurst.m
- lz_complexity
- sampEn.m
- t_SQI.m

The **Example** contains the files to test the software :
- the dataset : _signal_test_ in *.txt* or *.mat* extension
- the matlab live file : `test_ecgScorer_toolbox.mlx`

Please note the following points:
* All algorithms must be used with ECGs as standing vectors or matrices with leads columnwise arranged (temporal dimension in lines)

* We publish the software as it is and do not guarantee proper performance. Nevertheless, we highly acknowledge feedback. Use the issues functionality in github.

## Reference
[1] F. Kuetche, N. Alexendre, N. E. Pascal, and S. Thierry, ‘Simple, efficient, and generalized ECG signal quality assessment method for telemedicine applications’, Inform. Med. Unlocked, vol. 42, p. 101375, 2023, doi: [10.1016/j.imu.2023.101375](doi.org/10.1016/j.imu.2023.101375)
