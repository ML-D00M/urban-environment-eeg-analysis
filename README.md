# Green Urban Environments Enhance Brain-to-Brain Synchrony: An Electroencephalography (EEG) Study

This repository contains the analysis scripts for our study on neural synchronization in response to various urban environments, including both EEG analysis and statistical analysis of the results.

## Project Overview

Our research investigates how different urban environments (highways, boulevards, and parks) affect neural synchronization across individuals. We used EEG recordings and Inter-Subject Correlation (ISC) analysis to measure shared brain responses to urban environment videos, followed by statistical analysis to interpret the results.

## Repository Structure

- `EEG_preprocessing_and_ISC_analysis_scripts/`: Contains MATLAB scripts for EEG and ISC analysis
  - `step1_extractTriggers.m`: Extract triggers and align them with EEG recordings
  - `step2_segmentation.m`: Segment EEG data based on video stimuli
  - `step3_TemporalAlignment.m`: Align EEG data temporally
  - `step4_preprocessing.m`: Preprocess EEG data
  - `step5_isceeg_with_time_windows.m`: Perform ISC analysis
  - `step6_ISC_in_frequency_bands.m`: Analyze ISC in different frequency bands
  - `step7_topoplots_in_delta.m`: Generate topoplots for delta band
  - Additional utility functions used in scripts
- `Statistical_analysis_scripts/`: Contains R scripts for statistical analysis
  - `Statistical_analysis_All_&_Frequency_bands.R`: Perform ANOVA and post-hoc tests on ISC results. Generate plots for ISC results
  - `Statistical_analysis_Questionnaire.R`: Analyze questionnaire data

## Getting Started

### Prerequisites

For EEG analysis:
- MATLAB (version R2017b or higher)
- EEGLAB toolbox 

For statistical analysis:
- R (version R2017b or higher)
- Required R packages: rstatix, ggplot2, ggpubr, dplyr, tidyr

### Running the Analysis

1. Clone this repository
2. Download the data (see below)
3. In all scripts, replace the <your-path-to-data> placeholder for your actual paths to data
4. Run the MATLAB scripts in numerical order (step1 through step7) in the `EEG_preprocessing_and_ISC_analysis_scripts/` folder
5. Run the R scripts in the `Statistical_analysis_scripts/` folder:
  - first, the `Statistical_analysis_All_&_Frequency_bands.R` script - on the resulted data from the previous step (EEG and ISC)
  - second, the `Statistical_analysis_Questionnaire.R` script on the questionnaire data (see below)

## Data

For access to the full dataset with fully anonymized EEG data and questionnaires data from our experiments, please use [this link](https://www.dropbox.com/scl/fo/05yzofv5m4jze5iz3lcxx/APqAeE_lEwok9cgd8TEqcXw?rlkey=y1mksvdmld28ohrudfyytfbjt&st=a2z61i54&dl=0).

### EEG Data Files

Each participant's EEG data consists of three files:
- `PilotXX.eeg`: This is the binary file that contains the raw EEG data. It stores the continuous voltage values recorded from each electrode over time.
- `PilotXX.vhdr`: This is the header file in plain text format. It contains metadata about the EEG recording, such as: recording parameters (e.g., sampling rate, number of channels, data format), channel labels and names, reference electrode information, recording start time and date, and links to the associated `.eeg` and `.vmrk` files.
- `PilotXX.vmrk`: This is the marker file in plain text format. It contains event markers or triggers that were recorded during the EEG session. These markers indicate when certain events occurred, such as the onset of a stimulus (e.g., start of a video), participant responses or actions, technical events (e.g., recording start/stop).

These three files are interrelated and must be used together to properly interpret the EEG data. EEGLAB toolbox `pop_loadbv` function from the `step1_extractTriggers.m` script handles these file formats to load the data.

### Questionnaire Data

`Questionnaire_responses_anonymized.xlsx` contains participants' responses to questions about each urban environment (video) they watched during the experiment as well as to general questions about their preferences and experiences in urban environments.

Structure:
- Sheets: Each sheet (tab) in the Excel file represents a different set questions:
  - `General`: participants' responses to general questions about their preferences in urban environments.
  - `h_msk`, `b_msk`, `p_msk`, `h_spb`, `b_spb`, `p_spb`: participants' responses to questions about urban environments they watched (highways, boulevards, parks).
  - `Questions`: summary of responses (for all 6 environments) for each question about urban environments.
  - `Index`, `Indicies`: technical sheets used during analysis.
- Rows: Each row represents a participant.
- Columns: Each column represents participants' answers to a specific question on Likert scales (1 to 7) indicating agreement.

## Project Link

[https://github.com/ML-D00M/urban-environment-eeg-analysis](https://github.com/ML-D00M/urban-environment-eeg-analysis)

## Acknowledgments

- Thanks to all participants in our study
