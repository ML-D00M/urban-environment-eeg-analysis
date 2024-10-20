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
- MATLAB (version X.X or higher)

For statistical analysis:
- R (version X.X.X or higher)
- Required R packages: rstatix, ggplot2, ggpubr, dplyr, tidyr

### Running the Analysis

1. Clone this repository
2. Clone the data (see below)
3. In all scripts, replace the <your-path-to-data> placeholder for your actual paths to data
4. Run the MATLAB scripts in numerical order (step1 through step7) in the `EEG_preprocessing_and_ISC_analysis_scripts/` folder
5. Run the R scripts in the `Statistical_analysis_scripts/` folder:
  - first, the `Statistical_analysis_All_&_Frequency_bands.R` script - on the resulted data from the previous step (EEG and ISC)
  - second, the `Statistical_analysis_Questionnaire.R` script on the questionnaire data (see below)

## Data

For access to the full dataset with fully anonymized EEG data and questionnaires data from our experiments, please use [this link](add_link_to_the_data).

## License

This project is licensed under the [MIT License] - see the [LICENSE.md](LICENSE.md) file for details.

## Citation

If you use this code or data in your research, please cite our paper:

[Citation information for our paper]

## Contact

[Name] - [Email]

Project Link: [https://github.com/ML-D00M/urban-environment-eeg-analysis](https://github.com/ML-D00M/urban-environment-eeg-analysis)

## Acknowledgments

- [Add acknowledgments]
- Thanks to all participants in our study
