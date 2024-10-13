# Green Urban Environments Enhance Brain-to-Brain Synchrony: An Electroencephalography (EEG) Study

This repository contains the analysis scripts for our study on neural synchronization in response to various urban environments, including both EEG analysis and statistical analysis of the results.

## Project Overview

Our research investigates how different urban environments (highways, boulevards, and parks) affect neural synchronization across individuals. We used EEG recordings and Inter-Subject Correlation (ISC) analysis to measure shared brain responses to urban environment videos, followed by statistical analysis to interpret the results.

## Repository Structure

- `eeg_analysis/`: Contains MATLAB scripts for EEG and ISC analysis
  - `step1_extractTriggers.m`: Extract triggers and align them with EEG recordings
  - `step2_segmentation.m`: Segment EEG data based on video stimuli
  - `step3_TemporalAlignment.m`: Align EEG data temporally
  - `step4_preprocessing.m`: Preprocess EEG data
  - `step5_isceeg_with_time_windows.m`: Perform ISC analysis
  - `step6_ISC_in_frequency_bands.m`: Analyze ISC in different frequency bands
  - `step7_topoplots_in_delta.m`: Generate topoplots for delta band
- `statistical_analysis/`: Contains R scripts for statistical analysis
  - `anova_and_post_hoc.R`: Perform ANOVA and post-hoc tests on ISC results
  - `plot_results.R`: Generate plots for ISC results
  - `questionnaire_analysis.R`: Analyze questionnaire data
- `data/`: Contains example data (or instructions for accessing full dataset)
- `results/`: Contains output from analysis scripts

## Getting Started

### Prerequisites

For EEG analysis:
- MATLAB (version X.X or higher)
- [Add other necessary MATLAB toolboxes]

For statistical analysis:
- R (version X.X.X or higher)
- Required R packages: rstatix, ggplot2, ggpubr, dplyr, tidyr

### Running the Analysis

1. Clone this repository
2. [Add Instructions for accessing or simulating data]
3. Run the MATLAB scripts in numerical order (step1 through step7) in the `eeg_analysis/` folder
4. Run the R scripts in the `statistical_analysis/` folder

## Data

Due to privacy concerns, the full dataset is not included in this repository. A sample dataset is provided in the `data/` folder for testing purposes. For access to the full dataset, please contact [your contact information].

## License

This project is licensed under the [MIT License] - see the [LICENSE.md](LICENSE.md) file for details.

## Citation

If you use this code or data in your research, please cite our paper:

[Citation information for our paper]

## Contact

[Name] - [Email]

Project Link: [https://github.com/yourusername/urban-environment-eeg-analysis](https://github.com/yourusername/urban-environment-eeg-analysis)

## Acknowledgments

- [Add acknowledgments you want to include]
- Thanks to all participants in our study
