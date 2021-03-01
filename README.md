# Heart-rate-variability-and-neuromuscular-performance

This repository provides the scripts used in "Heart Rate Variability and Neuromuscular Performance Responses following Resistance Exercise" (Huick, 2021). The code is generalized for further usage. It processes force plate, dynamometry and heart rate data to calculate heart rate variability and neuromuscular performance metrics. Together with the settings in [Only_Script_You_Need_To_Care_About.m](https://github.com/lepremiere/heart-rate-variability-and-neuromuscular-performance/blob/main/Only_Script_You_Need_To_Care_About.m) and file labels that contain specific information blocks, the script sorts the variables following study logic and present subjects. Information about the correct labeling as well as the codebook for calculated metrics can be found in the supplements section of the reproduction documentation in the project repository (https://osf.io/43hnv/). Heart rate variability and neuromuscular performance markers are then correlated and depending on export settings, data, figures, and tables finally saved to an output folder.

Prerequisites:
-   MATLAB R2020b or later
    -   Signal Processing Toolbox
    -   Statistics and Machine Learning Toolbox
    -   Parrallel Computing Toolbox

Installation:
-   No installation required. Just download all files and place them in your working directory.
    
    **Note**: The main scripts will call upon "Functions" folder. This folder needs to be within the same folder as main scripts.

Usage:
-   Open [Only_Script_You_Need_To_Care_About.m](https://github.com/lepremiere/heart-rate-variability-and-neuromuscular-performance/blob/main/Only_Script_You_Need_To_Care_About.m). The main script consists of two sections, the first inherits settings that are parsed into the main functions in the second section.
-   Execute the script
    - A GUI opens to select the input files (e.g., [data of the study]()
    - Files are imported by main function [`ImportData()`](https://github.com/lepremiere/heart-rate-variability-and-neuromuscular-performance/blob/main/Functions/ImportData.m)
    - The imported data is then processed by [`DataProcessing()`](https://github.com/lepremiere/heart-rate-variability-and-neuromuscular-performance/blob/main/Functions/DataProcessing.m)
    - The processesed data is then analyzed by [`DataAnalysis()`](https://github.com/lepremiere/heart-rate-variability-and-neuromuscular-performance/blob/main/Functions/DataAnalysis.m)
    - After selecting the output folder via GUI, the results finally exported by [`Export()`](https://github.com/lepremiere/heart-rate-variability-and-neuromuscular-performance/blob/main/Functions/Export.m)

Settings:
-   Study logic
    - `study_design` is an array of character vectors. The input file labels are compared to this array to identify files of interest and to sort them following the study design. 
    - `output_variables` is an cell array with 4 cell array in it. These cells hold character vectors with the variable names for every variable in every test.
    - `baseline` is a character vector that contains the phases which should be considered as baseline to normalize variables to it.
-   Heart rate variability 
    - `artefact_recognition` is a scalar. If it is 0, heart rate data will not be searched for artefacts, if 1 is selected, artefacts are identified by a   moving median ± `artefact_threshold`. With the option 3, a [custom recognition algorithm](https://osf.io/z78jx/) algorithmn is applied. All artefacts are replaced by cubic spline interpolation at 10 Hz. 
    - `detrending` is a scalar. If it is 0, heart rate data is not detrended. If it is 1, a 3rd order polynomial is removed from the data. If it is 2, the data is filtered by a zero-phase low-pass Butterworth filter with an cut-off frequency of 0.035 Hz.    -
-   Analysis options
    - `includedPhases` is a cell array containing character vectors. These include the study phases that indicate the intervention and recovery or simply the phases that are later correlated.
    - `correlation_Type` is a 2 x 1 integer vector between 1 and 3. It indicates which correlation tpye should be used for analysis of the reproducibility (first) and the association (second). 1 = Pearson, 2 = Intraclass correlation , 3 = ANCOVA with the options of  "parallel lines" ([repeated measures correlation](https://doi.org/10.3389/fpsyg.2017.00456)), or "separate lines" between subjects specified in `ANCOVA_type` in the first row 
    - `ICC_type` is a character vector that determines the method of ICC that is used. Options: '1-1', '1-k', '2-1', '2-k', '3-1', '3-k'
    - `ANCOVA_type` is a cell array that contains 2 character vectors. They determine the type of ANCOVA that is used for the analysis of the reproducibility (first) and the association (second). Options: 'parallel lines' ([repeated measures correlation](https://doi.org/10.3389/fpsyg.2017.00456)) or 'separate lines'.
    - `NULL` is a 2 x 1 array. It holds 2 numbers from 0 to 1 that represent the correlation coefficient the results should be tested against. First number to test reproducibility results against, second to test association results.
    - `mode` is a cell array with 2 character vectors. These vectors determine if it should be tested if the results are lower than `NULL` ('inferiority'), greater than `NULL` ('superiority'), or just different from it ('difference').
    - `alpha` is a scalar that determines the confidence level.
    - `HR_variable` is a scalar that determines the single heart rate variable that the performance metrics should be compared to. The number refers to the position of the desired variable in `output_variables` (e.g., 4 = RMSSD).
- Export options
    - `export_options` is a 4 x 1 logical vector. Each row enables a certain part of the output. First row: enables the export of importable-, analysis- and individual data in MATLAB formt (e.g., for re-import). Second row: enables the export of summary, statistics and data tables in Excel format. Third row: enables the export of graphical illustrations of the statistics for each correlated metric (time consuming, 5 min.). Fourth row: enables the export of individual and group effect size plots.
    - `disp_vars` is a cell array that contains 4  integer vectors of maximal length 4. Each vector determines the varibles that should be shown on effect size plots. These integers refer to the `output_variables` vectors. The rows of `disp_vars` match the rows of `output_variables`.
    - `x_labels` is a cell array containing n character vectors. These vectors are used as labels for effect size plots and must therefore equal in number compared to the study phases that are analyzed by the choice of `includedPhases`.
    - `show_baseline` is a logical that determines if the baseline will also be shown on effect size plots.