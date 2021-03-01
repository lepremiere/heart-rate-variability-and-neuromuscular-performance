# Heart-rate-variability-and-neuromuscular-performance

This repository provides the scripts used in "Heart Rate Variability and Neuromuscular Performance Responses following Resistance Exercise" (Huick, 2021). The code is generalized for further usage. It processes force plate, dynamometry and heart rate data to calculate heart rate variability and neuromuscular performance metrics. Together with the settings in "Only_Script_You_Need_To_Care_About.m" and file labels that contain specific information blocks, the script sorts the variables following study logic and present subjects. Information about the correct labeling as well as the codebook for calculated metrics can be found in the supplements section of the reproduction documentation in the project repository (https://osf.io/43hnv/). Heart rate variability and neuromuscular performance markers are then correlated and depending on export settings, data, figures, and tables finally saved to an output folder.

Prerequisites:
-   MATLAB R2020b or later
    -   Signal Processing Toolbox
    -   Statistics and Machine Learning Toolbox
    -   Parrallel Computing Toolbox

Installation:
-   No installation required. Just download all files and place them in your working directory.
    
    **Note**: The main scripts will call upon "Functions" folder. This folder needs to be within the same folder as main scripts.

Usage:
-   Open [Only_Script_You_Need_To_Care_About.m](https://github.com/lepremiere/heart-rate-variability-and-neuromuscular-performance/blob/main/Only_Script_You_Need_To_Care_About.m)  The main script consists of two sections, the first inherits settings that are parsed into the main functions
-   Manipulate settings to match study design and methods      
-   Execute the script
    - A GUI opens. Select the input files
    - Files are imported by main function [`ImportData()`](https://github.com/lepremiere/heart-rate-variability-and-neuromuscular-performance/blob/main/Functions/ImportData.m)
    - The imported data is then processed by [`DataProcessing()`](https://github.com/lepremiere/heart-rate-variability-and-neuromuscular-performance/blob/main/Functions/DataProcessing.m)
    - The processesed data is then analyzed by [`DataAnalysis()`](https://github.com/lepremiere/heart-rate-variability-and-neuromuscular-performance/blob/main/Functions/DataAnalysis.m)
    - After selecting the output folder via GUI, the results finally exported by [`Export()`](https://github.com/lepremiere/heart-rate-variability-and-neuromuscular-performance/blob/main/Functions/Export.m)

