%% General Comment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                                                %                                                                                                                           
% This script processes force plate, dynamometry, and heart rate data and analyzes the correlation between          %         
% performance metrics and heart rate variability. It resembles the entire process from raw data to results of the   % 
% investigation "Relationship between Heart Rate Variability and Neuromuscular Performance Responses following      % 
% Resistance Exercise". Detailed information on input file labeling, properties, and content as well as on          % 
% processing and calculations can be found in the supplement section of the reproduction documentation              %  
% (https://osf.io/43hnv/).                                                                                          %
%                                                                                                                   %
% The first section sets necessary options to run the main functions in the later section. When running the script, %
% a GUI opens to choose files that should be analyzed. If it was choosen to export the results, a second GUI opens  % 
% to select the output foulder. Export may take several minutes.                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User Input/ Necessary Options
% Specifying phases of the study design. Strings need to include specific information:
% NumberOfTestCycle_StudyPhase_NumberOfPhaseRepetition
% Note: Phases can be altered in name and/or number
study_design = ['1_BASE_1',"1_BASE_2","1_BASE_3","1_BASE_4","1_BASE_5","1_BASE_6","1_BASE_7",...                              
                "1_INT_1",...                                                       
                "1_REC_1","1_REC_2","1_REC_3"];                                                          

% Specifying names of calculated variables for HRV (1. line), CMJ (2. line), SQV (3. line) and DYNO (4. line) 
% Note: It is possible to change variable names, but DO NOT alter number of variables, unless specific function output is adjusted,
%       Names must be unique, units must be included and in brackets '[]',
%       Capital letters will be displayed as acronym in the graphs section
output_variables = {{'HR [1/min]','RR [ms]','SDNN [ms]','RMSSD [ms]','TP [ms�]','VLF [ms�]','LF [ms�]','HF [ms�]','LF/HF [arb.]','SD1 [ms]', 'SD2 [ms]','SD2/SD1 [arb.]','Artifacts [%]'},...
                    {'Jump Height Airtime [cm]','Jump Height Impulse [cm]','Net Impulse [Ns]','Peak Force [N]', 'rel. Peak Force [N/kg]','Peak Power [W]','TakeOff Velocity [m/s]','Peak Velocity [m/s]','Peak RFD [N/s]','Time to Peak RFD [ms]'},...
                    {'Mean Propulsive Velocity Squat [m/s]', 'Mean Velocity Squat [m/s]', 'Peak Velocity Squat [m/s]', 'Peak Power Squat [W]', 'Peak RFD Squat [N/s]'},...
                    {'Maximum Voluntary Torque [Nm]','Superimposed Twitch Torque [Nm]','Resting Twitch Torque [Nm]','Peak RTD Rest [Nm/s]','Voluntary Activation [%]'}};       
                
% Specifies the data subset, indicated by phase prefix, from which baseline values are calculated
% Note: Can be altered to match internal study logic
baseline    = {'BASE','1_INT_1'};                               

% HRV data processing options
artefact_recognition = 2;                               % 0 = no recognition, 1 = recognition with moving median +- artefactThreshold, 2 = custom recognition algorithm
artefact_threshold   = 250;                             % Threshold in ms 
detrending           = 1;                               % 0 = no detrending, 1 = detrending by removing 3rd order polynomial, 2 = zero phase butterworth filter cut-off frequency: 0.035 Hz

% Analysis options
includedPhases      = {'INT', 'REC'};                   % Phases from study_design that will be analysed. Number of measurements in these phases must be greater than 1

correlation_Type    = [   3  ,...                       % Options Reproducibility:  1 = Pearson correlation, 2 = Intraclass correlation, 3 = ANCOVA 
                          3      ];                     % Options Association:      1 = Pearson correlation, 3 = ANCOVA 
                      
    ICC_type        = '3-1';                            % Options: '1-1', '1-k', '2-1', '2-k', '3-1', '3-k'     
    
    ANCOVA_type     = {'parallel lines',...             % Options Reproducibility:  'parrallel lines' - similar to repeated measures correlation, 'separate lines' - individual fit
                       'parallel lines'};               % Options Association:      'parrallel lines' - similar to repeated measures correlation, 'separate lines' - individual fit
                   
NULL                = [   0.875  ,...                   % Reproducibility Null-Hypothesis , Range 0 to 1 (absolute) 
                          0.875     ];                  % Association Null-Hypothesis ,     Range 0 to 1 (absolute)
                      
mode                = {'inferiority',...                % Line one for reproducibility, Line two for association
                       'inferiority'};                  % Options: 'inferiority', tests if r is smaller than NULL   (H0 = r > NULL),          
                                                        %          'superiority', tests if r is greater than NULL   (H0 = r < NULL)    
                                                        %          'difference',  tests if r is different from NULL (H0 = r == NULL)
                                                        
alpha               = 0.1;                              % Confidence level, Type 1 Error 

HR_variable         = 4;                                % Selects the HR variable association is tested for. 
                                                        % Options refer to 'output_variables' index of the first row, e.g. 4 = RMSSD (ms)
% Export options
export_options      = [   1    ,...                     % 1 = exports importable-, analysis-, and individual data files, 0 otherwise
                          1    ,...                     % 1 = exports summary, statistic, and data tables, 0 otherwise
                          1    ,...                     % 1 = exports graphical illustration of statistics for each metric, 0 otherwise; Time consuming
                          1       ];                    % 1 = exports effect size plots for each individual, 0 otherwise
                      
% Plot options, refers to "output_variables"
% Maximum size for each array: 4
disp_vars           = {[4, 8    ],...                   % Determines the HRV variables to be plotted for each individual. E.g.: 4  = "RMSSD (ms)"
                       [15, 19, 22],...                 % Determines the CMJ variables to be plotted for each individual. E.g.: 15 = "Jump Height Impulse (cm)"
                       [24, 26, 27],...                 % Determines the SQV variables to be plotted for each individual. E.g.: 24 = "Mean Propulsive Velocity (m/s)"
                       [29, 31, 33]};                   % Determines the ITT variables to be plotted for each individual. E.g.: 29 = "Maximum Voluntary Torque (Nm)"
                   
xlabels = {'PRE';'POST24';'POST48';'POST72'};...        % Sets the xlabels for individual plots, must match number of phases that will be included due to 'includedPhases',
show_baseline = 0;                                      % Shows the baseline in individual and group mean effect size plot.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main functions
% Do not alter following commands, unless you know what you do.
clc
timerVal = tic;
addpath(string(pwd) + '\Functions');                    % Adds the folder \Functions to the working folder to access distinct functions

data                                                       = ImportData();
[individualData, groupData, analysisData]                  = DataProcessing(data, study_design, output_variables, artefact_recognition, artefact_threshold, detrending, baseline);
[reproducibilityResults, associationResults, corrMatrices] = DataAnalysis(analysisData, includedPhases, HR_variable, correlation_Type, ICC_type, ANCOVA_type, alpha, NULL, mode);
Export(export_options);

timerVal = toc(timerVal);
fprintf('\n###################################\n\nTotal Duration: \t%-2.2f seconds\n\n', timerVal);
clearvars -except groupData individualData analysisData reproducibilityResults associationResults corrMatrices data
