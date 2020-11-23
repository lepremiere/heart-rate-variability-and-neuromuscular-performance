%% General Comment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                            
%      
%
%
%
%
%
%
%
% Requires external function 'uipickfiles.m' (Copyright (c) 2007, Douglas M. Schwarz)                               % 
% https://www.mathworks.com/matlabcentral/fileexchange/10867-uipickfiles-uigetfile-on-steroids                      %                                                                                                                           
%                                                                                                                   %         
% More than 2 files must be selected.                                                                               % 
% If CHECK functions have been used before and either 'HRVdata, CMJdata, ITTdata, SQVdata' exists in workspace      % 
% this data will be processed and input prompt surpressed.                                                          %
%                                                                                                                   %     
% Dynamometry data must be a struct including fields '.Torque' with data stored in '.Torque.values'                 % 
% and sampling frequency in '.Torque.interval' and '.Keyboard' with time of press in '.Keyboard.times'              %  
% Every keypress must indicate a stimulation. Superimposed stimulations need to occur first.                        %  
% Two stimulations per repetition required (superimposed, resting). Single files need to include all                % 
% stimulations of the test session. For further detail compare example file ('ITT_1_EXAMPLE_1.mat')                 % 
%                                                                                                                   % 
% Counter-Movement data must be ASCII table with vertical ground reaction forces stored in variables '.z1, .z2,     %  
% .z3, .z4' and time in the last column. For further detail compare example file ('CMJ_1_EXAMPLE_1_1.csv')          % 
% Note: German number format ('1,4') has to be changed to international notation ('1.4'). Use 'convertForcePlate.m' %  
%                                                                                                                   % 
% Heart rate variability data must be RR-intervals stored as n-by-1 vector in a textfile (.txt).                    %  
% RR-intervals must be in 'ms'. For further detail compare example file ('HRV_1_EXAMPLE_1.csv')                     % 
%                                                                                                                   % 
% Files have to be named in following pattern (underscore as seperator):                                            % 
% Subject_NumberOfTestCycles_StudyPhase_NumberOfPhaseRepetitions_(NumberOfTrial(only CMJ and SQV))                  % 
%                                                                                                                   % 
% Example:  XY_1_BASELINE_1_1, XY_1_BASELINE_1_2,..., XY_1_INTERVENTION_1_1,..., XY_2_RECOVERY_2_6 for CMJ and SQV  % 
%           XY_1_BASELINE_1,   XY_1_BASELINE_1,...,   XY_1_INTERVENTION_1,...,   XY_2_RECOVERY_2   for ITT and HRV  %
%                                                                                                                   % 
% Information structure must fall within following logic:                                                           % 
% - Pieces of information must be seperated by underscores ('_')                                                    %     
% - Number of phases and trials must be less or equal to 9                                                          % 
% - Phase names must fall in alphabetical order (e.g. Baseline, Intervention, Recovery)                             % 
% - Case sensitive                                                                                                  % 
%                                                                                                                   %
% External function 'uipickfiles.m' used (Copyright (c) 2007, Douglas M. Schwarz)                                   % 
% External function 'corrdist.m' used (Copyright (c) 2012, Xu Cui)                                                  % 
%                                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User Input
% Specifying phases of the study design. Strings need to include specific information:
% NumberOfTestCycles_StudyPhase_NumberOfPhaseRepetitions
% Note: Phases can be altered in name and/or number
study_design = ['1_BASE_1',"1_BASE_2","1_BASE_3","1_BASE_4","1_BASE_5",...                              
                "1_INT_1",...                                                       
                "1_REC_1","1_REC_2","1_REC_3",...                                                              
                "1_WASH_1","1_WASH_2","1_WASH_3",...        
                "2_INT_1",...
                "2_REC_1","2_REC_2","2_REC_3"];
            
% Specifies the data subset, indicated by phase prefix, from which baseline values are calculated
% Note: Can be altered to match internal study logic
baseline    = 'BASE';                               

% Specifying names of calculated variables for HRV (1. line), CMJ (2. line), SQV (3. line) and DYNO (4. line) 
% Note: It is possible to change variable names, but DO NOT alter number of variables, unless specific function output is adjusted,
%       Names must be unique, 
%       Capital letters will be displayed as acronym in the graphs section
output_variables = {{'HR (1/min)','RR (ms)','SDNN (ms)','RMSSD (ms)','TP (ms²)','VLF (ms²)','LF (ms²)','HF (ms²)','LF/HF (arb.)','SD1 (ms)', 'SD2 (ms)','SD2/SD1 (arb.)','Artifacts (%)'},...
                    {'Body Weight (kg)','Jump Height Airtime (cm)','Jump Height Impulse (cm)','Net Impulse (Ns)','Peak Force (N)', 'rel. Peak Force (N/kg)','Peak Power (W)','TakeOff Velocity (m/s)','Peak Velocity(m/s)','Peak RFD (N/s)','Time to Peak RFD (ms)'},...
                    {'Total Weight (kg)','Mean Propulsive Velocity (m/s)','Intensity MPV (%1RM)','Mean Velocity (m/s)','Intensity MV (%1RM)', 'Peak Velocity (m/s)','Intensity PV (%1RM)','Peak Power Squat (W)'},...
                    {'Maximum Voluntary Contraction (Nm)','Superimposed Twitch (Nm)','Resting Twitch (Nm)','Voluntary Activation (%)'}};     

% HRV data processing options
artefact_recognition = 1;                               % 0 = no recognition, 1 = recognition with moving median +- artefactThreshold
artefact_threshold   = 250;                             % Threshold in ms 
detrending           = 1;                               % 0 = no detrending, 1 = detrending by removing 3rd order polynomial, 2 = zero phase butterworth filter cut-off frequency: 0.035 Hz

% Analysis options
includedPhases      = {'INT', 'REC'};                   % Phases from study_design that will be analysed. Number of measurement in these phases must be greater than 1

correlation_Type    = [   3  ,...                       % Options Reproducibility:  1 = Pearson correlation, 2 = Intraclass correlation, 3 = ANCOVA 
                          3      ];                     % Options Association:      1 = Pearson correlation, 3 = ANCOVA 
    ICC_type        = '3-1';                            % Options: '1-1', '1-k', '2-1', '2-k', '3-1', '3-k'     
    ANCOVA_type     = {'parallel lines',...             % Options Reproducibility:  'parrallel lines' - similar to repeated measures correlation, 'seperate lines' - individual fit
                       'parallel lines'};               % Options Association:      'parrallel lines' - similar to repeated measures correlation, 'seperate lines' - individual fit
alpha               = 0.1;                             % Alpha level, Type 1 Error rate 
NULL                = [   0.875  ,...                   % Reproducibility Null-Hypothesis , Range 0 to 1 (absolute), Inferiority test. Tests if measurement is smaller than NULL
                          0.875     ];                  % Association Null-Hypothesis ,     Range 0 to 1 (absolute), Inferiority test. Tests if measurement is smaller than NULL
side                = 2;                                % Options: 1 = one-sided test, 2 = two-sided test. For one-sided test, the NULL has to be in the proper direction

HR_variable         = 4;                                % Selects the HR variable, to test association to. Options refer to 'outputVariables' index of the first row, e.g. 4 = RMSSD (ms)

% Export options
export_options      = [   1    ,...                     % 1 = exports groupData and statistics tables to desired folder, 0 otherwise
                          1    ,...                     % 1 = exports graphical illustration of statistics for each metric, 0 otherwise
                          1       ];                    % 1 = exports plots effect sizes for each individual, 0 otherwise
% Plot options                      
disp_vars           = {[4         ],...                 % Determines the HRV variables to be plotted for each individual
                       [15, 20, 23],...                 % Determines the CMJ variables to be plotted for each individual
                       [26, 30, 32],...                 % Determines the SQV variables to be plotted for each individual
                       [33, 35, 36]};                   % Determines the ITT variables to be plotted for each individual
                   
xlabels = {'PRE';'POST24';'POST48';'POST72';...         % Sets the xlabels for individual plots, must match number of phases that will be included due to 'includedPhases',
           'PRE_2';'POST24_2';'POST48_2';'POST72_2';};  % Must match number of measurement repetitions (E.g.: in this study, measurements were conducted 2 times: 'PRE_2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main functions
% Do not alter following commands, unless you know what you do.
clc
timerVal = tic;
addpath(string(pwd) + '\Functions');                    % Adds the folder \Functions to the working folder to access distinct functions
 
[groupData, individualData] = DataProcessing(study_design, output_variables, artefact_recognition, artefact_threshold, detrending);
transformedData = DataTransformation(groupData, baseline);
[reproducibilityResults, associationResults] = DataAnalysis(transformedData, includedPhases, HR_variable, correlation_Type, ICC_type, ANCOVA_type, alpha, NULL, side);
Export(export_options);

timerVal = toc(timerVal);
fprintf('\n###################################\n\nTotal Duration: \t%-2.2f seconds\n\n', timerVal);
clearvars -except groupData individualData transformedData reproducibilityResults associationResults 
