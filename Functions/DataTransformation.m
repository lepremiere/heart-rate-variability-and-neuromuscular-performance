function analysisData = DataTransformation (groupData, baseline_name)
% This function creates transformed versions of the table "groupData" w.r.t. a baseline.
%
% USAGE
% E.g.:
% transformedData = DataTransformation (groupData, 'BASE')
%
% INPUT
% groupData:                Table that 

%
% ouput_variables:          String array that determines how the calculated variables will 
%                           be called. Number of inputs must match number of calculated
%                           variables returned by functions "CounterMovementJump.m",
%                           "BarVelocity.m" "Dynamometry.m", "HRV_Analysis.m".
%
% artefact_recognition:     Scalar that determines if artefacts should be identified and 
%                           removed during "HRV_Analysis.m". Options: 0 = no recognition,
%                           1 = recognition with moving median channel
%
% artefact_threshold:       Scalar that determines the width of the moving median channel for 
%                           artefact recognition in ms. 
%                           
% detrending:               Scalar that determines if HRV data should be detrended during
%                           "HRV_Analysis.m". Options: 0 = no detrending, 1 = detrending
%                           by removing 3rd order polynomial, 2 = zero phase butterworth 
%                           filter cut-off frequency: 0.035 Hz
%
% OUTPUT
% groupData:                Output table with dimension depending on input. Number of columns
%                           equal number of elements in "study_design" + 2. Number of Rows 
%                           depends on number of participants and variables to be calculated.
%                           Table includes every variable for every participant over the 
%                           distinct events. Missing values are filled with NaNs. If a file
%                           for one of the four variable groups is found, the table includes
%                           every metric of that analysis for each participant. This table 
%                           includes only the best trial of each variable for each event.
%
% individualData:           Table (5xparticipants) that stores structs for each variable group
%                           (CMJ, HRV, SQV, DYN) for every participant and a summary table of
%                           best trials for each event.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 09.July.2020
% lepremiere
%---------------------------------------------------------------------------------------------------

%% Calculating baseline mean and standard deviation
% Baseline is calculated on all columns that contain the user input 'baselineName'. e.g 'BASE' 
% Calculates mean and std for each variable and participant with respect to the individual baseline
study_design    = groupData.Properties.VariableNames;
baseLine        = find(contains(study_design, baseline_name) == 1);
                                               
baselineMean    = nanmean(groupData{:, baseLine},2);
baselineSD      = nanstd(groupData{:, baseLine}')';

unique_identifier = numel(unique(groupData{:, 3}));

%% Calculates group mean and standard deveiation for each variable
% Mean
ABS_mean = varfun(@nanmean,groupData(:, 3:end),'GroupingVariables','Variable');                       % Calulates the mean for each variable
ABS_mean = array2table(ABS_mean{:, 2:end},...
           'VariableNames',['GroupCount',groupData.Properties.VariableNames(4:end)],...      % Converts results to table with similar Var-/RowNames
           'RowNames'     ,ABS_mean{:,1});
ABS_mean = [groupData(1:height(ABS_mean),2:3),...
            ABS_mean(groupData{1:unique_identifier, 3},:)];                                             % Sorts table to match input

% SD
ABS_sd = varfun(@nanstd,groupData(:, 3:end),'GroupingVariables','Variable');                          % Calulates the standard deviation for each variable
ABS_sd = array2table(ABS_sd{:, 2:end},...
                                'VariableNames',['GroupCount',groupData.Properties.VariableNames(4:end)],...    
                                'RowNames'     ,ABS_sd{:,1});
ABS_sd = [groupData(1:height(ABS_sd),2:3),...
          ABS_sd(groupData{1:unique_identifier, 3},:)];

%% Calulates effect size glass' d (ES = (x - mean)/std) for each cell and the mean and std for each variable
ES_individual = [groupData(:,1:3),...
                 array2table((groupData{:, 4:end}-baselineMean)./baselineSD,...               
                 'VariableNames',groupData.Properties.VariableNames(4:end),...    
                 'RowNames'     ,groupData.Properties.RowNames)];
% GroupMean
ES_mean = varfun(@nanmean,ES_individual(:, 3:end),'GroupingVariables','Variable');
ES_mean = array2table(ES_mean{:,2:end},...
          'VariableNames',['GroupCount',groupData.Properties.VariableNames(4:end)],...
          'RowNames',ES_mean{:,1});
ES_mean = [groupData(1:height(ES_mean),2:3),...
           ES_mean(ES_individual{1:unique_identifier, 3},:)];

% SD
ES_sd = varfun(@nanstd,ES_individual(:, 3:end),'GroupingVariables','Variable');
ES_sd = array2table(ES_sd{:,2:end},...
        'VariableNames',['GroupCount',groupData.Properties.VariableNames(4:end)],...
        'RowNames',ES_sd{:,1});
ES_sd = [groupData(1:height(ES_sd),2:3),...   
        ES_sd(groupData{1:unique_identifier, 3},:)];

%% Calulates percentage of baseline (perc = (x /mean)*100) for each cell and the mean and std for each variable
percentage_individual = [groupData(:,1:3),...                                                                   
                         array2table((groupData{:, 4:end}./baselineMean)*100,...
                         'VariableNames',groupData.Properties.VariableNames(4:end),...   
                         'RowNames',groupData.Properties.RowNames)];
% GroupMean
percentage_mean = varfun(@nanmean,percentage_individual(:,3:end),'GroupingVariables','Variable');
percentage_mean = array2table(percentage_mean{:,2:end},...
                  'VariableNames',['GroupCount',groupData.Properties.VariableNames(4:end)],...
                  'RowNames',percentage_mean{:,1});
percentage_mean = [groupData(1:height(percentage_mean),2:3),...
                   percentage_mean(percentage_individual{1:unique_identifier, 3},:)];

% SD
percentage_sd = varfun(@nanstd,percentage_individual(:,3:end),'GroupingVariables','Variable');
percentage_sd = array2table(percentage_sd{:,2:end},...
                'VariableNames',['GroupCount',groupData.Properties.VariableNames(4:end)],...
                'RowNames',percentage_sd{:,1});
percentage_sd = [groupData(1:height(percentage_mean),2:3),...
                percentage_sd(groupData{1:unique_identifier, 3},:)];

%% Delta-Values

deltaVarNames = join([study_design(4:end-1);study_design(5:end)]');

% Absolutes
delta_ABS_individual     =  groupData{:, 5:end} - groupData{:, 4:end-1};
delta_ABS_individual     = [ groupData(:, 1:3),array2table(delta_ABS_individual,'VariableNames',deltaVarNames)];

% GroupMean   
delta_ABS_mean                                   = varfun(@nanmean,delta_ABS_individual(:, 3:end),'GroupingVariables','Variable');
delta_ABS_mean.Properties.VariableNames(3:end)   = delta_ABS_individual.Properties.VariableNames(4:end);
delta_ABS_mean                                   =  array2table(delta_ABS_mean{:,2:end},'VariableNames',delta_ABS_mean.Properties.VariableNames(2:end),'RowNames',delta_ABS_mean{:,1});
delta_ABS_mean                                   = [groupData(1:height(delta_ABS_mean),2:3),... 
                                                    delta_ABS_mean(delta_ABS_individual{1:unique_identifier,3},:)];

% SD
delta_ABS_sd                                     = varfun(@nanstd,delta_ABS_individual(:, 3:end),'GroupingVariables','Variable');
delta_ABS_sd.Properties.VariableNames(3:end)     = delta_ABS_individual.Properties.VariableNames(4:end);
delta_ABS_sd                                     = array2table(delta_ABS_sd{:,2:end},'VariableNames',delta_ABS_sd.Properties.VariableNames(2:end),'RowNames',delta_ABS_sd{:,1});
delta_ABS_sd                                     = [groupData(1:height(delta_ABS_sd),2:3),... 
                                                    delta_ABS_sd(groupData{1:unique_identifier,3},:)];

% Percentages
delta_PERC_individual     =  percentage_individual{:,5:end} - percentage_individual{:,4:end-1};
delta_PERC_individual     = [percentage_individual(:,1:3),array2table(delta_PERC_individual,'VariableNames',deltaVarNames)];

% GroupMean   
delta_PERC_mean                                   = varfun(@nanmean,delta_PERC_individual(:,3:end),'GroupingVariables','Variable');
delta_PERC_mean.Properties.VariableNames(3:end)   = delta_PERC_individual.Properties.VariableNames(4:end);
delta_PERC_mean                                   = array2table(delta_PERC_mean{:,2:end},'VariableNames',delta_PERC_mean.Properties.VariableNames(2:end),'RowNames',delta_PERC_mean{:,1});
delta_PERC_mean                                   = [groupData(1:height(delta_PERC_mean),2:3),... 
                                                    delta_PERC_mean(delta_PERC_individual{1:unique_identifier,3},:)];

% SD
delta_PERC_sd                                     = varfun(@nanstd,delta_PERC_individual(:,3:end),'GroupingVariables','Variable');
delta_PERC_sd.Properties.VariableNames(3:end)     = delta_PERC_individual.Properties.VariableNames(4:end);
delta_PERC_sd                                     = array2table(delta_PERC_sd{:,2:end},'VariableNames',delta_PERC_sd.Properties.VariableNames(2:end),'RowNames',delta_PERC_sd{:,1});
delta_PERC_sd                                     = [groupData(1:height(delta_PERC_sd),2:3),... 
                                                    delta_PERC_sd(groupData{1:unique_identifier,3},:)];

% Effect sizes
delta_ES_individual     =  ES_individual{:,5:end} - ES_individual{:,4:end-1};
delta_ES_individual     = [ES_individual(:,1:3),array2table(delta_ES_individual,'VariableNames',deltaVarNames)];

% GroupMean   
delta_ES_mean                                   = varfun(@nanmean,delta_ES_individual(:,3:end),'GroupingVariables','Variable');
delta_ES_mean.Properties.VariableNames(3:end)   = delta_ES_individual.Properties.VariableNames(4:end);
delta_ES_mean                                   = array2table(delta_ES_mean{:,2:end},'VariableNames',delta_ES_mean.Properties.VariableNames(2:end),'RowNames',delta_ES_mean{:,1});
delta_ES_mean                                   = [groupData(1:height(delta_ES_mean),2:3),... 
                                                  delta_ES_mean(delta_ES_individual{1:unique_identifier,3},:)];

% SD
delta_ES_sd                                     = varfun(@nanstd,delta_ES_individual(:,3:end),'GroupingVariables','Variable');
delta_ES_sd.Properties.VariableNames(3:end)     = delta_ES_individual.Properties.VariableNames(4:end);
delta_ES_sd                                     = array2table(delta_ES_sd{:,2:end},'VariableNames',delta_ES_sd.Properties.VariableNames(2:end),'RowNames',delta_ES_sd{:,1});
delta_ES_sd                                     = [groupData(1:height(delta_ES_sd),2:3),... 
                                                  delta_ES_sd(groupData{1:unique_identifier,3},:)];



%% Puts output together
Results         = { groupData,...
                    percentage_individual,...
                    ES_individual,...
                    delta_ABS_individual,...
                    delta_PERC_individual,...
                    delta_ES_individual,...
                    ABS_mean,...
                    ABS_sd,...
                    percentage_mean,...
                    percentage_sd,...
                    ES_mean,...
                    ES_sd,...
                    delta_ABS_mean,...
                    delta_ABS_sd,...
                    delta_PERC_mean,...
                    delta_PERC_sd,...
                    delta_ES_mean,...
                    delta_ES_sd}';
                
analysisData = array2table(Results,'RowNames',{'Absolute_individual',...
                                                  'Percentage_individual',...
                                                  'ES_individual',...
                                                  'delta_ABS_individual',...
                                                  'delta_PERC_individual',...
                                                  'delta_ES_individual',...
                                                  'ABS_mean',...
                                                  'ABS_sd',...
                                                  'Percentage_mean',....
                                                  'Percentage_sd'...
                                                  'ES_mean',...
                                                  'ES_sd',...
                                                  'delta_ABS_mean',...
                                                  'delta_ABS_sd',...
                                                  'delta_PERC_mean',...
                                                  'delta_PERC_sd',...
                                                  'delta_ES_mean',...
                                                  'delta_ES_sd'});
end