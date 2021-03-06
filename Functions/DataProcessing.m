function [individualData, groupData, analysisData]...
            = DataProcessing(data, study_design, output_variables, artefact_recognition, artefact_threshold, detrending, baseline)
% This function processes raw data from force plates, dynamometer and hrv recordings and calculates
% variables of interest. Additionally, processed data will be sorted w.r.t name tags of the files 
% and the study design.
%
% USAGE
% E.g.:
% [groupData, individualData] = DataProcessing(study_design, output_variables, 1, 250, 1)
%
% INPUT
% data:                     Contains a struct with fields including the files for each 
%                           test (CMJ, SQV, DYNO, HRV) in cell arrays. The last field contains
%                           a cell array with participant identifiers. "data" is the output
%                           of the function "ImportData.m".
% study_design:             String array that determines which files will be combined to 
%                           distinct events. E.g.: {"1_BASE_1, 1_BASE_2,..."}
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
% individualData:           Table (n x participants) that stores structs for each variable group
%                           (CMJ, HRV, SQV, DYN) for every participant and a summary table of
%                           best trials for each event.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 09.July.2020
% lepremiere
%---------------------------------------------------------------------------------------------------

tic
fprintf('\nProcessing: \t\t')

%% Counter-Movement Jump (CMJ)
% Calculating the jump qualities. Function 'CounterMovementJump'. Creates 'Output_CMJ' with all results and files from input
try
    %% Getting orientation within input 
    % Converting the filenames of CMJ files from char to string 
    CMJ = data.CMJ;
    files = string(CMJ(:,1:2));

    % Identifying positions of seperator underscore in filenames.
    % Exceptional case for only 1 CMJ file is given
    if(length(files(:,2))==1)
        pos{1}  = strfind(files(:,2),'_'); 
    else
        pos     = strfind(files(:,2),'_');
    end

    % Extracting information blocks between seperator underscore
    for i = 1:length(files(:,2))
        patterns(i,1)   = extractBetween(files(i,2),1,pos{i,1}(end)-1);
        patterns(i,2)   = extractBetween(files(i,2),1,pos{i,1}(1)-1);
        temp_chars(i,1) = extractBetween(files(i,2),pos{i,1}(1)+1,pos{i,1}(end)-1);
    end

    % Assigning unique blocks regarding information 
    testOccasions   = unique(patterns(:,1));
    subject         = unique(patterns(:,2));
    temp_chars      = unique(temp_chars);

    % Identifiying the participant of the test occasions by finding name prefix in test occasion. 
    % Note: Test occasions include name prefix as well as study phase and number of test occasion
    for i = 1:length(subject)
        testedParticipant{i} = contains(testOccasions, subject(i)+'_'+temp_chars);
    end

    % Converting logic (1,0) output of function above to array indecies (1,2,3...)
    for i = 1:length(testedParticipant)
        participantIndex{i} = find(testedParticipant{:,i}==1);
    end

    % Identifying distinct study phases of the subjects
    for i = 1:length(testOccasions)
        trial_number{i} = (count(files(:,2),testOccasions(i)));
    end

    %Converting logic (1,0) output of function above to array indecies (1,2,3...)
    for i = 1:length(trial_number)
        orientation{i} = find(trial_number{:,i}==1);
    end

    %% Looping throug jump data
    % Looping through subjects
    for subjectNumber = 1:length(participantIndex)

        subjectArray            = participantIndex{subjectNumber};     % Selecting the arrays of test occasions that belong to participant with regarding subject number
        subjectTestOccasions	= testOccasions(subjectArray);         % Selecting only the test occasions of the subject

        % Looping through subjects test occasions. n = number of test oaccasion
        for n = 1:numel(subjectArray)

            % Looping through subject's trials on specific test occasion. m = number of trial
            for m = 1:length(orientation{subjectArray(n)})          

                  % Actual analysis of ground reaction forces from raw data 
                  % Further infromation can be found in function 'CounterMovementJump'
                  Result(m,:) = CounterMovementJump(CMJ{orientation{subjectArray(n)}(m),3});
            end
            
            % Gathering filenames for result description with prefix CMJ
            for i = 1:length(orientation{subjectArray(n)})
                 RowNames{i} = ['CMJ_',files{orientation{subjectArray(n)}(i),2}];
            end

            Result                      = array2table(Result);              % Creating table from all trials
            Height                      = height(Result);
            Width                       = width(Result);
            Result.Properties.RowNames  = RowNames;                         % Assigning rownames from filenames
            Result                      = sortrows(Result,3,'descend');     % Sorting rows for third column. Preset = 'netImpulse' (3)
            
            % Calculating mean, standard deviation and finding max for each column
            for i = 1:Width
                 Result{Height+1,i} = mean(Result{1:Height,i},1);
                 Result{Height+2,i} = std(Result{1:Height,i},1)/Result{Height+1,i}*100;
                 Result{Height+3,i} = max(Result{1:Height,i});
            end
            
            % Declaring the variable names with units to columns of table
            Result.Properties.VariableNames     =   output_variables{2};
            Result.Properties.RowNames(end-2)   =   {'Mean'};
            Result.Properties.RowNames(end-1)   =   {'CV'};
            Result.Properties.RowNames(end)     =   {'Best'};
            
            % Saving final table to cell array to later join into struct
            Output_CMJ{n} = Result;

            clearvars Result RowNames 
        end
        
        % Creating table with only best trials over test occasions
        for i = 1:numel(subjectTestOccasions)
            Best_Trials(i,:)                        = (Output_CMJ{i}(end,:));
            Best_Trials.Properties.RowNames(end)    = string(['CMJ_',subjectTestOccasions{i}]);
        end

        Height  = height(Best_Trials);
        Width   = width(Best_Trials);
        
        % Calculating mean and standard deviation for best trials
        for i = 1:Width
            Best_Trials{Height+1,i} = mean(Best_Trials{1:Height,i},1);
            Best_Trials{Height+2,i} = std(Best_Trials{1:Height,i},1)/Best_Trials{Height+1,i}*100;
        end

        Best_Trials.Properties.RowNames(end-1)  = {'Mean'}; 
        Best_Trials.Properties.RowNames(end)    = {'CV'};
        
        % Adding prefix 'CMJ' to filename. DELETABLE if named correctly in the first place
        for i = 1:length(subjectTestOccasions)
            subjectTestOccasions(i) = string(['CMJ_',subjectTestOccasions{i}]);
        end
        
        % Additional field names for best trials, fields and raw data
        subjectTestOccasions{end+1} = 'Best_Trials';
        subjectTestOccasions{end+1} = 'Fields';
        subjectTestOccasions{end+1} = 'Data';
        
        % Adding best trials and fields to cell array with results
        Output_CMJ{end+1} = Best_Trials;
        Output_CMJ{end+1} = subjectTestOccasions;
        
        % Getting raw data of participant and attaching it to cell array
        sets                         = vertcat(orientation{subjectArray});
        Output_CMJ{end+1}            = CMJ(sets,:);
        
        % Saving struct of subject in placeholder 
        Output_CMJ1{subjectNumber,:} = cell2struct(Output_CMJ,subjectTestOccasions,2);

        clearvars Output_CMJ Best_Trials

    end
    
    % Creating struct of all participant and trials 
    Output_CMJ = cell2struct(Output_CMJ1,subject,1);
    clearvars -except CMJ SQV DYNO HRV input Output_CMJ Output_SQV ...
                      Output_DYNO Output_HRV arl brl crl artefact_threshold ...
                      variables subjectPrefix study_design output_variables ...
                      artefact_recognition detrending saving data baseline

catch

% Error message in case something above did not work. Does not stop script
    if(exist('CMJ','var') == 0)
        fprintf('\nNo CMJ files found (3)\n')
    else
        fprintf('\nERROR CMJ (4)')
    end

end

%% Squat Velocity (SQV)
% Calculating the jump qualities. Function 'CounterMovementJump'. Creates 'Output_CMJ' with all results and files from input
try    
    %% Getting orientation within input 
    % Converting the filenames of CMJ files from char to string 
    clearvars files
    SQV = data.SQV;
    files = string(SQV(:,1:2));

    % Identifying positions of seperator underscore in filenames.
    % Exceptional case for only 1 SQV file is given
    if(length(files(:,2))==1)
        pos{1}  = strfind(files(:,2),'_'); 
    else
        pos     = strfind(files(:,2),'_');
    end

    % Extracting information blocks between seperator underscore
    for i = 1:length(files(:,2))
        patterns(i,1) = extractBetween(files(i,2),1,pos{i,1}(end)-1);
        patterns(i,2) = extractBetween(files(i,2),1,pos{i,1}(1)-1);
        temp_chars(i,1) = extractBetween(files(i,2),pos{i,1}(1)+1,pos{i,1}(end)-1);
    end

    % Assigning unique blocks regarding information 
    testOccasions   = unique(patterns(:,1));
    subject         = unique(patterns(:,2));
    temp_chars = unique(temp_chars);

    % Identifiying the participant of the test occasions by finding name prefix in test occasion. 
    % Note: Test occasions include name prefix as well as study phase and number of test occasion
    for i = 1:length(subject)
        testedParticipant{i} = contains(testOccasions, subject(i)+'_'+temp_chars);
    end

    % Converting logic (1,0) output of function above to array indecies (1,2,3...)
    for i = 1:length(testedParticipant)
        participantIndex{i} = find(testedParticipant{:,i}==1);
    end

    % Identifying distinct study phases of the subjects
    for i = 1:length(testOccasions)
        trial_number{i} = (count(files(:,2),testOccasions(i)));
    end

    %Converting logic (1,0) output of function above to array indecies (1,2,3...)
    for i = 1:length(trial_number)
        orientation{i} = find(trial_number{:,i}==1);
    end

    %% Looping throug jump data
    % Looping through subjects
    for subjectNumber = 1:length(participantIndex)

        subjectArray            = participantIndex{subjectNumber};     % Selecting the arrays of test occasions that belong to participant with regarding subject number
        subjectTestOccasions	= testOccasions(subjectArray);         % Selecting only the test occasions of the subject

        % Looping through subjects test occasions. n = number of test oaccasion
        for n = 1:numel(subjectArray)

            % Looping through subject's trials on specific test occasion. m = number of trial
            for m = 1:length(orientation{subjectArray(n)})          
                  % Actual analysis of ground reaction forces from raw data 
                  % Further infromation can be found in function 'CounterMovementJump'
                  Result(m,:) = SquatVelocity(SQV{orientation{subjectArray(n)}(m),3});
            end
            
            % Gathering filenames for result description with prefix CMJ
            for i = 1:length(orientation{subjectArray(n)})
                 RowNames{i} = ['SQV_',files{orientation{subjectArray(n)}(i),2}];
            end

            Result                      = array2table(Result);              % Creating table from all trials
            Height                      = height(Result);
            Width                       = width(Result);
            Result.Properties.RowNames  = RowNames;                         % Assigning rownames from filenames
            Result                      = sortrows(Result,2,'descend');     % Sorting rows for third column. Preset = 'netImpulse' (3)
            Result.Properties.VariableNames = output_variables{3};          % Declaring the variable names with units to columns of table
            
            % Calculating mean, standard deviation and finding max for each column
            for i = 1:Width
                 Result{Height+1,i} = mean(Result{1:Height,i},1);
                 Result{Height+2,i} = std(Result{1:Height,i},1)/Result{Height+1,i}*100;             
                 if(contains(Result.Properties.VariableNames(i), "%1RM"))
                    Result{Height+3,i} = min(Result{1:Height,i});
                 else
                     Result{Height+3,i} = max(Result{1:Height,i});
                 end
            end
            
            Result.Properties.VariableNames     =   output_variables{3};
            Result.Properties.RowNames(end-2)   =   {'Mean'};
            Result.Properties.RowNames(end-1)   =   {'CV'};
            Result.Properties.RowNames(end)     =   {'Best'};
            
            % Saving final table to cell array to later join into struct
            Output_SQV{n} = Result;

            clearvars Result RowNames 

        end
        
        % Creating table with only best trials over test occasions
        for i = 1:numel(subjectTestOccasions)
            Best_Trials(i,:)                        = (Output_SQV{i}(end,:));
            Best_Trials.Properties.RowNames(end)    = string(['SQV_',subjectTestOccasions{i}]);
        end

        Height  = height(Best_Trials);
        Width   = width(Best_Trials);
        
        % Calculating mean and standard deviation for best trials
        for i = 1:Width
            Best_Trials{Height+1,i} = mean(Best_Trials{1:Height,i},1);
            Best_Trials{Height+2,i} = std(Best_Trials{1:Height,i},1)/Best_Trials{Height+1,i}*100;
        end

        Best_Trials.Properties.RowNames(end-1)  = {'Mean'}; 
        Best_Trials.Properties.RowNames(end)    = {'CV'};
        
        % Adding prefix 'CMJ' to filename. DELETABLE if named correctly in the first place
        for i = 1:length(subjectTestOccasions)
            subjectTestOccasions(i) = string(['SQV_',subjectTestOccasions{i}]);
        end
        
        % Additional field names for best trials, fields and raw data
        subjectTestOccasions{end+1} = 'Best_Trials';
        subjectTestOccasions{end+1} = 'Fields';
        subjectTestOccasions{end+1} = 'Data';
        
        % Adding best trials and fields to cell array with results
        Output_SQV{end+1} = Best_Trials;
        Output_SQV{end+1} = subjectTestOccasions;
        
        % Getting raw data of participant and attaching it to cell array
        sets                         = vertcat(orientation{subjectArray});
        Output_SQV{end+1}            = SQV(sets,:);
        
        % Saving struct of subject in placeholder 
        Output_SQV1{subjectNumber,:} = cell2struct(Output_SQV,subjectTestOccasions,2);

        clearvars Output_SQV Best_Trials
    end
    
    % Creating struct of all participant and trials 
    Output_SQV = cell2struct(Output_SQV1,subject,1);
    clearvars -except CMJ SQV DYNO HRV input Output_CMJ Output_SQV Output_DYNO Output_HRV...
                arl brl crl artefact_threshold variables subjectPrefix study_design...
                output_variables artefact_recognition detrending saving data baseline

catch

% Error message in case something above did not work. Does not stop script
    if(exist('SQV','var') == 0)
        fprintf('\nNo SQV files found (3)\n')
    else
        fprintf('\nERROR SQV (4)')
    end

end
%% Maximum Voluntary Contraction (MVC), Voluntary Activation (VA)
% Extracting maximum voluntary, superimposed and resting twitch torque, to calculate voluntary activation. Function 'Dynamometry'. Creates 'Output_DYNO' with all results and files from input
try
    %% Getting orientation within input 
    % Converting the filenames of DYNO files from char to string    
    clearvars files
    DYNO   = data.DYNO;
    files = string(DYNO(:,1:2));

    % Identifying positions of seperator underscore in filenames.
    % Exceptional case for only 1 DYNO file is given
    if(length(files(:,2)) == 1)
        pos{1}  = strfind(files(:,2),'_'); 
    else
        pos     = strfind(files(:,2),'_');
    end

    % Extracting information blocks between seperator underscore
    for i = 1:length(files(:,2))
        patterns(i,1) = extractBetween(files(i,2),1,pos{i,1}(end)+1);
        patterns(i,2) = extractBetween(files(i,2),1,pos{i,1}(1)-1);
        temp_chars(i,1) = extractBetween(files(i,2),pos{i,1}(1)+1,pos{i,1}(end)+1);
    end

    % Assigning unique blocks regarding information 
    testOccasions   = unique(patterns(:,1));
    subject         = unique(patterns(:,2));
    temp_chars = unique(temp_chars);

    % Identifiying the participant of the test occasions by finding name prefix in test occasion. 
    % Note: Test occasions include name prefix as well as study phase and number of test occasion
    for i = 1:length(subject)
        testedParticipant{i} = contains(testOccasions, subject(i)+'_'+temp_chars);
    end

    % Converting logic (1,0) output of function above to array indecies (1,2,3...)
    for i = 1:length(testedParticipant)
        participantIndex{i} = find(testedParticipant{:,i}==1);
    end

    % Identifying distinct study phases of the subjects
    for i = 1:length(testOccasions)
        trial_number{i} = (count(files(:,2),testOccasions(i)));
    end

    %Converting logic (1,0) output of function above to array indecies (1,2,3...)
    for i = 1:length(trial_number)
        orientation{i} = find(trial_number{:,i}==1);
    end
    
    %% Looping through data
    % Looping through subjects
    for subjectNumber = 1:length(participantIndex)

        subjectArray = participantIndex{subjectNumber};         % Selecting the arrays of test occasions that belong to participant with regarding subject number
        subjectTestOccasions = testOccasions(subjectArray);     % Selecting only the test occasions of the subject
        
        % Looping through subject's test occasions. n = number of test oaccasion
        for n = 1:numel(subjectTestOccasions)
            
            % Actual analysis of torque and keypresses from raw data 
            % Further infromation can be found in function 'Dynamometry'
            Results = Dynamomety(DYNO{orientation{subjectArray(n)},3},subjectTestOccasions,n);
            Results.Properties.VariableNames = output_variables{4};
            Output_DYNO{n} = Results;
        end

        % Creating table with only best over test occasions
        for i = 1:numel(subjectTestOccasions)
            Best_Trials(i,:) = Output_DYNO{i}(end,:);
            Best_Trials.Properties.RowNames(end)=(Output_DYNO{i}(1,1).Properties.RowNames);
        end

        Height = height(Best_Trials);
        Width = width(Best_Trials);
        
        % Calculating mean and coefficient of variation for each table column
        for i = 1:Width
            Best_Trials{Height+1,i} = mean(Best_Trials{1:Height,i},1);
            Best_Trials{Height+2,i} = std(Best_Trials{1:Height,i},1)/Best_Trials{Height+1,i}*100;
        end
        
        % Declaring variable names for output
        Best_Trials.Properties.RowNames(end-1) = {'Mean'};
        Best_Trials.Properties.RowNames(end) ={'CV'};

        subjectTestOccasions{end+1} = 'Best_Trials';
        subjectTestOccasions{end+1} = 'Fields';
        subjectTestOccasions{end+1} = 'Data';
        
        % Filling cell array with results and summary of test occasions
        Output_DYNO{end+1} = Best_Trials;
        Output_DYNO{end+1} = subjectTestOccasions(1:end);
        
        % Getting file orientations of subject 
        sets = vertcat(orientation{subjectArray});
        
        % Adding raw data to subjects output
        Output_DYNO{end+1} = DYNO(sets,:);
        Output_DYNO1{subjectNumber,:} = cell2struct(Output_DYNO,subjectTestOccasions,2);

        clearvars Output_DYNO Best_Trials
    end
    
    % Creating struct of all participant and trials  
    Output_DYNO = cell2struct(Output_DYNO1,subject,1); 
    clearvars -except CMJ SQV DYNO HRV input Output_CMJ Output_SQV Output_DYNO Output_HRV...
                    arl brl crl artefact_threshold variables subjectPrefix study_design...
                    output_variables artefact_recognition detrending saving data baseline
    
catch
    % Error message in case something above did not work. Does not stop script
    if(exist('DYNO','var') == 0)
        fprintf('\nNo DYNO files found (6)\n')
    else
        fprintf('\nERROR DYNO (7)')
    end
end

%% Heart Rate Variability (HRV)
% Identifying artifacts, detrending and calculating HRV-indices of RR-data. Function 'HRV_Analysis'. Creates 'Output_HRV' with all results and files from input
try
   %% Getting orientation within input 
    % Converting the filenames of HRV files from char to string    
    clearvars files
    HRV = data.HRV;
    files = string(data.HRV(:,1:2));

    % Identifying positions of seperator underscore in filenames.
    % Exceptional case for only 1 HRV file is given
    if(length(files(:,2)) == 1)
        pos{1}  = strfind(files(:,2),'_'); 
    else
        pos     = strfind(files(:,2),'_');
    end
    
    % Extracting information blocks between seperator underscore
    for i = 1:length(files(:,2))
        patterns(i,1) = extractBetween(files(i,2),1,pos{i,1}(end)+1);
        patterns(i,2) = extractBetween(files(i,2),1,pos{i,1}(1)-1);
        temp_chars(i,1) = extractBetween(files(i,2),pos{i,1}(1)+1,pos{i,1}(end)+1);
    end

    % Assigning unique blocks regarding information 
    testOccasions   = unique(patterns(:,1));
    subject         = unique(patterns(:,2));
    temp_chars = unique(temp_chars);

    % Identifiying the participant of the test occasions by finding name prefix in test occasion. 
    % Note: Test occasions include name prefix as well as study phase and number of test occasion
    for i = 1:length(subject)
        testedParticipant{i} = contains(patterns(:,1), subject(i)+'_'+temp_chars);
    end

    % Converting logic (1,0) output of function above to array indecies (1,2,3...)
    for i = 1:length(testedParticipant)
        participantIndex{i} = find(testedParticipant{:,i}==1);
    end

    % Identifying distinct study phases of the subjects
    for i = 1:length(testOccasions)
        trial_number{i} = (count(files(:,2),testOccasions(i)));
    end

    %Converting logic (1,0) output of function above to array indecies (1,2,3...)
    for i = 1:length(trial_number)
        orientation{i} = find(trial_number{:,i}==1);
    end
    
    %% Looping through data
    % Looping through subjects
    for subjectNumber = 1:length(participantIndex)

        subjectArray = participantIndex{subjectNumber};         % Selecting the arrays of test occasions that belong to participant with regarding subject number
        subjectTestOccasions = testOccasions(subjectArray);     % Selecting only the test occasions of the subject
        
        % Looping through subjects test occasions
        for n = 1:numel(subjectArray)
            
            % Actual analysis of torque and keypresses from raw data 
            % Further infromation can be found in function 'HRV'
            Result(n,:) = HRV_Analysis(HRV{subjectArray(n),3}, artefact_recognition, detrending, artefact_threshold);
        end
        
        % Assigning variable names with units to columns of result table
        Result = array2table(Result, 'VariableNames', output_variables{1});
        
        % Getting subjects test occasions                          
        for i = 1:length(subjectArray)
            RowNames{i} = patterns{subjectArray(i),1};
        end

        Height = height(Result);
        Width = width(Result);
        Result.Properties.RowNames = RowNames;                  % Assigning test occasions to row names
        
        % Calculating mean and coefficient of variation for each table column
        for i = 1:Width           
            Result{Height+1,i} = mean(Result{1:Height,i},1);
            Result{Height+2,i} = std(Result{1:Height,i},1)/Result{Height+1,i}*100;           
        end

        Result.Properties.RowNames(end-1) = {'Mean'};
        Result.Properties.RowNames(end) ={'CV'};
        
        % Getting indecies of subjects raw data and filenames
        index = cell2mat(orientation);
        
        % Adding subjects results, file names and raw data to a cell array
        Output_HRV{1,subjectNumber} = Result;
        Output_HRV{2,subjectNumber} =[HRV(index(subjectArray),2), HRV(index(subjectArray),3)];
        
        % Creating struct for the subject with all results and data
        Output_HRV1{subjectNumber,1} = cell2struct(Output_HRV(1:2,subjectNumber),{'HRV','Data'},1);

        clearvars Output_HRV Result RowNames
    end
    
    % Creating struct for all subjects and data
    Output_HRV = cell2struct(Output_HRV1,subject,1);
    
    clearvars -except CMJ SQV DYNO HRV input Output_CMJ Output_SQV Output_DYNO Output_HRV...
                    arl brl crl artefact_threshold variables subjectPrefix study_design...
                    output_variables artefact_recognition detrending saving data baseline
 
catch
    % Error message in case something above did not work. Does not stop script
    if(exist('HRV','var') == 0)
        fprintf('\nNo HRV files found (7)\n')
    else
        fprintf('\nERROR HRV (8)\n')
    end
end

%% Fill Tables
% Concatenating and summarizing function outputs. Filters data regarding test occasions declared in study design.
% Creates 'groupData' with all variables for each pearticipant and 'individualData' with all data and individual summary
try
    % Creates cell array with result structs, sorted by subjects and variables
    subjectPrefix   = data.Participants;
    variables       = fieldnames(data);
    variables       = variables(1:end-1);
    for i = 1:length(subjectPrefix)
        for j = 1:size(variables,1)
            dat = eval(['Output_',variables{j}]);                 
            try
                outputModel{j,i} = dat.(subjectPrefix{i});          % Gets subjects struct data
            catch
                outputModel{j,i} = [];                              % Empty cell if variable was not found for subject
            end
        end
    end
        
    Output       = array2table(outputModel,'VariableNames',subjectPrefix,'RowNames',string(variables));     % Creates table from cell array above 
    [vartype{1:length(study_design)}] = deal("double"); 

    var_list = {'HRV', 'CMJ', 'SQV', 'DYNO'};
    num = find(ismember(var_list,  variables));
    % Looping through participants found in table 'Output'
    for jj = 1:width(Output) 
        
        clearvars overPhasesHRV overPhasesCMJ overPhasesDYNO affiliationHelperHRV affiliationHRV affiliationHelperCMJ affiliationCMJ...
                  affiliationHelperDYNO affiliationDYNO matchStudydesign matchStudydesignSelection nanTable   
        
        % Looping through found variable groups
        for i = 1:height(Output)
            
            % Grid to navigate incoming file
            try
                % Checks if cell arrays row name of present input is according to variable manipulation below and additionally not empty 
                if(string(Output(i,1).Properties.RowNames) == 'HRV' && isempty(Output(i,1)) == 0)                   
                    num_HRV = find(var_list == "HRV");
                    nanTable_HRV = array2table(nan(length([output_variables{num_HRV}]),length(study_design)),...
                                   'VariableNames',study_design,...
                                   'RowNames',string([output_variables{num_HRV}]));
                    
                    % Flipping table with phases as colums and variable as rows for later group summary over phases
                    overPhasesHRV = array2table(table2array(Output(i,jj).(string(Output(i,jj).Properties.VariableNames)){1,1}.HRV(1:end-2,:)).',...
                                                'VariableNames',Output(i,jj).(string(Output(i,jj).Properties.VariableNames)){1,1}.HRV(1:end-2,:).Properties.RowNames,...
                                                'RowNames',Output(i,jj).(string(Output(i,jj).Properties.VariableNames)){1,1}.HRV(1:end-2,:).Properties.VariableNames);
                                       
                    % Compares each name of phases contained in subjects file description to phases of the study design.
                    % Creates matrix with indecies of matched strings where each row represents study design phases and columns the compared phases
                    affiliationHelperHRV = {};
                    affiliationHRV = [];
                    for ii = 1:length(study_design)                        
                        if(length(string(overPhasesHRV.Properties.VariableNames)) == 1)
                            affiliationHelperHRV{ii,1} = strfind(string(overPhasesHRV.Properties.VariableNames),study_design(ii));
                        else
                            affiliationHelperHRV(ii,:) = strfind(string(overPhasesHRV.Properties.VariableNames),study_design(ii));
                        end                        
                    end
                    
                    % Converting cell array with indecies to matrix with integers (logical values)
                    for ii = 1:size(affiliationHelperHRV,1)
                        for j = 1:size(affiliationHelperHRV,2)
                          affiliationHRV(ii,j) = ~isempty(affiliationHelperHRV{ii,j});
                        end
                    end
                    
                    % Getting indecies from true phase comparisons. Replacing placeholder with subject data in specific phases
                    % Note: If number of variables has changed, 'nanTable(XX:length(...)' needs to be adjusted
                    [matchStudydesign,matchSubject] = find(affiliationHRV == 1);
                    nanTable_HRV(:,matchStudydesign) = overPhasesHRV(:,matchSubject);
                
                elseif(exist('nanTable_HRV', 'var') == 0)
                    nanTable_HRV = [];
                end
                
            catch
            end
                        
            % Grid to navigate incoming file
            try
                % Checks if cell arrays row name of present input is according to variable manipulation below and additionally not empty 
                if(string(Output(i,1).Properties.RowNames) == 'CMJ' && isempty(Output(i,1)) == 0 )
                    
                    num_CMJ = find(var_list == "CMJ");
                    nanTable_CMJ = array2table(nan(length([output_variables{num_CMJ}]),length(study_design)),...
                                   'VariableNames',study_design,...
                                   'RowNames',string([output_variables{num_CMJ}]));
                    
                    
                    % Flipping table with phases as colums and variable as rows for later group summary over phases
                    overPhasesCMJ = array2table(table2array(Output(i,jj).(string(Output(i,jj).Properties.VariableNames)){1,1}.Best_Trials(1:end-2,:)).',...
                                                'VariableNames',Output(i,jj).(string(Output(i,jj).Properties.VariableNames)){1,1}.Best_Trials(1:end-2,:).Properties.RowNames,...
                                                'RowNames',Output(i,jj).(string(Output(i,jj).Properties.VariableNames)){1,1}.Best_Trials(1:end-2,:).Properties.VariableNames);
                   
                   % Compares each name of phases contained in subjects file description to phases of the study design.
                   % Creates matrix with logicals where each row represents study design phases and columns the compared phases       
                   affiliationHelperCMJ = {};
                   affiliationCMJ = [];
                   for ii = 1:length(study_design)                        
                        if(length(string(overPhasesCMJ.Properties.VariableNames)) == 1)
                            affiliationHelperCMJ{ii,1} = strfind(string(overPhasesCMJ.Properties.VariableNames),study_design(ii));
                        else
                            affiliationHelperCMJ(ii,:) = strfind(string(overPhasesCMJ.Properties.VariableNames),study_design(ii));
                        end                       
                    end
                    % Converting cell array with indecies to matrix with integers (logical values)
                    for ii = 1:size(affiliationHelperCMJ,1)
                        for j = 1:size(affiliationHelperCMJ,2)
                          affiliationCMJ(ii,j) = ~isempty(affiliationHelperCMJ{ii,j});
                        end
                    end

                    % Getting indecies from true phase comparisons. Replacing placeholder with subject data in specific phases
                    % Note: If number of variables has changed, 'nanTable(XX:length(...)' needs to be adjusted
                    [matchStudydesign,matchSubject] = find(affiliationCMJ == 1);
%                     nanTable_CMJ((14:13+length(overPhasesCMJ{:,1})),matchStudydesign) = overPhasesCMJ(:,matchSubject);
                    nanTable_CMJ(:,matchStudydesign) = overPhasesCMJ(:,matchSubject);
                elseif(exist('nanTable_CMJ', 'var') == 0)
                    nanTable_CMJ = [];
                end
                               
            catch
            end
            
             % Grid to navigate incoming file
            try
                % Checks if cell arrays row name of present input is according to variable manipulation below and additionally not empty 
                if(string(Output(i,1).Properties.RowNames) == 'SQV' && isempty(Output(i,1)) == 0 )                    
                    num_SQV = find(var_list == "SQV");
                    nanTable_SQV = array2table(nan(length([output_variables{num_SQV}]),length(study_design)),...
                                   'VariableNames',study_design,...
                                   'RowNames',string([output_variables{num_SQV}]));
                    
                    % Flipping table with phases as colums and variable as rows for later group summary over phases
                    overPhasesSQV = array2table(table2array(Output(i,jj).(string(Output(i,jj).Properties.VariableNames)){1,1}.Best_Trials(1:end-2,:)).',...
                                                'VariableNames',Output(i,jj).(string(Output(i,jj).Properties.VariableNames)){1,1}.Best_Trials(1:end-2,:).Properties.RowNames,...
                                                'RowNames',Output(i,jj).(string(Output(i,jj).Properties.VariableNames)){1,1}.Best_Trials(1:end-2,:).Properties.VariableNames);
                   
                   % Compares each name of phases contained in subjects file description to phases of the study design.
                   % Creates matrix with logicals where each row represents study design phases and columns the compared phases
                   affiliationHelperSQV = {};
                   affiliationSQV = [];
                   for ii = 1:length(study_design)
                        if(length(string(overPhasesSQV.Properties.VariableNames)) == 1)
                            affiliationHelperSQV{ii,1} = strfind(string(overPhasesSQV.Properties.VariableNames),study_design(ii));
                        else
                            affiliationHelperSQV(ii,:) = strfind(string(overPhasesSQV.Properties.VariableNames),study_design(ii));
                        end                       
                    end
                    % Converting cell array with indecies to matrix with integers (logical values)
                    for ii = 1:size(affiliationHelperSQV,1)
                        for j = 1:size(affiliationHelperSQV,2)
                          affiliationSQV(ii,j) = ~isempty(affiliationHelperSQV{ii,j});
                        end
                    end

                    % Getting indecies from true phase comparisons. Replacing placeholder with subject data in specific phases
                    % Note: If number of variables has changed, 'nanTable(XX:length(...)' needs to be adjusted
                    [matchStudydesign,matchSubject] = find(affiliationSQV == 1);
                    nanTable_SQV(:,matchStudydesign) = overPhasesSQV(:,matchSubject);
                
                elseif(exist('nanTable_SQV', 'var') == 0)
                    nanTable_SQV = [];
                end
                               
            catch
            end
             
            % Grid to navigate incoming file
            try
                % Checks if cell arrays row name of present input is according to variable manipulation below and additionally not empty 
                if(string(Output(i,1).Properties.RowNames) == 'DYNO' && isempty(Output(i,1)) == 0)                 
                    num_DYNO = find(var_list == "DYNO");
                    nanTable_DYNO = array2table(nan(length([output_variables{num_DYNO}]),length(study_design)),...
                                   'VariableNames',study_design,...
                                   'RowNames',string([output_variables{num_DYNO}]));
                    
                    % Flipping table with phases as colums and variable as rows for later group summary over phases
                    overPhasesDYNO = array2table(table2array(Output(i,jj).(string(Output(i,jj).Properties.VariableNames)){1,1}.Best_Trials(1:end-2,:)).',...
                                            'VariableNames',Output(i,jj).(string(Output(i,jj).Properties.VariableNames)){1,1}.Best_Trials(1:end-2,:).Properties.RowNames,...
                                            'RowNames',Output(i,jj).(string(Output(i,jj).Properties.VariableNames)){1,1}.Best_Trials(1:end-2,:).Properties.VariableNames);

                   % Compares each name of phases contained in subjects file description to phases of the study design.
                   % Creates matrix with logicals where each row represents study design phases and columns the compared phases   
                   affiliationHelperDYNO = {};
                   affiliationDYNO = [];
                    for ii = 1:length(study_design)                       
                        if(length(string(overPhasesDYNO.Properties.VariableNames)) == 1)
                            affiliationHelperDYNO{ii,1} = strfind(string(overPhasesDYNO.Properties.VariableNames),study_design(ii));
                        else
                            affiliationHelperDYNO(ii,:) = strfind(string(overPhasesDYNO.Properties.VariableNames),study_design(ii));
                        end
                    end

                    % Converting cell array with indecies to matrix with integers (logical values)
                    for ii = 1:size(affiliationHelperDYNO,1)
                        for j = 1:size(affiliationHelperDYNO,2)
                          affiliationDYNO(ii,j) = ~isempty(affiliationHelperDYNO{ii,j});
                        end
                    end

                    % Getting indecies from true phase comparisons. Replacing placeholder with subject data in specific phases
                    % Note: If number of variables has changed, 'nanTable(XX:length(...)' needs to be adjusted
                    [matchStudydesign,matchSubject] = find(affiliationDYNO == 1);
                    nanTable_DYNO(:,matchStudydesign) = overPhasesDYNO(:,matchSubject);
                    
                elseif(exist('nanTable_DYNO', 'var') == 0)
                    nanTable_DYNO = [];
                end
                
            catch
            end
            
            nanTable = [nanTable_HRV; nanTable_CMJ; nanTable_SQV; nanTable_DYNO];
            Output{length(variables)+1,jj} = {nanTable};          % Adding subjects summary table to original output struct
            
        end      
        Output.Properties.RowNames{end} = 'Summary';       
    end
        
    % Creating ghost table in size and properties specific to number of subjects
    groupData = table('Size',[length([output_variables{num}])*length(subjectPrefix),length(study_design)+2],...
                      'VariableTypes',['string','string',vartype{:}],...
                      'VariableNames',['Participant','Variable',study_design]); 
                  
    % Extracts summary tables  
    for i = 1:width(Output)       
           dats{i} = [string(repmat(Output.Properties.VariableNames{i},[height(Output.(i){end,1}),1])),...
                      string(Output.(i){end,1}.Properties.RowNames),table2array(Output.(i){end,1})];
    end

    % Create "Test" variable to enhance readability of the output table
    tests = [];
    var_list = {'HRV', 'CMJ', 'SQV', 'DYNO'};
    for i = 1:numel(string(variables))
       num_temp = find(ismember(var_list, string(variables(i))));
       if(i == 1)
            [tests{1:numel(output_variables{num_temp}), 1}] = deal(variables(i, :));
       else
            [tests{end +1:end + numel(output_variables{num_temp}), 1}] = deal(variables(i, :)); 
       end
        
    end
    
    Test = [];
    for i = 1:numel(subjectPrefix)
        Test = [Test; tests];
    end   
    
    % Concatenates summry tables and creates final output 'groupData',
    % as well as table 'individualData' with individual results and raw data 
    dats                = vertcat(dats{:,:});
    groupData{:,1:2}    = dats(:,1:2);
    groupData{:,3:end}  = str2double(dats(:,3:end));   
    groupData = addvars(groupData, string(Test), 'After', 'Participant', 'NewVariableNames', 'Test'); 
    individualData      = Output; 
    analysisData        = DataTransformation(groupData, baseline);
    
    t2 = toc;
    fprintf('%-2.2f seconds\n', t2);
    
catch
     fprintf('\nERROR filling tables (9)\n')
end

end



