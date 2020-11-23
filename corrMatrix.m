function [reproducibilityResults, associationResults] = corrMatrix()
addpath(string(pwd) + '\Functions');



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
alpha               = 0.05;                             % Alpha level, Type 1 Error rate 
NULL                = [   0.875  ,...                   % Reproducibility Null-Hypothesis , Range 0 to 1 (absolute), Inferiority test. Tests if measurement is smaller than NULL
                          0.875     ];                  % Association Null-Hypothesis ,     Range 0 to 1 (absolute), Inferiority test. Tests if measurement is smaller than NULL
side                = 2;                                % Options: 1 = one-sided test, 2 = two-sided test. For one-sided test, the NULL has to be in the proper direction

HR_variable         = 4;                                % Selects the HR variable, to test association to. Options refer to 'outputVariables' index of the first row, e.g. 4 = RMSSD (ms)


%%
[groupData, ~] = DataProcessing(study_design, output_variables, artefact_recognition, artefact_threshold, detrending);
transformedData = DataTransformation(groupData, baseline);

%%
tic
fprintf('\nAnalysis: \t\t\t')
err = {};
mute = [0, 0];
missingReps = [];

% Loop through 2 iterations (Absolute, Delta)
for m = [0,1]
    clearvars -except mute missingReps err includedPhases HR_variable correlation_Type ICCtype ANCOVA_type alpha NULL side transformedData Association Reproducibility m k 
    
    % Loop through 3 iterations (Raw, Percentage, Normalized units)
    for k = 1:3
        
        data                = transformedData.Results{k + 3*m};                                     % get data accoring to m and k
        variables           = string(transformedData.Results{'ABS_mean'}.Properties.RowNames);      % gather names of all variables
        participants        = unique(data.Participant);                                             % identifying participants
                                                       
        study_design        = string(data.Properties.VariableNames);                                % get study phases of table from column names
        targetPhases        = study_design(find(contains(study_design,includedPhases) == 1));       % only keep phases equal to "includedPhases"
        underscores         = strfind(targetPhases,'_');                                            % find delimeter of phase names
        analysisData        = data(:,targetPhases);                                                 % cut data with respect to containing "includedPhases"
        delta_Phases        = [];
        
        %% Special Case Deltas
        spaces = cell2mat(strfind(analysisData.Properties.VariableNames, ' '));                     % find all spaces in phase names, which indicates delta phases
        
        % sum of spaces only > 0 if m = 1, meaning delta data is beeing analyzed
        if(sum(spaces) > 0)     
            
            % loop through all included phases
            for i = 1:numel(targetPhases)
                tPhase_char      = char(targetPhases(i));                                           % conversion from string to char "XX_INT_1 XX_INT_2"
                delta_parts(1,:) = string(tPhase_char(1:spaces(i) - 1));                            % separating first part of phase names "XX_INT_1"
                delta_parts(2,:) = string(tPhase_char(spaces(i) + 1:end));                          % separating second part of phase names "XX_INT_2"
                
                % checks if 2 space-separated parts in phase name exist,  
                % meaning that the phase only includes "includedPhases" and
                % more than 2 components
                if(sum(contains(delta_parts, includedPhases)) == 2)               
                    delta_Phases = [delta_Phases, targetPhases(i)];                                 % concatenates all delta phases           
                end
            end           
            targetPhases = delta_Phases;                                                            % assigns delta phases to target phases            
        end
        % error check for insufficient number of recordings in target phases
        if(isempty(targetPhases))
           break; 
        end
        %% Getting data of interest excluding NaN, sorting and arranging it
        
        % Getting integer from phase name to identify repeated measures design
        for i = 1:numel(targetPhases)
            numRepetition(i,1) = str2num(targetPhases{i}(1:underscores{i}(1)-1));
        end
        
        % Getting number of distinct study repetitions
        num   = unique(numRepetition);
        
        % Getting data of interest for every study repetition
        for i = 1:numel(num)
                phaseIndeces{i} = find(numRepetition == num(i));                                    % Identifies the column/phase in data that matches the repeated measures
                dataOI{i}       = analysisData(:,targetPhases(phaseIndeces{i}));                    % Splits analysisData into subsets of measurement repetitions     
        end

        % Concatenating measurement repetitions for each participant
        % Looping through entire data
        for jj = 1:size(data, 1)
            temp = [];                                                                              
            
            % Looping through number of measurement repetitions
            for i = 1:size(dataOI, 2)
                temp(:, i) = dataOI{:, i}{jj,:};                                                    % Fetches data temporarily
            end
            tempData{jj, 1} = temp;                                                                 % Creates cells containing all valid measurement pairs for each participant
        end
            
        unsplitData = [dataOI{:,:}];                                                                % Creates table of dataOI that contains NaN 
        
        % Arranging dataOI in columns for either reproducibility or association analysis
        % Looping through all variables of initial dataset
        for jj = 1:numel(variables)
            
            % Looping through all participants
            for i = 1:numel(participants)                
                repetitive_data{jj,i} = [ones(size(tempData{jj + numel(variables) * (i-1) ,:}, 1), 1)*i,...         % Concatenating participant identifier and its respective measurements each in one column
                                         tempData{jj + numel(variables) * (i-1) ,:}];
                associative_data{jj,i} = [ones(size(tempData{jj + numel(variables) * (i-1) ,:}(:), 1), 1)*i,...     % Concatenating participant identifier and its respective measurements in two columns
                                         tempData{jj + numel(variables) * (i-1) ,:}(:)];
            end
            targetData{jj,1}        = vertcat(repetitive_data{jj,:});                               % Concatenating the measuremnts of all participants for repetetive data
            targetData{jj,2}        = vertcat(associative_data{jj,:});                              % Concatenating the measuremnts of all participants for associative data            
        end
        
        % Creating tables with respective variable names for later identification and notation
        % Looping through targetData
        for i = 1:length(targetData)
            
            % Control for special case comparing HR_variable to itself
            if(i == HR_variable)
                targetData{i, 3} = [targetData{HR_variable, 2}, targetData{i, 2}(:, 2)]; 
                targetData{i, 3} = array2table(targetData{i, 3}(sum(~isnan(targetData{i, 3}), 2) == 3, :),...        
                                        'VariableNames', ['Subjects', variables(HR_variable),variables(i) + '_2']);         % Creates table with association data and labels it
                targetData{i, 1} = array2table(targetData{i, 1}(sum(~isnan(targetData{i, 1}), 2) == 3, :),...
                                        'VariableNames', ['Subjects', variables(i),variables(i) + '_2']);                   % Creates table with reproducibility data and labels it 
            else
                targetData{i, 3} = [targetData{HR_variable, 2}, targetData{i, 2}(:, 2)];               
                targetData{i, 3} = array2table(targetData{i, 3}(sum(~isnan(targetData{i, 3}), 2) == 3, :),...
                                        'VariableNames', ['Subjects', variables(HR_variable),variables(i)]); 
                targetData{i, 1} = array2table(targetData{i, 1}(sum(~isnan(targetData{i, 1}), 2) == 3, :),...
                                        'VariableNames', ['Subjects', variables(i),variables(i) + '_2']);
            end
        end
        
        % Final tables
        reproducibilityData = table(variables,targetData(:,1));
        associationData     = table(variables, targetData(:,3));

        %% Reproducibility
        % Error catch function
        if(sum(~cellfun(@isempty, reproducibilityData{:, 2})) == 0)
            if(mute(1) == 0)
                err{1} = sprintf('ERROR: Insufficient data for analysis of reproducibility\n\n');
                mute(1) = 1;
            end
        else
            try
                % Looping through all variables
                for i = 1:size(reproducibilityData,1)

                    % Fetching current variable data to table M
                    M       = reproducibilityData{i,2}{:,:};                                           

                    % Differentiating between selected analysis method
                    switch correlation_Type(1)
                        case 1  % Pearson correlation
                            try
                                [r, p, LB, UB] = corrcoef(M{:,2:3}, 'Alpha', alpha);                            % Calling upon function corrcoef()
                                p_NULL = p_approx(r(1,2)*sign(r(1,2)), NULL(1), size(M, 1));                    % Aprroximating p value against NULL with function p_approx()
                                reproducibilityData{i,3:7}  = [r(1, 2), LB(1, 2), UB(1, 2), p(1,2), p_NULL];    % Concatenating results in cells
                            catch
                                missingReps = [missingReps, i];
                            end
                        case 2 % Intra-class correlation
                            try
                                [r, LB, UB, F, df1, df2, p] = ICC(M{:,2:3}, ICCtype, alpha, side, 0);           % Calling upon extern function ICC()
                                p_NULL = p_approx(r, NULL(1), size(M, 1));                                      % Aprroximating p value against NULL with function p_approx()
                                reproducibilityData{i,3:10}  = [r, LB, UB, F, df1, df2, p, p_NULL];             % Concatenating results in cells     
                            catch 
                                missingReps = [missingReps, i];
                            end
                        case 3 % ANCOVA   
                            try
                                names = M.Properties.VariableNames;                                             % Getting column names of table for grouping
                                [~, atab, ctab, ~] = aoctool(M{:,2}, M{:,3}, M{:,1}, alpha, ...                 % Calling upon ANCOVA with aoctool() with selected type
                                                             names{2}, names{3}, names{1},...                   % grouping variable = participants
                                                             'off', ANCOVA_type{1});
                                r = sqrt(sum([atab{3, 3}])/sum([[atab{3,3}], [atab{end, 3}]]))...               % Calculating correlation coefficient by ignoring subject as dummy variable
                                    *sign(ctab{find(strcmp(ctab(:, 1), 'Slope')), 2});                          % Only calculating on the explained variance of target variable
                                F = atab{3, end-1};                                                             % Getting F value for target variable
                                df = atab{end, 2};                                                              % Getting degrees of freedom for target variable
                                p = atab{3, end};                                                               % Getting p value for target value, tested against no effect (p0)
                                [LB, UB] = analyticalCI(r, df, alpha);                                  % Calculating confidence interval by calling upon analyticalCI(), uses Fisher's z transformation
                                p_NULL = p_approx(r*sign(r), NULL(1), df);                              % Approximating p value for specified NULL hypothesis by calling upon p_approx(), using probability density estimates
                                % Adding statistics tables to output
                                reproducibilityData{i,3:4}   = {array2table(atab(2:end,:), 'VariableNames', atab(1,:)),...
                                                                array2table(ctab(2:end,:), 'VariableNames', ctab(1,:))};
                                reproducibilityData{i,5:11}  = [r, LB, UB, F, df, p, p_NULL];                   % Concatenating results in cells   
                            catch 
                                missingReps = [missingReps, i];
                            end
                    end  
                end

                % Selecting names depending on selected correlation type
                switch correlation_Type(1)
                    case 1
                        varNames = {'Variable','ReproducibilityData','r','LB','UB','p(0)',['p(', num2str(NULL(1)), ')']};
                    case 2
                        varNames = {'Variable','ReproducibilityData','r','LB','UB','F','df1','df2','p(0)',['p(<', num2str(NULL(1)), ')']};
                    case 3
                        varNames = {'Variable','ReproducibilityData', 'AnovaTable', 'CovariationTable', ...
                                    'r','LB','UB','F','df','p(0)',['p(', num2str(NULL(1)), ')']};
                end

                % Creating output of for reproducibility results
                reproducibilityData.Properties.VariableNames    = varNames;                             
                Reproducibility{k, m + 1}                       = reproducibilityData;

            catch
                % Returns error message in case data has not the required size for the tests selected
                if(isempty(M) == 1)
                    err{2} = sprintf('ERROR: No repeated measures found for variables:\n\n');
                end
            end
        end
        %% Individual Association      
        % Error catch function
        if(sum(~cellfun(@isempty, associationData{:, 2})) == 0)
            if(mute(2) == 0)
                err{3} = sprintf('ERROR: Insufficient data for analysis of association\n\n');
                mute(2) = 1;
            end
        else
            try
                % Looping through all variables
                for i = 1:size(associationData,1)

                    % Fetching current variable data to table M
                    M       = associationData{i,2}{:,:};                                           

                    % Differentiating between selected analysis method
                    switch correlation_Type(1)
                        case 1  % Pearson correlation
                            try
                                [r, p, LB, UB] = corrcoef(M{:,2:3}, 'Alpha', alpha);                            % Calling upon function corrcoef()
                                p_NULL = p_approx(r(1,2)*sign(r(1,2)), NULL(1), size(M, 1));                    % Aprroximating p value against NULL with function p_approx()
                                associationData{i,3:7}  = [r(1, 2), LB(1, 2), UB(1, 2), p(1, 2), p_NULL];       % Concatenating results in cells
                            catch
                            end
                        case 2 % Intra-class correlation
                            try
                                [r, LB, UB, F, df1, df2, p] = ICC(M{:,2:3}, ICCtype, alpha, side, 0);           % Calling upon extern function ICC()
                                p_NULL = p_approx(r, NULL(1), size(M, 1));                                      % Aprroximating p value against NULL with function p_approx()
                                associationData{i,3:10}  = [r, LB, UB, F, df1, df2, p, p_NULL];                 % Concatenating results in cells      
                            catch
                            end
                        case 3 % ANCOVA 
                            try
                                names = M.Properties.VariableNames;                                             % Getting column names of table for grouping
                                [~, atab, ctab, ~] = aoctool(M{:,2}, M{:,3}, M{:,1}, alpha, ...                 % Calling upon ANCOVA with aoctool() with selected type
                                                             names{2}, names{3}, names{1},...                   % grouping variable = participants
                                                             'off', ANCOVA_type{1});
                                r = sqrt(sum([atab{3, 3}])/sum([[atab{3,3}], [atab{end, 3}]]))...               % Calculating correlation coefficient by ignoring subject as dummy variable
                                    *sign(ctab{find(strcmp(ctab(:, 1), 'Slope')), 2});                          % Only calculating on the explained variance of target variable
                                F = atab{3, end-1};                                                             % Getting F value for target variable
                                df = atab{end, 2};                                                              % Getting degrees of freedom for target variable
                                p = atab{3, end};                                                               % Getting p value for target value, tested against no effect (p0)
                                [LB, UB] = analyticalCI(r, df, alpha);                                  % Calculating confidence interval by calling upon analyticalCI(), uses Fisher's z transformation
                                p_NULL = p_approx(r*sign(r), NULL(2), df);                              % Approximating p value for specified NULL hypothesis by calling upon p_approx(), using probability density estimates
                                % Adding statistics tables to output
                                associationData{i,3:4}   = {array2table(atab(2:end,:), 'VariableNames', atab(1,:)),...
                                                                array2table(ctab(2:end,:), 'VariableNames', ctab(1,:))};
                                associationData{i,5:11}  = [r, LB, UB, F, df, p, p_NULL];                   % Concatenating results in cells       
                            catch 
                                aa = 1;
                            end
                    end  
                end

                % Selecting names depending on selected correlation type
                switch correlation_Type(1)
                    case 1
                        varNames = {'Variable','AssociationData','r','LB','UB','p(0)',['p(', num2str(NULL(1)), ')']};
                    case 2
                        varNames = {'Variable','AssociationData','r','LB','UB','F','df1','df2','p(0)',['p(<', num2str(NULL(1)), ')']};
                    case 3
                        varNames = {'Variable','AssociationData', 'AnovaTable', 'CovariationTable', ...
                                    'r','LB','UB','F','df','p(0)',['p(', num2str(NULL(1)), ')']};
                end

                % Creating output of for association results
                associationData.Properties.VariableNames    = varNames;                             
                Association{k, m + 1}                       = associationData;

            catch
                % Returns error message in case data has not the required size for the tests selected
                if(numel(unique(M.Subjects)) < 2 )
                    err{4} = sprintf('ERROR: Not enough participants for ANCOVA (association)\n\n');
                end
            end
        end
    end
end

%% Creating output
% Error catch function if analysis has not been conducted
try
    reproducibilityResults = table(Reproducibility(:,1), Reproducibility(:,2),...
                                   'VariableNames', {'Absolute', 'Delta'},...
                                   'RowNames',      {'Raw', 'Percentage', 'Effect Size'});
catch
    reproducibilityResults = [];
end

try
    associationResults     = table(Association(:,1), Association(:,2),...
                                   'VariableNames', {'Absolute', 'Delta'},...
                                   'RowNames', {'Raw', 'Percentage', 'Effect Size'});
catch
    associationResults = [];
end

t5 = toc;
fprintf('%4.2f seconds\n\n', t5);
for e = find(~cellfun(@isempty, err) == 1)
    if(e == 2)
        fprintf(err{2});
        fprintf('%s \n\n', num2str(unique(missingReps)));
    else
        fprintf(err{e});
    end
end

end