function [reproducibilityResults, associationResults] = DataAnalysis(transformedData, includedPhases, HR_variable, correlation_Type, ICCtype, ANCOVA_type, alpha, NULL, side)
% This function statistically analyzes the transformed data w.r.t. the input parameters.
%
% USAGE
% E.g.:
% [reproducibilityResults, associationResults] = DataAnalysis(transformedData, includedPhases, 4, 3, 'C-1', 'parallel', 0.05, [0.9, 0.9], 1)
%
% INPUT
% transformedData:          Output table (18x1) from function 'DataProcessing'.
% includedPhases:           String array that determines which columns of the input table 
%                           will be analyzed. Must at least partly match phases declared 
%                           in 'study_design' used for 'DataProcessing'. Case sensitive.
%                           E.g.: {'INT'}, study_design = {...,'1_INT_1',...}
%
% HR_variable:              Scalar that determines which variable will be correlated to all
%                           variables found in the input table for analysis of association.
%                           Refers to the variable names declared in 'outputVariables' used in
%                           'DataProcessing'.
%
% correlation_Type:         Scalar that determines which analysis method will be used. 
%                           Options: 1 = Pearson, 2 = Intraclass correlation, 3 = Analysis of
%                           covariance.
%
% ICC_Type:                 If correlation_Type = 2, string that determines the ICC to be used.
%                           Options: '1-1', 'A-1', 'C-1', '1-k', 'A-k', 'C-k'.
%                           For further information consult function 'ICC'.
%
% ANCOVA_type:              String that determines the analysis of covariance type to be used.
%                           Options: 'separate line', creates individual fit
%                           for every participant, 'parallel lines', creates one fit for all
%                           participants with individual intercepts.
%
% alpha:                    Double that determines the alpha error rate for significance tests.
%
% NULL:                     Double (2x1) that determines the null hypotheses to be tested.
%                           Refers to correlation coefficient r. First value
%                           for reproducibility, second for association.
%
% side:                     Scalar that sets if hypotheses will be tested one- or two-sided.
%                           Options: 1 = one-sided, 2 = two-sided.
%
% OUTPUT
% reproducibilityResults:   Output table (3x2) that stores reproducibility analysis results for 
%                           raw data, percentages and effect sizes in rows with results for 
%                           absolute and delta values in columns. Each field contains the 
%                           statistics for the specific combination of row and column, as well
%                           as analyzed data, and statistic tables from ANCOVA if chosen.
%
% associationResults:       Output table (3x2) that stores association analysis results for 
%                           raw data, percentages and effect sizes in rows with results for 
%                           absolute and delta values in columns. Each field contains the 
%                           statistics for the specific combination of row and column, as well
%                           as analyzed data, and statistic tables from ANCOVA if chosen.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 09.July.2020
% lepremiere
%---------------------------------------------------------------------------------------------------
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
