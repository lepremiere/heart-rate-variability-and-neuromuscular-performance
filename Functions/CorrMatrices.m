function corrMatrices = CorrMatrices(CA, includedPhases, type)
% This function cross-correlates all variables from test category "HRV" with every
% other variable available from tables in a cell array. Repeated masures correlation 
% is utilized for every participant. 
%
% USAGE
% E.g.:
% corrMatrices = CorrMatrices(transformedData, {'INT', 'REC'}, 'parallel lines')
%
% INPUT
% CA:       Input must be a cell array (6 x 1) containing tables with
%           results of function "DataTransformation". The tables must
%           include variables: "Participant", "Test", "Variable", and the
%           measurements taken. E.g.: "1_BASE_1". Number of measurements
%           and participants must be greater than 2. Additionally, the
%           tables must include the test "HRV" and at least one other.
%
% type:     String input to chose options. "parallel lines" for fitting
%           individual y-intercepts but common slopes for each participant
%           to determine correlation coefficient. "seperate lines" for fitting 
%           individual y-intercepts as well as slopes for each participant.
%
% OUTPUT
% Results:  Cell array (6 x 1) of tables with the correlation coefficients
%           for each "HRV" variable with every non-"HRV" variable and each
%           input table.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 09.July.2020
% lepremiere
%---------------------------------------------------------------------------------------------------

% Entry condition. At least tests "HRV" and one other must be included 
if(numel(unique(CA{1,1}{:}.Test)) > 1)

    % Picking all individual absolute and delta tables from transformed
    % data output
    datasets = {CA{1,1},...
                CA{4,1},...
                CA{2,1},...
                CA{5,1},...
                CA{3,1},...
                CA{6,1}};
    
    % Looping through the tables
    for k = 1:length(datasets)

        data = datasets{k}{:};
        new_data = [];
        [h, w] = size(data{:,4:end});
        
        num_study_design = sum(contains(data.Properties.VariableNames(4:end), includedPhases)); % Identyfing number of measurements that match "includedPhases"
        num_participants = size(unique(data.Participant),1);                                    % Identyfing number of participants
        num_variables = size(unique(data.Variable),1);                                          % Identyfing number of variables
        num_hrv = sum(data.Test == 'HRV')/num_participants;                                     % Identyfing number of variables of test "HRV"
        participant_indexes = [];
        out = [];
        
        % Preparing array with participant identefiers for later processing
        for i = 1:size(unique(data.Participant),1)
            participant_indexes = [participant_indexes; ones(num_study_design, 1)*i]; 
        end
        
        % Extracting values of each participant and test occasion and
        % arranging it to matrix
        for i = 1:num_variables 
            indexes = data.Variable == data.Variable(i);                                        % Extracting all data of iterated variable
            A = data(indexes, [1, 3:end]);
            data_storage{i,1} = ...                                                             % Extracting the measurements that match "includedPhases"
                A(:, [true,true, contains(A.Properties.VariableNames(3:end), includedPhases)]);
            
            % Rearanging results to matrix 
            for j = 1:num_participants  
                new_data = [new_data; [data_storage{i}{j, 3:end}']];   
            end
            out = [out, new_data];
            new_data = [];
        end

        output = array2table([participant_indexes, out], 'VariableNames', {'Participant', data.Variable{1:num_variables}});
        
        % Looping through all "HRV" variables
        for i = 1:num_hrv
            % Looping through all non-"HRV" variables
            for j = num_hrv:num_variables
               if(i ~= j)
                    data_share = [output(:,1), output(:,1+i), output(:,1+j)];                   % Getting participant identifier and data pair for correlation
                    nan_excludes = sum(~isnan(data_share{:,2:end}),2) > 1;                      % Excluding rows that contain NaN
                    M = data_share(nan_excludes, :);
                    
                    % If only one participant can be identified, regular
                    % Pearson correlation is performed. Else analysis of
                    % covariation.
                    if(size(unique(M.Participant),1) <= 1)
                        r = corrcoef(M{:,2}, M{:,3});                                           % Pearson correlation
                        r = r(2,1);
                    else
                        names = M.Properties.VariableNames;                                     % Getting column names of table for grouping
                        [~, atab, ctab, ~] = aoctool(M{:,2}, M{:,3}, M{:,1}, 0.05, ...          % Calling upon ANCOVA with aoctool() with selected type
                                                     names{2}, names{3}, names{1},...           % grouping variable = participants
                                                     'off', type);
                        r = sqrt(sum([atab{3, 3}])/sum([[atab{3,3}], [atab{end, 3}]]))...       % Calculating correlation coefficient by ignoring subject as dummy variable
                            *sign(ctab{find(strcmp(ctab(:, 1), 'Slope')), 2});                  % Only calculating on the explained variance of target variable
                    end
                    R(i,j) = r;    
               end  
            end 
        end
        
        % Rearranging the correlation result matrix and creating table
        % "num_hrv-+1" excludes variable "Num Artefacts"
        matrix = round(R(1:num_hrv-1,num_hrv+1:end), 3);
        corr_matrix = array2table(matrix,...
                    'VariableNames', data.Variable(num_hrv+1:num_variables),...
                    'RowNames', data.Variable(1:num_hrv-1));
        
        % Saving temporary final result in cell array
        final_output{k,1} = corr_matrix;
    end
    
    % Creating output cell array that matches input 
    corrMatrices = array2table(final_output,...
        'VariableNames', {'CorrelationMatrix'},...
        'RowNames', {CA.Properties.RowNames{[1,4,2,5,3,6]}});

else
    % Error catch function
    corrMatrices = [];
    fprintf('\nERROR: Insufficient Data\n')
end
end

