function Export(export_options)
% This function exports the results produced by "DataAnalysis.m" and
% "CorrMatrices.m". It is fitted towards that purpose and can hardly be
% used under different circumstances.
%
% USAGE
% E.g.:
% Export(export_options)
%
% INPUT
% The required variables will be called from the workspace.
%
% OUTPUT
% The function outputs figures consisting of the calculated effect sizes of
% each individual and a group mean, graphical illustrations of the
% statistics as well as correlation matrices. Additionally, tables of a
% variaty of data and the workspace variables "individualData" as well as
% "analysisData" will be saved. All outputs will be saved to the folder
% chosen by the ui-window in the beginning.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 19.January.2021
% lepremiere
%---------------------------------------------------------------------------------------------------

% Required variables will be drawn from workspace if any output is wanted
if(sum(export_options) > 0)
    tic
    fprintf('Export: \t\t\t')

    selpath                 = uigetdir;                                 % Choosing output folder
    groupData               = evalin('base', 'groupData');
    individualData          = evalin('base', 'individualData');
    reproducibilityResults  = evalin('base', 'reproducibilityResults');
    associationResults      = evalin('base', 'associationResults');
    corrMatrices            = evalin('base', 'corrMatrices');
    analysisData            = evalin('base', 'analysisData');
    correlation_Type        = evalin('base', 'correlation_Type');
    ANCOVA_type             = evalin('base', 'ANCOVA_type');
    xlabels                 = evalin('base', 'xlabels');
    disp_vars               = evalin('base', 'disp_vars');
    output_variables        = evalin('base', 'output_variables');
    show_baseline           = evalin('base', 'show_baseline');
    baseline                = evalin('base', 'baseline');
    includedPhases          = evalin('base', 'includedPhases');
    study_design            = evalin('base', 'study_design');
    data                    = evalin('base', 'data');
    alpha                   = evalin('base', 'alpha');
    
    % Error catch function
    if(selpath == 0)
        fprintf('\nNo folder selected\n');
        return
    end

end

% Setting name prefixes of different purposes. Iinitializes invisible figure
namePrefix = {'Reproducibility', 'Association'};
fig = figure('visible', 'off', 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1], 'color', 'w'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section generates the subfigures of FIGURE 3 in the paper. The subfigures were then cut    %
% together to produce FIGURE 3.                                                                   %
% FOLDER: "\Figures\Statistics\2_Association RMSSD\"                                              %
% FILES:  "15_Association_RMSSD_Jump Height Impulse.png",                                         %
%         "24_Association_RMSSD_Mean Propulsive Velocity Squat.png",                              %
%         "31_Association_RMSSD_Resting Twitch.png"                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exporting graphical illustration for the statistics
if(export_options(3) == 1) 
    % Initial condition if reproducibility data is present
    if(size(reproducibilityResults, 1) > 0)
        try
            % Folders for reproducibility data
            dirs  = selpath + "\Figures\Statistics\1_Reproducibility\";
            ind = 1;
            mkdir(dirs);
            
            % Iterative plotting of graphical illustration 
            for i = 1:height(reproducibilityResults.Absolute{1,1})
                try
                    % Assigning cell array with all results of interest
                    M = {reproducibilityResults.Absolute{1,1}{i,2}{:,:},...
                         reproducibilityResults.Delta{1,1}{i,2}{:,:},...
                         reproducibilityResults.Absolute{2,1}{i,2}{:,:},...
                         reproducibilityResults.Delta{2,1}{i,2}{:,:},...
                         reproducibilityResults.Absolute{3,1}{i,2}{:,:},...
                         reproducibilityResults.Delta{3,1}{i,2}{:,:}};
                    
                    % Assigning cell array with t-tables for analysis of
                    % covariation
                    if(correlation_Type(1) == 3) 
                        idx = [5, 9, 11]; % Setting indeces for t-tables beeing present
                        T = {reproducibilityResults.Absolute{1,1}{i,4}{:,:},...
                             reproducibilityResults.Delta{1,1}{i,4}{:,:},...
                             reproducibilityResults.Absolute{2,1}{i,4}{:,:},...
                             reproducibilityResults.Delta{2,1}{i,4}{:,:},...
                             reproducibilityResults.Absolute{3,1}{i,4}{:,:},...
                             reproducibilityResults.Delta{3,1}{i,4}{:,:}};
                    elseif(correlation_Type(1) == 2)
                        idx = [3, 9, 11];   % Setting indeces for t-tables not beeing present
                        T = [];
                    else
                        idx = [3, 6, 8];    % Setting indeces for t-tables not beeing present
                        T = [];
                    end
                    
                    % Assigning cell array with the statistics of interest
                    stats   = {reproducibilityResults.Absolute{1,1}{i, idx};...
                               reproducibilityResults.Delta{1,1}{i, idx};...
                               reproducibilityResults.Absolute{2,1}{i, idx};...
                               reproducibilityResults.Delta{2,1}{i, idx};...
                               reproducibilityResults.Absolute{3,1}{i, idx};...
                               reproducibilityResults.Delta{3,1}{i, idx}};

                    % Getting the null hypothesis specified in the
                    % options secton of "Only_Script_You_Need_To_Care_About.m"
                    h0      = reproducibilityResults.Absolute{1,1}.Properties.VariableNames{idx(end)};
                    idx1    = find(h0 == '(',1,'last');
                    idx2    = find(h0 == ')',1,'first');
                    h0      = ['p', h0(idx1:idx2)];
                                        
                    % Plotting the statistical results
                    plot_stats_LM(fig, M, T, stats, h0, namePrefix{1}, correlation_Type(1), ANCOVA_type{1}, 1);
                    drawnow
                    % Getting the image and saving it 
                    imName = [num2str(i), '_', namePrefix{1}, '_',...
                              strrep(M{1}.Properties.VariableNames{2}(1:find(M{1}.Properties.VariableNames{2} == '[') - 2), '/', '_'),...
                              '.png'];
                    print(gcf,  fullfile(dirs, imName),'-dpng', '-r300');
                    ind = ind + 1;
                catch
                end
            end
        catch
        end
    else
        % Suppressing error message if no reproducibility study design is
        % found
        reps = strfind(study_design, '_');
        reps = vertcat(reps{:});
        reps = reps(:,1);
        
        for i = 1:length(reps)
            num = study_design{i}(reps(i)-1);
        end
        
        if(numel(unique(num)) > 1)
            fprintf('\nNo reproducibility data found!\n');
        end
    end

    % Initial condition if association data is present
    if(size(associationResults, 1) > 0)
        try
            % Folders for reproducibility data
            var = associationResults.Absolute{1,1}.AssociationData{1,1}.Properties.VariableNames{2};
            dirs  = selpath + "\Figures\Statistics\2_Association " + var(1:strfind(var, '[')-2);
            ind = 1;
            mkdir(dirs);
            
            % Iterative plotting of graphical illustration 
            for i = 1:height(associationResults.Absolute{1,1})
                try
                    % Assigning cell array with all results of interest
                    M = {associationResults.Absolute{1,1}{i,2}{:,:},...
                         associationResults.Delta{1,1}{i,2}{:,:},...
                         associationResults.Absolute{2,1}{i,2}{:,:},...
                         associationResults.Delta{2,1}{i,2}{:,:},...
                         associationResults.Absolute{3,1}{i,2}{:,:},...
                         associationResults.Delta{3,1}{i,2}{:,:}};
                     
                    % Assigning cell array with t-tables for analysis of
                    % covariation
                    if(correlation_Type(2) == 3) 
                        idx = [5, 9, 11]; % Setting indeces for t-tables beeing present
                        T = {associationResults.Absolute{1,1}{i,4}{:,:},...
                             associationResults.Delta{1,1}{i,4}{:,:},...
                             associationResults.Absolute{2,1}{i,4}{:,:},...
                             associationResults.Delta{2,1}{i,4}{:,:},...
                             associationResults.Absolute{3,1}{i,4}{:,:},...
                             associationResults.Delta{3,1}{i,4}{:,:}};
                    elseif(correlation_Type(2) == 2)
                        idx = [3, 9, 11];   % Setting indeces for t-tables not beeing present
                        T = [];
                    else
                        idx = [3, 6, 8];    % Setting indeces for t-tables not beeing present
                        T = [];
                    end
                    
                    % Assigning cell array with the statistics of interest
                    stats   = {associationResults.Absolute{1,1}{i, idx};...
                               associationResults.Delta{1,1}{i, idx};...
                               associationResults.Absolute{2,1}{i, idx};...
                               associationResults.Delta{2,1}{i, idx};...
                               associationResults.Absolute{3,1}{i, idx};...
                               associationResults.Delta{3,1}{i, idx}};
                           
                    % Getting the null hypothesis specified in the
                    % options secton of "Only_Script_You_Need_To_Care_About.m"
                    h0      = associationResults.Absolute{1,1}.Properties.VariableNames{idx(end)};
                    idx1    = find(h0 == '(',1,'last');
                    idx2    = find(h0 == ')',1,'first');
                    h0      = ['p', h0(idx1:idx2)];
                    
                    % Plotting the statistical results
                    plot_stats_LM(fig, M, T, stats, h0, namePrefix{2}, correlation_Type(2), ANCOVA_type{2}, 2);
                    drawnow
                    % Getting the image and saving it
                    imName = [num2str(i), '_', namePrefix{2},'_',...
                                 [strrep(M{1}.Properties.VariableNames{2}(1:find(M{1}.Properties.VariableNames{2} == '[') - 2), '/', '_'),'_'...
                                  strrep(M{1}.Properties.VariableNames{3}(1:find(M{1}.Properties.VariableNames{3} == '[') - 2), '/', '_')],...  
                                 '.png'];
                    print(fullfile(dirs, imName),'-dpng', '-r300');
                catch
                end                    
            end

        % Folder for association matrices
        dirs  = selpath + "\Figures\Statistics\3_Association Matrices\";
        ind = 1;
        mkdir(dirs);
        
        % Iterative plotting and saving of all possible correlation matrices
        for i = 1:numel(corrMatrices)
            plot_corr_matrix(fig, corrMatrices{i, :} {:});
            drawnow
            imName = num2str(i) + "_Association Matrix_" + corrMatrices.Properties.RowNames{i} + ".png";
            print(gcf, fullfile(selpath + "\Figures\Statistics\3_Association Matrices\", imName),'-dpng', '-r300');
        end
        catch
        end
    else
        fprintf('\nERROR: No association data found!\n');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section generates the subfigures of FIGURE 2 in the paper. The subfigures were then cut    %
% together to produce FIGURE 2.                                                                   %
% FOLDER: "\Figures\Individuals\                                                                  %
% FILES: "A.png", "B.png", "C.png", "D.png"                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exporting graphical illustration of effect sizes for each individual and group mean
if(export_options(4) == 1)
    try
        % Setting name prefixes of different purposes. Iinitializes
        % invisible figure and creating output folder
        fig = figure('visible', 'off', 'Units', 'Normalized',...
                     'OuterPosition', [0, 0, 0.5, 1], 'color', 'w');
        participants = unique(analysisData.Results{1, 1}.Participant);
        dirs  = selpath + "\Figures\Individuals\";
        mkdir(dirs);
        
        % Iterative plotting and saving for each individual
        for i = 1:numel(participants)  
            plot_individuals(fig, analysisData, participants(i), output_variables,...
                             disp_vars, xlabels, includedPhases, show_baseline, baseline, alpha);

            imName = [participants{i}, '.png'];
            print(gcf, fullfile(dirs, imName),'-dpng', '-r300');
        end
            
        % Plotting and saving of group mean   
        plot_individuals(fig, analysisData, 'mean', output_variables,...
                        disp_vars, xlabels, includedPhases, show_baseline, baseline);
        imName = ['Group Mean', '.png'];
        print(gcf, fullfile(selpath + "\Figures\", imName),'-dpng', '-r300');
    catch
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exporting data 
if(export_options(1) == 1)
    try
        % Creating output folder
        dirs = selpath + "\Data\";
        mkdir(dirs);
        
        % Exporting "individualData" and "analysisData" from workspace
        save(selpath + "\Data\" + "individualData.mat", 'individualData');
        save(selpath + "\Data\" + "analysisData.mat",   'analysisData');
        save(selpath + "\Data\" + "importableData.mat", 'data');
    catch 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section generates the tables that contain statistics presented in TABLE 2 of the paper.    %
% TABLE 2 contains the values of following excel spreadsheets:                                    %
% FOLDER: "\Tables\2_AnalysisData\"                                                               %
% FILES: "7_ABS_mean.xlsx", "8_ABS_sd.xlsx", "11_ES_mean.xlsx", "12_ES_sd.xlsx"                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exporting tables
if(export_options(2) == 1)
    try
        % Creating output folder
        dirs = selpath + "\Tables\";
        mkdir(dirs);
        
        % Exporting groupData
        varNames = groupData.Properties.VariableNames;
        outputTable = table2cell(groupData);
        idx = ismissing(groupData);
        outputTable(idx) = {' '};
        writecell(vertcat(varNames, outputTable), fullfile(dirs, 'groupData.xlsx'));  
        
        % Creating subfolder for data tables
        dirs  = selpath + "\Tables\2_AnalysisData\";
        mkdir(dirs);
        
        % Exporting all tables available in "analysisData"
        for k = 1:size(analysisData, 1)
            name = analysisData.Properties.RowNames(k);
            varNames = analysisData.Results{k, 1}.Properties.VariableNames;
            outputTable = table2cell(analysisData.Results{k, 1});
            idx = ismissing(analysisData.Results{k, 1});
            outputTable(idx) = {' '};
            writecell(vertcat(varNames, outputTable),...
                      fullfile(dirs, [num2str(k), '_', name{:}, '.xlsx'])); 
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section generates the tables that contain statistics presented in TABLE 3 of the paper.    %
% TABLE 3 contains the values of following excel spreadsheets:                                    %
% FOLDER: "\Tables\3_Statistics\2_Association\"                                                   %
% FILES: "1_Absolute_Raw.xlsx", "2_Absolute_Percentage.xlsx", "3_Absolute_Effect Size.xlsx"       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dr = {"1_Reproducibility\", "2_Association\"};       % Output folder names for statistics
        data = {reproducibilityResults, associationResults}; 
        
        % Looping through reproducibility and association data
        for n = 1:2
            count = 1;
            % Selecting and creating folder of interest
            dirs  = selpath + "\Tables\3_Statistics\" + dr{n}; 
            mkdir(dirs);                                      
            
            % Looping through statistics for raw and delta statistics
            for k = 1:2
                % Looping through raw, percental and effect size statistics
                for j = 1:3
                    % Exporting the statistic table of the iteration
                    try
                        name = data{1, n}.Properties.RowNames{j};
                        name2 = data{1, n}.Properties.VariableNames{k};
                        idx = find(strcmp(data{1, n}.(k){j}.Properties.VariableNames, 'r') == 1);
                        varNames = data{1, n}.(k){j}.Properties.VariableNames([1, idx:end]);
                        outputTable = table2cell(data{1, n}.(k){j}(:, [1, idx:end]));
                        idx = ismissing(data{1, n}.(k){j}(:, [1, idx:end]));
                        outputTable(idx) = {' '};
                        writecell(vertcat(varNames, outputTable), ...
                                  fullfile(dirs, [num2str(count), '_', name2, '_', name, '.xlsx'])); 
                        count = count+1;
                    catch
                        count = count +1;
                    end
                end
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

        % Exporting Correlation Matrices
        dirs = selpath + "\Tables\3_Statistics\3_Association Matrices";
        mkdir(dirs);
        
        for k = 1:size(corrMatrices, 1)
            name = corrMatrices.Properties.RowNames(k);
            varNames = corrMatrices.CorrelationMatrix{k, 1}.Properties.VariableNames;
            rowNames = corrMatrices.CorrelationMatrix{k, 1}.Properties.RowNames;
            outputTable = table2cell(corrMatrices.CorrelationMatrix{k, 1});
            idx = ismissing(corrMatrices.CorrelationMatrix{k, 1});
            outputTable(idx) = {' '};
            writecell(horzcat([{'Variable'}; rowNames], vertcat(varNames, outputTable)),...
                                fullfile(dirs, [num2str(k), '_', name{:}, '.xlsx'])); 
        end
    catch 
    end
end
if(sum(export_options) > 0)
    t = toc;
    fprintf('%-2.2f seconds\n', t);
end
end


