function Export(export_options)

    if(sum(export_options) > 0)
        tic
        fprintf('Export: \t\t\t')
        
        selpath                 = uigetdir;
        groupData               = evalin('base', 'groupData');
        reproducibilityResults  = evalin('base', 'reproducibilityResults');
        associationResults      = evalin('base', 'associationResults');
        transformedData         = evalin('base', 'transformedData');
        correlation_Type        = evalin('base', 'correlation_Type');
        ANCOVA_type             = evalin('base', 'ANCOVA_type');
        xlabels                 = evalin('base', 'xlabels');
        disp_vars               = evalin('base', 'disp_vars');
        study_design            = evalin('base', ' study_design');
        
        if(selpath == 0)
            fprintf('\nNo folder selected\n');
            return
        end

    end
        % Params
        namePrefix = {'Reproducibility', 'Association'};
        fig = figure('visible', 'off', 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1], 'color', 'w'); 
        
        if(export_options(2) == 1)
            % Reproducibility data
            if(size(reproducibilityResults, 1) > 0)
                % Folders for reproducibility data
                dirs  = selpath + "\Figures\Statistics\1_Reproducibility\";
                ind = 1;
                % Recursive plotting
                mkdir(dirs);
                for i = 1:height(reproducibilityResults.Absolute{1,1})
                    try
                        M = {reproducibilityResults.Absolute{1,1}{i,2}{:,:},...
                             reproducibilityResults.Delta{1,1}{i,2}{:,:},...
                             reproducibilityResults.Absolute{2,1}{i,2}{:,:},...
                             reproducibilityResults.Delta{2,1}{i,2}{:,:},...
                             reproducibilityResults.Absolute{3,1}{i,2}{:,:},...
                             reproducibilityResults.Delta{3,1}{i,2}{:,:}};

                        if(correlation_Type(1) == 3) 
                            idx = [5, 10, 11];
                            T = {reproducibilityResults.Absolute{1,1}{i,4}{:,:},...
                                 reproducibilityResults.Delta{1,1}{i,4}{:,:},...
                                 reproducibilityResults.Absolute{2,1}{i,4}{:,:},...
                                 reproducibilityResults.Delta{2,1}{i,4}{:,:},...
                                 reproducibilityResults.Absolute{3,1}{i,4}{:,:},...
                                 reproducibilityResults.Delta{3,1}{i,4}{:,:}};
                        else
                            idx = [3, 6, 7];
                            T = [];
                        end
                            stats   = {reproducibilityResults.Absolute{1,1}{i, idx};...
                                       reproducibilityResults.Delta{1,1}{i, idx};...
                                       reproducibilityResults.Absolute{2,1}{i, idx};...
                                       reproducibilityResults.Delta{2,1}{i, idx};...
                                       reproducibilityResults.Absolute{3,1}{i, idx};...
                                       reproducibilityResults.Delta{3,1}{i, idx}};

                            h0 = reproducibilityResults.Absolute{1,1}.Properties.VariableNames{end};

                        if(ind == 1)
                           plot_stats_LM(fig, M, T, stats, h0, namePrefix{1}, correlation_Type(1), ANCOVA_type{1}, 1);
                           im = getframe(gcf);
                           ind = ind + 1;
                        end

                        plot_stats_LM(fig, M, T, stats, h0, namePrefix{1}, correlation_Type(1), ANCOVA_type{1}, 1);
                        im = getframe(gcf);
                        imName = [num2str(i), '_', namePrefix{1}, '_',...
                                  strrep(M{1}.Properties.VariableNames{2}(1:find(M{1}.Properties.VariableNames{2} == '(') - 2), '/', '_'),...
                                  '.png'];
                        imwrite(im.cdata, fullfile(dirs, imName));
                        ind = ind + 1;
                    catch
                    end
                end


            else
                fprintf('\nERROR: No reproducibility data found!\n');
            end

    %         % Association data
            if(size(associationResults, 1) > 0)
                % Folders for reproducibility data
                dirs  = selpath + "\Figures\Statistics\2_Association\";
                ind = 1;
                % Recursive plotting
                mkdir(dirs);
                for i = 1:height(associationResults.Absolute{1,1})
                    try
                        M = {associationResults.Absolute{1,1}{i,2}{:,:},...
                             associationResults.Delta{1,1}{i,2}{:,:},...
                             associationResults.Absolute{2,1}{i,2}{:,:},...
                             associationResults.Delta{2,1}{i,2}{:,:},...
                             associationResults.Absolute{3,1}{i,2}{:,:},...
                             associationResults.Delta{3,1}{i,2}{:,:}};

                        if(correlation_Type(2) == 3) 
                            idx = [5, 10, 11];
                            T = {associationResults.Absolute{1,1}{i,4}{:,:},...
                                 associationResults.Delta{1,1}{i,4}{:,:},...
                                 associationResults.Absolute{2,1}{i,4}{:,:},...
                                 associationResults.Delta{2,1}{i,4}{:,:},...
                                 associationResults.Absolute{3,1}{i,4}{:,:},...
                                 associationResults.Delta{3,1}{i,4}{:,:}};
                        else
                            idx = [3, 6, 7];
                            T = [];
                        end
                        stats   = {associationResults.Absolute{1,1}{i, idx};...
                                   associationResults.Delta{1,1}{i, idx};...
                                   associationResults.Absolute{2,1}{i, idx};...
                                   associationResults.Delta{2,1}{i, idx};...
                                   associationResults.Absolute{3,1}{i, idx};...
                                   associationResults.Delta{3,1}{i, idx}};

                        h0 = associationResults.Absolute{1,1}.Properties.VariableNames{end};

                        if(ind == 1)
                           plot_stats_LM(fig, M, T, stats, h0, namePrefix{2}, correlation_Type(2), ANCOVA_type{2}, 2);
                           im = getframe(gcf);
                           ind = ind + 1;
                        end

                        plot_stats_LM(fig, M, T, stats, h0, namePrefix{2}, correlation_Type(2), ANCOVA_type{2}, 2);
                        im = getframe(gcf);
                        imName = [num2str(i), '_', namePrefix{2},'_',...
                                     [strrep(M{1}.Properties.VariableNames{2}(1:find(M{1}.Properties.VariableNames{2} == '(') - 2), '/', '_'),'_'...
                                      strrep(M{1}.Properties.VariableNames{3}(1:find(M{1}.Properties.VariableNames{3} == '(') - 2), '/', '_')],...  
                                     '.png'];
                        imwrite(im.cdata, fullfile(dirs, imName));
                    catch
                    end
                end


            else
                fprintf('\nERROR: No association data found!\n');
            end
        end
    if(export_options(3) == 1)

        fig = figure('visible', 'off', 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.5, 1], 'color', 'w');
        participants = unique(transformedData.Results{1, 1}.Participant);
        
        dirs  = selpath + "\Figures\Individuals\";
        mkdir(dirs);

        for i = 1:numel(participants)  

            plot_individuals(fig, transformedData, participants(i), disp_vars, xlabels);
            im = getframe(gcf);
            imName = [participants{i}, '.png'];
            imwrite(im.cdata, fullfile(dirs, imName));
            
        end
    end
    
    if(export_options(1) == 1)
        
        dirs  = selpath + "\Tables\";
        mkdir(dirs);
        
        varNames = groupData.Properties.VariableNames;
        outputTable = table2cell(groupData);
        idx = ismissing(groupData);
        outputTable(idx) = {' '};
        
        writecell(vertcat(varNames, outputTable), fullfile(dirs, 'groupData.xlsx'));  
        
        dirs  = selpath + "\Tables\Data\";
        mkdir(dirs);
        
        for k = 1:size(transformedData, 1)
            
            name = transformedData.Properties.RowNames(k);
            varNames = transformedData.Results{k, 1}.Properties.VariableNames;
            rowNames = [{' '}; transformedData.Results{k, 1}.Properties.RowNames];
            outputTable = table2cell(transformedData.Results{k, 1});
            idx = ismissing(transformedData.Results{k, 1});
            outputTable(idx) = {' '};
            
            if(k > 6)
                writecell(horzcat(rowNames, vertcat(varNames, outputTable)), fullfile(dirs, [num2str(k), '_', name{:}, '.xlsx'])); 
            else
                writecell(vertcat(varNames, outputTable), fullfile(dirs, [num2str(k), '_', name{:}, '.xlsx'])); 
            end
        end
        
        dr = {"1_Reproducibility\", "2_Association\"};
        data = {reproducibilityResults, associationResults};
                
        for n = 1:2
            count = 1;
            dirs  = selpath + "\Tables\Statistics\" + dr{n};
            mkdir(dirs);
            
            for k = 1:2
                for j = 1:3
                    try
                        name = data{1, n}.Properties.RowNames{j};
                        name2 = data{1, n}.Properties.VariableNames{k};
                        idx = find(strcmp(data{1, n}.(k){j}.Properties.VariableNames, 'r') == 1);
                        varNames = data{1, n}.(k){j}.Properties.VariableNames([1, idx:end]);
                        outputTable = table2cell(data{1, n}.(k){j}(:, [1, idx:end]));
                        idx = ismissing(data{1, n}.(k){j});
                        outputTable(idx) = {' '};

                        writecell(vertcat(varNames, outputTable), fullfile(dirs, [num2str(count), '_', name2, '_', name, '.xlsx'])); 
                        count = count+1;
                    catch
                        count = count +1;
                    end
                end
            end
        end
        
        
        
    end
    if(sum(export_options) > 0)
        t = toc;
        fprintf('%-2.2f seconds\n', t);
    end
end


