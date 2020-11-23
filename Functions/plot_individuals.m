function plot_individuals(fig, transformedData, participants, disp_vars, xlabels)

    
    var_names = transformedData.Results{3, 1}.Variable([disp_vars{:}]);

    data = transformedData.Results{3, 1}(transformedData.Results{3, 1}.Participant == participants,:);
    mu   = data([disp_vars{:}],:);

    swc = 0.5;

    includedPhases = {'INT', 'REC'};
    study_design   = string(mu.Properties.VariableNames);                                % get study phases of table from column names
    targetPhases   = study_design(find(contains(study_design, includedPhases) == 1));    % only keep phases equal to "includedPhases"
    underscores    = strfind(targetPhases,'_');                                          % find delimeter of phase names
    analysisData   = {mu(:,targetPhases)};                        % cut data with respect to containing "includedPhases"

    titles =  {'Counter-Movement Jump', 'Bar Velocitiy', 'Dynamometry'};

    set(0,'CurrentFigure',fig);
    clf;
    line_types = {'-d', '--d', '-.d', ':d'};
    colors = [{[0, 0.4470, 0.7410]}, {[1, 0.5940, 0.1250]}, {[0.9500, 0.3250, 0.0980]}];
    %%
    spaces = cell2mat(strfind(analysisData{1, 1}.Properties.VariableNames, ' '));        % find all spaces in phase names, which indicates delta phases

    % Getting integer from phase name to identify repeated measures design
    for i = 1:numel(targetPhases)
        numRepetition(i,1) = str2num(targetPhases{i}(1:underscores{i}(1)-1));
    end

    % Getting number of distinct study repetitions
    num   = unique(numRepetition);

    % Getting data of interest for every study repetition
    for i = 1:numel(num)
        for j = 1:size(analysisData, 1)
            phaseIndeces{i} = find(numRepetition == num(i));                             % Identifies the column/phase in data that matches the repeated measures
            dataOI{j, i}    = analysisData{j,:}(:,[targetPhases(phaseIndeces{i})]);        % Splits analysisData into subsets of measurement repetitions
        end
    end
    %%
    m = size(dataOI{1, 1}, 2);
    n = size(dataOI, 2);
    l = size(dataOI{1, 1}, 1);
    k = size(disp_vars, 2);

    for i = 1:size(var_names, 1)
        limit = find(var_names{i} == '(');
        idx = isstrprop(var_names{i}(1:limit(1) - 1),'upper');
        legend_names{i,:} = var_names{i}(idx); 
    end

    mu = [];
    sd = [];
    for i = 1:n
        mu = [mu, dataOI{1, i}{:,:}, nan(l, 1)];
%         sd = [sd, dataOI{2, i}{:,:}, nan(l, 1)];
    end

    for i = 1:k
        disp_idx{i,:} = [1:size(disp_vars{i}, 2)];
        if(i > 1)
            disp_idx{i,:} = disp_idx{i - 1, :}(end) + [1:size(disp_vars{i}, 2)];
        end
    end

    for i = 1:n
       swc_idx(i, :) = [1:m];
       if(i > 1)
          swc_idx(i, :) = swc_idx(i-1, end) + [1:m] + 1;  
       end
    end

    for ii = 1:k-1
        plt = subplot(k-1, 1, ii);
        hold on

        for j = 1:n
            fill([swc_idx(j, 1), swc_idx(j, 1), swc_idx(j,end), swc_idx(j,end)],...
                 [-swc swc swc -swc],[0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeColor','none');
        end

        for i = [0.7, 1.1, 1.7, 2.5, 4.5]
            for j = 1:n
               xline(swc_idx(j, 1), '--k', 'LineWidth', 2);
               plot(swc_idx(j,:), ones(m,1)*i,  'LineStyle', '--', 'Color', [0.8 0.8 0.8]);
               plot(swc_idx(j,:), ones(m,1)*-i, 'LineStyle', '--', 'Color', [0.8 0.8 0.8]);
            end
        end

        for i = 1:size(disp_idx{1}, 2)
            h(i) = plot(mu(disp_idx{1}(i), :), line_types{i},...
                            'Color', [0 0 0],'LineWidth', 1.5,'MarkerSize', 4);
        end

        for i = 1:size(disp_idx{ii+1}, 2)
            g(i) = plot(mu(disp_idx{ii+1}(i), :), line_types{i},...
                            'Color', colors{ii},'LineWidth', 1.5,'MarkerSize', 4);
        end

    %     max_y = 
        title(titles(ii), 'FontSize', 12);
        xticks(sort(reshape(swc_idx, [], 1)));
        xticklabels(xlabels);
        xlim([0.5, m*n + 1 + 0.5 ]);
        yticks(round(yticks));
        ylim(plt.YLim * 1.0);
        set(gca,'FontWeight','bold');
        legend([h, g], legend_names{[disp_idx{1}, disp_idx{ii + 1}]}, 'Location', 'north');
        set(gca,'LineWidth',2','TickDir','out','TickLength',[.005 .005]);
        ylabel('n.u.','FontSize',8,'FontWeight','bold','FontSize',10);
        hold off
    end


end
