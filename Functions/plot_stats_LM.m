function plot_stats_LM(fig, M, T, stats, h0, namePrefix, correlation_Type, ANCOVA_type, type)
% This function plots a graphical illustration of the statistical results for 
% reproducibility and association results from "DataAnalysis".
%
% USAGE
% E.g.:
% plot_stats_LM(fig, M, T, stats, 'p(0.875)', 'Association', 3, 'paralell lines', 2) 
%
% INPUT
% fig:               Figure handle that determines the target to be plotted on. 
% M:                 Cell array (1 x 6) that contains tables of the data to be
%                    correlated as well as a participant identifier for raw values,
%                    percentages, and effect sizes, as well as delta values for each one.
%                    Organization: raw, delta raw, perc., delta perc., es, delta es.
% T:                 Cell array (1 x 6) that contains t-tables for each table in M.
% stats:             Cell array (6 x 1) that contains the correlation coefficient, 
%                    the p value tested against 0, and the p value tested against 
%                    the value specified by "h0" for each table in M.
% namePrefix:        String or character vector that will be used to label the figure.
% correlation_Type:  Scalar that determines which correlation will be performed.
%                    Options: 1 = Pearson, 2 = ICC, 3 = ANCOVA (repeated measures correlation)
% ANCOVA_type:       String or character vector that specifies the type of ANCOVA to use.
%                    Options: 'parallel lines', separate lines'
% type:              Scalar that determines which of reproducibility (1) or association (2) 
%                    is analyzed. This changes the titles of the plot.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 20.January.2021
% lepremiere
%---------------------------------------------------------------------------------------------------

set(0,'CurrentFigure',fig)
clf;
m = length(M);

% Label specifications
mode        = {'', '\Delta', '', '\Delta', '', '\Delta'};
version     = {'Raw', 'Raw', 'Percental', 'Percental', 'Normalized', 'Normalized'};

idx         = isstrprop(M{1}.Properties.VariableNames{2}(1:find(M{1}.Properties.VariableNames{2} == '[') - 1),'upper');
variable    = M{1}.Properties.VariableNames{2}(idx);
idx         = [find(M{1}.Properties.VariableNames{2} == '['),...
               find(M{1}.Properties.VariableNames{2} == ']')];
unit        = M{1}.Properties.VariableNames{2}(idx(1):idx(2));
label       = [variable, ' ', (unit)];
labels      = {label, label,...
              [variable , ' [%]'],   [variable , ' [%]'],...
              [variable, ' [n.u.]'], [variable, ' [n.u.]']};

% Special label specifiation for association analysis
if(namePrefix == "Association")
    idx         = isstrprop(M{1}.Properties.VariableNames{3}(1:find(M{1}.Properties.VariableNames{3} == '[') - 1),'upper');
    variable2    = M{1}.Properties.VariableNames{3}(idx);
    idx         = [find(M{1}.Properties.VariableNames{3} == '['),...
                   find(M{1}.Properties.VariableNames{3} == ']')];
    unit        = M{1}.Properties.VariableNames{3}(idx(1):idx(2));
    label       = [variable2, ' ', (unit)];
    labels2      = {label, label,...
                  [variable2 , ' [%]'],   [variable2 , ' [%]'],...
                  [variable2, ' [n.u.]'], [variable2, ' [n.u.]']};
end

% Looping through M
for i = 1:m
    data            = M{i};                 
    participants    = unique(data(:,1));    % Getting individuals present
    
    % Splitting data into separate datasets for each individual
    for j = 1:numel(participants) 
        participant_data{j} = data(data.Subjects(:, 1) == participants.Subjects(j), 2:3);
    end
    
    % Choosing which of the total 6 subplots to plot onto
    subplot(m/2, m/3, i);
    hold on
    
    % Plotting participant data iteratively to get distinct colors
    for j = 1:numel(participants)
      h = plot(participant_data{:, j}{:, 1}, participant_data{:, j}{:, 2},...
          'o', 'MarkerSize', 3, 'LineWidth', 3); 
      colors(j,:) = get(h, 'Color');
    end
    
    % Plot title. Depending on the type of analysis (reproducibility vs. association)
    if(type == 2)
        tit = horzcat(M{i}.Properties.VariableNames{2}(1:find(M{1}.Properties.VariableNames{2} == '[') - 2),'/',... 
              M{i}.Properties.VariableNames{3}(1:find(M{1}.Properties.VariableNames{3} == '[') - 2));
    else
        tit = M{i}.Properties.VariableNames{2}(1:find(M{1}.Properties.VariableNames{2} == '[') - 2);
    end

    if(correlation_Type < 3)
        f = fitlm(data(:, 2:3));                                % Fitting a linear model
        h = plot(f);                                            % Plotting the linear model
        delete(h(1));                                           % Deleting the black data points      
        
        % Appereance behavior
        set(h(2),'Color','k','LineWidth',2);
        set(h(3),'Color','k','LineWidth',1,'LineStyle','--');
        set(h(4),'Color','k','LineWidth',1,'LineStyle','--'); 
        wid = range(xlim)*1.5;
        xlim([mean(xlim)-wid/2, mean(xlim)+wid/2]);
        wid = range(ylim)*1.5;
        ylim([mean(ylim)-wid/2, mean(ylim)+wid/2]);
        axis square;
        l = legend('Location','eastoutside');
        title({[mode{i}, version{i}, ' ', tit, ' ', namePrefix], ' '},'Interpreter','tex');
        % Text box for stats
        pos = l.Position;
        a = annotation('textbox',pos,'String',...
                       {string(['f(x)          = ', num2str(round(f.Coefficients{2,1},2)),'x + ', num2str(round(f.Coefficients{1,1},2))]),...    
                        string(['RMSE          = ', num2str(round(f.RMSE,2))]),...
                        string(['R²            = ', num2str(round(f.Rsquared.Ordinary,3))]),...
                        string(['r             = ', num2str(round(stats{i}(1), 3))]),...
                        string(['p(r=0.000)    = ', num2str(round(stats{i}(2), 4))]),...
                        string([h0,          ' = ' , num2str(round(stats{i}(3), 4))])},...
                        'FontSize', 10, 'EdgeColor', 'k', 'LineWidth', 2,...
                        'BackgroundColor', [1 1 1], 'FontName', 'FixedWidth',...
                        'FontWeight', 'bold', 'FitBoxToText', 'on');
        set(l, 'visible', 'off');   

    else
        
        % Getting the slopes and intercepts of the covariation tables
        coT             = T{i};
        intercepts      = [coT.Estimate{1:1 + numel(participants), :}]';
        slopes          = [coT.Estimate{1 + numel(participants)+1:end, :}]';
        
        % Plotting individual data with previously specified colors
        for j = 1:numel(participants)
            
            % Distinguishing between parallel or separate lines
            if(ANCOVA_type == "parallel lines")
                z = 0;
            else
                z = slopes(j+1);
            end
            
            % Specifying final intercept and slope to be plotted as well as standard error 
            % for the individual 
            y = (intercepts(1) + intercepts(j+1)) + (slopes(1) + z)*participant_data{:, j}{:, 1};
            h = plot(participant_data{:, j}{:, 1}, y,'-', 'LineWidth', 2, 'Color', colors(j,:));
            SE(j) = nanmean((y - participant_data{:, j}{:, 2}).^2);
        end

        RMSE = sqrt(mean(SE));                              % Calculating RMSE for the group
        % Appereance behavior
        wid = range(xlim)*1.5;
        xlim([mean(xlim)-wid/2, mean(xlim)+wid/2]);
        wid = range(ylim)*1.5;
        ylim([mean(ylim)-wid/2, mean(ylim)+wid/2]);
        axis square;
        l = legend('Location','eastoutside');
        title({[mode{i}, version{i}, ' ',  tit, ' ', namePrefix], ' '});
        % Text box for stats
        pos = l.Position;
        a = annotation('textbox',pos,'String',...
                       {string(['<intercept>   = ' , num2str(round(intercepts(1), 3))]),...
                        string(['<slope>       = ' , num2str(round(slopes(1), 3))]),...
                        string(['RMSE          = ' , num2str(round(RMSE, 3))]),...
                        string(['r             = ' , num2str(round(stats{i}(1), 3))]),...
                        string(['p(r=0.000)    = ' , num2str(round(stats{i}(2), 4))]),...
                        string([h0,          ' = ' , num2str(round(stats{i}(3), 4))])},...
                        'FontSize', 10, 'EdgeColor', 'k', 'LineWidth', 2,...
                        'BackgroundColor', [1 1 1], 'FontName', 'FixedWidth',...
                        'FontWeight', 'bold', 'FitBoxToText', 'on');
        set(l, 'visible', 'off');        
    end
    
    % Appereance behavior
    set(gca,'FontWeight','bold');
    set(gca,'LineWidth',2','TickDir','out','TickLength',[.005 .005]); 
    xlabel([mode{i},labels{i}],'Interpreter','tex');
    
    % Distinguishing labels
    if(namePrefix == "Association")
        ylabel([mode{i},labels2{i}],'Interpreter','tex');
    else
        ylabel([mode{i},labels{i}, '\_2'],'Interpreter','tex');
    end
    
    % Adding a line with slope 1 for reproducibility analysis
    switch namePrefix
        case 'Reproducibility'
            lower_lim = min([xlim, ylim]);
            upper_lim = max([xlim, ylim]);
            xlim([lower_lim, upper_lim]);
            ylim([lower_lim, upper_lim]);
            plot(xlim, ylim, 'Color', [0.01 0.01 0.01 0.25]);  
    end
    hold off
end
end