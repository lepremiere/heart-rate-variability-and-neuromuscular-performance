function plot_individuals(fig, analysisData, participants, output_variables,...
                          disp_vars, xlabels, includedPhases, show_baseline, baseline, alpha)
% This function calculates qualities of maximum voluntary contractions
% paired with superimposed and resting electrical nerve stimulations.
%
% USAGE
% E.g.:
% Results = Dynamomety(dyno_data, {'XY_INT_1', 'XY_INT_2'}, 1)
%
% INPUT
% data:             
%
% test_occasions:   Must be a string array (n x 1) that contains names to
%                   map the trials of the data to. E.g. {'XY_INT_1', 'XY_INT_2'}. Trials 
%                   will be named accordingly 'XX_INT_1_1', 'XY_INT_1_2'.
%
% n:                Points at the position in test_occasions, which should
%                   be used for name mapping.
%
% OUTPUT
% Results:          Table (m + 3 x 4) containing every trial plus mean,
%                   coefficient of variation, and the best trial for 4 metrics: Maximum
%                   voluntary contraction torque, superimposed twitch torque, resting twitch
%                   torque and voluntary activation as percentage.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 09.July.2020
% lepremiere
%---------------------------------------------------------------------------------------------------

type = participants;
err_bar = false;

% Assigning data
if(type == "mean")
    data = analysisData.Results{11, 1};
    data2 = analysisData.Results{12, 1};
else
    data = analysisData.Results{3, 1};
    data2 = analysisData.Results{3, 1};
    data = data(data.Participant == participants,:);
end

% Checking for available tests in data. HRV, CMJ, SQV, ITT
for i = 1:length(output_variables)
    av_tests(i) = sum(contains(data.Variable,[output_variables{i}])) > 0;
end

if(show_baseline == 1)
    includedPhases = [baseline, includedPhases];
    xlabels = ['BASELINE'; xlabels];
end

% Filtering variables that should be displayed for availability
disp_vars = disp_vars(av_tests);
ind = [output_variables{:}];
[~,inds] = max(data.Variable == ind([disp_vars{:}]));
var_names = data.Variable(inds);

mu   = data(inds,:);                                                                 % Final data
sig = data2(inds,:);
swc = 0.5;                                                                             % Smallest worthwhile change. Width of the gray shaded area in the plot                                                  
study_design   = string(mu.Properties.VariableNames);                                % get study phases of table from column names
targetPhases   = study_design(find(contains(study_design, includedPhases) == 1));    % only keep phases equal to "includedPhases"
underscores    = strfind(targetPhases,'_');                                          % find delimeter of phase names
analysisData_mu   = {mu(:,targetPhases)};                                            % cut data with respect to containing "includedPhases"
analysisData_sig   = {sig(:,targetPhases)};

% Appereance behavior
titles =  {'Counter-Movement Jump', 'Bar Velocitiy', 'Dynamometry'};
titles = titles(av_tests(2:end));
set(0,'CurrentFigure',fig); clf;
line_types = {'-d', '--d', '-.d', ':d'};
colors = [{[0, 0.4470, 0.7410]}, {[1, 0.5940, 0.1250]}, {[0.9500, 0.3250, 0.0980]}];
colors = colors(av_tests(2:end));

spaces = cell2mat(strfind(analysisData_mu{1, 1}.Properties.VariableNames, ' '));     % find all spaces in phase names, which indicates delta phases

% Getting integer from phase name to identify repeated measures design
for i = 1:numel(targetPhases)
    numRepetition(i,1) = str2num(targetPhases{i}(1:underscores{i}(1)-1));
end

% Getting number of distinct study repetitions
num   = unique(numRepetition);

% Getting data of interest for every study repetition
for i = 1:numel(num)
    for j = 1:size(analysisData_mu, 1)
        phaseIndeces{i} = find(numRepetition == num(i));                             % Identifies the column/phase in data that matches the repeated measures
        dataOI_mu{j, i}    = analysisData_mu{j,:}(:,[targetPhases(phaseIndeces{i})]);% Splits analysisData into subsets of measurement repetitions
        dataOI_sig{j, i} = analysisData_sig{j,:}(:,[targetPhases(phaseIndeces{i})]);
    end
end
%%
m = size(dataOI_mu{1, 1}, 2);
n = size(dataOI_mu, 2);
l = size(dataOI_mu{1, 1}, 1);
k = size(disp_vars, 2);

for i = 1:size(var_names, 1)
    limit = find(var_names{i} == '[');
    idx = isstrprop(var_names{i}(1:limit(1) - 1),'upper');
    legend_names{i,:} = var_names{i}(idx); 
end

mu = [];
sd = [];
for i = 1:n
    mu = [mu, dataOI_mu{1, i}{:,:}, nan(l, 1)];
     sd = [sd, dataOI_sig{1, i}{:,:}, nan(l, 1)];
end
disp_idx = {};
for i = 1:k
    if(sum(contains(data.Variable, ind(disp_vars{i})),1) > 0)
        disp_idx = [disp_idx; [1:size(disp_vars{i}, 2)]];
        if(i > 1)
            disp_idx{end,:} =  disp_idx{end,:} + max([disp_idx{1:i-1, :}]);
        end
    end
end

for i = 1:n
   if(show_baseline == 1 && i < 2 )
       temp_ticks = sum(contains(dataOI_mu{i}.Properties.VariableNames, 'BASE'))/2+1;
   elseif(i == 1)
      temp_ticks = []; 
   end
    swc_idx(i, :) = phaseIndeces{i}(~contains(dataOI_mu{i}.Properties.VariableNames, 'BASE'));
   
   if(i > 1)
      swc_idx(i, :) = swc_idx(i, :)+1;  
   end
end

if(length(includedPhases) > 2)
    num_subplots = 3;
else
    num_subplots = 1;
end

for ii = 1:length(disp_idx)-1        
    plt = subplot(3, 1, ii);
    hold on

    for j = 1:n
        fill([swc_idx(j,1),  swc_idx(j,1), swc_idx(j,end), swc_idx(j,end)],...
             [-swc swc swc -swc],[0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeColor','none');
    end

    for i = [1.6449, 1.96, 2.5758]
        for j = 1:n
           xline(swc_idx(j, 1), '--k', 'LineWidth', 2);
           plot(swc_idx(j,:), ones(length(swc_idx(j,:)),1)*i,  'LineStyle', '--', 'Color', [0.8 0.8 0.8], 'LineWidth', 1.15);
           plot(swc_idx(j,:), ones(length(swc_idx(j,:)),1)*-i, 'LineStyle', '--', 'Color', [0.8 0.8 0.8], 'LineWidth', 1.15);
        end
    end

    if(type == "mean" && err_bar == true)
        for i = 1:size(disp_idx{1}, 2)
            h(i) = errorbar(mu(disp_idx{1}(i), :), sd(disp_idx{1}(i), :), line_types{i},...
                             'Color', [0 0 0],'LineWidth', 1.5,'MarkerSize', 4);
        end

        for i = 1:size(disp_idx{ii+1}, 2)
            g(i) = errorbar(mu(disp_idx{ii+1}(i), :), sd(disp_idx{ii+1}(i), :), line_types{i},...
                            'Color', colors{ii},'LineWidth', 1.5,'MarkerSize', 4);                   
        end
    else
        for i = 1:size(disp_idx{1}, 2)
            h(i) = plot(mu(disp_idx{1}(i), :), line_types{i},...
                            'Color', [0 0 0],'LineWidth', 1.5,'MarkerSize', 4);
        end

        for i = 1:size(disp_idx{ii+1}, 2)
            g(i) = plot(mu(disp_idx{ii+1}(i), :), line_types{i},...
                            'Color', colors{ii},'LineWidth', 1.5,'MarkerSize', 4);                   
        end
    end
    
    title(titles(ii), 'FontSize', 12);
    xticks(sort([temp_ticks; reshape(swc_idx, [], 1)]));
    xticklabels(xlabels);
    if(numel(xlabels) > 8)
        ax = gca;
        ax.XTickLabelRotation = -45;
    end
    xlim([0.5, m*n + 1 + 0.5]);
    ylim([min([plt.YLim(1) * 1.1, -6]), max([plt.YLim(2) * 1.1, 6])]);
    yticks(round(yticks));
    set(gca,'FontWeight','bold');

    if(n <= 1 || show_baseline == 1)
        loc = 'northeastoutside';
        if(range(xlim) < 6)
            pbaspect([1.25 1 1]);
        elseif(range(xlim) < 10)
            pbaspect([2 1 1]);
        else
            pbaspect([3 1 1]);
        end
        xlim([0.5, max(max(swc_idx))+0.5]);
    else
        loc = 'north';
        xlim([0.5, max(max(swc_idx))+0.5]);
    end
    legend([h, g], legend_names{[disp_idx{1}, disp_idx{ii+1}]}, 'Location', loc);
    set(gca,'LineWidth',2','TickDir','out','TickLength',[.005 .005]);
    ylabel('n.u.','FontSize',8,'FontWeight','bold','FontSize',10);
    hold off

end
end
