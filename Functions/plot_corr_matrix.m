function plot_corr_matrix(fig, input)  
% This function plots a heatmap of a correlation matrix on a specified figure.
%
% USAGE
% E.g.:
% plot_corr_matrix(fig_handle, CorrMatrices)  
%
% INPUT
% fig:  Figure handle that determines the target to be plotted on. 
% k:    Matrix (n x m) that holds values for the heatmap.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 20.January.2021
% lepremiere
%---------------------------------------------------------------------------------------------------

fig;
clf
% Plotting the heatmap
h = heatmap(input.Properties.VariableNames,...
    input.Properties.RowNames,...
    input{:,:},...
    'Colormap', parula,...
    'FontSize', 10);

% Appereance behavior
caxis(h, [-1,1]);
ylabel('Heart Rate Variability');
xlabel('Neuromuscular Performance');
end