function [threshold, two_means] = get_threshold(X, verbose, output_plot)
% This function estimates a threshold that can be used alongside a moving mean 
% on RR-data to identify outlier. The estimate is based on identifying two groups
% inside the data by a k-means algorigthm, fitting gamma distributions to those 
% groups and determining if those two groups differ too much. The estimate
% is then choosen either to cut off the outlier group or including both
% groups.
%
% USAGE
% E.g.:
% [threshold, two_means] = get_threshold(X, 1, figure_handle)
%
% INPUT
% X:            RR-intervals stored in an array (n x 1).
% verbose:      Options: 1 = graphical illustration will be plotted.
%               0 = otherwise.
% output_plot:  If verbose equals 1, a figure handle must be passed to
%               to which the graphical illustration will be plotted. 
%
% OUTPUT
% threshold:    Threshold in ms.
% two_means:    Logical that indicates wether two groups have been
%               identified.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 19.January.2021
% lepremiere
%---------------------------------------------------------------------------------------------------

y = X;                                              % RR-intervals   
movm = movingmean(X, 30);                           % Moving mean
y_diff = abs(y -  movm);                            % Difference vector 
[means, idx] = k_means(y_diff);                     % Returns the means of the identified groups and their indeces
x = [-max(abs(y_diff))*1.1:max(abs(y_diff))*1.1];   % Vector over which distribution function will be evaluated

[~,idx1] = min(abs(means));                         % Indeces of group 1
[~,idx2] = max(abs(means));                         % Indeces of group 2

% Fitting gamma distributions to the groups and the overall data
[mu1,  ci_left1, ci_right1, gauss1] = fit_dist(x, y_diff(idx==idx1), "Gamma",      0.05);
[mu2,  ci_left2, ci_right2, gauss2] = fit_dist(x, y_diff(idx==idx2), "Gamma",      0.05);
[mu3,  ci_left3, ci_right3, gauss3] = fit_dist(x, y_diff,            "Gamma",      0.05);

% Determining the threshold. If the difference of the upper limit of the
% confidence interval of the group with the lower mean and the lower limit
% of the group with the higher mean is smaller than a certain margin, the
% higher group will be excluded by the threshold.
if(ci_right1 - ci_left2 < -mu1/3)
    threshold =  ci_right1;
    two_means = 1;
else       
    threshold = ci_right3;
    two_means = 0;
end

% Plot
if(verbose == 1)
    output_plot;
    cla
    hold on
    
    % Plotting group data
    scatter(y_diff(idx==idx1), zeros(length(y_diff(idx==idx1)), 1), 'xg'); % Group 1
    scatter(y_diff(idx==idx2), zeros(length(y_diff(idx==idx2)), 1), 'xr'); % Group 2
    
    % Plotting means
    scatter(mu1, 0, 'og', 'filled');    
    scatter(mu2, 0, 'or', 'filled');    
    scatter(mu3, 0, 'om', 'filled');    
    
    % Plotting distribution functions
    plot(x, gauss1, 'g','LineWidth', 2);
    plot(x, gauss2, 'r','LineWidth', 2);
    plot(x, gauss3, 'm','LineWidth', 2);
    
    % Orientation lines
    xline(ci_left1, '--g');
    xline(ci_right1, '--g');
    xline(ci_left2, '--r');
    xline(ci_right2, '--r');
    xline(ci_left3, '--m');
    xline(ci_right3, '--m');
    h = xline(threshold, 'LineWidth', 2);
    
    % Apereance behavior
    legend(h, 'Threshold');
    xlim([0 max([ci_right2, ci_right3, threshold])*1.1]);
    ylim([0, 1.5]);
    title('Artefact Recognition');
    xlabel('Successive Differences [ms]');
    ylabel('Probability Density Function [p/max(p)]');
    ax = gca;
    ax.LineWidth = 2;
    ax.FontWeight = 'bold';
    hold off
end
end



