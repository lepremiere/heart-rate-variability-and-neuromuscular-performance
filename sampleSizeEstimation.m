 function sampleSizeEstimation(varargin)
% This function calculates the confidence intervals for repeated measures
% correlations with varying degrees of freedom. Additionally, a target
% sample size can be included to see which correlation coefficient will
% have which confidence interval at a given degree of freedom. Confidence
% intervals are calculated with Fishers's Z-transformation. The output is a
% figure containing those information for analysis of reproducibility and
% validity.
%
% USAGE
% E.g.:
% sampleSizeEstimation('n', 4, 'rep_measures', 2, 'num_cycles', 2, 'alpha', 0.1)
%
% INPUT
% n:            Determines the sample size, or how many participants you are
%               going to include in your study. Must be greater than 1.
% n_sims:       Determines the maximal sample size for which this function
%               will calculate the metrics. Must be greater than 1.
% rep_measures: Determines the measurement repetitions you conduct. Must be
%               greater or equal to 3.
% num_cycles:   Indicates how many cycles of the repeated measurement you
%               are conducting. Must be greater or equal to 1.
%               Reproducibility will only be obtained if it is greater than 1.
% alpha:        Sets the alpha error for the confidence intervals (E.g.
%               alpha = 0.1 --> 90% confidence interval). Must be greater
%               than 0 and smaller than 1.
%
% OUTPUT
% This function only outputs a figure with two plots. Left: reliability,
% right: validity. The plots have the correlation [r] on the x-axis, the
% corresponding confidence intervals [r] on the y-axis and the degrees of
% freedom [df] as color. Further, black solid lines are included on the
% common interpretation categories for the strength of correlation. Lastly,
% a black dashed line indicates how confidence intervals would fall with
% your study set-up at a given correlation.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 15.October.2020
% lepremiere
%---------------------------------------------------------------------------------------------------

%% Settings
addpath([pwd, '\Functions\']);
fig = figure('units','normalized','outerposition',[0 0 1 1]);
set(fig, 'color', 'w');
titles = {'Reliabity', 'Validity'};
clf
tic;
% Options
default_n = 4;             	% Number of participants
default_n_sims = 0;         % Upper limit of degrees of freedom for which should be simulated
default_rep_measures = 4;   % Number of repeated measures within one cycle
default_num_cycles = 1;     % Number of cycles that the repeated measures are conducted
default_alpha = 0.1;        % Alpha error level for confidence interval calculation
default_resolution = 1;     % Resolution of the increments to 'n_sims'

% Creating input parser and setting default values
p = inputParser;
addParameter(p, 'n',                default_n,              @(x) x > 0);
addParameter(p, 'n_sims',           default_n_sims,         @(x) x > 1);
addParameter(p, 'rep_measures',     default_rep_measures,   @(x) x >= 3);
addParameter(p, 'num_cycles',       default_num_cycles,     @(x) x >= 1);
addParameter(p, 'alpha',            default_alpha,          @(x) x > 0 && x < 1);
addParameter(p, 'resolution',       default_resolution,     @(x) x == 0.1 || x == 0.5 || floor(x) == x);
parse(p, varargin{:});
ip = p.Results;

%% Assigning the varargin to variables
n = ip.n;
rep_measures = ip.rep_measures;
num_cycles = ip.num_cycles;
alpha = ip.alpha;
n_sims = ip.n_sims;
resolution = ip.resolution;

% Either taking input as 'n_sims' or participants times measurements times
% cycles times 125% 
if n_sims == 0
    n_sims = n*rep_measures*num_cycles*1.25;
else
    n_sims = (n*(rep_measures*n_sims-1)-1)*1.25;
end

r_OI =      linspace(0, 1, 100);                % Correlation coefficients of interest. Determines the resolution of the x-axis
N =         [4:resolution:n_sims];              % Array of number of participants over which will be simulated
r_limits =  [0, 0.2, 0.5, 0.75,  0.9, 0.99;...  % Interpretation limits of correlation coefficients. 
             0, 0.1, 0.3,  0.5,  0.7, 0.9];     % Adapted from Hopkins (2015) (sportsci.org/2015/ValidRely.htm).

% Calculates the degrees of freedom for the analysis of reproducibility
% (top) and validity. Formula used: n(k-1)-1 from Bakdash (2017) (doi:10.3389/fpsyg.2017.00456)
% With n = number of participants, k = number of repeated measures.
% Additionally, k is multiplied with the number of cycles. 
df =        [(n * (rep_measures * (num_cycles - 1) - 1) - 1) ,    (n * (rep_measures * (num_cycles - 1) - 1) - 1);...
              n * (rep_measures *  num_cycles - 1) - 1,            n *(rep_measures  *  num_cycles - 1) - 1];
% Looping       
for m = [1,2]
    
    subplot(1, 2, m);
    hold on

    % Plotting the interpreation categories (r_limits) as grid
    % Plotted as 3D cube, so that these lines are later on top
    for k = [1, -1]
        for i = m
            for j = 1:+size(r_limits(1, 1:end), 2)
                plot3([0, 1, 1, 0],     [r_limits(i, j) * k, r_limits(i, j) * k, r_limits(i, j) * k, r_limits(i, j) * k], [min(N), min(N), max(N), max(N)], 'k');
                plot3([r_limits(i, j),   r_limits(i, j),     r_limits(i, j),     r_limits(i, j)],    [0, 1*k, 1*k, 0],    [min(N), min(N), max(N), max(N)], 'k');
            end
        end
    end

    % Appereance behavior
    xlabel('Correlation [r]', 'FontWeight', 'bold', 'FontSize', 14);
    ylabel([num2str((1-alpha)*100),'% Confidence Interval [r]'], 'FontWeight', 'bold', 'FontSize', 14);
    set(gca,'LineWidth', 2','TickDir','in','TickLength',[.005 .005]);   
    set(gca,'FontWeight','bold');  
    xlim([-0.003 1.0]);
    box on

    % Plotting the correlaiton coefficients of interest as orientation
    plot(r_OI, r_OI, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2);
    
    % Calculating the confidence interval for the corresponding correlation
    % with function 'analyticalCI.m' which utilizes Fisher's
    % Z-transformation for each number of participants in 'N'.
    for j = 1:numel(N)
        for i = 1:numel(r_OI)
            [LB, UB] = analyticalCI(r_OI(i), N(j), alpha);
            r_CI(i, :) = [LB, UB];
        end
        CI{j, :} = r_CI;
    end

    % Gathering the confidence interval bounds as grids for later plotting
    % Additionally, plotting the dashed line for the study design declared
    % in settings.
    LLB = [];
    UUB = [];
    for i = 1:numel(N)
        LLB = [LLB, CI{i}(:, 1)];
        UUB = [UUB, CI{i}(:, 2)];
        if(N(i) == df(m))
            h = plot3(r_OI', LLB(:,end), ones(size(LLB, 1), 1)*N(end), '--k', 'LineWidth', 2); 
                plot3(r_OI', UUB(:,end), ones(size(UUB, 1), 1)*N(end), '--k', 'LineWidth', 2); 
        end
    end
    
    % Transforming the vectors 'N' and 'r' to matrices for surface plots
    n = ones(size(r_OI, 2), size(N, 2)) .* N;
    r = ones(size(r_OI, 2), size(N, 2)) .* r_OI';
    
    % Plotting the confidence intervals as surface
    s1 = surf(r, LLB, n, 'EdgeColor', 'none', 'FaceAlpha', 1);
    s2 = surf(r, UUB, n, 'EdgeColor', 'none', 'FaceAlpha', 1);
    
    % Appereance behavior
    colormap('hsv');
    c = colorbar;
    c.Label.String = "df";
    c.Label.FontSize = 14;
    c.Label.Rotation = 90;
    view([0 90])
    title(titles{m}, 'FontWeight', 'bold', 'FontSize', 16);
    
    % Error catch function, if sample configuration is not sufficient for
    % the simulation
    try 
        legend(h, ['df = ', num2str(df(m))], 'Location', 'southeast')
    catch
        legend( 'df = 0, insufficient sample', 'Location', 'southeast')
    end
end
t = toc 
