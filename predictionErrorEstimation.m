function [C,RMSE] = predictionErrorEstimation(varargin)
% This function uses Cholasky factorization to simulate the root mean squared 
% error of two normally distributed variables with a correlation from 0 to 1.
%
% USAGE
% E.g.:
% [C, RMSE] = predictionErrorEstimation('mu', 0, 'sig', 1,...
%                                       'nSims', 1e3, 'n', 1e4, 'precision', 0.01,...
%                                       'error_request', 0.5, 'verbose', 1);
%
% INPUT
% mu:               Determines the mean that is used to generate random datasets.  
% sig:              Determines the standard deviation for the datasets.
% nSims:            Sets the number of simulations per r.
% n:                Sets the length of the random data sets.
% precision:        This number represents the step size for r.
% error_request:    Determines the error rate that should be found along
%                   the correlation coefficients.
% verbose:          Shows plot if 1. 
%
% OUTPUT
%
% C:                Vector containing the mean correlation that were
%                   calculated for each step in r.
% RMSE:             Vector conatining the mean root mean square error for
%                   each step in r.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 15.October.2020
% lepremiere
%---------------------------------------------------------------------------------------------------

% Options
default_mu = 0;                 % Mean of artificial data samples
default_sig = 1;                % Standard deviation of artificial data samples
default_nSims = 1e4;            % Number of simulations per r
default_n = 1e4;                % Number data points in data samples
default_precision = 0.001;      % Step size for r
default_error_request = 0.5;    % Error rate that should looked for
default_verbose = 1;            % Shows plots if 1

% Creating input parser and setting default values
p = inputParser;
addParameter(p, 'mu',               default_mu);
addParameter(p, 'sig',              default_sig);
addParameter(p, 'nSims',            default_nSims,          @(x) x > 0);
addParameter(p, 'n',                default_n,              @(x) x > 3);
addParameter(p, 'precision',        default_precision,      @(x) x <= 0.1);
addParameter(p, 'error_request',    default_error_request);
addParameter(p, 'verbose',          default_verbose,        @(x) x == 1 || x == 0);
parse(p, varargin{:});
ip = p.Results;

%% Assigning the varargin to variables
mu = ip.mu;
sig = ip.sig;
nSims = ip.nSims;
n = ip.n;
precision = ip.precision;
error_request = ip.error_request;
verbose = ip.verbose;

r = linspace(0, 1-precision, 1/precision);

% Looping through correlation coefficients r
parfor m = 1:length(r) 
    disp([num2str(m),'/',num2str(length(r))])
    rMSE = [];
    Cor = [];
    
    % Looping through number of simulations
    for i = 1:nSims
        Dis = sig*randn(n,2);               % Target distribution
        Rwish = [1, r(m);r(m),1];           % Desired correlation matrix
        L = chol(Rwish);                    % Cholasky factorization
        M = Dis*L+mu;                       % Output distribution
        x = M(:,1); y = M(:,2);
        c = corrcoef(x,y);                  % Correlation
        Cor(i) = c(2,1);                    % Getting corr coef from matrix
        rMSE(i) = sqrt(mean((x - y).^2));   % Root mean squared error
    end
    
    RMSE(m) = mean(rMSE);
    C(m) = mean(Cor);  
end

%% Plot
if(verbose == 1)
    figure(1);
    clf;
    plot(C, RMSE, '-', 'LineWidth', 4);
    set(gca, 'FontWeight', 'bold', 'FontSize', 12, 'LineWidth', 2);
    xlabel('Correlation [r]', 'FontWeight', 'bold', 'FontSize', 14);
    ylabel('Root Mean Squared Error [SD]', 'FontWeight', 'bold', 'FontSize', 14);

    idx = find(RMSE <= error_request, 1, 'first'); % Find error rquest
    if(idx > 0)
        yline(RMSE(idx), '--k', 'Label',['r', char(8776), ' ', num2str(r(idx)), ' @ 0.5*SD'], 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 2, 'LabelHorizontalAlignment', 'center');
        xline(r(idx), '--k', 'LineWidth', 2);
    else
        disp('Requested error rate not found!');
    end
end
end
