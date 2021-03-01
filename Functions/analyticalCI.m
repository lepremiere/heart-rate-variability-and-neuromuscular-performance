function [LB, UB] = analyticalCI(r, n, alpha)
% This function calculates the confidence interval for a given
% correlation coefficient, sample size and alpha level using
% Fisher's z-transformation.
%
% USAGE
% E.g.:
% [LB, UB] = analyticalCI(0.8, 43, 0.05)
%
% INPUT
% r:        Correlation coefficient for which the confidence interval should
%           be calculated, limits -1:1
% n:        Number of observations from which r is estimated
% alpha:    Maximum type 1 error rate allowed
%
% OUTPUT
% [LB, UB]: Lower and upper bound of the confidence interval
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 09.July.2020
% lepremiere
%---------------------------------------------------------------------------------------------------

% Fisher's z-transformation
z_critical  = norminv(abs(1 - alpha/2));               % Getting critical z values according to alpha level
z_prime     = 0.5 * log((1+r) / (1-r));                % r to z transformation
se          = 1 / sqrt(n-3);                           % Calculating standard error
LB          = tanh(z_prime - z_critical * se);         % Calculates z value for lower bound and transforms it via tanh
UB          = tanh(z_prime + z_critical * se);         % Calculates z value for upper bound and transforms it via tanh
    
end

    
    