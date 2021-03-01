function [mu, ci_left, ci_right, pdf] = fit_dist(x, y, mode, alpha)
% This function calculates a probability distribution function for a given
% array of values.
%
% USAGE
% E.g.:
% [mu, ci_left, ci_right, pdf] = fit_dist(x, y, 'Gaussian', 0.05)
%
% INPUT
% x:            Array (n x 1) over which the probability distribution 
%               function will be evaluated.
% y:            Array (n x 1) which contains the values that the
%               distribution function will be fitted to.
% mode:         String or character vector that defines the type of 
%               distribution that will be used. Options: "Gaussian", 
%               "Gamma".
% alpha:        Scalar that sets the alpha error rate for calculating
%               the limits of the confidence intervals.
%
% OUTPUT
% mu:           Mean (Integral(p) = 0.5) of the distribution.
% ci_left:      Lower limit (Integral(p) = "alpha") of the confidence 
%               interval.
% ci_right:     Upper limit (Integral(p) = 1-"alpha") of the confidence 
%               interval.
% pdf:          Array (n x 1) that contains the normalized probability 
%               distribution function.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 20.January.2021
% lepremiere
%---------------------------------------------------------------------------------------------------

switch mode
    case 'Gaussian'
        [mu, sig]   = normfit(y);                           % Fitting normal distribution
        pdf         = normpdf(x,mu,sig);                    % Numerical distribution function
        cpdf        = cumsum(pdf);                          % Cummulative distribution function

        ci_left     = x(find(cpdf >= alpha, 1, 'first'));   % Lower limit of the confidence interval
        ci_right    = x(find(cpdf >= 1-alpha, 1, 'first')); % Upper limit of the confidence interval
        pdf         = pdf./max(pdf);                        % Normalizing distribution function

    case 'Gamma'
        x1 = [0:max(y)*5];                                  % Changing x to match special requirements
        params = gamfit(y);                                 % Estimating parameters 
        a = params(1);
        b = params(2);
        pdf = gampdf(x1,a,b);                               % Numerical distribution function
        
        % Special case if the mean of the function is
        % to close to zero, and therefore the distribuition
        % changes from skewed bell-shaped curve to exponential
        % function
        if(pdf(1) == max(pdf))                              
            cpdf = cumsum(1 - (1 - pdf));                    
        else
            cpdf = cumsum(pdf);
        end

        mu = x1(find(cpdf >= 0.5, 1, 'first'));             % Mean of distribution function
        ci_left = x1(find(cpdf >= alpha, 1, 'first'));      % Lower limit of the confidence interval
        ci_right = x1(find(cpdf >= 1-alpha, 1, 'first'));   % Upper limit of the confidence interval 
        pdf = gampdf(x,a,b);                                % Numerical distribution function
        pdf = pdf./max(pdf);                                % Normalizing distribution function
end
end