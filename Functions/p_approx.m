function p = p_approx(r, r0, n)

% This function approximates the probability for a correlation coefficient
% beeing inferior to a certain true population effect.
%
% USAGE
% E.g.:
% p = p_approx(0.77, 0.8, 42)
%
% INPUT
% r:    Correlation coefficient that should be tested for 
%       inferiority compared to r0. 
%
% r0:   The Null-hypothesis to be tested for.
%
% n:    Number of observation in determining r.
%
% OUTPUT
% p:    Probability that r comes from a population with 
%       true effect r0.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 09.July.2020
% lepremiere
%---------------------------------------------------------------------------------------------------

precision   = 1e4;                                          % Determines the precision for the approximation
r_interval  = -1-1/precision:1/precision:1+(1/precision);   % Array for which the probability density function should be calculated
r0_pdf      = corrdist(r_interval, r0, n);                  % Probability density function for r0
r0_cum      = cumtrapz(r0_pdf)/trapz(r0_pdf);               % Cummulative density function for r0
p           =  r0_cum(find(r_interval >= r, 1));            % Cummulative probability at r
if(r > r0)
  p = 1 - r0_cum(find(r_interval >= r, 1));
end

end