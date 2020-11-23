function y = corrdist(r, ro, n) 
%This function computes the probability density function for the
%correlation coefficient of a bivariate random variable.
%
% USAGES
% y = corrdist(r, ro, n)
%
% INPUT
% r:    Vector of possible correlation random variables, i.e. the values at
%       which the pdf is evaluated at.
%
% ro:   The given (true) correlation coefficient, i.e. the population
%       correlation coefficient. length(ro) > 1 supported.
%
% n:    The number of samples in the correlated data. Only length(n) = 1
%       supported.
% 
% OUTPUT
% y:    The probability density function for r, given ro, for n data
%       samples of a bivariate normal distribution.
%
%-----------------------------------------------------------------------
% Latest Edit: 11.June.2012
% Joshua D Carmichael
% josh.carmichael@gmail.com
%
% Original Author: Xu Cui, Stanford University (retrieved 11.June.2012)
%-----------------------------------------------------------------------
%accept vectorized inputs.
if(length(ro)> 1.5)
    r   = repmat(r(:),1,length(ro));
    ro  = repmat(ro(:)', length(r),1);
end
if( n < 120 )
    
    y = (n-2) * gamma(n-1) * ((1-ro.^2).^((n-1)/2)).* (1-r.^2).^((n-4)/2);
    y = y./ (sqrt(2*pi) * gamma(n-1/2) * (1-ro.*r).^(n-3/2));
    y = y.* (1+ 1/4*(ro.*r+1)/(2*n-3) + 9/16*(ro.*r+1).^2 / (2*n-1)/(2*n+1));
    
else
    
    y = (n-2) * (1-ro.^2)^((n-1)/2) * (1-r.^2).^((n-4)/2);
    y = y./ (sqrt(2*pi) * (1-ro.*r).^(n-1/2)) * n.^(-1/2);
    y = y.* (1+ 1/4*(ro.*r+1)/(2*n-1) + 9/16*(ro.*r+1).^2 / (2*n-1)/(2*n+1));
    
end
y(r>1)              = 0;
y(r<-1)             = 0;
y(~isfinite(y))     = 0;
return;