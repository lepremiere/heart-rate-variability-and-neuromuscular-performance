function movm = movingmean(x, k)  
% This function calculates a moving mean over the k/2 next and previous
% non-NaN values. Therefore, its temporal coherence is impaired.
%
% USAGE
% E.g.:
% movm = movingmean(x, 30) 
%
% INPUT
% x:    Array (n x 1) for which the moving mean is calculated. 
% k:    Scalar that determines the window width.
%
% OUTPUT
% movm: Array (n x 1) that contains the moving mean.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 20.January.2021
% lepremiere
%---------------------------------------------------------------------------------------------------

n = length(x);
movm = nan(n,1);

for i = 1:n
    if(i == 1)
        idx = find(~isnan(x(i:end)), k/2, 'first');     % Finding the first k/2 non-NaN values
        movm(i) = mean(x(idx));
    elseif(i == n)
        idx = find(~isnan(x(1:i)), k/2, 'last');        % Finding the last k/2 non-NaN values
        movm(i) = mean(x(idx));
    else
        idx = find(~isnan(x(i:end)), k/2, 'first');
        idx2 = find(~isnan(x(1:i)), k/2, 'last');
        movm(i) = mean(x([idx2;i-1+idx]));    
    end            
end
end