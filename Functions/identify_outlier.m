function [inds] = identify_outlier(RR)
% This function resembles the outlier identification algorithm by 
% Lipponen & Taravainen (2019). https://doi.org/10.1080/03091902.2019.1640306.
% Number in paranthesis indicate the formula of the paper it is referring
% to. E.g.: (1), (2).
% Additionally, MATLABs native function "isoutlier" is applied for further
% enhancement of artefact recognition.
%
% USAGE
% E.g.:
% [inds] = identify_outlier(RR, verbose)
%
% INPUT
% RR:       RR-intervals stored in an array (n x 1).
%
% OUTPUT
% inds:     Logical array (n x 1) that indicates which values
%           of the input array is NOT corrupted.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 20.January.2021
% lepremiere
%---------------------------------------------------------------------------------------------------

time = cumsum(RR);  
n = length(RR);

% Parameters suggested by authors
alpha = 5.2;    
c1 = 0.13;
c2 = 0.17;
k = 45;                                % Window width 

dRRs = [0; RR(2:end) - RR(1:end-1)];   % Delta values. (1)

for i = 1:n
    % Calculating threshold 1 (Th1). (2)
    if(i < length(dRRs) - k)
        if(i < k+1)
            Th1(i,:) = alpha * iqr(dRRs(1:i+k), 'all');
        else
            Th1(i,:) = alpha * iqr(dRRs(i-k:i+k), 'all');
        end
    else
        Th1(i,:) = alpha * iqr(dRRs(i-k:end), 'all');
    end
    
    % Difference vector between RR(i) and 11-beat median interval (3)
    if(i < length(RR) - 5)
        if(i < 6)
            mRRs(i,:) = RR(i) - median(RR(1:i+5));
        else
            mRRs(i,:) = RR(i) - median(RR(i-5:i+5));
        end
    else
        mRRs(i,:) = RR(i) - median(RR(i-5:end));
    end  
    
    % Special case for mRR < 0. (5)
    if(mRRs(i) < 0)
        mRRs(i) = 2 * mRRs(i);
    end  
end

% Calculating threshold 2 (Th2). (6)
for i = 1:n
    if(i < n - k)
        if(i < k+1)
            Th2(i,:) = alpha * iqr(mRRs(1:i+k), 'all');
        else
            Th2(i,:) = alpha * iqr(mRRs(i-k:i+k), 'all');
        end
    else
        Th2(i,:) = alpha * iqr(mRRs(i-k:end), 'all');
    end   
end

dRR = dRRs./Th1;    % (3)
mRR = mRRs./Th2;    % (7)

S11 = dRR;          % (8)
S21 = dRR;          % (11)

% Calculating subpsace S12. (9)
for i = 1:length(RR)
    if(i < 2 && i < length(RR) - 1)      
        S12(i,:) = dRR(i+1);      
    elseif(i > length(RR) - 1)
        S12(i,:) = dRR(i-1);      
    else
        if(dRR(i) > 0 )
            S12(i,:) = max([dRR(i-1), dRR(i+1)]);      
        else
            S12(i,:) = min([dRR(i-1), dRR(i+1)]);     
        end 
    end
end

% Calculating subpspace S22. (12)
for i = 1:length(RR)
    if(i < 2 && i < length(RR) - 2)      
        S22(i,:) = dRR(i+2);      
    elseif(i > length(RR) - 2)
        S22(i,:) = dRR(i-1);      
    else
        if(dRR(i) >= 0 )
            S22(i,:) = min([dRR(i-1), dRR(i+2)]);      
        else
            S22(i,:) = max([dRR(i-1), dRR(i+2)]);     
        end 
    end
end

Eqs = [];   % Helper variable to gather results of the equations

for i = 1:n
    % Creating logicals during the first step of the decision algorithm.
    % Figure 1, top panel.
    if(abs(dRR(i)) > 1)
        [Eq1, Eq2]       = equations(S11(i), S12(i), c1, c2);
        if(Eq1 == true || Eq2 == true)
             Eqs(i,1) = true;
        else
            Eqs(i,1) = false;
        end
    else
        Eqs(i,1) = false;
    end
   
    % Creating logicals during the first step of the decision algorithm.
    % Figure 1, second panel.
    if((abs(dRR(i)) > 1 || abs(mRR(i)) > 3) && Eqs(i,1) == false && i < n-2)
        [Eq3, Eq4, Eq5] = equations2(dRR(i:i+2), mRR(i));
        if(Eq3 == true || Eq4 == true || Eq5 == true)
            Eqs(i,1) = true;
        else
            Eqs(i,1) = false;
        end
    end
end

RR(Eqs == true) = deal(NaN);        % Removing identified outlier by algorithm of Lipponen et al.
RR(isoutlier(RR)==1) = deal(NaN);   % Removing outlier by MATLABs native function "isoutlier()"
inds = ~isnan(RR);                  % Returning the indeces of not-corrupted RR-intervals

% Debug plot
verbose = 0;
if(verbose == 1)
    figure(1)
    clf
    subplot(3,1,1)
    hold on
    plot(time, RR)
    scatter(time(isoutlier(RR)), RR(isoutlier(RR)));
    scatter(time(Eqs==1), RR(Eqs==1));

    subplot(3,1,2)
    hold on
    plot(time, abs(dRR))
    yline(1);
    ylim([0, 1.2]);

    subplot(3,1,3)
    hold on
    plot(time, abs(mRR))
    yline(1);
    ylim([0, 1.2]);

    figure(2)
    clf
    subplot(2,1,1)
    hold on
    scatter(S11, S12)
    xline(0);
    yline(0);
    
    subplot(2,1,2)
    hold on
    scatter(S21, S22)
    xline(0);
    yline(0);
end
end

% These functions are the mathematical function stated out in figure 1 of
% of the paper and are labeled equally.
function [Eq1, Eq2] = equations(S11, S12, c1, c2)
    if(S11 > 1 && S12 < -c1 * S11 + c2)
        Eq1 = true;
    else
        Eq1 = false;
    end

    if(S11 < -1 && S12 > c1 * S11 - c2)
        Eq2 = true;
    else
        Eq2 = false;
    end
end

function [Eq3, Eq4, Eq5] = equations2(dRR, mRR)
    if(sign(dRR(1))*dRR(2) < -1)
        Eq3 = true;
    else
        Eq3 = false;
    end
    if(abs(mRR) > 3)
        Eq4 = true;
    else
        Eq4 = false;
    end
    if(sign(dRR(1)) * dRR(3) < -1)
        Eq5 = true;
    else
        Eq5 = false;
    end
end

% This function will not be used, as it is part of a decision process how
% to replace outliers. In this study all outliers were replaced by spline
% interpolation.
function [Eq6, Eq7] = equations3(RR, medRR, Th2)
    if(abs(RR(1)/2 - medRR) < Th2)
        Eq6 = true;
    else
        Eq6 = false;
    end
    
    if(abs(RR(1)+RR(2) - medRR) < Th2)
        Eq7 = true;
    else
        Eq7 = false;
    end    
end
