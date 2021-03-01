function Results = HRV_Analysis(rawdata, AC, DT, threshold)
% This function calculates heart rate variability indices of RR-intervals.
%
% USAGE
% E.g.:
% Results = HRV_Analysis(HR_data, 1, 1, 200)
%
% INPUT
% rawdata:      Must be a table (n x 1) containing only RR-intervals in ms.
%
% AC:           0 = no artefact correction will be performed. 
%               1 = artefacts will be detected by falling out of the mean +- a threshold. 
%               2 = artefacts will be recognized by three mechanisms.
%                   The artefact recognition algorithmn by Lipponen & Tarvainen 
%                   (2019, https://doi.org/10.1080/03091902.2019.1640306) was rebuild and
%                   is applied. Additionally, Matlab's 'isoutlier()' function is applied.
%                   Lastly, a k-means-based algorithm was utilized to fit two gamma
%                   distributions which provide the final threshold value used alongside
%                   the a moving mean. (Default)
%
% DT:           0 = no detrending will be performed.
%               1 = data will be detrended with a third order polynomial.
%               2 = data will be detrended with a zero-phase Butterworth filter,
%                   with a cut-off frequency at 0.035 Hz.
%
% threshold:    Determines the threshold for artefact correction. In ms.
%
% OUTPUT
% Results:      Array (1 x 13) containing: Heart rate, mean RR-interval,
%               standard deviation of RR-intervals, root mean square of successive
%               differences, total power, very low frequency power, 
%               low frequency power, high frequency power, low/high frequency
%               power ratio, standard deviation along line of identity, standard 
%               deviation rectangular to line of identity, ratio of SD along 
%               and rectangular to line of identity, artefacts as percentage
%               of total sample points.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 02.January.2021
% lepremiere
%---------------------------------------------------------------------------------------------------

%% Setup
% Note: Input (data) must be 1-by-n table with RR-intervals 
if(istable(rawdata) == 1)
    data    = table2array(rawdata);
else
    data = rawdata;
end

N   = length(data);
Fs  = 1000;                                                                 % Sampling frequency of HR monitor
Fi  = 10;                                                                   % Interpolation frequency for equidistant data
raw_time_idx = cumsum(data)./Fs;
time_idx = raw_time_idx;

% Artifact recognition. Replaces hits with NaN
if(AC > 0)
    if(AC == 1)
        movm = movingmedian(data, 30);                                      % Moving median
        idx = sqrt((data - movm).^2) >= threshold;                          % Identifying outlier
        data(idx) = deal(NaN);                                              % Replacing outlier
        
    elseif(AC == 2)
        inds = identify_outlier(data);                                      % Identifying outlier by custom algorithm
        tmp_data = data;
        tmp_data(inds == 0) = deal(NaN);
        movm = movingmean(tmp_data, 30);                                    % Calculating moving mean without first sweep of outliers
        [threshold, two_means] = get_threshold(data(inds), 0);              % Getting threshold for final sweep
        idx = sqrt((data - movm).^2) < threshold;                           % Getting indeces of values outside of moving mean  +- threshold
        idx2 = [1; sqrt((data(2:end) - data(1:end-1)).^2) < 2*threshold];   % Getting indeces of deltas bigger than  2 times the threshold
        idx(sum(idx + idx2 + inds) == 3) = deal(1);                         % Creating logicals for values that were not identified as outlier  
        data(idx==0) = deal(NaN);                                           % Excluding all values that were identified as outlier
    else
        error('Invalid AC option');
    end

    % Replacing NaNs with cubic spline interpolation
    x       = raw_time_idx(1):1/Fi:raw_time_idx(end)+1;                     % Time of interpolation
    data_spline = spline(time_idx(~isnan(data)), data(~isnan(data)), x);    % Cubic spline interpolation
    artefacts = sum(isnan(data));                                           % Counting artefacts
    
    for i = 1:length(data)
        if(isnan(data(i)) == 1)
            k = find(x >= time_idx(i), 1, 'first');                         % Finding matching time indeces
            if(range([data_spline(k) movm(i)]) < threshold)                 % Replacing NaN with interpolation if not ouside boundaries
                data(i) = data_spline(k);
            else
                data(i) = NaN;                                              % NaN if interpolation not possible
            end
        end
    end
    
    artefacts   = artefacts + sum(isnan(data));                             % Counting artefacts
    time        = time_idx(~isnan(data));                                   % New time indeces for frequency analysis
    x           = time(1):1/Fi:time(end);                                   % New time vector for interpolation
    data        = data(~isnan(data));                                       % Removing all NaNs from array
   
else
    artefacts = 0;
end

%% Calculating HR and mean RR before detrending
HR      = 60/(nanmean(data)/Fs);
RR = nanmean(data);

%% Detrending data 
% DT = 1 applies 3rd order polynomial detrend. DT = 2 applies zero phase butterworth filter 
switch DT
    case 1    
        data        = detrend(data, 3);
    case 2   
        fc          = 0.035;                        % Cutt-off frequency
        fs          = abs(1/(mean(data)/Fs));       % Average sampling (RR-interval) frequency
        [B,A]       = butter(3, fc/(fs/2)*0.8022);  % Filter coefficients
        data        = data - filtfilt(B ,A, data);  % Zero phase filtering        
end

%% Calculating time domain and non-linear indices
% Note: Non-linear indices are calculated with RMSSD since it is equivalent to SDSD (std(SSD)) if signal is stationary
SSD     = data(2:end) - data(1:end-1);      % Successive differences
rMSSD   = sqrt(nanmean(SSD.^2));            % Root mean square of successive RR intervals
sDNN    = nanstd(data);                     % Standard deviation of RR intervals

%% Non-linear
SD1 = sqrt(1/2 * var(data(1:end-1) - data(2:end))); % Standard deviation perpendicular to line of identity in poincare plot
SD2 = sqrt(1/2 * var(data(1:end-1) + data(2:end))); % Standard deviation along the line of identity in poincare plot

%% Calculating frequency domain indecies
data    = spline(time, data, x);
n       = length(data);         
[pxx,f] = pwelch(data./Fs, min([500*Fi,n]), min([500*Fi,n])/2, 500*Fi^2, Fi);   % Returns two-sided Welch's power density estimate for frequency vector f. 
                                                                                % Inputs: x, segmentation window, number of overlap, number of discrete fourier transform (DFT) points, sampling frequency 
% Calculating power in the specific frequency bands and LF/HF ratio (TP 0 - 0.4, VLF 0 - 0.04, LF 0.04 -  0.15, HF 0.15 - 0.4)
p(1) = trapz(f(find(f > 0.000 & f <= 0.4))  ,pxx(find(f > 0.000 & f <= 0.4)) .*1000^2);
p(2) = trapz(f(find(f > 0.000 & f <= 0.04)) ,pxx(find(f > 0.000 & f <= 0.04)).*1000^2);
p(3) = trapz(f(find(f > 0.04  & f <= 0.15)) ,pxx(find(f > 0.04  & f <= 0.15)).*1000^2);
p(4) = trapz(f(find(f > 0.15  & f <= 0.4))  ,pxx(find(f > 0.15  & f <= 0.4)) .*1000^2);
p(5) = p(3)/p(4);

%% Concatenates results 
Results = [floor(HR), floor(RR), sDNN, rMSSD, p(1), p(2), p(3), p(4), p(5), SD1, SD2, SD2/SD1, (artefacts/N)*100];

end

