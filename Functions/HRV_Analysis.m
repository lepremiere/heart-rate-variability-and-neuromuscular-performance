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
% AC:           Scalar with two options: 0 = no artefact correction will be
%               performed. 1 = artefacts will be detected by falling out of 
%               the median +- a threshold.
%
% DT:           Scalar with three options: 0 = no detrending will be performed.
%               1 = data will be detrended with a third order polynomial.
%               2 = data will be detrended with a zero-phase Butterworth filter,
%               with a cut-off frequency at 0.035 Hz.
%
% threshold:    Determines the threshold for artefact correction.
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
% Latest Edit: 09.July.2020
% lepremiere
%---------------------------------------------------------------------------------------------------

%% Setup
% Note: Input (data) must be 1-by-n table with RR-intervals 
data    = table2array(rawdata);
Fs      = 1000;                         % Sampling frequency of HR monitor. Preset = 1000 (1/1000s)
Fi      = 4;                            % Interpolation frequency for equidistant data
time    = cumsum(data)./Fs;             % Time of RR-recording
x       = time(1):1/Fi:time(end);       % Time of interpolation
movm    = movmedian(data, 10);          % Moving median. Artifact recognition identifies relying on this value. Window width = 10

%% Artifact recognition. Replaces hits with NaN
if(AC == 1)
    for i = 1:length(data)
        if(range([data(i) , movm(i)]) >= threshold)
            data(i) = NaN;
        else
            data(i) = data(i);
        end
    end   
    artefacts = sum(isnan(data))/length(data)*100;      % Calculating artifact ratio in percent
    data_spline = spline(time, data, x);                % Cubic spline interpolation

    % Replacing NaNs with cubic spline interpolation
    for i = 1:length(data)
        if(isnan(data(i)) && i < length(data))
            data(i) = data_spline(find(x >= time(i), 1, 'first'));
        end
        if(i == 1 || i == length(data))
            data(i) = nanmean(data);                    % Taking mean if 1st value must be replaced. Due to high starting variability of spline method 
        end
    end
else
    artefacts = 0;
end

%% Calculating HR and mean RR before detrending
HR      = 60/(mean(data)/Fs);
RR      = mean(data);

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

SD1 = sqrt(1/2 * var(data(1:end-1) - data(2:end)));
SD2 = sqrt(1/2 * var(data(1:end-1) + data(2:end)));
% SD1     = sqrt(0.5*rMSSD^2);                % Standard deviation perpendicular to line of identity in poincare plot
% SD2     = sqrt(2*sDNN^2 - 0.5*rMSSD^2);     % Standard deviation along the line of identity in poincare plot

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
Results = [floor(HR), floor(RR), sDNN, rMSSD, p(1), p(2), p(3), p(4), p(5), SD1, SD2, SD2/SD1, artefacts];

end
