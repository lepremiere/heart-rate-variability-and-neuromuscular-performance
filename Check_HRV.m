function HRVdata = Check_HRV(varargin)
% This function calculates heart rate variability measures for a nx1 array of RR-intervals
% stored in a ".txt" file. Results will be visualized and files can be selected for 
% further processing.
%
% USAGE
% E.g.:
% HRVdata = Check_HRV() -> Uses default parameters
% HRVdata = Check_HRV('artefact_recognition', 1, 'artefact_threshold', 200, 'detrending', 1) -> Sets parameters
%
% INPUT
% artefact_recognition: 0 = no recognition, 1 = recognition with moving median +- artefact_threshold
% artefact_threshold:   Determines the recognition threshold. In ms.
% detrending:           0 = no detrending, 1 = detrending with 3rd order polynomial,
%                       2 = zero phase butterworth filter cut-off frequency: 0.035 Hz
%
% OUTPUT
% HRVdata:  Cell array (3xn) with filenames, -paths, and data for each selected file will be
%           returned. Press button next to accept file or discard to drop it.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 01.September.2020
% lepremiere
%---------------------------------------------------------------------------------------------------

% HRV options
default_artefact_recognition = 1;     % 0 = no recognition, 1 = recognition with moving median +- artefactThreshold
default_artefact_threshold   = 200;   % Threshold in ms
default_detrending           = 1;     % 0 = no detrending, 1 = detrending with 3rd order polynomial, 2 = zero phase butterworth filter cut-off frequency: 0.035 Hz
default_median_length        = 14;    % Determines how many points will be used for calculating the median

% Creating input parser and setting default values
p = inputParser;
addParameter(p, 'artefact_recognition',  default_artefact_recognition,   @(x) x == 1 || x == 0);
addParameter(p, 'artefact_threshold',    default_artefact_threshold,     @(x) x >= 1);
addParameter(p, 'detrending',            default_detrending,             @(x) x >= 0 && x < 3);
addParameter(p, 'median_length',         default_median_length,          @(x) x > 0);
parse(p, varargin{:});
ip = p.Results;

% Accessing required functions and selecting files
warning off;
addpath(string(pwd) + '\Functions');    % Includes functions folder for access
global F ii;                            % Global to interact with Callbacks of the figure
files   = loadFiles();                  % Getting files via GUI. Requires external function 'uipickfiles.m'
F       = [];
h       = figure(1);

%% Looping through files
for ii = 1:length(files(:,1))

    clearvars -except files plots artefact_recognition detrending  ii artefact_threshold input median_length F h ip % Cleaning up after each iteration

    rawdata = table2array(files{ii,3});         % Getting raw data
    data    = rawdata;
    Fs      = 1000;                             % Sampling frequency of HR monitor
    Fi      = 4;                                % Interpolation frequency for equidistant data
    time    = cumsum(data)./Fs;                 % Time of RR-recording
    x       = time(1):1/Fi:time(end);           % Time of interpolation
    movm    = movmean(data, ip.median_length);   % Moving median. Artifact recognition identifies relying on this value. Window width = 10

    % Artifact recognition. Replaces hits with NaN
    if(ip.artefact_recognition == 1)
        for i = 1:length(data)
            if(range([data(i) , movm(i)]) >= ip.artefact_threshold)
                data(i) = NaN;
            else
                data(i) = data(i);
            end
        end

        artefacts = sum(isnan(data))/length(data)*100;        % Calculating artifact ratio in percent
        data_spline = spline(time, data, x);                  % Cubic spline interpolation

        % Replacing NaNs with cubic spline interpolation
        for i = 1:length(data)
            if(isnan(data(i)))
                if(i == 1 || i == length(data))
                    data(i) = nanmean(data);                  % Taking mean if 1st value must be replaced. High starting variability of spline method 
                else
                    data(i) = data_spline(find(x >= time(i), 1 ,'first'));
                end
            end
        end
        data_ac = data;
    else
        data_ac = rawdata;
        artefacts = 0;
    end

    % Calculating HR and mean RR before detrending
    HR      = 60/(mean(data)/Fs);
    mean_RR = mean(data);

    % Detrending data 
    % DT = 1 applies 3rd order polynomial detrend. DT = 2 applies zero phase butterworth filter 
    switch ip.detrending

        case 1    
        data = detrend(data,3);
        trend = data_ac - data;

        case 2   
        fc          = 0.035;                        % Cutt-off frequency
        fs          = abs(1/(mean(data)/Fs));       % Average sampling (RR-interval) frequency
        [B,A]       = butter(3,fc/(fs/2)*0.8022);   % Filter coefficients
        data        = data - filtfilt(B,A,data);    % Zero phase filtering
        trend       = mean(data_ac) + filtfilt(B,A,data);
    end

    % Calculating time domain and non-linear indices
    % Note: Non-linear indices are calculated with RMSSD since it is equivalent to SDSD (std(SSD)) if signal is stationary
    SSD     = data(2:end) - data(1:end-1);      % Successive differences
    rMSSD   = sqrt(nanmean(SSD.^2));            % Root mean square of successive RR intervals
    sDNN    = std(data);                        % Standard deviation of RR intervals
    
    %% Non-Linear indices
    adj_RR = [data(1:end-1), data(2:end)];
    SD1 = sqrt(1/2 * var(adj_RR(:,1) - adj_RR(:,2))); % Standard deviation perpendicular to line of identity in poincare plot
    SD2 = sqrt(1/2 * var(adj_RR(:,1) + adj_RR(:,2))); % Standard deviation along the line of identity in poincare plot

    %% Calculating frequency domain indecies
    data    = spline(time,data,x);
    n       = length(data);         
    [pxx,f] = pwelch(data./Fs,min([500*Fi,n]),min([500*Fi,n])/2,500*Fi^2,Fi);       % Returns twosided Welch's power density estimate for frequency vector f. 
                                                                                    % Inputs: x, segmentation window, number of overlap, number of discrete fourier transform (DFT) points, sampling frequency 
    % Calculating power in the specific frequency bands and LF/HF ratio (TP 0 - 0.4, VLF 0 - 0.04, LF 0.04 -  0.15, HF 0.15 - 0.4)
    p(1) = trapz(f(find(f > 0.000 & f <= 0.4))  ,pxx(find(f > 0.000 & f <= 0.4)) .*1000^2);
    p(2) = trapz(f(find(f > 0.000 & f <= 0.04)) ,pxx(find(f > 0.000 & f <= 0.04)).*1000^2);
    p(3) = trapz(f(find(f > 0.04  & f <= 0.15)) ,pxx(find(f > 0.04  & f <= 0.15)).*1000^2);
    p(4) = trapz(f(find(f > 0.15  & f <= 0.4))  ,pxx(find(f > 0.15  & f <= 0.4)) .*1000^2);
    p(5) = p(3)/p(4);

    %% Concatenates results 
    Results = [floor(HR),floor(mean_RR), sDNN, rMSSD, p(1), p(2), p(3), p(4), p(5), SD1, SD2, SD2/SD1, artefacts];

    %% Plots
    legend_names = {'Raw data', 'Median (10)', 'Corrected data', 'Median + threshold ', '3rd Order polynomial'};  
    clf
    leg = [];                                                   % Legend selection helper
    
    % Raw data
    subplot(2,4,[1:4]);
    hold on
    h(1) = plot(rawdata, '-x');                                 % Raw data             
    leg = [leg, 1];
    xlabel('Time (s)')
    ylabel('RR (ms)')
    title(files(ii,2), 'Interpreter', 'none')
    hold off

    % Artefact recognition
    if(ip.artefact_recognition == 1)
    hold on
    h(2) = plot(movm,'--');                                     % Moving median
    h(3) = plot(data_ac,'-o');                                  % Artefact corrected data
    h(4) = plot(movm - ip.artefact_threshold,'k');              % Moving median threshold boundaries
    plot(movm + ip.artefact_threshold,'k');
    title(['File: ', files{ii,2}], 'Interpreter', 'none');
    hold off
    leg = [leg, 2, 3, 4];
    end

    % Detrending
    if(ip.detrending > 0)
        hold on
        h(5) = plot(trend,'-r');                                % Trend that has been removed
        leg = [leg, 5];
        if(ip.detrending == 2)
           legend_names{end} = '0.035Hz Butterworth filter';
        hold off 
        end
    end
    
    % Legend
    names = {'Raw data', 'Median (10)', 'Corrected data', 'Median + threshold '};
    legend(h(leg), legend_names(leg), 'Location', 'southeast');
    
    %% Poincare plot
        
    mu =  mean(adj_RR, 1)';                                     % Get the center of the adjacent RR data
    
    % Projection of a point (x) onto a line (g) can be represented as:
    % P_g(x) = r0 + lambda * u, with r0 and u beeing local vectors of g
    r0 = [adj_RR(1,1); adj_RR(1,1)];                            % A vector on the line of identity
    u0 = r0 - [adj_RR(end,1); adj_RR(end,1)];                   % Another vector on the line of identity
    v = adj_RR';                                                % Adjacent RR data that will be projected onto line of identity
    lambda0 = ((v - r0)'*u0)./(norm(u0)^2);                     % Vector containing lambda for every element in v
    P0 = r0 + lambda0'.*u0;                                     % Projections of adjacent RR data onto line of identity
    
    D = (P0(:,1) - v(:,1))/norm(P0(:,1) - v(:,1));              % Normalized difference vector for a single data point, to get a vector
                                                                % perpendicular to the line of identity
    
    vert_line = max(vecnorm(P0-v));                             % Getting length for the line orthogonal to line of identity
    D_line = [mu + D * vert_line, mu - D * vert_line];          % Defining start and end of the line orthogonal to line of identity

    SD1_line = sort([mu , mu + D * SD1], 1);                             % Defining start and end of the SD1 line
    SD2_line = [mu , mu + abs(u0/norm(u0) * SD2)];              % Defining start and end of the SD2 line
    
    % Plot
    subplot(2,4,7)
    scatter(adj_RR(:,1), adj_RR(:,2), '.');                     % Adjacent RR data
    tit = sprintf('Poincare Plot');
    title(tit);
    hold on
    plot(adj_RR(:,1), adj_RR(:,1), '-k', 'LineWidth', 2);       % LOI
    plot(D_line(:, 1), D_line(:,2), '-k', 'LineWidth', 2);      % Line perpendicular to LOI
    k(1) = plot(SD1_line(1, :), SD1_line(2, :), '-g', 'LineWidth', 2);  % SD1 line
    k(2) = plot(SD2_line(1, :), SD2_line(2, :), '-r', 'LineWidth', 2);  % SD2 line
    
    angle = [D/norm(D), u0/norm(u0)];                           % Getting the angle between basis and LOI
    transform = [1 0; 0 1] * inv(angle);                        % Creating transform matrix for this angle
    theta = [0 : 0.01 : 2*pi];                                  % Array of angles for which the ellipse will be plotted
    X = mu  + ([SD1 * cos(theta) ; SD2 * sin(theta) ]' *transform )';   % Defining the ellipse coordinates
    plot(X(1,:), X(2,:), 'k', 'LineWidth', 2);                  % Ellipse
    legend([k(1), k(2)], {['SD1: ', num2str(round(SD1, 1))], ['SD2: ', num2str(round(SD2,1))],}, 'Location', 'southeast');
    xlim([min(data), max(data)]);
    ylim([min(data), max(data)]);
    axis square;
    grid on;

    %% Frequency domains
    subplot(2,4,[5:6])
    hold on
    plot(f,pxx,'k','HandleVisibility','off');
    area(f(find(f>=0.0 & f<=0.04)),pxx(find(f>=0.00 & f<=0.04)),'FaceColor',[0.7 0.7 0.7]);             % Area for VLF
    area(f(find(f>0.04 & f<=0.15)),pxx(find(f>0.04 & f<=0.15)),'FaceColor',[1 0 0],'FaceAlpha',0.4);    % Area for LF
    area(f(find(f>0.15 & f<=0.4)),pxx(find(f>0.15 & f<=0.4)),'FaceColor',[0 0.9 0],'FaceAlpha',0.4);    % Area for HF
    xlabel('Frequency (Hz)');
    ylabel('PSD (s�/Hz)');
    title('FFT spectrum');
    xline(0.04);
    xline(0.15);
    xline(0.4);
    ylim([0 max(pxx(30:end))*1.1]);
    xlim([0 0.5]);
    legend('VLF (0-0.04 Hz)','LF (0.04-0.15 Hz)','HF (0.15-0.4 Hz)');
    hold off

    %% Countinuous wavelet function
    subplot(2,4,8)
    helperPlotScalogram3d(data,Fi);                             % Continuous wavelet function. MatLab function
    ylim([0 0.6]);
    view(290,45);
    grid on;
    
    %% Control loop
    c = uicontrol('Position',[10 80 75 50],'String', 'Accept',  'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'g', 'Callback','uiresume(gcf);');                        % Callback for button accept. uiresume stops script temporarily
    d = uicontrol('Position',[10 5 75 50], 'String', 'Discard', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'r', 'Callback','global F ii; F = [F,ii];uiresume()');    % Callback for button discard. Saves discarded file for later removement
    uiwait();
    if(~ishandle(h) == 1)
       error('Figure has been closed')     
    end

end

files(F,:)=[];      % Deletes all input files that were discarded
clear('global')
HRVdata = files;    % Returns accepted set of files
end
