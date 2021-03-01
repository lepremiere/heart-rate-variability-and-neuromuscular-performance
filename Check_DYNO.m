function DYNOdata = Check_DYNO(varargin)
% This function calculates heart rate variability measures for a nx1 array of RR-intervals
% stored in a ".txt" file. Results will be visualized and files can be selected for 
% further processing.
%
% USAGE
% E.g.:
% ITTdata = Check_ITT('key', 'h')
%
% INPUT
% key:      Character that defines which button press indicates a stimulation. Default: "g".
%
% OUTPUT
% ITTdata:  Cell array (3xn) with filenames, -paths, and data for each selected file will be
%           returned. Press button next to accept file or discard to drop it.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 01.September.2020
% lepremiere
%---------------------------------------------------------------------------------------------------

% Creating input parser and setting default values
default_key = 'g';
p = inputParser;
addParameter(p, 'key',  default_key,   @(x) ischar(x));
parse(p, varargin{:});
ip = p.Results;

% Accessing required functions and selecting files
addpath(string(pwd) + '\Functions');
global F f i ii                          % Global variables to interact with Callback of figure
files                  = loadFiles();    % Getting files via GUI. Requires external function 'uipickfiles.m'
F                      = [];
screeningWindow        = 5000;           % Window width on data to be plotted. In ms
h                      = figure(1);

%% Looping through files
for ii = 1:length(files(:,1))
    f = [];
    data        = files{ii,3};                                    % Getting raw data
    Fs          = 1/data.Torque.interval;                         % Getting sampling frequency
    Fc          = 20;                                             % Cut-off frequency
    [B,A]       = butter(3,(Fc/(Fs/2)/0.8022),'low');             % Calculating coefficients of 3rd order lowpass butterworth filter
    torque      = filtfilt(B,A,data.Torque.values);               % Filtering torque with zero phase butterworth filter
    key_presses = data.Keyboard.codes(:,1) == double(ip.key);
    stim_time   = round((fix(data.Keyboard.times(key_presses)*Fs)),0); % Getting the indices of stimulations indicated by keyboardpresses
    time        = [1:length(torque)]./Fs';                        % Getting time in ms or n/Fs    
    max_stim_w  = 250;                                            % Setting the window following stimulation to analyse peak to peak torque in
    min_stim_w  = 50;                                             % Setting the window following stimulation to analyse peak to peak torque in
    rfdWindow   = 10;
    
    % Error catch function if keypresses do not contained desired key
    if(length(stim_time) == 0)
        msg = sprintf('No keypress "%s" found \nAvailable keys: %s', ip.key, native2unicode(unique(data.Keyboard.codes)));
        error(msg);
    end
    
    % Looping through stimulation pairs (superimposed, resting) 
    for i = 1:floor(length(stim_time)/2)

        index   = stim_time(2*i-1);
        index_2 = stim_time(2*i);

        twitchTorqueMVC     = range([min(torque(index:index + min_stim_w)),...                  % Peak to peak torque within superimposed stimulation + window 
                                     max(torque(index:index + max_stim_w))]);    
        twitchTorqueRest    = range([min(torque(index_2:index_2 + min_stim_w)),...              % Peak to peak torque within resting stimulation + window
                                     max(torque(index_2:index_2 + max_stim_w))]); 
        MVC                 = range([max(torque(index - 5000:index)),...                        % Peak to peak torque within stimulation - window*20 (~5s) 
                                     min(torque(index_2 - min_stim_w:index_2))]);
        VA                  = (1 - (twitchTorqueMVC*(torque(index)/MVC)/twitchTorqueRest))*100; % VA = (1 - (superimposed twitch*(Tb/MVC)*resting twitch^-1)) *100.
        rest_contraction    = torque(index_2:index_2 + max_stim_w*2);
        
        % Rate of force development of the resting twitch
        for j = 1:length(rest_contraction) - rfdWindow
           rfd(j) = (rest_contraction(j+rfdWindow) - rest_contraction(j))/rfdWindow*Fs;
        end
        [peakRFD_rest, ind_rfd] = max(rfd);  
        ind_rfd = ind_rfd + index_2;

        %% Plot
        clf
        hold on
        index   = stim_time(2*i-1) - screeningWindow;
        index_2 = stim_time(2*i)   + screeningWindow/2;
        h = plot(time(index:min([index_2, length(time)])),...  % Plots data of the specific stimulation
                      torque(index:index_2), 'LineWidth', 2);
        
        % Plotting markers that indicate the data points that are used for calculation
        for j = [0, 1]
            [~, ind]  =  min(torque(stim_time(2*i-j):stim_time(2*i-j) + min_stim_w));
            h(1) = plot(time(ind + stim_time(2*i-j)), torque(stim_time(2*i-j) + ind), 'xg', 'MarkerSize', 8, 'LineWidth', 2);
            [~, ind]  = max(torque(stim_time(2*i-j):stim_time(2*i-j) + max_stim_w));
            h(2) = plot(time(ind + stim_time(2*i-j)), torque(stim_time(2*i-j) + ind), 'xr', 'MarkerSize', 8, 'LineWidth', 2);
            h(3) = plot(time(ind_rfd), torque(ind_rfd), 'xc', 'MarkerSize', 8, 'LineWidth', 2);
            if(j == 1)
                [~, ind]  = max(torque(stim_time(2*i-j)-2000:stim_time(2*i-j)));
                h(4) = plot(time(stim_time(2*i-j)-2000+ind), torque(stim_time(2*i-j)-2000+ind), 'xb', 'MarkerSize', 8, 'LineWidth', 2);
            end
        end
        
        % Plot behavior
        xlim([index, index_2]/Fs);
        ylim([min(torque(index:index_2))-50, max(torque(index:index_2))*1.1]);
        titl = sprintf('File: %s \n MVT: %.0f Nm, Superimposed: %.0f Nm, Resting: %0.f Nm, Peak RTD: %.0f Nm/s VA: %.0f %% ',...
            [char(files{ii,2}), '_', num2str(i)], MVC, twitchTorqueMVC, twitchTorqueRest, peakRFD_rest, VA );
        title(titl, 'Interpreter', 'none');
        xlabel('Time [s]');
        ylabel('Torque [Nm]');
        xline(time(stim_time(2*i-1)),'','Superimposed','LabelHorizontalAlignment','center','LabelVerticalAlignment','mid','FontWeight','bold');     % Sets vertical line to first (superimposed) stimulation
        xline(time(stim_time(2*i)),'','Resting','LabelHorizontalAlignment','center','LabelVerticalAlignment','mid','FontWeight','bold');            % Sets vertical line to second (resting) stimulation
        legend(h(1:4), {'Min', 'Max', 'PeakRFD', 'PeakTorque'})
        ax = gca;
        ax.FontWeight = 'bold';
        ax.LineWidth = 2;
        hold off
        
        %% Control loop
        c = uicontrol('Position',[10 80 75 50],'String', 'Accept',  'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'g', 'Callback','uiresume(gcf);');                        % Callback for button accept. uiresume stops script temporarily
        d = uicontrol('Position',[10 5 75 50], 'String', 'Discard', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'r', 'Callback','global f i; f = [f,i];uiresume()');      % Callback for button discard. Saves discarded file for later removement
        uiwait();

        if(~ishandle(h) == 1)
           error('Figure has been closed')     
        end
    end
    % Filters Data
    data.Keyboard.times([f+f-1, f+f]) = [];         % Deletes all discarded keyboard presses/stimulations pairs
    files{ii,3} = data;                             % Replaces input file with cleaned data
    
    if(isempty(data.Keyboard.times))
        F = [F,ii];
    end
end
%% Create outputs
files(F,:) = [];                                    % Deletes all input files that do not conatain any selected stimulation
clear('global')
DYNOdata = files;                                    % Returns cleaned files
end
