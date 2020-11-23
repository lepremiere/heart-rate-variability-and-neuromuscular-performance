 function CMJdata = Check_CMJ()
% This function calculates jump qualities for a counter-movement jump from 
% a table containing vertical forces and time. Results will be visualized
% and files can be selected for further processing.
%
% USAGE
% E.g.:
% CMJdata = Check_CMJ()
%
% Selected files must be ".csv" of force plate data, with vertical ground
% reaction forces stored in variables "z1-z4" and time in last column.
%
% OUTPUT
% CMJdata:  Cell array (3xn) with filenames, -paths, and data for each selected file will be
%           returned. Press button next to accept file or delete to drop it.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 31.August.2020
% lepremiere
%---------------------------------------------------------------------------------------------------
warning off;
addpath(string(pwd) + '\Functions');    % Include folder "Functions" to access required functions
global f ii                             % Initializes global variables to work with during animation
f = []; 
files = loadFiles();                    % Getting files via GUI. Requires external function 'uipickfiles.m'
ax     = figure(1);                     % Initialize figure

%% Looping through files
for ii = 1:length(files(:,1))
    %% Determining total vertical ground reaction force (VGRF)
    % Note: Force plate data must be table with VGRF stored in columns named 'z1-z4' and time stored in last column. Time must be consistent
    data    = files{ii,3};
    Fz      = data.z1 + data.z2 + data.z3 + data.z4;

    % Determines the the sampling frequency 'Fs'
    Fs = 1/(data.(data.Properties.VariableNames{end})(11) - data.(data.Properties.VariableNames{end})(10));

    %% Filtering the data
    Fc          = 50;                                   % Cutt-off frequency
    [b,a]       = butter(3,(Fc/(Fs/2)/0.8022),'low');   % Calculating coefficients of second order lowpass butterworth filter. Adjusted by 0.8022
    Fz_filtered = filtfilt(b,a,Fz);                     % Filtering VGRF with zero phase lowpass butterworth filter 

    %% Adjusting VGRF in case force plate was relocated
    if (mean(Fz_filtered) < 50)
        minimum     = sortrows(Fz_filtered).*-1;        
        relocation  = mean(minimum(1:500));
        Fz_filtered = Fz_filtered + relocation;
    end

    %% Finding resting VGRF and calculating body weight
    % Note: Force plate data has to include a 3 second resting period
    interval = 100;     % Window width in ms

    % Calculating mean and standard deviation for n-windows for the first 2000 ms
    for i = 1:2000 - interval
        if(Fz_filtered(i,1) > 10)
            q(i,1) = std(Fz_filtered(i:i + interval));
            q(i,2) = mean(Fz_filtered(i:i + interval));
        else
            q(i,1) = inf;                       
            q(i,2) = inf;
        end
    end

    % Finding window with smallest variability (std), assuming to be the best baseline and calculating bodyweight from its mean VGRF
    % Also, sets variability as base line noise
    [row ~]   = find(q(:,1) == min(q(:,1)));
    BW        = mean(Fz_filtered(1:1500)./9.81);                      % m = F/a
    noise     = Fz - Fz_filtered;                 
    threshold = range(noise(row:row + interval));   % Sets variability as base line noise

    %% Determining start, takeoff, landing and airtime of CMJ
    % Getting orientation in the jump
    subzero     = find(Fz_filtered<10);                                                                          % Getting indices of no contact to force plate
    takeoff     = find(Fz_filtered(row:subzero(1)) < threshold,1,'first')+row;                                  % Determines takeoff as first VGRF smaller than baseline noise. +row accounts for subzero offset
    landing     = find(Fz_filtered(takeoff:end) > threshold,1,'first');                                         % Determines landing as first VGRF greater than baseline noise after takeoff
    airtime     = landing*1/Fs;                                                                                 % Landing ontains number of sampling points for the jump
    posStart    = find(Fz_filtered(row:subzero(1)-1) > Fz_filtered(row) + threshold,1,'first') + row - 1000;    % Determines start of movement by finding the first VGRF greater than baseline VGRF + noise. 1000 ms offset 
    negStart    = find(Fz_filtered(row:subzero(1)-1) < Fz_filtered(row) - threshold,1,'first') + row - 1000;    % Determines start of movement by finding the first VGRF greater than baseline VGRF - noise. 1000 ms offset
    start       = min([posStart,negStart]);                                                                     % Finds true start

    %% Calculatng net force and net impulse 
    netForce     = Fz_filtered - BW*9.81 ;                  % F(net) = VGRF - BW*g
    netImpulse   = trapz((netForce(start:takeoff)))*1/Fs;   % Integrating netForce with respect to time

    %% Calculating variables of interest
    jumpHeightAirtime   = (1/8*9.81*airtime^2)*100;                                             % h = 1/8*g*t^2. In cm (100)
    velocity            = cumtrapz((netForce)*1/Fs)/BW;                                         % v = delta(p)/m
    peakVelocity        = max(velocity(start:takeoff));                
    takeoffVelocity     = velocity(takeoff);
    jumpHeightImpulse   = (takeoffVelocity^2/2*9.81);                                       % h = v^2/2g. In cm (100)
    Power               = Fz_filtered.*velocity;                                                % W = F*v    
    peakPower           = max(Power(start:takeoff));
    peakForce           = max(Fz_filtered(start:takeoff));
    relpeakForce        = max(Fz_filtered(start:takeoff))/BW;

    turningPoint        = find(Fz_filtered == min(Fz_filtered(start:takeoff-250)),1,'first');   % Determines the start of increasing VGRF after counter-movement
    jump                = Fz_filtered(turningPoint:takeoff);
    rfd_interval        = 50;                                                                   % Determines the window width in ms for finding maximum rate of force development

    % Calculates the RFDs during jump with respect to time
    for i = 1:length(jump)-rfd_interval
       rFD(i) = ((jump(i+rfd_interval)-jump(i))/rfd_interval)*Fs;
    end

    % Finds max RFD and time to max RFD
    RFD         = max(rFD);
    [col1 row1] = find(rFD == RFD);                                                                                                                                         
    time2maxRFD = row1;

    %% Plot
    clf;
    [ax, h(1), h(2)] = plotyy([1:size(Fz_filtered, 1)], netForce, [1:size(Fz_filtered, 1)], velocity); % Plots filtered VGRF and velocity
    
    % Aligns y-axes origins
    h(2).Color = [1 0 0];
    ylim2 = get(ax(2),'Ylim');
    ratio = ylim2(1)/ylim2(2);
    ylim1 = get(ax(1),'Ylim');
    set(ax(1),'Ylim',[ylim1(2)*ratio ylim1(2)]);
    
    hold on 
    h1 = area([start:takeoff],netForce(start:takeoff)); % Plots impulse as filled area

    % Setting labels
    titl = sprintf('File: %s \nJumpheight: %.2f / %.2f cm, Peak Power: %.2f W', string(files(ii,2)), jumpHeightImpulse, jumpHeightAirtime, peakPower);
    title(titl, 'Interpreter', 'none');
    ylabel(ax(1), 'Force (N)');
    ylabel(ax(2), 'Velocity (m/s)');
    xlabel('Time (ms)');
    legend([h(1), h(2), h1], {'VGRF', 'Velocity', 'Impulse'});
    
    % Line markers for different events
    xline(row,'','Baseline','LabelHorizontalAlignment','center','FontWeight','bold');                                           % Baseline start
    xline(row + interval,'','LabelHorizontalAlignment','center','FontWeight','bold');                                           % Baseline end
    xline(start,'.',{'Start'},'LabelHorizontalAlignment','center','FontWeight','bold');                                         % Start of counter-movement
    xline(takeoff,'',{'Take off'},'LabelHorizontalAlignment','center','FontWeight','bold');                                     % Take off of platform
    xline(subzero(1) + landing,'',{'Landing'},'LabelHorizontalAlignment','center','FontWeight','bold');                         % Landing on platform
    xline(turningPoint,'',{'Turning Point'},'LabelHorizontalAlignment','center','FontWeight','bold');                           % Begin of breaking phase
    xline(find(Fz_filtered == peakForce,1,'first'),'',{'Peak Force'},'LabelHorizontalAlignment','center','FontWeight','bold');  % Peak Force
    xline(turningPoint + row1 + rfd_interval/2,'',{'Peak RFD'},'LabelHorizontalAlignment','center','FontWeight','bold');        % Peak rate of force developement
    xline(find(Power == peakPower),'',{'Peak Power'},'LabelHorizontalAlignment','center','FontWeight','bold');                  % Peak power

    % Callbacks to control loop and figure behavior. Selecting files.
    c = uicontrol('Position',[10 80 75 50],'String', 'Accept',  'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'g', 'Callback','uiresume(gcf);');
    d = uicontrol('Position',[10 5 75 50], 'String', 'Discard', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'r', 'Callback','global f ii; f = [f,ii];uiresume()');
    uiwait();
    if(~ishandle(ax) == 1)
       error('Figure has been closed')     
    end
end

% Returning selected files
files(f,:) = [];
clear('global')
CMJdata = files;
 end

