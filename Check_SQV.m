 function SQVdata = Check_SQV()
% This function calculates velocities for a squat from a table containing 
% vertical forces and time. Results will be visualized and files can be 
% selected for further processing.
%
% USAGE
% SQVdata = Check_SQV()
%
% Selected files must be ".csv" of force plate data, with vertical ground
% reaction forces stored in variables "z1-z4" and time in last column.
%
% OUTPUT
% SQVdata:  Cell array (3xn) with filenames, -paths, and data for each selected file will be
%           returned. Press button next to accept file or delete to drop it.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 11.September.2020
% lepremiere
%---------------------------------------------------------------------------------------------------
warning off
addpath(string(pwd) + '\Functions');
global f ii
f     = [];
files = loadFiles();            % Getting files via GUI. Requires external function 'uipickfiles.m'
h     = figure(1);

%% Looping through files
for ii = 1:size(files(:,1), 1)
    %% Determining total vertical ground reaction force (VGRF)
    % Note: Force plate data must be table with VGRF stored in columns named 'z1-z4' and time stored in last column. Time must be consistent
    data = files{ii,3};
    Fz = data.z1 + data.z2 + data.z3 + data.z4;

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
    % Note: Force plate data has to include a 3 second rest period
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
    noise     = Fz - Fz_filtered;                 
    threshold = range(noise(row:row + interval));   % Sets variability as base line noise

    %% Determining start of the movement                           
    posStart    = find(Fz_filtered(row:end) > Fz_filtered(row) + 2*threshold,1,'first') + row ;    % Determines start of movement by finding the first VGRF greater than baseline VGRF + 2*noise. 1000 ms offset 
    negStart    = find(Fz_filtered(row:end) < Fz_filtered(row) - 2*threshold,1,'first') + row ;    % Determines start of movement by finding the first VGRF greater than baseline VGRF - 2*noise. 1000 ms offset
    start       = min([posStart, negStart]);   
    BW          = mean(Fz_filtered(1:start))/9.81;                                                 % m = F/a

    %% Calculating variables of interest
    netForce        = Fz_filtered - BW*9.81 ;                  % F(net) = VGRF - BW*g
    acceleration    = Fz_filtered/BW-9.81;                     % a = F/m                    
    velocity        = cumtrapz(acceleration)/Fs;               % v = Integral(a)*dt                               
    [peak_velocity, index_pv]   = max(velocity);
    turningPoint    = find(velocity == min(velocity(1:find(velocity == peak_velocity))));
    startCon        = find(velocity(turningPoint:end) > 0 ,1, 'first') + turningPoint;
    endCon          = find(velocity == max(velocity(startCon:find(velocity == peak_velocity))));
    mpv_ind         = find(acceleration(startCon:endCon) > 0) + startCon;

    MPV = mean(velocity(mpv_ind));                             % Mean propulsive velocity
    MV  = mean(velocity(startCon:endCon));                     % Mean velocity
    PV  = peak_velocity;                                       % Peak velocity

    % Formulas from http://dx.doi.org/10.1055/s-0043-102933
    Intensity_MPV = -5.961*MPV^2 - 50.71*MPV + 117;
    Intensity_MV  = -12.87*MV^2  - 46.31*MV  + 116.3;
    Intensity_PV  = -10.85*PV^2  - 25.1*PV   + 130.3;

    power       = Fz_filtered.*velocity;                       % W = F*v
    peak_power  = max(power);

    turningPoint = find(Fz_filtered == min(Fz_filtered(start:endCon)),1,'first');   % Determines the start of increasing VGRF after counter-movement
    lift         = Fz_filtered(turningPoint:index_pv);
    rfd_interval = 50;                                                              % Sets the window width in ms for finding maximum rate of force development

    % Calculates the RFDs during lift with respect to time
    for i = 1:length(lift)-rfd_interval
       rFD(i) = ((lift(i + rfd_interval) - lift(i))/rfd_interval)*Fs;
    end

    % Finds max RFD and time to max RFD
    RFD         = max(rFD);
    [col1 row1] = find(rFD == RFD);                                                                                                                                         
    time2maxRFD = row1;
    %% Plot

    clf;
    h = plot(Fz_filtered);                                     % Force data
    hold on;
    h1 = area(mpv_ind, Fz_filtered(mpv_ind), BW*9.81);         % Propulsion area
    yl = ylim;
    
    % Title
    titl = sprintf('File: %s \nMPV: %.2f m/s, 1-RM: %.0f, PV: %.2f m/s, 1-RM: %.0f, PRFD: %.0f N/s', string(files(ii,2)), MPV, Intensity_MPV, PV, Intensity_PV, RFD);
    title(titl, 'Interpreter', 'none');   
    ylabel('Force (N)');
    xlabel('Time (ms)');

    % Line markers for different events
    xline(start,'.',{'Start'},'LabelHorizontalAlignment','center','FontWeight','bold');                                     % Start of counter-movement
    xline(turningPoint,'',{'Turning Point'},'LabelHorizontalAlignment','center','FontWeight','bold');                       % Begin of breaking phase
    xline(startCon,'',{'Start Concentric'},'LabelHorizontalAlignment','center','FontWeight','bold');                        % Start of concentric phase of the lift
    xline(find(velocity == peak_velocity),'',{'Peak Velocity'},'LabelHorizontalAlignment','center','FontWeight','bold');    % Peak velocity
    xline(find(power == peak_power),'',{'Peak Power'},'LabelHorizontalAlignment','center','FontWeight','bold');             % Peak power

    yyaxis right;
    v = plot(velocity, 'r');                                   % Velocity
    
    % Aligning the two y-axes
    yl = [(yl(1) - BW*9.81)/BW*9.81 (yl(2) - BW*9.81)/BW*9.81]./(2*9.81);
    ylim (yl);
    xl = [1 length(Fz_filtered)];
    xlim (xl);
    yline(0);

    % Legend
    legend([h, v, h1], {'VGRF', 'Velocity', 'Propulsion'}, 'Location', 'southwest');

    % Callbacks to control loop and figure behavior. Selecting files.
    c = uicontrol('Position',[10 80 75 50],'String', 'Accept',  'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'g','Callback','uiresume(gcf);');
    d = uicontrol('Position',[10 5 75 50], 'String', 'Discard', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'r', 'Callback','global f ii; f = [f,ii];uiresume()');
    uiwait();
    if(~ishandle(h) == 1)
       error('Figure has been closed')     
    end
end

% Returning selected files
files(f,:) = [];
clear('global')
SQVdata = files;
end

