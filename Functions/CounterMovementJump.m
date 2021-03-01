function Results = CounterMovementJump(data)

% This function calculates jump qualities for a counter-movement jump from 
% a table containing vertical forces and time.
%
% USAGE
% E.g.:
% Results = CounterMovementJump(jump_data)
%
% INPUT
% data:     Data needs to be a table (dimensions: m > 3000, n > 5),
%           including variables z1-z4 with vertical ground reaction forces
%           and time stored in the last column. Additionally, the data must
%           contain a non-contact phase to determine liftoff. A 3 second 
%           resting period in the beginnging of each recording is required.
%
% OUTPUT
% Results:  Array (1x11) with outputs containing: body weight, jump height
%           calculated from airtime, jump height calculated from net impulse, peak
%           force, peak force normalized to bodyweight, peak power, take off
%           velocity, peak velocity, peak rate of force development, time to peak
%           rate of force development
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 09.July.2020
% lepremiere
%---------------------------------------------------------------------------------------------------

%% Determining total vertical ground reaction force (VGRF)
% Note: Force plate data must be table with VGRF stored in columns with names containing "z" and time stored in last column. Time must be consistent
inds    = contains(data.Properties.VariableNames, 'z'); % Finding columns with z-forces
Fz      = sum(data{:,inds}, 2);

% Determines the the sampling frequency 'Fs'
Fs = 1/(data.(data.Properties.VariableNames{end})(11) - data.(data.Properties.VariableNames{end})(10));

%% Filtering the data
Fc          = 50;                                   % Cutt-off frequency
[b,a]       = butter(3,(Fc/(Fs/2)/0.8022),'low');   % Calculating coefficients of second order lowpass butterworth filter. Adjusted by 0.8022
Fz_filtered = filtfilt(b,a,Fz);                     % Filtering VGRF with zero phase lowpass butterworth filter 

%% Adjusting VGRF in case force plate was relocated
if (mean(Fz_filtered) < 250)
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
BW        = mean(Fz_filtered(1:2000)./9.81);    % m = F/a  
threshold = range(Fz(row:row + interval));      % Sets variability as base line noise

%% Determining start, takeoff, landing and airtime of CMJ
% Getting orientation in the jump
subzero     = find(Fz_filtered < 0);                                                                          % Getting indices of no contact to force plate
takeoff     = find(Fz_filtered(row:subzero(1)) < threshold,1,'first')+row;                                    % Determines takeoff as first VGRF smaller than baseline noise. +row accounts for subzero offset
landing     = find(Fz_filtered(takeoff:end) > threshold,1,'first');                                           % Determines landing as first VGRF greater than baseline noise after takeoff
airtime     = landing * 1/Fs;                                                                                 % Landing obtains number of sampling points for the jump. By multiplying with sampling frequency, airtime is calculated
posStart    = find(Fz_filtered(row:subzero(1) - 1) > Fz_filtered(row) + threshold*2,1,'first') + row - 500;   % Determines start of movement by finding the first VGRF greater than baseline VGRF + noise. 500 ms offset 
negStart    = find(Fz_filtered(row:subzero(1) - 1) < Fz_filtered(row) - threshold*2,1,'first') + row - 500;   % Determines start of movement by finding the first VGRF greater than baseline VGRF - noise. 500 ms offset
start       = min([posStart,negStart]);                                                                       % Finds true start

%% Calculatng net force and net impulse 
netForce     = Fz_filtered - BW*9.81 ;                  % F(net) = VGRF - BW*g
netImpulse   = trapz((netForce(start:takeoff)))*1/Fs;   % Integrating netForce with respect to time

%% Calculating variables of interest
jumpHeightAirtime   = (1/8*9.81*airtime^2)*100;                                             % h = 1/8*g*t^2. In cm (100)
velocity            = cumtrapz((netForce)*1/Fs)/BW;                                         % v = delta(p)/m             
takeoffVelocity     = velocity(takeoff);
jumpHeightImpulse   = (takeoffVelocity^2/2*9.81);                                           % h = v^2/2g. In cm (100)
Power               = Fz_filtered.*velocity;                                                % W = F*v 
peakVelocity        = max (velocity(start:takeoff)); 
peakPower           = max(Power(start:takeoff));
peakForce           = max(Fz_filtered(start:takeoff));
relpeakForce        = max(Fz_filtered(start:takeoff))/BW;

turningPoint        = find(Fz_filtered == min(Fz_filtered(start:takeoff-250)),1,'first');   % Determines the start of increasing VGRF after counter-movement
jump                = Fz_filtered(turningPoint:takeoff);
rfd_interval        = 50;                                                                   % Sets the window width in ms for finding maximum rate of force development

% Calculates the RFDs during jump with respect to time
for i = 1:length(jump)-rfd_interval
    if(jump(i) > BW*9.81)
        rFD(i) = ((jump(i+rfd_interval)-jump(i))/rfd_interval)*Fs;
    else
        rFD(i) = 0;
    end
end

% Finds max RFD and time to max RFD
RFD         = max(rFD);
[col1 row1] = find(rFD == RFD);                                                                                                                                         
time2maxRFD = row1;

% Concatenates all variables to ouput 'Results'
Results = [jumpHeightAirtime, jumpHeightImpulse, netImpulse, peakForce, relpeakForce, peakPower, takeoffVelocity, peakVelocity, RFD, time2maxRFD];   

end
