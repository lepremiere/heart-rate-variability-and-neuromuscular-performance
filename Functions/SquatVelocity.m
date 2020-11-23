function Results = SquatVelocity(data)

% This function calculates jump qualities for a counter-movement jump from 
% a table containing vertical forces and time.
%
% USAGE
% E.g.:
% Results = SquatVelocityp(jump_data)
%
% INPUT
% data:     Data needs to be a table (dimensions: m > 3000, n > 5),
%           including variables z1-z4 with vertical ground reaction forces
%           and time stored in the last column. A 3 second resting period in
%           the beginnging of each recording is required.
%
% OUTPUT
% Results:  Array (1x8) with outputs containing: total weight, mean
%           propulsive velocity, mean velocity, peak velocity, intensity 
%           estimates for each velocity and peak power
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 31.August.2020
% lepremiere
%---------------------------------------------------------------------------------------------------

%% Determining total vertical ground reaction force (VGRF)
% Note: Force plate data must be table with VGRF stored in columns named 'z1-z4' and time stored in last column. Time must be consistent
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

%% Determining start, takeoff, landing and airtime of CMJ                            
posStart    = find(Fz_filtered(row:end) > Fz_filtered(row) + 2*threshold,1,'first') + row ;    % Determines start of movement by finding the first VGRF greater than baseline VGRF + 2*noise. 1000 ms offset 
negStart    = find(Fz_filtered(row:end) < Fz_filtered(row) - 2*threshold,1,'first') + row ;    % Determines start of movement by finding the first VGRF greater than baseline VGRF - 2*noise. 1000 ms offset
start       = min([posStart, negStart]);   
BW          = mean(Fz_filtered(1:start))/9.81;                                                 % m = F/a

%% Calculatng net force and net impulse 
netForce     = Fz_filtered - BW*9.81 ;                  % F(net) = VGRF - BW*g
% netImpulse   = trapz((netForce(start:takeoff)))*1/Fs;   % Integrating netForce with respect to time

%% Calculating variables of interest
acceleration    = Fz_filtered/BW-9.81;                                         
velocity        = cumtrapz((netForce)*1/Fs)/BW;                                        
peak_velocity   = max(velocity);
turningPoint    = find(velocity == min(velocity(1:find(velocity == peak_velocity))));
startCon        = find(velocity(turningPoint:end) > 0 ,1, 'first') + turningPoint;
endCon          = find(velocity == max(velocity(startCon:find(velocity == peak_velocity))));
mpv_ind         = find(acceleration(startCon:endCon) > 0) + startCon;

MPV = mean(velocity(mpv_ind));
MV  = mean(velocity(startCon:endCon));
PV  = peak_velocity;

Intensity_MPV = -5.961*MPV^2 - 50.71*MPV + 117;
Intensity_MV  = -12.87*MV^2  - 46.31*MV  + 116.3;
Intensity_PV  = -10.85*PV^2  - 25.1*PV   + 130.3;

power       = Fz_filtered.*velocity;
peak_power  = max(power);

turningPoint = find(Fz_filtered == min(Fz_filtered(start:endCon)),1,'first');   % Determines the start of increasing VGRF after counter-movement
jump         = Fz_filtered(turningPoint:end);
rfd_interval = 50;                                                           % Sets the window width in ms for finding maximum rate of force development

% Calculates the RFDs during jump with respect to time
for i = 1:length(jump)-rfd_interval
   rFD(i) = ((jump(i + rfd_interval) - jump(i))/rfd_interval)*Fs;
end

% Finds max RFD and time to max RFD
RFD         = max(rFD);
[col1 row1] = find(rFD == RFD);                                                                                                                                         
time2maxRFD = row1;
% Concatenates all variables to ouput 'Results'
Results = [BW, MPV, Intensity_MPV, MV, Intensity_MV, PV, Intensity_PV, peak_power];   

end
