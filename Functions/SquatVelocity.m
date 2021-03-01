function Results = SquatVelocity(data)
% This function calculates lift qualities for a loaded backsquat from 
% a table containing vertical forces and time.
%
% USAGE
% E.g.:
% Results = SquatVelocityp(data)
%
% INPUT
% data:     Data needs to be a table (dimensions: m > 3000, n >= 5),
%           including variables containing "z" with vertical ground reaction forces
%           and time stored in the last column. A 3 second resting period in
%           the beginnging of each recording is required.
%
% OUTPUT
% Results:  Array (1x8) with outputs containing: total weight, mean
%           propulsive velocity, mean velocity, peak velocity, intensity 
%           estimates for each velocity and peak power
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 19.January.2021
% lepremiere
%---------------------------------------------------------------------------------------------------

%% Determining total vertical ground reaction force (VGRF)
% Note: Force plate data must be table with VGRF stored in columns with names containing
% 'z' and time stored in the last column. Time must be consistent.
inds = contains(data.Properties.VariableNames, 'z'); % Finding columns with z-forces
Fz   = sum(data{:,inds}, 2);

% Cutting VGRF to the first peak. This reduces drift in case bar oscilation is present
if(std(Fz(1:300)) > 10)
    [~, ind1] = min(Fz(1:300));
    [~, ind2] = max(Fz(1:300));  
    idx = max([ind1, ind2]);
    Fz = Fz(idx:end);
end

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
interval = 1000;     % Window width in ms

% Calculating mean and standard deviation for n-windows for the first 3000 ms
for i = 1:3000 - interval
    if(Fz_filtered(i,1) > 10)
        q(i,1) = std(Fz_filtered(i:i + interval));
        q(i,2) = mean(Fz_filtered(i:i + interval));
    else
        q(i,1) = inf;                       
        q(i,2) = inf;
    end
end

% Finding window with smallest variability (std), assuming to be the best baseline and calculating bodyweight from its mean VGRF
[row, ~]   = find(q(:,1) == min(q(:,1)));                 
threshold  = q(row,1);                       % Sets variability as base line noise

%% Determining start of the movement                           
posStart    = find(Fz_filtered(row+interval:end) > q(row,2) + 3*threshold,1,'first') + row + interval - 1;    % Determines start of movement by finding the first VGRF greater than baseline VGRF + 3*noise.  
negStart    = find(Fz_filtered(row+interval:end) < q(row,2) - 3*threshold,1,'first') + row + interval - 1;    % Determines start of movement by finding the first VGRF greater than baseline VGRF - 3*noise. 
start       = min([posStart, negStart])-250;                                                                  % Sets the true start. 250 ms offset
BW          = find_BW(Fz_filtered(1:start), Fs)/9.81;                                                         % Determines the total weight of participant and load. m = F/a                                                                         

%% Calculating variables of interest
netForce        = Fz_filtered - BW*9.81;                  % F(net) = VGRF - BW*g
acceleration    = [netForce/BW];                          % a = F/m                          
velocity        = cumtrapz(acceleration)/Fs;              % v = Integral(a)*dt                               
[peak_velocity, index_pv] = max(velocity);
turningPoint    = find(velocity == min(velocity(1:find(velocity == peak_velocity))));
startCon        = find(velocity(turningPoint:end) > 0 ,1, 'first') + turningPoint;
endCon          = find(velocity == max(velocity(startCon:find(velocity == peak_velocity))));
mpv_ind         = find(acceleration(startCon:endCon) > 0) + startCon;

MPV = mean(velocity(mpv_ind));                             % Mean propulsive velocity
MV  = mean(velocity(startCon:endCon));                     % Mean velocity
PV  = peak_velocity;                                       % Peak velocity

% Formulas from http://dx.doi.org/10.1055/s-0043-102933
% Intensity_MPV = -5.961*MPV^2 - 50.71*MPV + 117;
% Intensity_MV  = -12.87*MV^2  - 46.31*MV  + 116.3;
% Intensity_PV  = -10.85*PV^2  - 25.1*PV   + 130.3;

power       = Fz_filtered.*velocity;                       % W = F*v
peak_power  = max(power);

turningPoint = find(Fz_filtered == min(Fz_filtered(start:endCon)),1,'first');   % Determines the start of increasing VGRF after counter-movement
lift         = Fz_filtered(turningPoint:index_pv);
rfd_interval = 50;                                                              % Sets the window width in ms for finding maximum rate of force development

% Calculates the RFDs during lift with respect to time
for i = 1:length(lift)-rfd_interval
   if(lift(i) > BW*9.81)
       rFD(i) = ((lift(i + rfd_interval) - lift(i))/rfd_interval)*Fs;
   else
       rFD(i) = 0;
   end
end

% Finds max RFD and time to max RFD
RFD         = max(rFD);
[~, row1]   = find(rFD == RFD);                                                                                                                                         
time2maxRFD = row1;

% Concatenates all variables to ouput 'Results'
Results = [MPV, MV, PV, peak_power, RFD];   

end
