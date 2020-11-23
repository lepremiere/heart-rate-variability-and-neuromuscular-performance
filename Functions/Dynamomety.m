function Results = Dynamomety(data, test_occasions, n)

% This function calculates qualities of maximum voluntary contractions
% paired with superimposed and resting electrical nerve stimulations.
%
% USAGE
% E.g.:
% Results = Dynamomety(dyno_data, {'XY_INT_1', 'XY_INT_2'}, 1)
%
% INPUT
% data:             Data must be a struct containing the fields 'Torque' with subfields 
%                   'Torque.values' (n x 1) where the acutal torque data is stored and 
%                   'Torque.interval' (scalar) which indicates the sampling frequency. 
%                   Further, it must contain field 'Keyboard' with subfield 
%                   'Keyboard.times' (n x 1) that indicates the time of key
%                   presses. Each key press is interpreted as electrical nerve
%                   stimulation. Always two stimulations per contraction expected.
%                   First stimualtion interpreted as 'superimposed', second as
%                   'resting' twitch. Both peak to peak torques are calculated with
%                   respect to the torque right before resting stimulation.
%
% test_occasions:   Must be a string array (n x 1) that contains names to
%                   map the trials of the data to. E.g. {'XY_INT_1', 'XY_INT_2'}. Trials 
%                   will be named accordingly 'XX_INT_1_1', 'XY_INT_1_2'.
%
% n:                Points at the position in test_occasions, which should
%                   be used for name mapping.
%
% OUTPUT
% Results:          Table (m + 3 x 4) containing every trial plus mean,
%                   coefficient of variation, and the best trial for 4 metrics: Maximum
%                   voluntary contraction torque, superimposed twitch torque, resting twitch
%                   torque and voluntary activation as percentage.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 09.July.2020
% lepremiere
%---------------------------------------------------------------------------------------------------

%% Filtering torque data 
% Note: Dynamometer data must be a struct including fields '.Torque' with data stored in '.Torque.values' and sampling frequency in '.Torque.interval' and '.Keyboard' with time of press in '.Keyboard.times' 
% Every keypress must indicate stimulation. Superimposed stimulations need to occur first. Two stimulations per repetition required (superimposed, resting).

Fs    = 1/data.Torque.interval;                         % Getting sampling frequency
Fc    = 20;                                             % Cut-off frequency
[B,A] = butter(3,(Fc/(Fs/2)/0.8022), 'low');            % Calculating coefficients of 3rd order lowpass butterworth filter

Torque      = filtfilt(B, A, data.Torque.values);       % Filtering torque with zero phase butterworth filter

%% Identifying stimulations and calculating peak to peak torque
Stim_number            = numel(data.Keyboard.times);                % Counting number of keypresses. Note: Stimulations must be induced by kepresses. No other keypresses allowed.
Stim_time              = round((fix(data.Keyboard.times*Fs)), 0);   % Getting the indices of stimulations
maxStimulationWindow   = 250;                                       % Setting the window following stimulation to analyse peak to peak torque in
minStimulationWindow   = 50;                                        % Setting the window following stimulation to analyse peak to peak torque in

% Looping through stim times and allocating superimposed and resting stimulation torque as well as maximum voluntary torque (MVC) and voluntary activation (VA)
% Note: Every 1st stimulation needs to be superimposed and every 2nd the resting stimulation. No extra keypresses allowed
for i = 1:length(Stim_time)
    
    % Every second stimulation will be assigned to resting condition
    switch rem(i,2)
        
    case  1
        
        twitchTorqueMVC(floor(i/2)+1,1) = range([min(Torque(Stim_time(i):Stim_time(i) + minStimulationWindow)),...              % Peak to peak torque within superimposed stimulation + window 
                                                 max(Torque(Stim_time(i):Stim_time(i) + maxStimulationWindow))]);    
        if( i <= length(Stim_time) - 1)
            
            MVC(floor(i/2)+1,1)         = range([max(Torque(Stim_time(i) - maxStimulationWindow*20:Stim_time(i))),...           % Peak to peak torque within stimulation - window*20 (~5s) 
                                                 min(Torque(Stim_time(i + 1):Stim_time(i + 1) + minStimulationWindow))]);
        end
        
    case  0
        
        twitchTorqueRest(i/2,1)         = range([min(Torque(Stim_time(i):Stim_time(i) + minStimulationWindow)),...              % Peak to peak torque within resting stimulation + window
                                                 max(Torque(Stim_time(i):Stim_time(i) + maxStimulationWindow))]); 
                                            
        VA(i/2,1)                       = (1 - (twitchTorqueMVC(i/2)*(Torque(Stim_time(i-1))/MVC(i/2))/twitchTorqueRest(i/2)))*100;   % VA = (1 - (superimposed twitch*(Tb/MVC)*resting twitch^-1)) *100.
                                                                                                                                      % Tb = torque at stiumaltion. Accounts for not stimulation at true maximum torque
    end
end

Results = array2table([MVC twitchTorqueMVC twitchTorqueRest VA ]);      

% Assigning row names for all trials with respect to test occasion
for i = 1:Stim_number/2
    RowNames{i,1} = [test_occasions{n}, '_', num2str(i)];
end

Results.Properties.RowNames = RowNames;
Results                     = sortrows(Results, 1, 'descend');          % Sorting results for specific column. Preset = 1, torque
Height                      = height(Results);

% Calculating mean and standard deviation of result columns
for i = 1:width(Results)
    
    Results{Height+1,i} = mean(Results{1:Height,i},1);
    Results{Height+2,i} = std(Results{1:Height,i},1)/Results{Height+1,i}*100;
    
    % Gets best result of columns
    % Note: Flipped evaluation logic for the superimposed twitch
    if(i == 2)
        Results{Height+3,i} = min(Results{1:Height,i});
    else
        Results{Height+3,i} = max(Results{1:Height,i}); 
    end
end

Results.Properties.RowNames(end-2:end)  = {'Mean','CV','Best'};

end