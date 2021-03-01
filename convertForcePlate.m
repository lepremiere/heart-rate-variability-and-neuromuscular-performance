function convertForcePlate()
% This function converts ASCII with comma as decimal indicator to dot.
%
% USAGE
% E.g.:
% convertForcePlate()
%
% Selected files must be ".csv".
% If the selected output folder is the original location of the files, they will be replaced. 
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 11.September.2020
% lepremiere
%---------------------------------------------------------------------------------------------------
addpath(string(pwd) + '\Functions'); 

% Loads files, finds file name and sets output folder 
files   = loadFiles ();
st      = strfind(files{1,1},'\');
folder  = extractBetween(files{1,1},1,st(end));
selpath = uigetdir();
n = length(files(:,1));

    % Replaces comma with dot for decimal seperation and writes new tables to
    % output folder
    parfor i = 1:n
        data = convertGermanBullshit((files{i,1}));
        dir = join(['\',files{i,2},'.csv']);
        dir = char(dir);
        dir = dir(find(~isspace(dir)));
        writetable(data,[selpath,dir]);
        fprintf('File: %.0f/%.0f \n',i, n);
    end
end

function newTable = convertGermanBullshit(oldTable)

% Set input data type to 'char' to equalize data type and be able to apply
% strrep
opts                       = detectImportOptions(oldTable);
opts                       = setvartype(opts,'char');
opts.PreserveVariableNames = true;

% import and replacement of targets. doesnt work for tables, therefore
% conversion to char array then from str to double and back from array to
% table
x           = readtable(oldTable,opts);
varNames    = x.Properties.VariableNames;
x           = str2double(strrep(table2array(x),',','.'));
newTable    = array2table(x,'VariableNames',varNames);
end


%% Get files
function Output = loadFiles()
    file = uipickfiles;
    file = string(file);
    %% Load Files

    for i = 1:length(file) 
        [~,name,ext] = fileparts(file(i));
        Output{i,1} = file(i);
        Output{i,2} = name; 
        Output{i,3} = readtable(file(i),'ReadVariableNames',true, 'PreserveVariableNames', true);
    end
end