function Output  = loadFiles()

% This function loads data using an user interface.
%
% USAGE
% E.g.:
% Output  = loadFiles()
%
% INPUT
% Function opens user interface from which inputs can be selected.
%
% OUTPUT
% Output:   Returns cell array (n x 3) with filepath, filename 
%           and actual data for all selected files. 
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 09.July.2020
% lepremiere
%---------------------------------------------------------------------------------------------------


%% Load files
file = uipickfiles;     % Gets filepath from uipickfiles function.
                        % Function 'loadfiles' requires 'uipickfiles.m' (Copyright (c) 2007, Douglas M. Schwarz) --> https://www.mathworks.com/matlabcentral/fileexchange/10867-uipickfiles-uigetfile-on-steroids
file = string(file);    

%% Sort files
% Seperates name and extension of filepath and creates 'Ouput' with filepath, filename and the actual data 
for i = 1:length(file)
    
    [~,name,ext] = fileparts(file(i));
    Output{i,1}  = file(i);
    Output{i,2}  = name; 
    
    % File type specific loading methods
    switch ext 
        case '.csv'
            w           = importdata(file(i), ',', 1);
            Output{i,3} = array2table(w.data, 'VariableNames', w.colheaders);
        case '.txt'
            Output{i,3} = array2table(load(file(i)));
        case '.mat'
            Output{i,3} = load(file(i));   
    end
    
end

end
