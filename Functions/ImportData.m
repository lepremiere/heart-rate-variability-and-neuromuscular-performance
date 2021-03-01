function data = ImportData()
% This function imports files containing counter-movement jump, loaded
% back squat, dynamometry and heart rate data.
%
% USAGE
% E.g.:
% data = ImportData()  
%
% INPUT
% Files will be selected via GUI. 
%
% OUTPUT
% data:     Struct (2:n x 1) with fields containing the files for the specific
%           test and the last field containing subject identifier found within
%           the files. The size of the struct depends on how many tests are
%           identified within the files. Only valid files will be processed.
%
%---------------------------------------------------------------------------------------------------
% Latest Edit: 26.January.2021
% lepremiere
%---------------------------------------------------------------------------------------------------
     
%% Getting data files

addpath(string(pwd) + '\Functions');                    % Adds the folder \Functions to the working folder to access distinct functions
% Loading target files with GUI explorer. 
try
    tic
    fprintf('Importing Data: \t')
    warning('off'); 
    input = loadFiles();
                
catch 
    if(input{1,1} == '0')
        error('ERROR no files selected (1)')
    end 
    
    error('ERROR getting data files (1)')
    
end

%% Extracting filenames and sorting for data type
% Possible datatypes: .csv = ForcePlate(CMJ); .mat = Dynamometer(MVC, VA); .txt = RR-Recording(HRV)
% Specification to required appereance and orginization of data can be found in the specific section 
try   
    files       = string(input(:, 1:2));
    
    % Exceptional case for only 1 file is given
    if(length(files(:,2)) == 1)
        number = {strfind(files(:, 2),'_')};           % Identifying information blocks of filename seperated by underscore
    else
        number = strfind(files(:, 2),'_');           
    end
    
    tests   = [];
    err         = 0;
    errCounter  = [];
    
    % Seperating files to an array of same data types
    for i = 1:length(files(:,1))        
        err(i) = 0;
        [~,name,ext] = fileparts(files(i,1));       % Extracting appendix of filename

        if (ext == '.csv' && ~contains(name, 'CMJ') == 0)                      
            CMJ(i,:) = input(i,:);
            CMJ(i, 2) = {CMJ{i, 2}{:}(1:number{i}(end) - 1)};
        elseif (ext == '.csv' && ~contains(name,  'SQV') == 0)  
            SQV(i, :) = input(i, :);
            SQV(i, 2) = {SQV{i, 2}{1}(1:number{i}(end) - 1)};
        elseif (ext == '.mat' && ~contains(name,  'DYNO') == 0)  
            DYNO(i,:) = input(i,:);
            DYNO(i, 2) = {DYNO{i, 2}{1}(1:number{i}(end) - 1)};
        elseif (ext == '.txt' && ~contains(name,  'HRV') == 0)  
            HRV(i,:) = input(i,:);
            HRV(i, 2) = {HRV{i, 2}{1}(1:number{i}(end) - 1)};
        else 
            err(i)             = 1;                 % Recognizing wrong file types
            errCounter(end+1)  = i;                 % Counting wrong file types
        end
        
        if(err(i) == 0)
           subjectPrefix{i,1} = files{i,2}(1:number{i}(1)-1);  % Extracting name prefix of file
        else
           subjectPrefix{i,1} = 'zzz999';           % Adds participant 'zzz999' as placeholder for wrong data types
        end        
    end
    
    subjectPrefix = unique(subjectPrefix);          % Identifying distinct name prefixes
    
    % Excludes participant 'zzz999' if wrong data types are present
    if(contains([subjectPrefix{:}],'zzz999'))
        subjectPrefix = subjectPrefix(1:end-1);
    else
        subjectPrefix = subjectPrefix;
    end
    
    % Empty array cells from sorting function above were removed.
    % Additionally, found methods will be saved and returned in 'variables'
    data = {};
    try
        intersection    = find(~cellfun('isempty', HRV(:,1)) == 1);
        HRV             = HRV(intersection,:);
        tests           = [tests;"HRV"];
        n               = length(data);
        data{n+1}       = HRV;
    catch
    end
    try
        intersection    = find(~cellfun('isempty', CMJ(:,1)) == 1);
        CMJ             = CMJ(intersection,:);
        tests           = [tests; "CMJ"];
        n               = length(data);
        data{n+1}       = CMJ;
    catch
    end
    try
        intersection    = find(~cellfun('isempty', SQV(:,1)) == 1);
        SQV             = SQV(intersection,:);
        tests           = [tests; "SQV"];
        n               = length(data);
        data{n+1}       = SQV;
    catch
    end
    try
        intersection    = find(~cellfun('isempty', DYNO(:,1)) == 1);
        DYNO            = DYNO(intersection,:);
        tests           = [tests; "DYNO"];
        n               = length(data);
        data{n+1}       = DYNO;
    catch
    end
        
    if(sum(err) > 0)
        fprintf('\nWrong file type. File ignored: %s',files(errCounter,1)); % Returning wrong files
        fprintf('\n\nOnly .csv, .mat and .txt allowed\n');                  % Returns allowed data types
    end
    
    data = cell2struct([data, {subjectPrefix}], [tests; "Participants"], 2);
    
catch    
    if(sum(err) > 0)
        error('ERROR no files with adequate data type found (.csv, .mat, .txt) (3)')
    end
    
    error('ERROR sorting files. More than one file required (4)')   
end
t1 = toc;
fprintf('%4.2f seconds\n', t1);
end