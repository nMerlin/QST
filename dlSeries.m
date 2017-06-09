function [ dataStruct ] = dlSeries( varargin )
%DLSERIES Batch processing of diode laser current series
%   Detailed explanation goes here

% Optional input arguments
verbose = 0;
quiet = 'notquiet';
if nargin > 0
    for i = 1:nargin
        eval([varargin{i} '=1;']);
    end
end
if verbose == 0
    quiet = 'quiet';
end

%%% Variables
LOStruct = struct('filename',{},'number',{});
dataStruct = struct('filename',{});

%%% Create data overview
dispstat('','init',quiet);
dispstat('Checking filenames ...','timestamp','keepthis',quiet);
rawDataContents = dir('raw-data');

% Loop over files
for name = {rawDataContents.name}
    % Ignore unnecessary files
    filename = cell2mat(name);
    if strcmp(filename,'.') || strcmp(filename,'..')
        continue
    end
    if ~isempty(regexpi(filename,'.raw.','match'))
        continue
    end
    
    % Fetch number of measurement
    numberToken = regexpi(filename,'^([0123456789]*)-','tokens');
    number = str2double(cell2mat(numberToken{1}));
    
    % Distinguish between LOonly and other files
    if ~isempty(regexpi(filename,'LOonly','match'))
        LOStruct(number).filename = filename;
        LOStruct(number).number = number;
    else
        dataStruct(number).filename = filename;
        
        % Get LO power
        powerToken = regexpi(filename,'-([0123456789.]*)mW','tokens');
        dataStruct(number).powerLO = str2double(cell2mat(powerToken{1}));
    end
end % Loop over files
LOnumbers = cell2mat({LOStruct.number});

for iStruct = 1:length(dataStruct)
    % Process only data files
     filenameSIG = dataStruct(iStruct).filename;
    if isempty(filenameSIG)
        continue
    end
    dispstat(['Processing ',filenameSIG,' ...'], ...
        'timestamp','keepthis',quiet);
    
    %%% Load data
    try
        filename = regexprep(filenameSIG,'.raw','.mat');
        load(filename);
    catch
        LOnumber = max(LOnumbers(LOnumbers<=iStruct));
        filenameLO = LOStruct(LOnumber).filename;
        dispstat(['preparePhAvData number ',num2str(iStruct),' of ', ...
            num2str(length(dataStruct))]);
        X = preparePhAvData(filenameLO, filenameSIG);
    end % try
    
    % Compute and plot histograms
    filename = regexprep(filenameSIG,'.raw','');
    [nAv, width, dip, distTherm, distCoh] = plotPhAvData(X, filename);
    
    % Save quadratures and results
    save([filename,'.mat'],'X','nAv','width','dip','distTherm','distCoh');
end % iStruct

end % function

