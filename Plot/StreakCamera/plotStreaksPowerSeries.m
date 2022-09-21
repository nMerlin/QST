function plotStreaksPowerSeries(varargin)
% for a series of streak measurements located in folder 'raw-data' as .dat
% files.
%   Input Arguments:
%       'Plottype','lin': linear time- and intensity scale for integrated plot
%       'Plottype','log': logarithmic intensity scale for integrated plot
%       'Subtract','yes': subtract background data

%% Validate and parse input arguments
parser = inputParser;
defaultSubtract = 'yes'; 
addParameter(parser,'Subtract',defaultSubtract);
defaultPlottype = 'lin'; 
addParameter(parser,'Plottype',defaultPlottype);
defaultIntegrationArea = 'full'; 
addParameter(parser,'IntegrationArea',defaultIntegrationArea);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[integrationArea,plottype,subtract] = c{:};

%% Create data overview
dataStruct = struct('filename',{},'Power',{},'decaytime',{},'decaytimeError',{}, 'Max', {}, 'Sum', {});
dataStructBackground = struct('filename',{},'number',{});

rawDataContents = dir('raw-data');

%% Signal-files
for name = {rawDataContents.name}
    % Loop only over Signal files
    filename = cell2mat(name);
    if not(isempty(regexpi(filename,'background.img','match')))...
            || isempty(regexpi(filename,'.img','match'))
        continue
    end
    
    % Fetch number of measurement
    numberToken = regexpi(filename,'^([0123456789]*)-','tokens');
    number = str2double(cell2mat(numberToken{1}));
    dataStruct(number).filename = filename;
    
     % Fetch excitation power
    %powerToken = regexpi(filename,'-([0123456789.]*)mW','tokens'); %mira
    %powerToken = regexpi(filename,'MIRA-([0123456789.]*)mW','tokens'); %mira
    powerToken = regexpi(filename,'filter-([0123456789.]*)mW','tokens'); %mira
    power = str2double(cell2mat(powerToken{1}));
    dataStruct(number).Power = power;
    
end

%% Background-files
for name = {rawDataContents.name}
    % Loop only over background files
    filename = cell2mat(name);
    if isempty(regexpi(filename,'background.img','match'))
        continue
    end
    
    % Fetch number of measurement
    numberToken = regexpi(filename,'^([0123456789]*)-','tokens');
    number = str2double(cell2mat(numberToken{1}));
    dataStructBackground(number).filename = filename;
    dataStructBackground(number).number = number;
   
end
BGnumbers = cell2mat({dataStructBackground.number});

%% process the data
for number = 1:size(dataStruct,2)

    filenameSIG = dataStruct(number).filename;
    if isempty(filenameSIG)
        continue
    end    
   
   
    %find adequate BG-file, if this option is chosen
    if strcmp(subtract, 'yes')
        BGnumber = min(BGnumbers(BGnumbers>=number)); %background was measured after signal
        filenameBG = dataStructBackground(BGnumber).filename;
    else
        filenameBG = filenameSIG;
    end

    [decaytime, decaytimeError, Max, Sum] = plotStreak( filenameSIG, filenameBG,...
        'Plottype',plottype,'Subtract',subtract,'IntegrationArea',integrationArea );
     
    dataStruct(number).decaytime = decaytime; 
    dataStruct(number).decayError = decaytimeError;
    dataStruct(number).Max = Max;
    dataStruct(number).Sum = Sum;
end

%% write them in excel table
T = struct2table(dataStruct);
writetable(T,['powerSeries-subtract-' subtract '-IntArea-' integrationArea '.xls']);

end