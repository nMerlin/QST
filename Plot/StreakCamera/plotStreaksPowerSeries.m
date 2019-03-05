function plotStreaksPowerSeries(varargin)
% for a series of streak measurements located in folder 'raw-data' as .dat
% files.
%   Input Arguments:
%       'Plottype','lin': linear time- and intensity scale for integrated plot
%       'Plottype','log': logarithmic intensity scale for integrated plot

%% Validate and parse input arguments
parser = inputParser;
defaultPlottype = 'lin'; 
addParameter(parser,'Plottype',defaultPlottype);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[plottype] = c{:};

%% Create data overview
dataStruct = struct('filename',{},'Power',{},'decaytime',{});
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
    powerToken = regexpi(filename,'-([0123456789.]*)mW','tokens'); %mira
    %powerToken = regexpi(filename,'mW-([0123456789]*)mW','tokens'); %opo
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
   
    %find adequate BG-file
    BGnumber = min(BGnumbers(BGnumbers>=number)); %background was measured after signal
    filenameBG = dataStructBackground(BGnumber).filename;

    [decaytime, ~] = plotStreak( filenameSIG, filenameBG,'Plottype',plottype );
     
    dataStruct(number).decaytime = decaytime; 
       
end

%% write them in excel table
T = struct2table(dataStruct);
writetable(T,'powerSeries.xls');

end