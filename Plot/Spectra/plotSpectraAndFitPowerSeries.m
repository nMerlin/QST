function plotSpectraAndFitPowerSeries(varargin)
% this script makes spectrum plots for a series of spectrum measurements.
% The data should be in folder 'raw-data'.
% With 'Subtract', you can choose whether there are background measurements
% that should be subtracted. 

%% Validate and parse input arguments
parser = inputParser;
defaultSubtract = 'yes'; %
addParameter(parser,'Subtract',defaultSubtract);
parse(parser,varargin{:});
c = struct2cell(parser.Results);
[subtract] = c{:};

%% Create data overview
dataStruct = struct('filename',{},'number',{}, 'Power',{}, 'Max', {}, 'Integrated',{},...
    'peak', {}, 'FWHM', {}, 'Q', {});
dataStructBackground = struct('filename',{},'number',{});

rawDataContents = dir('raw-data');

%% Signal-files
for name = {rawDataContents.name}
    % Loop only over Signal files
    filename = cell2mat(name);
    if not(isempty(regexpi(filename,'background_1.txt','match')))...
            || isempty(regexpi(filename,'.txt','match'))
        continue
    end
    
    % Fetch number of measurement
    numberToken = regexpi(filename,'^([0123456789]*)-','tokens');
    number = str2double(cell2mat(numberToken{1}));
    dataStruct(number).filename = filename;
    
    % Fetch excitation power
    powerToken = regexpi(filename,'-([0123456789.]*)mW','tokens');
    power = str2double(cell2mat(powerToken{1}));
    dataStruct(number).Power = power;
    
end

%% Background-files
for name = {rawDataContents.name}
    % Loop only over background files
    filename = cell2mat(name);
    if isempty(regexpi(filename,'background_1.txt','match'))
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
        %[Max, peak, FWHM, Q] = 
        [Max,integratedInt, peak, FWHM, Q] = plotSpectrumAndFit( filenameSIG, filenameBG,...
            'Subtract','yes','Interpolate','yes','Fit','yes','Save','yes'); %%change this!
    else
        [Max,integratedInt, peak, FWHM, Q] = plotSpectrumAndFit( filenameSIG, filenameSIG, 'Subtract','no');
    end
    dataStruct(number).Max = Max;
    dataStruct(number).Integrated = integratedInt;
       
end

%% write them in excel table
T = struct2table(dataStruct);
writetable(T,'powerSeries.xls');
end